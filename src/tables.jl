# Tables.jl interfaces and DataAPI.jl metadata for TFI's result types (MomentArray, simulation
# and error panels, IRF results), and constructors that restore them from any Tables.jl source,
# e.g. a parquet file read with Parquet2.
#
# Maintainer notes
# - Every result type is a column table. Building the columns copies no data: data columns are
#   views of the stored arrays, and index columns (`__row`, `__run`, `__period`) are lazy
#   `mappedarray`s. IRF tables span several stored arrays (one per shock and statistic); their
#   columns are `mappedarray`s that look values up in those arrays.
# - Reserved column names: `__row` (MomentArray row names), `__run` and `__period` (panel
#   indices), `__shock` and `__stat` (IRF keys). Parquet has no Symbol type, so names are
#   written as Strings and turned back into Symbols when a table is read.
# - Table-level metadata (DataAPI.jl) has two parts. Structure: "tfi.kind" (the type that wrote
#   the table), "tfi.rowkey" (MomentArray's row-name column), "tfi.endog" (a SimulationPanel's
#   endogenous state names, comma-separated). Provenance: "simulate.<field>" for each field of a
#   panel's SimulationProvenance. Parquet2 writes this metadata into the file and exposes it on
#   reading.
# - Restoring checks "tfi.kind" when the source has it, so a table of one kind cannot be read as
#   another. Sources without metadata (CSV files, NamedTuples) restore too; a SimulationPanel then
#   needs its endogenous state names as an argument.

const ROW_KEY = :__row
const RUN_KEY = :__run
const PERIOD_KEY = :__period
const SHOCK_KEY = :__shock
const STAT_KEY = :__stat
const PROVENANCE_PREFIX = "simulate."

## Metadata helpers

# table-level metadata of any source, as Strings; empty if the source has none
function _source_metadata(table)
    DataAPI.metadatasupport(typeof(table)).read || return Dict{String,String}()
    Dict{String,String}(string(k) => string(DataAPI.metadata(table, k)) for k in DataAPI.metadatakeys(table))
end

# the source's metadata, after checking that it was written from a `kind` (if it says)
function _checked_metadata(table, kind::AbstractString)
    meta = _source_metadata(table)
    written = get(meta, "tfi.kind", kind)
    written == kind || throw(ArgumentError("table holds a $written, not a $kind"))
    meta
end

# a tuple of names as one metadata string and back
_encode_names(names) = join(names, ",")
_decode_names(s::AbstractString) = Tuple(Symbol(n) for n in eachsplit(s, ","; keepempty = false))

# provenance as metadata entries; fields that are `nothing` are left out
function _provenance_metadata(p::SimulationProvenance)
    meta = Dict{String,String}()
    for f in fieldnames(SimulationProvenance)
        v = getfield(p, f)
        isnothing(v) || (meta[PROVENANCE_PREFIX * string(f)] = string(v))
    end
    meta
end
_provenance_metadata(::Nothing) = Dict{String,String}()

# provenance from metadata entries; `nothing` if the source has none
function _provenance_from_metadata(meta)
    entry(f) = get(meta, PROVENANCE_PREFIX * string(f), nothing)
    isnothing(entry(:solution)) && return nothing
    SimulationProvenance(; solution = entry(:solution), tfi_version = entry(:tfi_version),
                         seed = _parse_or_nothing(Int, entry(:seed)), burn = parse(Int, entry(:burn)),
                         dynamic_align = Symbol(entry(:dynamic_align)), exog_start = entry(:exog_start),
                         endog_start = entry(:endog_start))
end
_parse_or_nothing(::Type, ::Nothing) = nothing
_parse_or_nothing(T::Type, s::AbstractString) = parse(T, s)

# the columns of a table other than the reserved ones, as a NamedTuple
function _value_columns(cols, reserved)
    names = Tuple(n for n in Tables.columnnames(cols) if !(n in reserved))
    NamedTuple{names}(map(n -> Tables.getcolumn(cols, n), names))
end

## MomentArray

Tables.istable(::Type{<:MomentArray}) = true
Tables.columnaccess(::Type{<:MomentArray}) = true
# `__row` (the row names as Strings), then a view of each data column
Tables.columns(m::MomentArray{ROWS,COLS}) where {ROWS,COLS} =
    NamedTuple{(ROW_KEY, COLS...)}((mappedarray(String, SVector(ROWS)),
                                    ntuple(j -> view(m.data, :, j), Val(length(COLS)))...))

tfi_metadata(::MomentArray) = Dict("tfi.kind" => "MomentArray", "tfi.rowkey" => string(ROW_KEY))

# row names, column names and data matrix of a table written from a MomentArray
function _moments_from_table(table, rowkey)
    meta = _checked_metadata(table, "MomentArray")
    rowkey = something(rowkey, Symbol(get(meta, "tfi.rowkey", string(ROW_KEY))))
    cols = Tables.columns(table)
    rowkey in Tables.columnnames(cols) || throw(ArgumentError("table has no row-name column $rowkey"))
    rows = Tuple(Symbol(r) for r in Tables.getcolumn(cols, rowkey))
    values = _value_columns(cols, (rowkey,))
    data = Matrix{mapreduce(eltype, promote_type, values)}(undef, length(rows), length(values))
    foreach((j, v) -> copyto!(view(data, :, j), v), 1:length(values), values)
    rows, keys(values), data
end

"""
    MomentArray(table; rowkey = nothing)
    MomentArray{ROWS,COLS}(table; rowkey = nothing)
    MomentArray{ROWS,COLS,T}(table; rowkey = nothing)
    MomentArray(table, rows, cols)

Restore a `MomentArray` from a table written from one, e.g. a parquet file. A table written
elsewhere needs a column of row names; `rowkey` names it. With `ROWS` and `COLS`, the table's
names are checked against them and the result's type is known at compile time; with `T`, its
element type too.
"""
function MomentArray(table; rowkey = nothing)
    rows, cols, data = _moments_from_table(table, rowkey)
    MomentArray{rows, cols}(data)
end

function MomentArray{ROWS,COLS}(table; rowkey = nothing) where {ROWS,COLS}
    rows, cols, data = _moments_from_table(table, rowkey)
    (rows == ROWS && cols == COLS) ||
        throw(ArgumentError("table has rows $rows and columns $cols, expected $ROWS and $COLS"))
    MomentArray{ROWS,COLS}(data)
end

MomentArray{ROWS,COLS,T}(data::AbstractMatrix) where {ROWS,COLS,T} =
    MomentArray{ROWS,COLS}(convert(AbstractMatrix{T}, data))
function MomentArray{ROWS,COLS,T}(table; rowkey = nothing) where {ROWS,COLS,T}
    cols = _moment_columns(table, rowkey, MomentArray{ROWS,COLS})
    data = Matrix{T}(undef, length(ROWS), length(COLS))
    foreach((j, v) -> copyto!(view(data, :, j), v), 1:length(COLS), cols)
    MomentArray{ROWS,COLS}(data)
end

"""
    MomentArray{ROWS,COLS,T,M}(table; rowkey = nothing)

Restore into the matrix type `M`, e.g. `SMatrix{3,2,Float64}`. A static matrix is filled
straight from the table's columns, without an intermediate `Matrix`.
"""
MomentArray{ROWS,COLS,T,M}(table; rowkey = nothing) where {ROWS,COLS,T,M<:AbstractMatrix{T}} =
    MomentArray{ROWS,COLS}(M(MomentArray{ROWS,COLS,T}(table; rowkey).data))
function MomentArray{ROWS,COLS,T,M}(table; rowkey = nothing) where {ROWS,COLS,T,M<:StaticMatrix{<:Any,<:Any,T}}
    cols = _moment_columns(table, rowkey, MomentArray{ROWS,COLS})
    R = length(ROWS)
    MomentArray{ROWS,COLS}(M(ntuple(k -> convert(T, cols[(k - 1) ÷ R + 1][(k - 1) % R + 1]), Val(R * length(COLS)))))
end

# the value columns of a table written from MomentArray{ROWS,COLS}, each read once (a parquet
# source decodes a column on every read), after checking the kind and the row and column names
# in place, without collecting the table's names. The names come as a type so that the
# columns' types are inferred.
function _moment_columns(table, rowkey, ::Type{MomentArray{ROWS,COLS}}) where {ROWS,COLS}
    kind = _metadata_entry(table, "tfi.kind", "MomentArray")
    kind == "MomentArray" || throw(ArgumentError("table holds a $kind, not a MomentArray"))
    rowkey = something(rowkey, Symbol(_metadata_entry(table, "tfi.rowkey", ROW_KEY)))
    cols = Tables.columns(table)
    _colnames_match(Tables.columnnames(cols), rowkey, COLS) ||
        throw(ArgumentError("table columns $(Tuple(Tables.columnnames(cols))) are not $rowkey and $COLS"))
    rownames = Tables.getcolumn(cols, rowkey)
    (length(rownames) == length(ROWS) && all(i -> Symbol(rownames[i]) === ROWS[i], eachindex(ROWS))) ||
        throw(ArgumentError("table rows $(Tuple(rownames)) are not $ROWS"))
    map(c -> Tables.getcolumn(cols, c), COLS)
end

# one metadata entry of any source; `default` if the source has no metadata or not this key
_metadata_entry(table, key, default) =
    DataAPI.metadatasupport(typeof(table)).read ? DataAPI.metadata(table, key, default) : default

# the names are `COLS` in order, plus `rowkey` anywhere
function _colnames_match(names, rowkey, COLS)
    length(names) == length(COLS) + 1 || return false
    j = 0
    for n in names
        n === rowkey && continue
        j += 1
        (j <= length(COLS) && n === COLS[j]) || return false
    end
    true
end

## Simulation and error panels

const _Panel = Union{SimulationPanel, SimulationErrorPanel}
Tables.istable(::Type{<:_Panel}) = true
Tables.columnaccess(::Type{<:_Panel}) = true
Tables.columns(sp::_Panel) = _panel_columns(sp.data)

# `__run` and `__period`, then each field as a vector. Table rows follow the array's memory
# order (periods within runs); the index columns map the array's CartesianIndices.
function _panel_columns(data::StructArray)
    ndims(data) <= 2 || throw(ArgumentError("panels with more than two dimensions are not tables"))
    idx = vec(CartesianIndices((size(data, 1), size(data, 2))))
    merge(NamedTuple{(RUN_KEY, PERIOD_KEY)}((mappedarray(_run_index, idx), mappedarray(_period_index, idx))),
          map(vec, components(data)))
end
_period_index(I::CartesianIndex) = I[1]
_run_index(I::CartesianIndex) = I[2]

tfi_metadata(sp::SimulationPanel) =
    merge(_provenance_metadata(sp.provenance),
          Dict("tfi.kind" => "SimulationPanel", "tfi.endog" => _encode_names(endog_states(sp))))
tfi_metadata(ep::SimulationErrorPanel) =
    merge(_provenance_metadata(ep.provenance), Dict("tfi.kind" => "SimulationErrorPanel"))

SimulationPanel(table) = _simulation_panel(table, nothing)
SimulationPanel(table, endog_names) = _simulation_panel(table, Tuple(endog_names))

function SimulationErrorPanel(table)
    meta = _checked_metadata(table, "SimulationErrorPanel")
    SimulationErrorPanel(_panel_data(table); provenance = _provenance_from_metadata(meta))
end

# both SimulationPanel table constructors: endogenous names from the argument, else the metadata
function _simulation_panel(table, endog)
    meta = _checked_metadata(table, "SimulationPanel")
    SimulationPanel(_panel_data(table), _endog_names(meta, endog); provenance = _provenance_from_metadata(meta))
end
_endog_names(meta, endog::Tuple) = endog
function _endog_names(meta, ::Nothing)
    haskey(meta, "tfi.endog") || throw(ArgumentError(
        "the table does not record the endogenous state names; supply them: SimulationPanel(table, endog_names)"))
    _decode_names(meta["tfi.endog"])
end

# periods × runs StructArray from a table with `__run` and `__period` columns
function _panel_data(table)
    cols = Tables.columns(table)
    names = Tables.columnnames(cols)
    (RUN_KEY in names && PERIOD_KEY in names) ||
        throw(ArgumentError("panel table needs $RUN_KEY and $PERIOD_KEY columns"))
    _unflatten_panel(Tables.getcolumn(cols, RUN_KEY), Tables.getcolumn(cols, PERIOD_KEY),
                     _value_columns(cols, (RUN_KEY, PERIOD_KEY)))
end

# values placed by their (run, period), so table rows may come in any order
function _unflatten_panel(run::AbstractVector, period::AbstractVector, values::NamedTuple)
    T, R = maximum(period), maximum(run)
    length(run) == T * R || throw(ArgumentError("panel table has $(length(run)) rows, expected $T periods × $R runs"))
    linear = LinearIndices((T, R))
    positions = map((r, t) -> linear[t, r], run, period)
    allunique(positions) || throw(ArgumentError("panel table has duplicate ($RUN_KEY, $PERIOD_KEY) rows"))
    StructArray(map(v -> (a = Matrix{eltype(v)}(undef, T, R); a[positions] = v; a), values))
end

## IRF results

Tables.istable(::Type{<:IRFResult}) = true
Tables.columnaccess(::Type{<:IRFResult}) = true
# `__shock`, `__stat` (mean, std, then other statistics), `__period`, then each variable. All
# lines have the same length, so the columns map the CartesianIndices of
# (period, shock, statistic) into the stored lines: rows by period within shock within statistic.
function Tables.columns(irf::IRFResult)
    lines = merge((mean = irf.means, std = irf.stdevs), irf.other_stats)
    shocks, stats = map(String, propertynames(irf.means)), map(String, keys(lines))
    blocks = Tuple(getproperty(getproperty(lines, Symbol(st)), Symbol(sh)) for st in stats for sh in shocks)
    T = length(first(blocks))
    all(b -> length(b) == T, blocks) || throw(ArgumentError("IRF lines differ in length"))
    idx = vec(CartesianIndices((T, length(shocks), length(stats))))
    block(I) = (I[3] - 1) * length(shocks) + I[2]
    vars = propertynames(first(blocks))
    values = map(v -> (cols = map(b -> getproperty(b, v), blocks); mappedarray(I -> cols[block(I)][I[1]], idx)),
                 NamedTuple{vars}(vars))
    merge(NamedTuple{(SHOCK_KEY, STAT_KEY, PERIOD_KEY)}((mappedarray(I -> shocks[I[2]], idx),
                                                        mappedarray(I -> stats[I[3]], idx),
                                                        mappedarray(_period_index, idx))), values)
end
tfi_metadata(::IRFResult) = Dict("tfi.kind" => "IRFResult")

"""
    irf_paths(irf::IRFResult)

The raw IRF draws as a table: `__shock`, `__run`, `__period`, then the saved path variables.
Empty if no paths were saved. Restore with `IRFResult(stats_table, paths_table)`.
"""
irf_paths(irf::IRFResult) = IRFPaths(irf)

# the paths table of an IRFResult, a separate type so that it gets its own Tables.jl methods
struct IRFPaths{I<:IRFResult}
    irf::I
end
Tables.istable(::Type{<:IRFPaths}) = true
Tables.columnaccess(::Type{<:IRFPaths}) = true
tfi_metadata(::IRFPaths) = Dict("tfi.kind" => "IRFPaths")

# All shocks' path panels have the same size, so the columns map the CartesianIndices of
# (period, run, shock) into the stored panels: rows by period within run within shock.
function Tables.columns(p::IRFPaths)
    saved = Tuple((String(shock), sa) for (shock, sa) in pairs(p.irf.paths) if sa isa StructArray)
    isempty(saved) && return NamedTuple{(SHOCK_KEY, RUN_KEY, PERIOD_KEY)}((String[], Int64[], Int64[]))
    shocks, panels = map(first, saved), map(last, saved)
    T, N = size(first(panels), 1), size(first(panels), 2)
    all(sa -> size(sa) == (T, N), panels) || throw(ArgumentError("IRF path panels differ in size"))
    idx = vec(CartesianIndices((T, N, length(panels))))
    vars = propertynames(first(panels))
    values = map(v -> (cols = map(sa -> getproperty(sa, v), panels); mappedarray(I -> cols[I[3]][I[1], I[2]], idx)),
                 NamedTuple{vars}(vars))
    merge(NamedTuple{(SHOCK_KEY, RUN_KEY, PERIOD_KEY)}((mappedarray(I -> shocks[I[3]], idx),
                                                       mappedarray(_run_index, idx),
                                                       mappedarray(_period_index, idx))), values)
end

IRFResult(::IsTable, stats_table) = IRFResult(IsTable(), stats_table, nothing)

function IRFResult(::IsTable, stats_table, paths_table)
    _checked_metadata(stats_table, "IRFResult")
    cols = Tables.columns(stats_table)
    order, groups = _row_groups(Tables.getcolumn(cols, STAT_KEY), Tables.getcolumn(cols, SHOCK_KEY))
    stats = Tuple(unique(Symbol(stat) for (stat, _) in order))
    shocks = Tuple(unique(Symbol(shock) for (_, shock) in order))
    period = Tables.getcolumn(cols, PERIOD_KEY)
    values = _value_columns(cols, (SHOCK_KEY, STAT_KEY, PERIOD_KEY))
    function line(stat, shock)
        rows = groups[(String(stat), String(shock))]
        rows = rows[sortperm(view(period, rows))]
        StructArray(map(v -> v[rows], values))
    end
    lines = NamedTuple{stats}(map(st -> NamedTuple{shocks}(map(sh -> line(st, sh), shocks)), stats))
    other = Base.structdiff(lines, NamedTuple{(:mean, :std)})
    IRFResult(lines.mean, lines.std, other, _irf_paths(paths_table, shocks))
end

_irf_paths(::Nothing, shocks) = NamedTuple{shocks}(map(_ -> (;), shocks))
function _irf_paths(paths_table, shocks)
    _checked_metadata(paths_table, "IRFPaths")
    cols = Tables.columns(paths_table)
    _, groups = _row_groups(Tables.getcolumn(cols, SHOCK_KEY))
    run, period = Tables.getcolumn(cols, RUN_KEY), Tables.getcolumn(cols, PERIOD_KEY)
    values = _value_columns(cols, (SHOCK_KEY, RUN_KEY, PERIOD_KEY))
    function shock_paths(shock)
        rows = get(groups, (String(shock),), nothing)
        isnothing(rows) ? (;) : _unflatten_panel(view(run, rows), view(period, rows), map(v -> view(v, rows), values))
    end
    NamedTuple{shocks}(map(shock_paths, shocks))
end

# row indices of each distinct combination of the key columns' values, in one pass; also the
# combinations in order of first appearance
function _row_groups(keycols::Vararg{AbstractVector,N}) where {N}
    groups = Dict{NTuple{N,String},Vector{Int}}()
    order = NTuple{N,String}[]
    for i in eachindex(first(keycols))
        key = map(c -> String(c[i]), keycols)
        rows = get!(groups, key) do
            push!(order, key)
            Int[]
        end
        push!(rows, i)
    end
    order, groups
end

## DataAPI.jl metadata (read-only) for all TFI tables

const _TFITable = Union{MomentArray, SimulationPanel, SimulationErrorPanel, IRFResult, IRFPaths}

DataAPI.metadatasupport(::Type{<:_TFITable}) = (read = true, write = false)
DataAPI.metadatakeys(x::_TFITable) = keys(tfi_metadata(x))
DataAPI.metadata(x::_TFITable; style::Bool = false) =
    style ? Dict(k => (v, :default) for (k, v) in tfi_metadata(x)) : tfi_metadata(x)
function DataAPI.metadata(x::_TFITable, key::AbstractString; style::Bool = false)
    v = tfi_metadata(x)[key]
    style ? (v, :default) : v
end
function DataAPI.metadata(x::_TFITable, key::AbstractString, default; style::Bool = false)
    v = get(tfi_metadata(x), key, default)
    style ? (v, :default) : v
end
