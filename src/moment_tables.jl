# Moments as a long table: the numbers that moments need, one row per economy and number, the same
# whether the results are in memory or saved with OutputVault. `MomentRows` reads moments from such
# a table, so `compare_economies` works on it unchanged.
#
# Maintainer notes
# - A row is (kind, group, variable, stat, period): kind "simulation" (group: the sample, period 0),
#   "impulse" (group: the shock) or "parameter" (group "parameters", stat "value", period 0). A leaf
#   moment needs one row, or two for an impulse moment whose transformation needs the basis path.
# - The OutputVault method of assemble_moments (in the extension) queries the same rows and must
#   return the same table; a test checks it.

const MOMENT_TABLE_COLUMNS = (:branch, :experdef, :variant, :kind, :group, :variable, :stat, :period, :value)

const MomentRow = NamedTuple{(:kind, :group, :variable, :stat, :period), Tuple{String, String, String, String, Int}}

# the rows a moment needs
_moment_rows(m::SimulationMoment) = MomentRow[("simulation", string(m.sample), string(m.var), string(m.stat), 0)]
function _moment_rows(m::ImpulseMoment)
    rows = MomentRow[("impulse", string(m.shock), string(m.var), string(m.stat), m.period)]
    if Transformations.requires_basis(transformation(m))
        push!(rows, ("impulse", string(m.basis), string(m.var), string(m.stat), m.period))
    end
    return rows
end
_moment_rows(m::ParameterMoment) = MomentRow[("parameter", "parameters", string(m.name), "value", 0)]
_moment_rows(m::CalculatedMoment) = unique(reduce(vcat, [_moment_rows(c) for c in values(m.inputs)]))

# a NamedTuple of moments, bare names read as SimulationMoments
_as_moments(moments::NamedTuple) = map(v -> v isa AbstractMoment ? v : SimulationMoment(v), moments)

# the distinct rows of all moments
_requested_rows(moments) = unique(reduce(vcat, [_moment_rows(m) for m in values(_as_moments(moments))]; init = MomentRow[]))

## Results in memory

"""
    assemble_moments(moments, economies::NamedTuple; branch = "memory", experdef = "default") -> NamedTuple
    assemble_moments(moments, runs)

The numbers that `moments` (a NamedTuple of [`AbstractMoment`](@ref)s or variable names) need, as
a long column table: `branch`, `experdef`, `variant`, `kind` (`"simulation"`, `"impulse"` or
`"parameter"`), `group` (the sample, the shock or `"parameters"`), `variable`, `stat`, `period`
(of an IRF; 0 otherwise) and `value`. Text columns are `String`s.

With `economies` (the results of each economy in memory, as [`compare_economies`](@ref) takes
them), every row gets `branch` and `experdef`, and the variant is the economy's name. With `runs`
(OutputVault `RunHandle`s or a `RunPattern`; needs `using OutputVault`) the numbers are read from
the saved runs, IRFs from the starting point `start` (keyword, default `:ss`). The same results
give the same rows either way. A number that a result lacks has no row.
"""
function assemble_moments(moments::NamedTuple, economies::NamedTuple; branch::AbstractString = "memory",
                          experdef::AbstractString = "default")
    rows = _requested_rows(moments)
    out = NamedTuple{MOMENT_TABLE_COLUMNS}((String[], String[], String[], String[], String[], String[], String[], Int[], Float64[]))
    for (name, source) in pairs(economies)
        for r in rows
            v = _row_value(source, r)
            if isnothing(v)
                continue
            end
            for (col, x) in zip(MOMENT_TABLE_COLUMNS, (branch, experdef, string(name), r.kind, r.group, r.variable, r.stat, r.period, v))
                push!(getproperty(out, col), x)
            end
        end
    end
    return out
end

# the value of a row in a source of results, or `nothing` if the source lacks it
function _row_value(r::EconomyResult, row)
    part = row.kind == "simulation" ? r.stats : row.kind == "impulse" ? r.irf : r.params
    return isnothing(part) ? nothing : _row_value(part, row)
end
function _row_value(s::SamplesTuple, row)
    if row.kind != "simulation" || !haskey(s, Symbol(row.group))
        return nothing
    end
    return _row_value(s[Symbol(row.group)], row)
end
function _row_value(m::MomentArray, row)
    var, stat = Symbol(row.variable), Symbol(row.stat)
    if !(var in _get_vars(m) && stat in _get_stats(m))
        return nothing
    end
    return Float64(getmoment(m, var, stat))
end
function _row_value(irf::IRFResult, row)
    if row.kind != "impulse"
        return nothing
    end
    lines = merge((mean = irf.means, std = irf.stdevs), irf.other_stats)
    stat, shock, var = Symbol(row.stat), Symbol(row.group), Symbol(row.variable)
    if !(stat in keys(lines) && shock in keys(lines[stat]))
        return nothing
    end
    line = lines[stat][shock]
    if !(var in propertynames(line) && 1 <= row.period <= length(line))
        return nothing
    end
    return Float64(getproperty(line, var)[row.period])
end
_row_value(params::NamedTuple, row) = row.kind == "parameter" && haskey(params, Symbol(row.variable)) ? Float64(params[Symbol(row.variable)]) : nothing
_row_value(source, row) = nothing

## Reading moments from a long table

"""
    MomentRows

The numbers of one economy from a long table of [`assemble_moments`](@ref), a source that
[`getmoment`](@ref) and [`compare_economies`](@ref) read every moment from.
"""
struct MomentRows
    values::Dict{MomentRow, Float64}
end

function momentvalues(d::MomentRows, m::Union{SimulationMoment, ImpulseMoment, ParameterMoment})
    values = map(_moment_rows(m)) do r
        if !haskey(d.values, r)
            throw(KeyError("the results have no $(r.kind) number $(r.group), $(r.variable), $(r.stat), period $(r.period)"))
        end
        d.values[r]
    end
    return first(values), length(values) == 2 ? values[2] : nothing
end

"""
    split_by_economy(table; order = nothing) -> NamedTuple

The long table of [`assemble_moments`](@ref) as one [`MomentRows`](@ref) per economy, keyed by
variant, the `economies` that [`compare_economies`](@ref) takes. Variants that repeat across
branches or experdefs are an error: select fewer runs. `order` lists the variants to keep, in that
order; without it, variants come in the order of the table.
"""
function split_by_economy(table; order = nothing)
    runs = unique(zip(table.branch, table.experdef, table.variant))
    variants = [r[3] for r in runs]
    if !allunique(variants)
        dup = unique(v for v in variants if count(==(v), variants) > 1)
        throw(ArgumentError("variants $(join(dup, ", ")) appear in several runs ($(join(("$(r[1])/res-$(r[2])-$(r[3])" for r in runs if r[3] in dup), ", "))); select runs with a narrower pattern"))
    end
    dicts = Dict(v => Dict{MomentRow, Float64}() for v in variants)
    for j in eachindex(table.value)
        row = MomentRow((table.kind[j], table.group[j], table.variable[j], table.stat[j], table.period[j]))
        d = dicts[table.variant[j]]
        if haskey(d, row)
            throw(ArgumentError("the table has two values of $row for variant $(table.variant[j])"))
        end
        d[row] = table.value[j]
    end
    names = isnothing(order) ? variants : String.(collect(order))
    for n in names
        if !haskey(dicts, n)
            dicts[n] = Dict{MomentRow, Float64}()
        end
    end
    return NamedTuple{Tuple(Symbol.(names))}(Tuple(MomentRows(dicts[n]) for n in names))
end
