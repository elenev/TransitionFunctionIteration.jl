# TFI's results with OutputVault: how TFI's types are stored (OutputVault does the I/O), the
# conveniences save_result, load_result and migrate_irf!, and moments read from saved runs.
#
# Storage
# - SimulationPanel, SimulationErrorPanel and MomentArray are Tables.jl tables with TFI's DataAPI
#   metadata, so OutputVault writes them as table files; from_storage restores them from the table.
# - An IRFResult is the folder (stats, paths): two tables, the statistics with the IRF's metadata.
# - ModelSolution, SteadyState and Interpolant are stored in their .jld2 files as their portable data
#   (to_portable, documented in src/portable.jl): one group per NamedTuple, one dataset per leaf,
#   plain HDF5-compatible arrays and strings that read without the code that wrote them. Files are
#   written uncompressed: h5py cannot read JLD2's compressed (chunked) datasets.
module TransitionFunctionIterationOutputVaultExt

using OutputVault
using OutputVault: JLD2, RunHandle, RunPattern, QueryShape, find_runs, query_results
import OutputVault: to_storage, from_storage, write_portable_jld2, read_portable_jld2
import TransitionFunctionIteration as TFI
using TransitionFunctionIteration: ModelSolution, SteadyState, MomentArray, SimulationPanel,
    SimulationErrorPanel, IRFResult, irf_paths, Interpolant, collect_parameters
using Tables

## Tables

from_storage(::Type{SimulationPanel}, data; context = nothing) = SimulationPanel(data)
from_storage(::Type{SimulationErrorPanel}, data; context = nothing) = SimulationErrorPanel(data)
from_storage(::Type{MomentArray}, data; context = nothing) = MomentArray(data)

# an IRF is a folder of two tables: statistics and raw paths. The statistics are a plain table with
# the IRF's metadata, since to_storage applies again to the folder's members.
to_storage(irf::IRFResult) = (stats = OutputVault.Formats.StoredTable(Tables.columntable(irf), TFI._source_metadata(irf)),
                              paths = irf_paths(irf))
function from_storage(::Type{IRFResult}, data; context = nothing)
    paths = get(data, :paths, nothing)
    return isnothing(paths) ? IRFResult(data.stats) : IRFResult(data.stats, paths)
end

## Portable .jld2 files: the portable data as groups and datasets

const PORTABLE = Union{ModelSolution, SteadyState, Interpolant}

write_portable_jld2(io, x::PORTABLE) = write_tree(io, "", TFI.to_portable(x))

# `context` is the model or model type of the load
read_portable_jld2(T::Type{<:PORTABLE}, io; context = nothing) = TFI.from_portable(T, read_group(io); model = context)

function write_tree(io, prefix, nt::NamedTuple)
    for (k, v) in pairs(nt)
        key = isempty(prefix) ? string(k) : "$prefix/$k"
        if v isa NamedTuple
            write_tree(io, key, v)
        else
            io[key] = v
        end
    end
end

read_group(g) = NamedTuple(Symbol(k) => (v = g[k]; v isa JLD2.Group ? read_group(v) : v) for k in keys(g))

## Saving and loading conveniences

function TFI.save_result(rf::RunHandle, result; start = :ss, dirty_repo_check::Symbol = :warn)
    has(p) = hasproperty(result, p) && !isnothing(getproperty(result, p))
    if has(:obj)
        solution = has(:stv) ? (; result.obj, result.stv) : (; result.obj)
        dir = save_solution(rf, solution; dirty_repo_check)
        # the new solution was left in staging (see the warning): keep the old run untouched
        if dir != OutputVault.run_path(rf)
            return dir
        end
        save_result_item(rf, "params/params", collect_parameters(result.obj))
    end
    # one item per panel: a NamedTuple of two panels of equal length would be one table
    for p in (:simseries, :errs)
        if has(p)
            save_result_item(rf, "simulation/$p", getproperty(result, p))
        end
    end
    if has(:errstats)
        save_result_item(rf, "errstats/errstats", result.errstats)
    end
    if has(:stats)
        save_result_item(rf, "stats", result.stats)
    end
    if has(:irf)
        save_result_item(rf, "irf/$start", result.irf)
    end
    if has(:recursions)
        for (k, v) in pairs(result.recursions)
            save_result_item(rf, "recursions/$k", v)
        end
    end
    return OutputVault.run_path(rf)
end

function load_solution(rf; model = nothing, native = false)
    if native
        return load_result_item(rf, "solution"; native = true)
    end
    parts = keys(load_result_item(rf, "solution"))
    T(p) = p === :obj ? ModelSolution : p === :stv ? SteadyState : Any
    return NamedTuple{parts}(map(p -> load_result_item(rf, "solution/$p", T(p); context = model), parts))
end

function load_simulation(rf; kw...)
    data = load_result_item(rf, "simulation")
    return (; simseries = SimulationPanel(data.simseries), errs = SimulationErrorPanel(data.errs))
end

function load_recursions(rf; kw...)
    names = keys(load_result_item(rf, "recursions"))
    return NamedTuple{names}(map(n -> load_result_item(rf, "recursions/$n", Interpolant), names))
end

# how each of TFI's items is loaded, in the order load_result(rf) loads them
const LOADERS = (
    solution   = load_solution,
    params     = (rf; kw...) -> load_result_item(rf, "params/params", MomentArray),
    simulation = load_simulation,
    errstats   = (rf; kw...) -> load_result_item(rf, "errstats/errstats", MomentArray),
    stats      = (rf; kw...) -> map(MomentArray, load_result_item(rf, "stats")),
    irf        = (rf; kw...) -> load_irfs(load_result_item(rf, "irf")),
    recursions = load_recursions,
)

function TFI.load_result(rf::RunHandle, item::Symbol; model = nothing, native::Bool = false)
    if !haskey(LOADERS, item)
        throw(ArgumentError("unknown item :$item; TFI's items are $(join(keys(LOADERS), ", "))"))
    end
    return LOADERS[item](rf; model, native)
end

function TFI.load_result(rf::RunHandle; model = nothing, native::Bool = false)
    held = Set(first(split(p, '/')) for p in OutputVault.list_run(rf).items.path)
    out = Pair{Symbol, Any}[]
    for item in keys(LOADERS)
        if !(string(item) in held)
            continue
        end
        value = TFI.load_result(rf, item; model, native)
        if item in (:solution, :simulation)
            append!(out, pairs(value))
        else
            push!(out, item => value)
        end
    end
    return NamedTuple(out)
end

# runs saved with OutputVault 0.3 hold one IRF from the steady state as irf/irf and irf/irf_paths
function irf_by_start(data::NamedTuple)
    if haskey(data, :irf)
        return (ss = (stats = data.irf, paths = get(data, :irf_paths, nothing)),)
    end
    return data
end
load_irfs(data) = map(d -> from_storage(IRFResult, d), irf_by_start(data))

function TFI.migrate_irf!(rf::RunHandle)
    items = OutputVault.list_run(rf).items.path
    moves = [p => "irf/ss/" * (p == "irf/irf" ? "stats" : "paths") for p in ("irf/irf", "irf/irf_paths") if p in items]
    for (from, to) in moves
        move_result_item(rf, from, to)
    end
    return moves
end

function TFI.migrate_irf!(pattern::RunPattern; move::Bool = false)
    for rf in find_runs(pattern)
        items = OutputVault.list_run(rf).items.path
        old = filter(in(items), ["irf/irf", "irf/irf_paths"])
        if isempty(old)
            continue
        end
        if move
            TFI.migrate_irf!(rf)
            println("moved: ", rf.branch, "/res-", rf.experdef, "-", rf.variant, ": ", join(old, ", "))
        else
            println("would move: ", rf.branch, "/res-", rf.experdef, "-", rf.variant, ": ", join(old, ", "))
        end
    end
    return nothing
end

## Moments from saved runs

# TFI's tables and where the fields of a MomentRow are in each table's query output: a column, or
# `nothing` for a fixed value (group "", period 0)
const SHAPES = (
    stats = (query = QueryShape("stats/{sample}"; keys = (:variable,), column = :stat),
             group = :sample, variable = :variable, stat = :stat, period = nothing),
    irf = (query = QueryShape("irf/{start}/stats"; keys = (:__shock, :__stat, :__period), column = :variable),
           group = :__shock, variable = :variable, stat = :__stat, period = :__period),
    params = (query = QueryShape("params/params"; keys = (:parameter,), column = :stat),
              group = nothing, variable = :parameter, stat = :stat, period = nothing))

function TFI.result_shape(kind::Symbol)
    if !haskey(SHAPES, kind)
        throw(ArgumentError("unknown shape :$kind; TFI's shapes are $(join(keys(SHAPES), ", "))"))
    end
    return SHAPES[kind].query
end

const Runs = Union{RunPattern, AbstractVector{<:RunHandle}, Tuple{Vararg{RunHandle}}}

function TFI.assemble_moments(moments::NamedTuple, runs::Runs; start = :ss)
    rows = TFI._requested_rows(moments)
    wanted = Set(rows)
    tables = [shape_table(runs, name, filter(r -> r.shape === name, rows), wanted, Symbol(start))
              for name in keys(SHAPES) if any(r -> r.shape === name, rows)]
    if isempty(tables)
        return Tables.columntable(TFI.MomentRecord[])
    end
    return reduce((a, b) -> map(vcat, a, b), tables)
end

# The numbers `rows` from the shape `name`, in the columns of the moment table. The query selects
# every combination of the rows' values, so only the numbers in `wanted` are kept.
function shape_table(runs, name, rows, wanted, start)
    spec = SHAPES[name]
    fields = (group = spec.group, variable = spec.variable, stat = spec.stat, period = spec.period)
    selections = (; (col => distinct(rows, f) for (f, col) in pairs(fields) if !isnothing(col))...)
    if :start in spec.query.placeholders
        selections = merge(selections, (; start))
    end
    t = query_results(runs, spec.query; selections...)
    n = length(t.value)
    group = isnothing(spec.group) ? fill("", n) : t[spec.group]
    period = isnothing(spec.period) ? zeros(Int, n) : t[spec.period]
    variable, stat = t[spec.variable], t[spec.stat]
    keep = [TFI.MomentRow((name, group[i], variable[i], stat[i], period[i])) in wanted for i in 1:n]
    return (branch = t.branch[keep], experdef = t.experdef[keep], variant = t.variant[keep],
            shape = fill(name, count(keep)), group = group[keep], variable = variable[keep], stat = stat[keep],
            period = period[keep], value = t.value[keep])
end

# the distinct values of field `f` of the rows, as a selection: text as Symbols, which match exactly
distinct(rows, f) = collect(Set(selector(getfield(r, f)) for r in rows))
selector(x::AbstractString) = Symbol(x)
selector(x) = x

run_order(runs::RunPattern) = find_runs(runs)
run_order(runs) = collect(RunHandle, runs)

function TFI.compare_economies(moments::NamedTuple, runs::Runs; start = :ss, kwargs...)
    moments = TFI._as_moments(moments)
    table = TFI.assemble_moments(moments, runs; start)
    order = run_order(runs)
    variants = [rf.variant for rf in order]
    if !allunique(variants)
        throw(ArgumentError("variants repeat across the runs ($(join(variants, ", "))); select runs with a narrower pattern"))
    end
    economies = TFI.split_by_economy(table; order = variants)
    return TFI.compare_economies(moments, economies; kwargs...)
end

end # module
