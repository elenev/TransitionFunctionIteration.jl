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

const SHAPES = (stats = QueryShape("stats/{sample}"; keys = (:variable,), column = :stat),
                irf = QueryShape("irf/{start}/stats"; keys = (:__shock, :__stat, :__period), column = :variable),
                params = QueryShape("params/params"; keys = (:parameter,), column = :stat))

function TFI.result_shape(kind::Symbol)
    if !haskey(SHAPES, kind)
        throw(ArgumentError("unknown shape :$kind; TFI's shapes are $(join(keys(SHAPES), ", "))"))
    end
    return SHAPES[kind]
end

const Runs = Union{RunPattern, AbstractVector{<:RunHandle}, Tuple{Vararg{RunHandle}}}

function TFI.assemble_moments(moments::NamedTuple, runs::Runs; start = :ss)
    rows = TFI._requested_rows(moments)
    wanted = Set(rows)
    out = TFI._empty_moment_table()
    sim = filter(r -> r.kind == "simulation", rows)
    if !isempty(sim)
        add_simulation!(out, runs, sim, wanted)
    end
    imp = filter(r -> r.kind == "impulse", rows)
    if !isempty(imp)
        add_impulse!(out, runs, imp, wanted, Symbol(start))
    end
    par = filter(r -> r.kind == "parameter", rows)
    if !isempty(par)
        add_parameter!(out, runs, par, wanted)
    end
    return out
end

# the distinct values of field `f` of the rows, as a query selection
selection(rows, f) = collect(Set(Symbol(getfield(r, f)) for r in rows))

# Each query selects every combination of the rows' values; only the rows a moment needs are kept.
function add_simulation!(out, runs, rows, wanted)
    t = query_results(runs, SHAPES.stats; sample = selection(rows, :group), variable = selection(rows, :variable),
                      stat = selection(rows, :stat))
    for i in eachindex(t.value)
        r = TFI.MomentRow(("simulation", t.sample[i], t.variable[i], t.stat[i], 0))
        if r in wanted
            TFI._push_moment!(out, t.branch[i], t.experdef[i], t.variant[i], r, t.value[i])
        end
    end
end

function add_impulse!(out, runs, rows, wanted, start)
    t = query_results(runs, SHAPES.irf; start, __shock = selection(rows, :group), __stat = selection(rows, :stat),
                      __period = collect(Set(r.period for r in rows)), variable = selection(rows, :variable))
    for i in eachindex(t.value)
        r = TFI.MomentRow(("impulse", t.__shock[i], t.variable[i], t.__stat[i], t.__period[i]))
        if r in wanted
            TFI._push_moment!(out, t.branch[i], t.experdef[i], t.variant[i], r, t.value[i])
        end
    end
end

function add_parameter!(out, runs, rows, wanted)
    t = query_results(runs, SHAPES.params; parameter = selection(rows, :variable), stat = :value)
    for i in eachindex(t.value)
        r = TFI.MomentRow(("parameter", "parameters", t.parameter[i], t.stat[i], 0))
        if r in wanted
            TFI._push_moment!(out, t.branch[i], t.experdef[i], t.variant[i], r, t.value[i])
        end
    end
end

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
