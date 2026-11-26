# TFI's results with OutputVault: how TFI's types are stored (OutputVault does the I/O), the
# conveniences save_result, load_result and migrate_irf!, and moments read from saved runs.
#
# Storage
# - SimulationPanel, SimulationErrorPanel and MomentArray are Tables.jl tables with TFI's DataAPI
#   metadata, so OutputVault writes them as table files; from_storage restores them from the table.
# - An IRFResult is the folder (stats, paths): two tables, the statistics with the IRF's metadata.
# - ModelSolution, SteadyState and Interpolant have a portable layout inside their .jld2 files
#   (write_jld2, read_jld2), plain HDF5-compatible arrays and metadata without Julia objects, so they
#   read without the code that wrote them. Layout (schema version 2):
#     meta/{schema_version, tfi_version, created, solution_id, model_type, singletons}
#     meta/model/... (model_metadata tree)
#     grid/<axis> (vectors; a singleton axis is a length-1 vector), grid/__axes (axis order)
#     exogenv/points/<var>, exogenv/points/__names (variable order), exogenv/mtrans, exogenv/ssidx
#     Pfct/<var>, Vfct/<var>, Ffct/<var>: arrays over the grid axes (Pfct in unconstrained space)
#     Tfct/<var>: arrays over the grid axes and next period's exogenous state
#     <fct>/__vars (variable order)
#     Pfct/__bounds/<var>, Tfct/__bounds/<var>: [lower, upper] for each variable, plus __names
#   Steady states: meta/{schema_version, tfi_version, model_type}, meta/model/..., and guess/, endog/,
#   exog/, values/ (one entry per name, plus __names for the order). Interpolants: meta/{schema_version,
#   tfi_version, created}, grid/ and fct/. The model is not stored: it comes from the `context` of
#   the load, a model or a model type (model_from_metadata).
#   The schema version changes when a change breaks readers of earlier files; the reader accepts
#   only the current version. Added entries do not change it (readers ignore unknown keys, such as
#   outputvault/kind and meta/provenance in files saved with OutputVault 0.3). Files are written
#   uncompressed: h5py cannot read JLD2's compressed (chunked) datasets.
module TransitionFunctionIterationOutputVaultExt

using OutputVault
using OutputVault: JLD2, RunHandle, RunPattern, QueryShape, find_runs, query_results
import OutputVault: to_storage, from_storage, write_jld2, read_jld2
import TransitionFunctionIteration as TFI
using TransitionFunctionIteration: AbstractModel, ModelSolution, SteadyState, MomentArray,
    SimulationPanel, SimulationErrorPanel, IRFResult, irf_paths, solution_id, Exogenv,
    StateSpaceGrid, Interpolant, unigrids, points, transition_matrix, ssidx, coefficients,
    varnames, dims, singletons, collect_parameters, model_metadata, model_from_metadata
using StaticArrays, StructArrays, Tables

const SCHEMA_VERSION = 2
const FCTS = (:Pfct, :Vfct, :Tfct, :Ffct)
const STV_PARTS = (:guess, :endog, :exog, :values)

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

## Portable .jld2 files: writing

function write_jld2(io, obj::ModelSolution)
    write_header(io)
    io["meta/solution_id"] = solution_id(obj)
    io["meta/model_type"] = TFI._model_type_name(obj.model)
    io["meta/singletons"] = collect(map(string, singletons(obj.ssg)))
    write_tree(io, "meta/model", model_metadata(obj.model))
    write_grid(io, obj.ssg)
    write_exogenv(io, obj.exogenv)
    for name in FCTS
        fct = getfield(obj, name)
        if !isnothing(fct)
            write_fct(io, string(name), fct)
        end
    end
    write_bounds(io, "Pfct/__bounds", obj.Pfct_bounds)
    write_bounds(io, "Tfct/__bounds", obj.Tfct_bounds)
end

function write_jld2(io, stv::SteadyState)
    io["meta/schema_version"] = SCHEMA_VERSION
    io["meta/tfi_version"] = TFI.tfi_version()
    io["meta/model_type"] = TFI._model_type_name(stv.model)
    write_tree(io, "meta/model", model_metadata(stv.model))
    for part in STV_PARTS
        write_namedtuple(io, string(part), getfield(stv, part))
    end
end

# an Interpolant over a state-space grid (e.g. recursions computed at a solution's policies): its
# grid and one array per variable, as in a solution file, so that it loads without the solution
function write_jld2(io, fct::Interpolant)
    write_header(io)
    write_grid(io, fct.ssg)
    write_fct(io, "fct", fct)
end

function write_header(io)
    io["meta/schema_version"] = SCHEMA_VERSION
    io["meta/tfi_version"] = TFI.tfi_version()
    io["meta/created"] = Libc.strftime("%Y-%m-%dT%H:%M:%S", time())
end

# nested NamedTuple as groups
function write_tree(io, prefix, nt::NamedTuple)
    for (k, v) in pairs(nt)
        if v isa NamedTuple
            write_tree(io, "$prefix/$k", v)
        else
            io["$prefix/$k"] = v
        end
    end
end

# a flat NamedTuple of numbers or numeric arrays as one entry per name, plus the order
function write_namedtuple(io, prefix, nt::NamedTuple)
    io["$prefix/__names"] = collect(map(string, keys(nt)))
    for (k, v) in pairs(nt)
        io["$prefix/$k"] = plain_value(v, k)
    end
end
plain_value(v::Number, k) = v
plain_value(v::AbstractArray{<:Number}, k) = collect(v)
plain_value(v, k) = throw(ArgumentError("cannot store $k::$(typeof(v)) as plain data"))

function write_grid(io, ssg)
    grids = unigrids(ssg)
    io["grid/__axes"] = collect(map(string, keys(grids)))
    for (k, g) in pairs(grids)
        io["grid/$k"] = g isa Number ? [g] : collect(g)
    end
end

function write_exogenv(io, exogenv)
    write_namedtuple(io, "exogenv/points", map(collect, StructArrays.components(points(exogenv))))
    io["exogenv/mtrans"] = Matrix(transition_matrix(exogenv))
    io["exogenv/ssidx"] = ssidx(exogenv)
end

# one array per variable over the full grid (singleton axes have length 1); a trailing axis for vector outputs
function write_fct(io, name, fct)
    coefs = coefficients(fct)
    io["$name/__vars"] = collect(map(string, varnames(fct)))
    for (i, var) in enumerate(varnames(fct))
        io["$name/$var"] = eltype(coefs) <: SMatrix ?
            permutedims(stack(c -> Vector(c[i, :]), coefs), (2:ndims(coefs)+1..., 1)) :
            map(c -> c[i], coefs)
    end
end

# bounds of each variable as [lower, upper] (±Inf where unbounded), in variable order
function write_bounds(io, prefix, nt::NamedTuple)
    io["$prefix/__names"] = collect(map(string, keys(nt)))
    for (k, b) in pairs(nt)
        io["$prefix/$k"] = bound_interval(b)
    end
end
bound_interval(::TFI.Unbounded) = [-Inf, Inf]
bound_interval(b::TFI.LowerBounded) = [Float64(b.lower), Inf]
bound_interval(b::TFI.UpperBounded) = [-Inf, Float64(b.upper)]
bound_interval(b::TFI.IntervalBounded) = [Float64(b.lower), Float64(b.upper)]

## Portable .jld2 files: reading

# the model of a load: a model instance (compared with the stored metadata) or a model type (rebuilt)
function file_model(io, context)
    if context isa AbstractModel
        check_model(io, context)
        return context
    elseif context isa Type && context <: AbstractModel
        return model_from_metadata(context, read_tree(io, "meta/model"))
    end
    throw(ArgumentError("a portable solution or steady state stores no model: load it with a model or a model type (`model` in load_result, `context` in load_result_item), or read the native copy"))
end

function read_jld2(::Type{<:ModelSolution}, io; context = nothing)
    mdl = file_model(io, context)
    if context isa AbstractModel
        check_exogenv(io, mdl)
    end
    return rebuild(io, mdl)
end

function read_jld2(::Type{<:SteadyState}, io; context = nothing)
    mdl = file_model(io, context)
    check_schema(io)
    parts = NamedTuple{STV_PARTS}(map(part -> read_namedtuple(io, string(part)), STV_PARTS))
    return SteadyState(; model = mdl, parts...)
end

function read_jld2(::Type{<:Interpolant}, io; context = nothing)
    check_schema(io)
    return read_fct(io, "fct", read_grid(io))
end

read_tree(io, prefix) = haskey(io, prefix) ? read_group(io[prefix]) : (;)
read_group(g) = NamedTuple(Symbol(k) => (v = g[k]; v isa JLD2.Group ? read_group(v) : v) for k in keys(g))

function read_namedtuple(io, prefix)
    names = io["$prefix/__names"]
    NamedTuple{Tuple(Symbol(n) for n in names)}(Tuple(io["$prefix/$n"] for n in names))
end

# warn if the supplied model's metadata differ from the stored ones
function check_model(io, mdl)
    differ = TFI._metadata_differences(read_tree(io, "meta/model"), model_metadata(mdl))
    if !isempty(differ)
        @warn "The supplied model's metadata differ from the stored model's" keys = first(differ, 10)
    end
end

# warn if the supplied model's exogenous environment differs from the stored one (which is used)
function check_exogenv(io, mdl)
    exogenv = try Exogenv(mdl) catch; nothing end
    if !isnothing(exogenv) && !isapprox(exogenv, read_exogenv(io))
        @warn "The supplied model's exogenous environment differs from the stored one (the stored one is used)"
    end
end

read_exogenv(io) = Exogenv(; pts = StructArray(read_namedtuple(io, "exogenv/points")),
                             mtrans = io["exogenv/mtrans"], ssidx = io["exogenv/ssidx"])

function check_schema(io)
    if io["meta/schema_version"] != SCHEMA_VERSION
        throw(ArgumentError("file schema $(io["meta/schema_version"]), this TransitionFunctionIteration reads $SCHEMA_VERSION"))
    end
end

function rebuild(io, mdl)
    check_schema(io)
    exogenv = read_exogenv(io)
    axes = Tuple(Symbol(a) for a in io["grid/__axes"])
    endog_axes = Tuple(a for a in axes if a != :exst)
    ssg = StateSpaceGrid(exogenv, NamedTuple{endog_axes}(map(a -> io["grid/$a"], endog_axes)))
    fcts = NamedTuple{FCTS}(map(name -> haskey(io, string(name)) ? read_fct(io, string(name), ssg) : nothing, FCTS))
    obj = ModelSolution(; model = mdl, exogenv, ssg, fcts.Pfct, fcts.Vfct, fcts.Tfct, fcts.Ffct,
                          Pfct_bounds = read_bounds(io, "Pfct/__bounds"), Tfct_bounds = read_bounds(io, "Tfct/__bounds"))
    if solution_id(obj) != io["meta/solution_id"]
        @warn "The restored solution's id differs from the stored one: the file may be corrupted"
    end
    return obj
end

# the grid of a file: the exogenous axis is the indices 1:n (as StateSpaceGrid(exogenv, ...) builds it)
function read_grid(io)
    axes = Tuple(Symbol(a) for a in io["grid/__axes"])
    grids = NamedTuple{axes}(map(a -> a == :exst ? (1:length(io["grid/exst"])) : io["grid/$a"], axes))
    StateSpaceGrid(grids)
end

read_bounds(io, prefix) = map(interval_bound, read_namedtuple(io, prefix))
function interval_bound((lo, hi))
    if isinf(lo) && isinf(hi)
        return TFI.Unbounded()
    elseif isinf(hi)
        return TFI.LowerBounded(lo)
    elseif isinf(lo)
        return TFI.UpperBounded(hi)
    end
    return TFI.IntervalBounded(lo, hi)
end

# interpolant from per-variable arrays; singleton axes dropped to match the grid's points
function read_fct(io, name, ssg)
    vars = Tuple(Symbol(v) for v in io["$name/__vars"])
    arrays = map(v -> io["$name/$v"], vars)
    gridsize = Tuple(values(dims(ssg)))
    keep = Tuple(i for (i, a) in enumerate(keys(dims(ssg))) if !(a in singletons(ssg)))
    pointsize = gridsize[collect(keep)]
    vector_valued = ndims(first(arrays)) == length(gridsize) + 1
    values_at(I) = NamedTuple{vars}(map(a -> vector_valued ? SVector{size(a, ndims(a))}(a[I, :]) : a[I], arrays))
    nodes = reshape([values_at(I) for I in CartesianIndices(gridsize)], pointsize)
    Interpolant(ssg, StructArray(nodes))
end

## Saving and loading conveniences

function TFI.save_result(rf::RunHandle, result; start = :ss, dirty_repo_check::Symbol = :warn)
    has(p) = hasproperty(result, p) && !isnothing(getproperty(result, p))
    if has(:obj)
        solution = has(:stv) ? (; result.obj, result.stv) : (; result.obj)
        save_solution(rf, solution; dirty_repo_check)
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

const ITEMS = (:solution, :params, :simulation, :errstats, :stats, :irf, :recursions)

function TFI.load_result(rf::RunHandle, item::Symbol; model = nothing, native::Bool = false)
    if item === :solution
        if native
            return load_result_item(rf, "solution"; native = true)
        end
        parts = keys(load_result_item(rf, "solution"))
        T(p) = p === :obj ? ModelSolution : p === :stv ? SteadyState : Any
        return NamedTuple{parts}(map(p -> load_result_item(rf, "solution/$p", T(p); context = model), parts))
    elseif item === :params
        return load_result_item(rf, "params/params", MomentArray)
    elseif item === :errstats
        return load_result_item(rf, "errstats/errstats", MomentArray)
    elseif item === :simulation
        data = load_result_item(rf, "simulation")
        return (; simseries = SimulationPanel(data.simseries), errs = SimulationErrorPanel(data.errs))
    elseif item === :stats
        return map(MomentArray, load_result_item(rf, "stats"))
    elseif item === :irf
        return load_irfs(load_result_item(rf, "irf"))
    elseif item === :recursions
        names = keys(load_result_item(rf, "recursions"))
        return NamedTuple{names}(map(n -> load_result_item(rf, "recursions/$n", Interpolant), names))
    end
    throw(ArgumentError("unknown item :$item; TFI's items are $(join(ITEMS, ", "))"))
end

function TFI.load_result(rf::RunHandle; model = nothing, native::Bool = false)
    held = Set(first(split(p, '/')) for p in OutputVault.list_run(rf).items.path)
    out = Pair{Symbol, Any}[]
    for item in ITEMS
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
    tables = Any[]
    sim = [r for r in rows if r.kind == "simulation"]
    if !isempty(sim)
        t = query_results(runs, SHAPES.stats; sample = unique(Symbol.(getfield.(sim, :group))),
                          variable = unique(Symbol.(getfield.(sim, :variable))), stat = unique(Symbol.(getfield.(sim, :stat))))
        push!(tables, long_rows(t, "simulation", t.sample, t.variable, t.stat, zeros(Int, length(t.value)), sim))
    end
    imp = [r for r in rows if r.kind == "impulse"]
    if !isempty(imp)
        t = query_results(runs, SHAPES.irf; start = Symbol(start), __shock = unique(Symbol.(getfield.(imp, :group))),
                          __stat = unique(Symbol.(getfield.(imp, :stat))), __period = unique(getfield.(imp, :period)),
                          variable = unique(Symbol.(getfield.(imp, :variable))))
        push!(tables, long_rows(t, "impulse", t.__shock, t.variable, t.__stat, Int.(t.__period), imp))
    end
    par = [r for r in rows if r.kind == "parameter"]
    if !isempty(par)
        t = query_results(runs, SHAPES.params; parameter = unique(Symbol.(getfield.(par, :variable))), stat = :value)
        push!(tables, long_rows(t, "parameter", fill("parameters", length(t.value)), t.parameter, t.stat, zeros(Int, length(t.value)), par))
    end
    return vcat_tables(tables)
end

# the query's rows in the columns of the moment table, keeping only the requested rows
function long_rows(t, kind, group, variable, stat, period, requested)
    wanted = Set(requested)
    keep = [i for i in eachindex(t.value) if TFI.MomentRow((kind, group[i], variable[i], stat[i], period[i])) in wanted]
    return NamedTuple{TFI.MOMENT_TABLE_COLUMNS}((String.(t.branch[keep]), String.(t.experdef[keep]), String.(t.variant[keep]),
        fill(kind, length(keep)), String.(group[keep]), String.(variable[keep]), String.(stat[keep]),
        Int.(period[keep]), Float64.(t.value[keep])))
end

function vcat_tables(tables)
    if isempty(tables)
        return NamedTuple{TFI.MOMENT_TABLE_COLUMNS}((String[], String[], String[], String[], String[], String[], String[], Int[], Float64[]))
    end
    return NamedTuple{TFI.MOMENT_TABLE_COLUMNS}(Tuple(reduce(vcat, [t[k] for t in tables]) for k in TFI.MOMENT_TABLE_COLUMNS))
end

run_order(runs::RunPattern) = find_runs(runs)
run_order(runs) = collect(RunHandle, runs)

function TFI.compare_economies(moments::NamedTuple, runs::Runs; start = :ss, kwargs...)
    moments = TFI._as_moments(moments)
    table = TFI.assemble_moments(moments, runs; start)
    order = run_order(runs)
    variants = [rf.variant for rf in order]
    if !allunique(variants)
        dup = unique(v for v in variants if count(==(v), variants) > 1)
        throw(ArgumentError("variants $(join(dup, ", ")) appear in several runs; select runs with a narrower pattern"))
    end
    economies = TFI.split_by_economy(table; order = variants)
    return TFI.compare_economies(moments, economies; kwargs...)
end

end # module
