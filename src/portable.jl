# Portable data: solutions, steady states and interpolants as nested NamedTuples of numbers, strings
# and numeric arrays, which any storage of plain data (HDF5, JSON, ...) can hold and other languages
# can read, and their reconstruction. No I/O here: the OutputVault extension maps the tree to a .jld2
# file, one group per NamedTuple.
#
# Layout (schema version 2):
#   meta/{schema_version, tfi_version, created, solution_id, model_type, singletons}
#   meta/model/... (model_metadata tree)
#   grid/<axis> (vectors; a singleton axis is a length-1 vector), grid/__axes (axis order)
#   exogenv/points/<var>, exogenv/points/__names (variable order), exogenv/mtrans, exogenv/ssidx
#   Pfct/<var>, Vfct/<var>, Ffct/<var>: arrays over the grid axes (Pfct in unconstrained space)
#   Tfct/<var>: arrays over the grid axes and next period's exogenous state
#   <fct>/__vars (variable order)
#   Pfct/__bounds/<var>, Tfct/__bounds/<var>: [lower, upper] for each variable, plus __names
# Steady states: meta/{schema_version, tfi_version, model_type}, meta/model/..., and guess/, endog/,
# exog/, values/ (one entry per name, plus __names for the order). Interpolants: meta/{schema_version,
# tfi_version, created}, grid/ and fct/. The model is not stored: it is passed to from_portable, as a
# model or a model type (model_from_metadata).
# The schema version changes when a change breaks readers of earlier data; the reader accepts only
# the current version. Added entries do not change it: readers ignore unknown keys.

const PORTABLE_SCHEMA = 2
const PORTABLE_FCTS = (:Pfct, :Vfct, :Tfct, :Ffct)
const STV_PARTS = (:guess, :endog, :exog, :values)

"""
    to_portable(x) -> NamedTuple

`x` (a [`ModelSolution`](@ref), [`SteadyState`](@ref) or [`Interpolant`](@ref)) as a nested
NamedTuple whose leaves are numbers, strings and numeric arrays, with no Julia types of TFI or of the
model: what other languages and later versions of the code can read. [`from_portable`](@ref)
rebuilds `x`. The model itself is not included, only its metadata ([`model_metadata`](@ref)).
"""
function to_portable(obj::ModelSolution)
    meta = (; portable_header()..., solution_id = solution_id(obj), model_type = _model_type_name(obj.model),
              singletons = collect(map(string, singletons(obj.ssg))), model = model_metadata(obj.model))
    fcts = Pair{Symbol, Any}[]
    for name in PORTABLE_FCTS
        fct = getfield(obj, name)
        if !isnothing(fct)
            push!(fcts, name => portable_fct(fct))
        end
    end
    fcts = Dict(fcts)
    fcts[:Pfct] = (; fcts[:Pfct]..., __bounds = portable_bounds(obj.Pfct_bounds))
    fcts[:Tfct] = (; get(fcts, :Tfct, (;))..., __bounds = portable_bounds(obj.Tfct_bounds))
    fcts = NamedTuple(k => fcts[k] for k in PORTABLE_FCTS if haskey(fcts, k))
    return (; meta, grid = portable_grid(obj.ssg), exogenv = portable_exogenv(obj.exogenv), fcts...)
end

function to_portable(stv::SteadyState)
    meta = (schema_version = PORTABLE_SCHEMA, tfi_version = tfi_version(), model_type = _model_type_name(stv.model),
            model = model_metadata(stv.model))
    parts = NamedTuple{STV_PARTS}(map(part -> portable_namedtuple(getfield(stv, part)), STV_PARTS))
    return (; meta, parts...)
end

# an Interpolant over a state-space grid (e.g. recursions computed at a solution's policies): its
# grid and one array per variable, as in a solution, so that it rebuilds without the solution
to_portable(fct::Interpolant) = (meta = portable_header(), grid = portable_grid(fct.ssg), fct = portable_fct(fct))

portable_header() = (schema_version = PORTABLE_SCHEMA, tfi_version = tfi_version(),
                     created = Libc.strftime("%Y-%m-%dT%H:%M:%S", time()))

# a flat NamedTuple of numbers or numeric arrays, with the order of its names
portable_namedtuple(nt::NamedTuple) = (; __names = collect(map(string, keys(nt))), NamedTuple{keys(nt)}(map(k -> plain_value(nt[k], k), keys(nt)))...)
plain_value(v::Number, k) = v
plain_value(v::AbstractArray{<:Number}, k) = collect(v)
plain_value(v, k) = throw(ArgumentError("cannot store $k::$(typeof(v)) as plain data"))

function portable_grid(ssg)
    grids = unigrids(ssg)
    return (; __axes = collect(map(string, keys(grids))), map(g -> g isa Number ? [g] : collect(g), grids)...)
end

portable_exogenv(exogenv) = (points = portable_namedtuple(map(collect, StructArrays.components(points(exogenv)))),
                             mtrans = Matrix(transition_matrix(exogenv)), ssidx = ssidx(exogenv))

# one array per variable over the full grid (singleton axes have length 1); a trailing axis for vector outputs
function portable_fct(fct)
    coefs = coefficients(fct)
    vars = varnames(fct)
    arrays = map(enumerate(vars)) do (i, var)
        eltype(coefs) <: SMatrix ?
            permutedims(stack(c -> Vector(c[i, :]), coefs), (2:ndims(coefs)+1..., 1)) :
            map(c -> c[i], coefs)
    end
    return (; __vars = collect(map(string, vars)), NamedTuple{Tuple(vars)}(Tuple(arrays))...)
end

# bounds of each variable as [lower, upper] (±Inf where unbounded), in variable order
portable_bounds(nt::NamedTuple) = (; __names = collect(map(string, keys(nt))), map(bound_interval, nt)...)
bound_interval(::Unbounded) = [-Inf, Inf]
bound_interval(b::LowerBounded) = [Float64(b.lower), Inf]
bound_interval(b::UpperBounded) = [-Inf, Float64(b.upper)]
bound_interval(b::IntervalBounded) = [Float64(b.lower), Float64(b.upper)]

## Rebuilding

"""
    from_portable(T, data; model = nothing)

Rebuild a `T` ([`ModelSolution`](@ref), [`SteadyState`](@ref) or [`Interpolant`](@ref)) from `data`
as [`to_portable`](@ref) returns it, or as read back from storage (nested NamedTuples with the same
names; unknown keys are ignored). For a solution or a steady state, `model` is the model (its
metadata are compared with the stored ones, with a warning if they differ) or a model type (the
model is rebuilt from the stored metadata, [`model_from_metadata`](@ref)). An interpolant needs no
model.
"""
function from_portable(::Type{<:ModelSolution}, data; model = nothing)
    mdl = portable_model(data, model)
    check_portable_schema(data)
    exogenv = portable_exogenv_from(data.exogenv)
    if model isa AbstractModel
        check_portable_exogenv(exogenv, mdl)
    end
    axes = Tuple(Symbol(a) for a in data.grid.__axes)
    endog_axes = Tuple(a for a in axes if a != :exst)
    ssg = StateSpaceGrid(exogenv, NamedTuple{endog_axes}(map(a -> data.grid[a], endog_axes)))
    has_fct(name) = haskey(data, name) && haskey(data[name], :__vars)
    fcts = NamedTuple{PORTABLE_FCTS}(map(name -> has_fct(name) ? fct_from_portable(data[name], ssg) : nothing, PORTABLE_FCTS))
    obj = ModelSolution(; model = mdl, exogenv, ssg, fcts.Pfct, fcts.Vfct, fcts.Tfct, fcts.Ffct,
                          Pfct_bounds = bounds_from_portable(data.Pfct.__bounds), Tfct_bounds = bounds_from_portable(data.Tfct.__bounds))
    if solution_id(obj) != data.meta.solution_id
        @warn "The restored solution's id differs from the stored one: the data may be corrupted"
    end
    return obj
end

function from_portable(::Type{<:SteadyState}, data; model = nothing)
    mdl = portable_model(data, model)
    check_portable_schema(data)
    parts = NamedTuple{STV_PARTS}(map(part -> namedtuple_from_portable(data[part]), STV_PARTS))
    return SteadyState(; model = mdl, parts...)
end

function from_portable(::Type{<:Interpolant}, data; model = nothing)
    check_portable_schema(data)
    return fct_from_portable(data.fct, grid_from_portable(data.grid))
end

# the model of a rebuild: a model instance (compared with the stored metadata) or a model type (rebuilt)
function portable_model(data, model)
    stored = haskey(data.meta, :model) ? data.meta.model : (;)
    if model isa AbstractModel
        differ = _metadata_differences(stored, model_metadata(model))
        if !isempty(differ)
            @warn "The supplied model's metadata differ from the stored model's" keys = first(differ, 10)
        end
        return model
    elseif model isa Type && model <: AbstractModel
        return model_from_metadata(model, stored)
    end
    throw(ArgumentError("portable data of a solution or steady state store no model: rebuild them with a model or a model type"))
end

function check_portable_schema(data)
    v = data.meta.schema_version
    if v != PORTABLE_SCHEMA
        throw(ArgumentError("portable schema $v, this TransitionFunctionIteration reads $PORTABLE_SCHEMA"))
    end
end

# warn if the supplied model's exogenous environment differs from the stored one (which is used)
function check_portable_exogenv(exogenv, mdl)
    own = try Exogenv(mdl) catch; nothing end
    if !isnothing(own) && !isapprox(own, exogenv)
        @warn "The supplied model's exogenous environment differs from the stored one (the stored one is used)"
    end
end

function namedtuple_from_portable(d)
    names = d.__names
    return NamedTuple{Tuple(Symbol(n) for n in names)}(Tuple(d[Symbol(n)] for n in names))
end

portable_exogenv_from(d) = Exogenv(; pts = StructArray(namedtuple_from_portable(d.points)), mtrans = d.mtrans, ssidx = d.ssidx)

# the grid of an interpolant: the exogenous axis is the indices 1:n (as StateSpaceGrid(exogenv, ...) builds it)
function grid_from_portable(d)
    axes = Tuple(Symbol(a) for a in d.__axes)
    StateSpaceGrid(NamedTuple{axes}(map(a -> a == :exst ? (1:length(d.exst)) : d[a], axes)))
end

bounds_from_portable(d) = map(interval_bound, namedtuple_from_portable(d))
function interval_bound((lo, hi))
    if isinf(lo) && isinf(hi)
        return Unbounded()
    elseif isinf(hi)
        return LowerBounded(lo)
    elseif isinf(lo)
        return UpperBounded(hi)
    end
    return IntervalBounded(lo, hi)
end

# interpolant from per-variable arrays; singleton axes dropped to match the grid's points
function fct_from_portable(d, ssg)
    vars = Tuple(Symbol(v) for v in d.__vars)
    arrays = map(v -> d[v], vars)
    gridsize = Tuple(values(dims(ssg)))
    keep = Tuple(i for (i, a) in enumerate(keys(dims(ssg))) if !(a in singletons(ssg)))
    pointsize = gridsize[collect(keep)]
    vector_valued = ndims(first(arrays)) == length(gridsize) + 1
    values_at(I) = NamedTuple{vars}(map(a -> vector_valued ? SVector{size(a, ndims(a))}(a[I, :]) : a[I], arrays))
    nodes = reshape([values_at(I) for I in CartesianIndices(gridsize)], pointsize)
    Interpolant(ssg, StructArray(nodes))
end
