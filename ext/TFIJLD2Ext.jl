# Solution files: plain HDF5-compatible arrays and metadata written with JLD2, plus an exact
# JLD2 copy of the model and the bounds. Layout (schema version 1):
#   meta/{schema_version, tfi_version, created, solution_id, model_type, singletons}
#   meta/model/... (model_metadata tree), meta/provenance/<key> (strings supplied by the caller)
#   grid/<axis> (vectors; a singleton axis is a length-1 vector), grid/__axes (axis order)
#   exogenv/points/<var>, exogenv/mtrans, exogenv/ssidx
#   Pfct/<var>, Vfct/<var>, Ffct/<var>: arrays over the grid axes (Pfct in unconstrained space)
#   Tfct/<var>: arrays over the grid axes and next period's exogenous state
#   <fct>/__vars (variable order), Pfct/__bounds and Tfct/__bounds (bounds as text)
#   exact/{model, Pfct_bounds, Tfct_bounds} (JLD2 objects)
# Files are written uncompressed: h5py cannot read JLD2's compressed (chunked) datasets.
module TFIJLD2Ext

using TransitionFunctionIteration, JLD2, StaticArrays, StructArrays
const TFI = TransitionFunctionIteration
import TransitionFunctionIteration: save_solution, load_solution

const SCHEMA_VERSION = 1
const FCTS = (:Pfct, :Vfct, :Tfct, :Ffct)

## Writing

function save_solution(path::AbstractString, obj::ModelSolution;
                       provenance::AbstractDict = Dict{String,String}())
    jldopen(path, "w") do io
        io["meta/schema_version"] = SCHEMA_VERSION
        io["meta/tfi_version"] = TFI.tfi_version()
        io["meta/created"] = Libc.strftime("%Y-%m-%dT%H:%M:%S", time())
        io["meta/solution_id"] = solution_id(obj)
        for (k, v) in provenance
            io["meta/provenance/$k"] = string(v)
        end
        io["meta/model_type"] = string(parentmodule(typeof(obj.model)), ".", nameof(typeof(obj.model)))
        io["meta/singletons"] = collect(map(string, TFI.singletons(obj.ssg)))
        write_tree(io, "meta/model", model_metadata(obj.model))
        write_grid(io, obj)
        write_exogenv(io, obj.exogenv)
        for name in FCTS
            fct = getfield(obj, name)
            isnothing(fct) || write_fct(io, string(name), fct)
        end
        io["Pfct/__bounds"] = [string(k, " = ", v) for (k, v) in pairs(obj.Pfct_bounds)]
        io["Tfct/__bounds"] = [string(k, " = ", v) for (k, v) in pairs(obj.Tfct_bounds)]
        io["exact/model"] = obj.model
        io["exact/Pfct_bounds"] = obj.Pfct_bounds
        io["exact/Tfct_bounds"] = obj.Tfct_bounds
    end
    path
end

"nested NamedTuple as groups"
function write_tree(io, prefix, nt::NamedTuple)
    for (k, v) in pairs(nt)
        v isa NamedTuple ? write_tree(io, "$prefix/$k", v) : (io["$prefix/$k"] = v)
    end
end

function write_grid(io, obj)
    grids = unigrids(obj.ssg)
    io["grid/__axes"] = collect(map(string, keys(grids)))
    for (k, g) in pairs(grids)
        io["grid/$k"] = g isa Number ? [g] : collect(g)
    end
end

function write_exogenv(io, exogenv)
    for (k, v) in pairs(StructArrays.components(points(exogenv)))
        io["exogenv/points/$k"] = collect(v)
    end
    io["exogenv/mtrans"] = Matrix(transition_matrix(exogenv))
    io["exogenv/ssidx"] = ssidx(exogenv)
end

"one array per variable over the full grid (singleton axes have length 1); a trailing axis for vector outputs"
function write_fct(io, name, fct)
    coefs = coefficients(fct)
    io["$name/__vars"] = collect(map(string, varnames(fct)))
    for (i, var) in enumerate(varnames(fct))
        io["$name/$var"] = eltype(coefs) <: SMatrix ?
            stack(c -> Vector(c[i, :]), coefs) |> a -> permutedims(a, (2:ndims(a)..., 1)) :
            map(c -> c[i], coefs)
    end
end

## Reading

load_solution(path::AbstractString) = jldopen(io -> rebuild(io, io["exact/model"]), path, "r")

function load_solution(path::AbstractString, mdl::AbstractModel)
    jldopen(path, "r") do io
        check_model(io, mdl)
        rebuild(io, mdl)
    end
end

function load_solution(path::AbstractString, ::Type{M}) where {M<:AbstractModel}
    jldopen(path, "r") do io
        rebuild(io, model_from_metadata(M, read_tree(io, "meta/model")))
    end
end

read_tree(io, prefix) = haskey(io, prefix) ? read_group(io[prefix]) : (;)
read_group(g) = NamedTuple(Symbol(k) => (v = g[k]; v isa JLD2.Group ? read_group(v) : v) for k in keys(g))

"warn if the supplied model's metadata or exogenous environment differ from the stored ones"
function check_model(io, mdl)
    stored = TFI.flatten_metadata(read_tree(io, "meta/model"), "")
    supplied = TFI.flatten_metadata(model_metadata(mdl), "")
    differ = sort!([k for k in union(keys(stored), keys(supplied)) if get(stored, k, nothing) != get(supplied, k, nothing)])
    isempty(differ) || @warn "The supplied model's metadata differ from the stored model's" keys = first(differ, 10)
    exogenv = try Exogenv(mdl) catch; nothing end
    isnothing(exogenv) || isapprox(exogenv, read_exogenv(io)) ||
        @warn "The supplied model's exogenous environment differs from the stored one (the stored one is used)"
end

function read_exogenv(io)
    g = io["exogenv/points"]
    names = Tuple(Symbol(k) for k in keys(g))
    pts = StructArray(NamedTuple{names}(map(k -> g[string(k)], names)))
    Exogenv(; pts, mtrans = io["exogenv/mtrans"], ssidx = io["exogenv/ssidx"])
end

function rebuild(io, mdl)
    io["meta/schema_version"] == SCHEMA_VERSION ||
        throw(ArgumentError("solution file schema $(io["meta/schema_version"]), this TFI reads $SCHEMA_VERSION"))
    exogenv = read_exogenv(io)
    axes = Tuple(Symbol(a) for a in io["grid/__axes"])
    endog_axes = Tuple(a for a in axes if a != :exst)
    ssg = StateSpaceGrid(exogenv, NamedTuple{endog_axes}(map(a -> io["grid/$a"], endog_axes)))
    fcts = NamedTuple{FCTS}(map(name -> haskey(io, string(name)) ? read_fct(io, string(name), ssg) : nothing, FCTS))
    obj = ModelSolution(; model = mdl, exogenv, ssg, fcts.Pfct, fcts.Vfct, fcts.Tfct, fcts.Ffct,
                          Pfct_bounds = io["exact/Pfct_bounds"], Tfct_bounds = io["exact/Tfct_bounds"])
    solution_id(obj) == io["meta/solution_id"] ||
        @warn "The restored solution's id differs from the stored one: the file may be corrupted"
    obj
end

"interpolant from per-variable arrays; singleton axes dropped to match the grid's points"
function read_fct(io, name, ssg)
    vars = Tuple(Symbol(v) for v in io["$name/__vars"])
    arrays = map(v -> io["$name/$v"], vars)
    gridsize = Tuple(values(dims(ssg)))
    keep = Tuple(i for (i, a) in enumerate(keys(dims(ssg))) if !(a in TFI.singletons(ssg)))
    pointsize = gridsize[collect(keep)]
    vector_valued = ndims(first(arrays)) == length(gridsize) + 1
    values_at(I) = NamedTuple{vars}(map(a -> vector_valued ? SVector{size(a, ndims(a))}(a[I, :]) : a[I], arrays))
    nodes = reshape([values_at(I) for I in CartesianIndices(gridsize)], pointsize)
    Interpolant(ssg, StructArray(nodes))
end

end
