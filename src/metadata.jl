# Model metadata, solution ids and the solution (de)serialization stubs.
#
# Plain metadata are values every format can store: reals, strings, symbols, small numeric
# arrays, and NamedTuples of these. Model metadata are stored with saved solutions (the JLD2
# extension), where `model_from_metadata` restores a model and a supplied model is checked
# against them. Simulations do not carry model metadata: they record the `solution_id` of the
# solution that produced them, which is also stored in the solution file.

const MAX_METADATA_ARRAY_LENGTH = 1_000
const MAX_METADATA_DEPTH = 4

"""
    model_metadata(mdl::AbstractModel) -> NamedTuple

Plain-data description of a model: stored with saved solutions and compared with a supplied
model when a solution is loaded.

Default: the model's fields, recursively (NamedTuples, dictionaries with Symbol or String keys,
structs; at most $MAX_METADATA_DEPTH levels), keeping reals, strings, symbols and numeric
arrays/tuples with at most $MAX_METADATA_ARRAY_LENGTH elements. Anything else (functions,
interpolants, large arrays) is left out. Overload it for models whose fields include derived
objects, and overload `model_from_metadata` to match.
"""
model_metadata(mdl::AbstractModel) = something(_plain_metadata(mdl, 0), (;))

"""
    model_from_metadata(::Type{M}, metadata::NamedTuple) -> M

Rebuild a model of type `M` from `model_metadata` output (used by `load_solution(path, M)`).
Default: `M(; metadata...)`, which works for flat `@kwdef` models.
"""
model_from_metadata(::Type{M}, metadata::NamedTuple) where {M<:AbstractModel} = M(; metadata...)

# plain values pass through, anything else becomes `nothing` (dropped)
_plain_metadata(x::Union{Real,AbstractString,Symbol}, depth) = x
_plain_metadata(x::Union{Array{<:Real},AbstractRange{<:Real},StaticArray{<:Tuple,<:Real}}, depth) =
    length(x) <= MAX_METADATA_ARRAY_LENGTH ? collect(x) : nothing
_plain_metadata(x::Tuple{Vararg{Real}}, depth) = length(x) <= MAX_METADATA_ARRAY_LENGTH ? collect(x) : nothing
_plain_metadata(x::NamedTuple, depth) = _plain_namedtuple(pairs(x), depth)
_plain_metadata(x::AbstractDict{<:Union{Symbol,AbstractString}}, depth) = _plain_namedtuple(x, depth)
_plain_metadata(::Union{AbstractArray,AbstractDict,Function,Module,Type,Ptr}, depth) = nothing
function _plain_metadata(x, depth)
    T = typeof(x)
    isstructtype(T) && fieldcount(T) > 0 || return nothing
    _plain_namedtuple((f => getfield(x, f) for f in fieldnames(T)), depth)
end

function _plain_namedtuple(prs, depth)
    depth >= MAX_METADATA_DEPTH && return nothing
    kept = Pair{Symbol,Any}[]
    for (k, v) in prs
        pv = _plain_metadata(v, depth + 1)
        pv === nothing || push!(kept, Symbol(k) => pv)
    end
    isempty(kept) ? nothing : NamedTuple(kept)
end

# nested metadata as String pairs with dotted keys (`prefix.field.subfield`), so that two models'
# metadata can be compared key by key (load_solution with a supplied model)
function flatten_metadata(nt::NamedTuple, prefix::AbstractString)
    out = Dict{String,String}()
    for (k, v) in pairs(nt)
        key = isempty(prefix) ? string(k) : string(prefix, ".", k)
        v isa NamedTuple ? merge!(out, flatten_metadata(v, key)) : (out[key] = string(v))
    end
    out
end

"TFI version string, recorded in simulation provenance and saved solutions"
tfi_version() = string(pkgversion(@__MODULE__))

"""
    solution_id(obj::ModelSolution) -> String

Identifier of a solution's content: a SHA-1 hash of the grid points and of the values of Pfct,
Vfct, Tfct and Ffct (if any) at the grid nodes. Solutions with identical values have identical
ids, also after `save_solution` and `load_solution`. Simulations record it in their provenance.
"""
function solution_id(obj::ModelSolution)
    values = Float64[]
    foreach(g -> _append_numbers!(values, g), unigrids(obj.ssg))
    for fct in (obj.Pfct, obj.Vfct, obj.Tfct, obj.Ffct)
        isnothing(fct) || _append_numbers!(values, coefficients(fct))
    end
    bytes2hex(SHA.sha1(reinterpret(UInt8, values)))
end

# every number in a (nested) container, in iteration order: grid axes, then the nodes' SVectors
# (Pfct, Vfct, Ffct) or SMatrices (Tfct)
_append_numbers!(out, x::Number) = push!(out, x)
_append_numbers!(out, x) = (foreach(y -> _append_numbers!(out, y), x); out)

"""
    save_solution(path, obj::ModelSolution; provenance = Dict{String,String}())

Write a solution to a single HDF5-compatible file: grids, exogenous environment and the values
of Pfct, Vfct, Tfct and Ffct (if any) at the grid nodes as plain arrays, metadata (schema and TFI
versions, `solution_id`, `model_metadata`), plus an exact JLD2 copy of the model. `provenance`
holds whatever else identifies the run that produced the solution (e.g. a git commit), as
strings. Requires JLD2 (`using JLD2`).
"""
function save_solution end

"""
    load_solution(path)                      # model restored from the exact JLD2 copy
    load_solution(path, mdl::AbstractModel)  # supplied model, compared with the stored metadata
    load_solution(path, M::Type)             # model rebuilt by model_from_metadata(M, metadata)

Restore a solution written by `save_solution`; warns if the restored solution's `solution_id`
differs from the stored one. Requires JLD2 (`using JLD2`).
"""
function load_solution end

function _jld2_hint(io, exc, argtypes, kwargs)
    if exc.f in (save_solution, load_solution) && isnothing(Base.get_extension(@__MODULE__, :TFIJLD2Ext))
        print(io, "\n$(exc.f) needs JLD2: add it to your environment and run `using JLD2`.")
    end
end
