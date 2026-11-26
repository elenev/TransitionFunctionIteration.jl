# Model metadata: a plain-data description of a model, stored with portable solution and
# steady-state files, compared with a supplied model when one is loaded, and used to rebuild a
# model when a model type is supplied instead.

const MAX_METADATA_ARRAY_LENGTH = 1_000
const MAX_METADATA_DEPTH = 4

"""
    model_metadata(mdl) -> NamedTuple

A plain-data description of a model: reals, strings, symbols, small numeric arrays and
NamedTuples of these. The default takes the model's fields recursively (NamedTuples, dictionaries
with Symbol or String keys and structs, at most $(MAX_METADATA_DEPTH) levels deep), keeping reals,
strings, symbols and numeric arrays or tuples of at most $(MAX_METADATA_ARRAY_LENGTH) elements;
functions, interpolants and large arrays are left out. Add a method, together with
[`model_from_metadata`](@ref), for a model whose fields include derived objects.
"""
model_metadata(mdl::AbstractModel) = something(_plain_metadata(mdl, 0), (;))

"""
    model_from_metadata(::Type{M}, metadata::NamedTuple) -> M

Rebuild a model of type `M` from [`model_metadata`](@ref) output, used when a portable solution or
steady state is loaded with a model type rather than a model. Default `M(; metadata...)`.
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
    if !(isstructtype(T) && fieldcount(T) > 0)
        return nothing
    end
    return _plain_namedtuple((f => getfield(x, f) for f in fieldnames(T)), depth)
end

function _plain_namedtuple(prs, depth)
    if depth >= MAX_METADATA_DEPTH
        return nothing
    end
    kept = Pair{Symbol,Any}[]
    for (k, v) in prs
        pv = _plain_metadata(v, depth + 1)
        if pv !== nothing
            push!(kept, Symbol(k) => pv)
        end
    end
    return isempty(kept) ? nothing : NamedTuple(kept)
end

# nested metadata as String pairs with dotted keys, so that two models' metadata compare key by key
function _flatten_metadata(nt::NamedTuple, prefix::AbstractString = "")
    out = Dict{String,String}()
    for (k, v) in pairs(nt)
        key = isempty(prefix) ? string(k) : string(prefix, ".", k)
        if v isa NamedTuple
            merge!(out, _flatten_metadata(v, key))
        else
            out[key] = string(v)
        end
    end
    return out
end

# the keys whose values differ between two models' metadata
_metadata_differences(a::NamedTuple, b::NamedTuple) = let fa = _flatten_metadata(a), fb = _flatten_metadata(b)
    sort!([k for k in union(keys(fa), keys(fb)) if get(fa, k, nothing) != get(fb, k, nothing)])
end

_model_type_name(mdl) = string(parentmodule(typeof(mdl)), ".", nameof(typeof(mdl)))
