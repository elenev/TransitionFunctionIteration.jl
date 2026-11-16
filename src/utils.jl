# Internal utilities

"""
    StaticStructArray(S::StaticArray{Sz, NT}) where {Sz, NT <: NamedTuple}
    Convert a StaticArray of NamedTuples into a StructArray of StaticArrays, preserving the original shape of each field.
    Avoids allocations by reconstructing StaticArrays directly from the original data.
"""
@generated function StaticStructArray(S::StaticArray{Sz, NT}) where {Sz, NT <: NamedTuple}
    keys = fieldnames(NT)
    L = length(S) # Total number of elements, works on the Type in StaticArrays
    
    field_exprs = map(keys) do k
        # Linear indexing (S[i]) works universally across all StaticArray shapes
        elements = [:(getproperty(S[$i], $(QuoteNode(k)))) for i in 1:L]
        # SArray{Sz} automatically reconstructs the original shape (vector, matrix, etc.)
        :(SArray{$Sz}($(elements...)))
    end
    
    return :(StructArray(NamedTuple{$keys}(tuple($(field_exprs...)))))
end

# Coverts a dynamically sized StructArray to a StructArray of  by dispatching on the length of the input array
function _to_static_core(sa::StructArray, ::Val{S}) where S
    cols = StructArrays.components(sa)
    
    # S is a tuple of dimensions. Tuple{S...} creates the required StaticArray size type.
    static_cols = map(c -> SArray{Tuple{S...}}(c), cols)
    
    return StructArray(static_cols)
end

"""
    StaticStructArray(sa::StructArray) -> StructArray of StaticArrays
    Convert a StructArray of dynamically sized arrays into a StructArray of StaticArrays, preserving the original shape of each field.
    This is useful for performance when the size of the arrays is known at compile time, as it allows for 
    better optimization and avoids dynamic dispatch.
"""
function StaticStructArray(sa::StructArray)
    S = size(sa)
    # Converts the runtime integer N into a compile-time type Val{N}
    return _to_static_core(sa, Val(S)) 
end



value(x::AbstractFloat) = x
value(x::ForwardDiff.Dual) = x.value  # Extract value from Dual numbers

function merge_structarrays(arrays::StructArray...)
    array_components = map(sa -> StructArrays.components(sa), arrays)
    return StructArray(merge(array_components...))
end


function flatten_structarray(sa::StructArray)
    # 1. Extract the underlying N x M component arrays
    comps = StructArrays.components(sa)
    
    # 2. Reinterpret each component from SVector{E, T} to a raw T array 
    # This automatically prepends the E dimension, resulting in E x N x M
    flat_comps = map(comps) do c
        T = eltype(eltype(c)) # Extract the scalar type (e.g., Float64)
        reinterpret(reshape, T, c)
    end
    
    # 3. Reassemble the StructArray around the new views
    return StructArray(flat_comps)
end

# Function maps

@generated function fmap(funcs::Tuple, args...)
    calls = [:(funcs[$i](args...)) for i in 1:fieldcount(funcs)]
    return :(($(calls...),))
end

function fmap(nt::NamedTuple{NAMES}, args...) where NAMES
    return NamedTuple{NAMES}(fmap(values(nt), args...))
end

# Generalization of `any(f, x)` for wide structs/NamedTuples
@generated function unrolled_any(f, x::T) where T
    N = fieldcount(T)
    N == 0 && return :(false)
    
    checks = [:(f(getfield(x, $i))) for i in 1:N]
    return reduce((a, b) -> :($a || $b), checks)
end


# Struct Array Select 
@inline select(sa::StructArray, vars::NTuple{N,Symbol} where N) = select(sa, Val(vars))

@generated function select(sa::StructArray, ::Val{NAMES}) where NAMES
    # NAMES is a tuple of symbols, e.g., (:K, :Y)
    # We construct an expression that accesses each field directly
    
    # 1. Get the components (the internal NamedTuple)
    # 2. Build the subsetting expression: (K = comps.K, Y = comps.Y)
    elements = [:(getproperty(comps, $(QuoteNode(s)))) for s in NAMES]
    
    return quote
        comps = StructArrays.components(sa)
        # Construct a new NamedTuple with the subset
        new_components = NamedTuple{$NAMES}(($(elements...),))
        # Re-wrap into a StructArray
        return StructArray(new_components)
    end
end

"""
    with_blas_threads(f, n)

Call `f()` with BLAS limited to `n` threads, then restore the previous number (also if `f` throws).
The BLAS thread count is global to the process.
"""
function with_blas_threads(f::F, n) where {F}
    original = BLAS.get_num_threads()
    BLAS.set_num_threads(n)
    try
        return f()
    finally
        BLAS.set_num_threads(original)
    end
end

"""
    threaded_foreach(f, iter; use_threads = true)

Call `f(x)` for every `x` in `iter`. With `use_threads` the calls run in parallel (`Threads.@threads`)
and BLAS is limited to one thread meanwhile, to avoid oversubscribing the cores; otherwise they run
in order and BLAS is left alone. `f` must be thread-safe and must not
depend on the order of the calls.
"""
function threaded_foreach(f::F, iter; use_threads::Bool = true) where {F}
    if use_threads
        with_blas_threads(1) do
            Threads.@threads for x in iter
                f(x)
            end
        end
    else
        for x in iter
            f(x)
        end
    end
    return nothing
end
