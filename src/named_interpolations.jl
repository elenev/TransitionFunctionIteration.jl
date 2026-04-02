"""
    StateSpaceGrid(; name=grid, ...)
    StateSpaceGrid(exogenv, endog_grids::NamedTuple)

Multidimensional state space grid with named axes.

Each keyword argument defines one axis. Scalar-valued fields are marked as
*singletons* — they are stored in the grid but excluded from interpolation axes.
The second constructor prepends the exogenous state indices from `exogenv` as the `:exst` axis.
"""
struct StateSpaceGrid{NAMES,S,N,T,D<:Tuple} <: AbstractStateSpaceGrid{T,N}
    data::D

    function StateSpaceGrid{NAMES,S}(data::D) where {NAMES,S,D}
        N = length(NAMES)
        element_types = map(eltype, fieldtypes(D))
        T = NamedTuple{NAMES, Tuple{element_types...}}
        return new{NAMES, S, N, T, D}(data)
    end
end

# Infer singletons from the data
@inline function StateSpaceGrid{NAMES}(data::D) where {NAMES, D<:Tuple }
    NT = NamedTuple{NAMES,D}
    S = singleton_keys(NT)
    return StateSpaceGrid{NAMES,S}(data)
end

# Convenience constructor from named tuple
@inline StateSpaceGrid(nt::NamedTuple{NAMES,D}) where {NAMES, D} = StateSpaceGrid{NAMES}(values(nt))

# Convenience constructor from kwargs
@inline function StateSpaceGrid(; kwargs...)
    nt = NamedTuple(kwargs)
    return StateSpaceGrid{keys(nt)}(values(nt))
end

"Return the names of scalar-valued (singleton) fields in `ssg`. Singletons are excluded from interpolation but appear in state NamedTuples."
singletons(::StateSpaceGrid{NAMES,S}) where {NAMES,S} = S

nonsingletons(::StateSpaceGrid{NAMES,S}) where {NAMES,S} = filter(n -> n ∉ S, NAMES)

"Return the NamedTuple of univariate grid vectors keyed by state variable name."
@inline unigrids(ssg::StateSpaceGrid{NAMES}) where {NAMES} = NamedTuple{NAMES}(ssg.data)

"Return the total number of grid points (product of all dimension lengths)."
@inline npt(ssg::AbstractStateSpaceGrid) = prod(dims(ssg))

"Return the number of grid points per dimension as a NamedTuple."
@inline dims(ssg::AbstractStateSpaceGrid) = map(length, unigrids(ssg))

"Return the named state dimensions of a grid, or the named output variables of an interpolant."
@inline varnames(ssg::AbstractStateSpaceGrid) = keys(unigrids(ssg))

_points(grids::NamedTuple{NAMES,T}) where {NAMES,T} = (NamedTuple{NAMES}(x) for x in Iterators.product(grids...))

"""
    points(ssg) -> Generator
    points(exogenv) -> StructArray

For a state space grid: lazy generator of NamedTuples over the full Cartesian product
(first dimension varies fastest). For an exogenous environment: the array of exogenous state vectors.
"""
points(ssg::AbstractStateSpaceGrid) = _points(unigrids(ssg))

# Array interface: Iteration
@inline Base.iterate(ssg::AbstractStateSpaceGrid) = iterate(points(ssg))
@inline Base.iterate(ssg::AbstractStateSpaceGrid, state) = iterate(points(ssg), state)
@inline Base.iterate(r::Iterators.Reverse{<:AbstractStateSpaceGrid}) = iterate(Iterators.reverse(points(r.itr)))
@inline Base.iterate(r::Iterators.Reverse{<:AbstractStateSpaceGrid}, state) = iterate(Iterators.reverse(points(r.itr)), state)

@inline Base.size(ssg::AbstractStateSpaceGrid) = values(dims(ssg))
@inline Base.length(ssg::AbstractStateSpaceGrid) = npt(ssg)
@inline Base.eltype(::Type{<:StateSpaceGrid{NAMES,S,N,T}}) where {NAMES,S,N,T} = T
@inline Base.isdone(ssg::AbstractStateSpaceGrid, state) = Base.isdone(points(ssg), state)

# Array interface: Indexing
@inline Base.IndexStyle(::Type{<:StateSpaceGrid}) = IndexCartesian()
@inline function Base.getindex(ssg::AbstractStateSpaceGrid{T,N}, I::Vararg{Int, N}) where {T,N}
    ugrids = unigrids(ssg)
    return NamedTuple{keys(ugrids)}(ntuple(i -> ugrids[i][I[i]], Val(N)))
end

@inline Base.:(==)(ssg1::AbstractStateSpaceGrid, ssg2::AbstractStateSpaceGrid) = (unigrids(ssg1) == unigrids(ssg2))

function _isapprox(v1::AbstractVector, v2::AbstractVector; kwargs...)
    if length(v1) != length(v2)
        return false
    end
    return all(isapprox(x, y; kwargs...) for (x, y) in zip(v1, v2))
end


@inline function Base.isapprox(ssg1::AbstractStateSpaceGrid, ssg2::AbstractStateSpaceGrid; kwargs...)
    ugrids1 = unigrids(ssg1)
    ugrids2 = unigrids(ssg2)

    return keys(ugrids1) == keys(ugrids2) &&
        all(Base.isapprox(v1, v2; kwargs...) for (v1, v2) in zip(values(ugrids1), values(ugrids2)))
end

"""
    Interpolant(ssg, returnvals::NamedTuple)
    Interpolant(ssg, returnvals::StructArray)

Gridded interpolation of multiple functions over a `StateSpaceGrid`.

The NamedTuple constructor fills every grid point with the same constant values.
The StructArray constructor uses per-point values; its shape must match `dims(ssg)`.

Uses `Gridded(Linear())` for endogenous dimensions and `NoInterp()` for `:exst`.
Extrapolation is `Flat()` for endogenous axes and `Throw()` for `:exst`.
"""
struct Interpolant{SSG, ITP<:Interpolations.AbstractInterpolation, NAMES} <: AbstractInterpolant
    ssg::SSG
    itp::ITP
    function Interpolant(ssg::SSG, itp::ITP, ::NamedTuple{NAMES,T}) where {SSG, ITP, NAMES,T}
        return new{SSG,ITP,NAMES}(ssg, itp)
    end
end

@inline varnames(::Interpolant{SSG,ITP,NAMES}) where {SSG, ITP, NAMES} = NAMES

# Non-singleton state variable names in the order itp expects them for positional evaluation.
function Interpolant(ssg::AbstractStateSpaceGrid, returnvals::NamedTuple)
    Y = StructArray(map(x -> deepcopy(returnvals), points(ssg)))
    return Interpolant(ssg, Y)
end

function Interpolant(ssg::AbstractStateSpaceGrid, returnvals::StructArray)
    first_values = values(first(returnvals))
    Nof = length(first_values)

    if eltype(first_values) <: AbstractArray
        # @debug eltype(first_values)
        # Output is vector-valued.
        # Must be the same length for each variable.
        all_lengths = map(length, first_values)
        if length(unique(all_lengths)) != 1
            throw(ArgumentError("All output vectors must have the same length; got lengths: $all_lengths"))
        end
        Nvals = length(first(first_values))
        #@debug "Output is vector-valued with length $Nvals; constructing SMatrix of size ($Nof, $Nvals) for each grid point."
        #@debug "Returnvals: " typeof(returnvals) size(returnvals)
        #@debug "First values: " first_values[1]
        #typed_reshape(v) = similar_type(v, Size(1, length(v)))(v)
        f_stack(x) = mapreduce(permutedims, vcat, x)
        Y = [f_stack(x) for x in returnvals]
        #@debug "Y[1...]" Y[1]
    else
        Y = [SVector(values(x)) for x in returnvals]
    end

    function interpolation_rule(statevar)
        if statevar==:exst
            return NoInterp()           # No interpolation in exog state
        else
            return Gridded(Linear())    # Linear interpolation for other dimensions
        end
    end

    function extrapolation_rule(statevar)
        if statevar==:exst
            return Throw()           # No extrapolation in exog state
        else
            return Flat()            # Constant extrapolation for other dimensions
        end
    end

    ugrids = Base.structdiff(unigrids(ssg), NamedTuple{singletons(ssg)})
    interp_rules = map(interpolation_rule, keys(ugrids))
    extrap_rules = map(extrapolation_rule, keys(ugrids))

    #@debug "Knots: " typeof(values(ugrids)) map(size, values(ugrids))
    #@debug "Values: " typeof(Y) size(Y) Y[1]
    #@debug "Interpolation rules: " interp_rules
    itp = extrapolate(interpolate(values(ugrids), Y, values(interp_rules)), values(extrap_rules) )
    example_returnvals = first(returnvals)
    return Interpolant(ssg, itp, example_returnvals)
end

@inline function _same_keys(fct::AbstractInterpolant, statevars)
    itp_keys = nonsingletons(statespacegrid(fct))
    return length(keys(statevars)) == length(itp_keys) &&
            all(∈(itp_keys), keys(statevars))
end

"""
    (fct::Interpolant)(state::NamedTuple) -> NamedTuple
    (fct::Interpolant)(; kwargs...) -> NamedTuple
    (fct::Interpolant)(pts) -> StructArray

Evaluate the interpolant at one state point or a collection of points.
Singleton dimensions in the state are silently dropped before evaluation.
Calling on a collection returns a StructArray.
"""
@inline function (fct::AbstractInterpolant)(; x...)
    return fct((;x...))
end

@inline function (fct::AbstractInterpolant)(statevars::NamedTuple)
    statevars = Base.structdiff(statevars, NamedTuple{singletons(statespacegrid(fct))})
    itp_keys = nonsingletons(statespacegrid(fct))
    if keys(statevars) === itp_keys
        return _evaluate(fct, statevars)
    elseif _same_keys(fct, statevars)
        return _evaluate(fct, NamedTuple{itp_keys}(statevars))
    else
        throw(ArgumentError("Interpolant expects state variables $(itp_keys), got $(keys(statevars))"))
    end
end

@inline function (fct::AbstractInterpolant)(pts::Union{AbstractArray,Base.Generator})
    return StructArray(_evaluate(fct,x) for x in pts)
end

@inline function _evaluate(fct, x)
    out = fct.itp(x...)
    return name_values(out, fct)
end

"Return the `StateSpaceGrid` underlying the interpolant."
@inline statespacegrid(fct::Interpolant) = fct.ssg

@inline coefficients(itp::Interpolations.GriddedInterpolation) = itp.coefs
@inline coefficients(itp::Interpolations.Extrapolation) = coefficients(itp.itp)

"Return interpolation coefficients shaped as an array matching `dims(statespacegrid(fct))`. Each element is a vector of length `nof(fct)`."
@inline coefficients(fct::Interpolant) = reshape(coefficients(fct.itp), values(dims(fct.ssg)))

@inline knots(itp::Interpolations.GriddedInterpolation) = itp.knots
@inline knots(itp::Interpolations.Extrapolation) = knots(itp.itp)

"Return the grid knots of the interpolant as a tuple of univariate vectors."
@inline knots(fct::Interpolant) = knots(fct.itp)

@inline npt(fct::AbstractInterpolant) = npt(statespacegrid(fct))

"Return the number of output functions (components) of the interpolant."
@inline nof(fct::AbstractInterpolant) = length(varnames(fct))

@inline elsize(fct::AbstractInterpolant) = (size ∘ first ∘ coefficients)(fct)

# 1. Helper to uniformly extract scalars from vectors and SVectors from matrices
@inline _extract_row(y::Real, i) = y
@inline _extract_row(y::AbstractVector, i) = @inbounds y[i]
@inline _extract_row(y::AbstractMatrix, i) = @inbounds SVector{size(y, 2)}(y[i, :])

# 2. The zero-allocation constructor
@inline function _pack_rows(::Type{<:NamedTuple{K}}, y) where {K}
    # Val(length(K)) forces the compiler to completely unroll the tuple construction
    vals = ntuple(i -> _extract_row(y, i), Val(length(K)))
    return NamedTuple{K}(vals)
end

"Wrap values `x` into a NamedTuple keyed by the output names of `fct`."
@inline name_values(y, fct::AbstractInterpolant) = _pack_rows(NamedTuple{varnames(fct)}, y)

@inline function Base.:(==)(fct1::AbstractInterpolant, fct2::AbstractInterpolant)
    return (
            (statespacegrid(fct1) == statespacegrid(fct2)) &&
            (varnames(fct1) == varnames(fct2)) &&
            all(coef1 == coef2 for (coef1, coef2) in zip(coefficients(fct1), coefficients(fct2)))
    )
end

@inline function Base.isapprox(fct1::AbstractInterpolant, fct2::AbstractInterpolant; kwargs...)
    return varnames(fct1) == varnames(fct2) &&
        Base.isapprox(statespacegrid(fct1), statespacegrid(fct2); kwargs...) &&
        all(Base.isapprox(coef1, coef2; kwargs...) for (coef1, coef2) in zip(coefficients(fct1), coefficients(fct2)))
end


"""
    Exogenv(; pts, mtrans, ssidx)

Exogenous environment for a discrete Markov chain.

- `pts`: StructArray of exogenous state vectors (one row per state)
- `mtrans`: transition probability matrix (rows sum to 1)
- `ssidx`: index of the steady-state (default starting) state
"""
@kwdef struct Exogenv{V<:StructArray,M} <: AbstractExogenv
    pts::V
    mtrans::M
    ssidx::Int
    function Exogenv(pts::V, mtrans::M, ssidx::Int) where {V<:StructArray,M}
        pts = StaticStructArray(pts) # Convert to SVector for performance
        if length(pts) <= 10
            mtrans = SMatrix{size(mtrans)...}(mtrans) # Convert to SMatrix for performance
            # not worth it for large exnpt
        end
        return new{typeof(pts),typeof(mtrans)}(pts, mtrans, ssidx)
    end
end

"Return the number of variables in each exogenous state vector."
@inline nex(exogenv::Exogenv) = length(exogenv.pts[1])

"Return the number of discrete exogenous states."
@inline exnpt(exogenv::Exogenv) = length(exogenv.pts)

@inline varnames(exogenv::Exogenv) = propertynames(exogenv.pts)

@inline points(exogenv::Exogenv) = exogenv.pts
@inline transition_matrix(exogenv::Exogenv) = exogenv.mtrans
@inline conditional_distro(exogenv::AbstractExogenv, state_idx) = view(transition_matrix(exogenv), state_idx, :)
#@inline conditional_distro(exogenv::Exogenv, state_idx) = transition_matrix(exogenv)[state_idx, :]


"Return the steady-state index of the exogenous environment."
@inline ssidx(exogenv::Exogenv) = exogenv.ssidx

"Return a StructArray of exogenous points given an AbstractArray of exogenous state integers."
function points(exogenv::AbstractExogenv, states::AbstractArray{<:Integer})
    return StructArray(points(exogenv)[states])
end

@inline function Base.:(==)(env1::Exogenv, env2::Exogenv)
    return (points(env1) == points(env2)) &&
        (transition_matrix(env1) == transition_matrix(env2)) &&
        (ssidx(env1) == ssidx(env2))
end
@inline function Base.isapprox(env1::Exogenv, env2::Exogenv; kwargs...)
    return exnpt(env1) == exnpt(env2) &&
        propertynames(points(env1)) == propertynames(points(env2)) &&
        all(Base.isapprox(p1, p2; kwargs...) for (p1, p2) in zip(components(points(env1)), components(points(env2)))) &&
        Base.isapprox(transition_matrix(env1), transition_matrix(env2); kwargs...) &&
        (ssidx(env1) == ssidx(env2))
end

function StateSpaceGrid(exogenv::Exogenv, endog_grids::NamedTuple)
    return StateSpaceGrid(; exst = eachindex(exogenv.pts), endog_grids...)
end
