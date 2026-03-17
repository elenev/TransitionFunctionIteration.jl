"Abstract base for multidimensional state space grids."
abstract type AbstractStateSpaceGrid{T, N} <: AbstractArray{T, N} end

@generated singleton_keys(::Type{NamedTuple{N,T}}) where {N,T} =
    tuple((N[i] for i in 1:length(N) if fieldtype(T, i) <: Number)...)

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

"Abstract base for interpolated functions over a state space."
abstract type AbstractInterpolant end

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
        @debug "Output is vector-valued with length $Nvals; constructing SMatrix of size ($Nof, $Nvals) for each grid point."
        Y = [SMatrix{Nof,Nvals}(reduce(vcat, x)) for x in returnvals]
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

"""
    (fct::Interpolant)(state::NamedTuple) -> NamedTuple
    (fct::Interpolant)(; kwargs...) -> NamedTuple
    (fct::Interpolant)(pts) -> StructArray

Evaluate the interpolant at one state point or a collection of points.
Singleton dimensions in the state are silently dropped before evaluation.
Calling on a collection returns a StructArray.
"""
@inline function (fct::Interpolant)(; x...)
    return fct((;x...))
end

@inline function (fct::Interpolant)(statevars::NamedTuple)
    statevars = Base.structdiff(statevars, NamedTuple{singletons(statespacegrid(fct))})
    return _evaluate(fct, statevars)
end

@inline function (fct::Interpolant)(pts::Union{AbstractArray,Base.Generator})
    return StructArray(_evaluate(fct,x) for x in pts)
end

@inline function _evaluate(fct, x)
    out = fct.itp(x...)
    return name_values(out, fct)
end

"Return the `StateSpaceGrid` underlying the interpolant."
@inline statespacegrid(fct::AbstractInterpolant) = fct.ssg

@inline coefficients(itp::Interpolations.GriddedInterpolation) = itp.coefs
@inline coefficients(itp::Interpolations.Extrapolation) = coefficients(itp.itp)

"Return interpolation coefficients shaped as an array matching `dims(statespacegrid(fct))`. Each element is a vector of length `nof(fct)`."
@inline coefficients(fct::AbstractInterpolant) = reshape(coefficients(fct.itp), values(dims(fct.ssg)))

@inline knots(itp::Interpolations.GriddedInterpolation) = itp.knots
@inline knots(itp::Interpolations.Extrapolation) = knots(itp.itp)

"Return the grid knots of the interpolant as a tuple of univariate vectors."
@inline knots(fct::AbstractInterpolant) = knots(fct.itp)

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


"Abstract base for exogenous environments (Markov chains)."
abstract type AbstractExogenv end

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

@inline varnames(exogenv::Exogenv) = varnames(eltype(exogenv.pts))

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
    return propertynames(points(env1)) == propertynames(points(env2)) && 
        all(Base.isapprox(p1, p2; kwargs...) for (p1, p2) in zip(components(points(env1)), components(points(env2)))) &&
        Base.isapprox(transition_matrix(env1), transition_matrix(env2); kwargs...) &&
        (ssidx(env1) == ssidx(env2))
end


@kwdef struct NextPeriod{E,C,A,V<:AbstractVector}
    expectations::E
    continuations::C
    addvars::A = (;)
    prnext::V  = SVector{1}(1.0)
end

E(x, np::NextPeriod) = dot(np.prnext, x)



### Transformations Into the Real Line 
abstract type Bounded end
Base.broadcastable(b::Bounded) = Ref(b)


struct Unbounded <: Bounded end
@inline Base.:(==)(::Unbounded, ::Unbounded) = ( true )

struct LowerBounded{L<:Real} <: Bounded
    lower::L
end
PositiveBounded() = LowerBounded(0.0)
@inline Base.:(==)(b1::LowerBounded, b2::LowerBounded) = (b1.lower == b2.lower)

struct UpperBounded{U<:Real} <: Bounded
    upper::U
end
NegativeBounded() = UpperBounded(0.0)
@inline Base.:(==)(b1::UpperBounded, b2::UpperBounded) = (b1.upper == b2.upper)

struct IntervalBounded{L<:Real,U<:Real} <: Bounded
    lower::L
    upper::U
end
IntervalBounded() = IntervalBounded(0.0, 1.0)
Base.:(==)(b1::IntervalBounded, b2::IntervalBounded) = ((b1.lower == b2.lower) && (b1.upper == b2.upper))

to_constrained(x, ::Unbounded) = x
to_constrained(x, b::LowerBounded) = exp(x) + b.lower
to_constrained(x, b::UpperBounded) = -exp(x) + b.upper
to_constrained(x, b::IntervalBounded) = (b.upper - b.lower) / (1 + exp(-x)) + b.lower

to_unconstrained(y, ::Unbounded) = y
to_unconstrained(y, b::LowerBounded) = log(y - b.lower)
to_unconstrained(y, b::UpperBounded) = log(b.upper - y)
to_unconstrained(y, b::IntervalBounded) = log((y - b.lower) / (b.upper - y))





"Abstract base for user-supplied economic model logic. Implement `calc_eqm` for each subtype."
abstract type AbstractModel end

"Abstract base for assembled solution objects."
abstract type AbstractModelSolution end

"""
    ModelSolution(; model, exogenv, ssg, Pfct, Vfct)

Assembled solution container for an economic model.

- `model`: user-supplied `AbstractModel` instance (economics logic)
- `exogenv`: `Exogenv` for the Markov exogenous process
- `ssg`: combined `StateSpaceGrid` (exogenous + endogenous axes)
- `Pfct`: policy function `Interpolant`
- `Vfct`: value function `Interpolant`
"""
@kwdef struct ModelSolution{M<:AbstractModel,
            E<:AbstractExogenv,
            S<:AbstractStateSpaceGrid,
            PF<:AbstractInterpolant,
            VF<:AbstractInterpolant,
            TF<:AbstractInterpolant,
            PB<:NamedTuple, 
            TB<:NamedTuple} <: AbstractModelSolution
    model::M
    exogenv::E
    ssg::S
    Pfct::PF
    Vfct::VF
    Tfct::TF
    Pfct_bounds::PB 
    Tfct_bounds::TB
end

function Base.:(==)(obj1::ModelSolution, obj2::ModelSolution)
    return (
        (obj1.model == obj2.model) &&
        (obj1.exogenv == obj2.exogenv) &&
        (obj1.ssg == obj2.ssg) &&
        (obj1.Pfct == obj2.Pfct) &&
        (obj1.Vfct == obj2.Vfct) &&
        (obj1.Tfct == obj2.Tfct) &&
        (obj1.Pfct_bounds == obj2.Pfct_bounds) &&
        (obj1.Tfct_bounds == obj2.Tfct_bounds)
    )
end

function Base.isapprox(obj1::ModelSolution, obj2::ModelSolution; kwargs...)
    return (
        (obj1.model == obj2.model) &&
        Base.isapprox(obj1.exogenv, obj2.exogenv; kwargs...) &&
        Base.isapprox(obj1.ssg, obj2.ssg; kwargs...) &&
        Base.isapprox(obj1.Pfct, obj2.Pfct; kwargs...) &&
        Base.isapprox(obj1.Vfct, obj2.Vfct; kwargs...) &&
        Base.isapprox(obj1.Tfct, obj2.Tfct; kwargs...) &&
        (obj1.Pfct_bounds == obj2.Pfct_bounds) &&
        (obj1.Tfct_bounds == obj2.Tfct_bounds)
    )
end


abstract type AbstractSteadyState end

@kwdef struct SteadyState{M<:AbstractModel, GNAMES, GTYPES, XNAMES, XTYPES, VNAMES, VTYPES, NNAMES, NTYPES} <: AbstractSteadyState
    model::M
    guess::NamedTuple{GNAMES,GTYPES}
    endog::NamedTuple{NNAMES,NTYPES}
    exog::NamedTuple{XNAMES,XTYPES}
    values::NamedTuple{VNAMES,VTYPES}
end

state(ss::SteadyState) = merge(ss.exog, ss.endog)

@inline function Base.:(==)(ss1::SteadyState, ss2::SteadyState)
    return (ss1.model == ss2.model) &&
        (ss1.guess == ss2.guess) &&
        (ss1.endog == ss2.endog) &&
        (ss1.exog == ss2.exog) &&
        (ss1.values == ss2.values)
end
@inline function Base.isapprox(ss1::SteadyState, ss2::SteadyState; kwargs...)
    return (ss1.model == ss2.model) &&
        keys(ss1.guess) == keys(ss2.guess) &&
        all(Base.isapprox(v1, v2; kwargs...) for (v1, v2) in zip(values(ss1.guess), values(ss2.guess))) &&
        keys(ss1.endog) == keys(ss2.endog) &&
        all(Base.isapprox(v1, v2; kwargs...) for (v1, v2) in zip(values(ss1.endog), values(ss2.endog))) &&
        keys(ss1.exog) == keys(ss2.exog) &&
        all(Base.isapprox(v1, v2; kwargs...) for (v1, v2) in zip(values(ss1.exog), values(ss2.exog))) &&
        keys(ss1.values) == keys(ss2.values) &&
        all(Base.isapprox(v1, v2; kwargs...) for (v1, v2) in zip(values(ss1.values), values(ss2.values)))
end


"Abstract dispatch tag for evaluation context."
abstract type EvaluationMode end

"Evaluation context: Solving for equilibrium policy functions."
struct SolveMode <: EvaluationMode end

"Evaluation context: Iterating (updating) the value function."
struct IterateMode <: EvaluationMode end

"Evaluation context: Simulating the model."
struct SimulateMode <: EvaluationMode end

"Abstract dispatch tag for solution method."
abstract type SolutionMethod end

"Solution method: Transition Function Iteration (TFI) of Elenev and Landvoigt."
struct TFIMethod <: SolutionMethod end

"Solution method: Simultaneous Transition and Policy Function Iteration (STPFI) of Cao, Luo, and Nie."
struct STPFIMethod <: SolutionMethod end

function StateSpaceGrid(exogenv::AbstractExogenv, endog_grids::NamedTuple)
    return StateSpaceGrid(; exst = eachindex(exogenv.pts), endog_grids...)
end



## Simulation Stats

struct MomentArray{T, R <: NamedTuple, C <: NamedTuple, M <: AbstractMatrix{T}} <: AbstractMatrix{T}
    data::M
    row_idx::R  # e.g., (exst = 1, K = 2)
    col_idx::C  # e.g., (mean = 1, std = 2)
end

function MomentArray(data::AbstractMatrix, rows::Tuple, cols::Tuple)
    # Lock the symbol tuples into NamedTuples mapped to their integer indices
    r_map = NamedTuple{rows}(ntuple(identity, length(rows)))
    c_map = NamedTuple{cols}(ntuple(identity, length(cols)))
    
    return MomentArray(data, r_map, c_map)
end


# --- AbstractArray Interface ---
Base.size(m::MomentArray) = size(m.data)
Base.IndexStyle(::Type{<:MomentArray}) = IndexLinear()

# Pass-through linear and cartesian indexing
@inline Base.getindex(m::MomentArray, i::Int) = getindex(m.data, i)
@inline Base.getindex(m::MomentArray, I::Vararg{Int, 2}) = getindex(m.data, I...)

# Pass-through mutation (optional, but highly recommended)
@inline Base.setindex!(m::MomentArray, v, i::Int) = setindex!(m.data, v, i)
@inline Base.setindex!(m::MomentArray, v, I::Vararg{Int, 2}) = setindex!(m.data, v, I...)

# ---  Beautiful REPL Printing ---
function Base.show(io::IO, ::MIME"text/plain", m::MomentArray)
    println(io, "MomentArray ($(size(m, 1)) × $(size(m, 2))):")
    pretty_table(
        io, 
        m.data; 
        column_labels = collect(keys(m.col_idx)), 
        row_labels = collect(keys(m.row_idx)),
    )
end

# This method handles "compact" printing (inside tuples/arrays/containers)
function Base.show(io::IO, m::MomentArray)
    print(io, "MomentArray{$(eltype(m))}($(size(m,1)) × $(size(m,2)))")
end

# --- 4. The getmoment Dispatch ---
# Base Case: Scalar / Scalar -> Number
@inline getmoment(m::MomentArray, var::Symbol, stat::Symbol) = 
    m.data[m.row_idx[var], m.col_idx[stat]]

_expand_colon(needle, haystack) = (needle isa Colon) ? haystack : needle

# Collection / Scalar -> NamedTuple
function getmoment(m::MomentArray, vars::Union{Tuple, AbstractVector, Colon}, stat::Symbol)
    vars = _expand_colon(vars, keys(m.row_idx))
    keys_tup = Tuple(vars)
    vals_tup = Tuple(getmoment(m, v, stat) for v in keys_tup)
    return map(v -> getmoment(m, v, stat), vars)
end

# Scalar / Collection -> NamedTuple
function getmoment(m::MomentArray, var::Symbol, stats::Union{Tuple, AbstractVector, Colon})
    stats = _expand_colon(stats, keys(m.col_idx))
    keys_tup = Tuple(stats)
    vals_tup = map(s -> getmoment(m, var, s), stats)
    return NamedTuple{keys_tup}(vals_tup)
end

# Collection / Collection -> MomentArray
@inline function getmoment(m::MomentArray, ::Val{vars}, ::Val{stats}) where {vars, stats}
    # These lookups now happen at COMPILE TIME because 'vars' is in the type header
    row_inds = SVector(map(v -> m.row_idx[v], vars))
    col_inds = SVector(map(s -> m.col_idx[s], stats))
    
    # The return type is now perfectly predictable
    return MomentArray(view(m.data, row_inds, col_inds), vars, stats)
end

# 3. The ergonomic wrapper 
@inline getmoment(m::MomentArray, vars::Tuple, stats::Tuple) = 
    getmoment(m, Val(vars), Val(stats))

## Stubs

# 1. Define the custom exception type
struct MethodNotImplementedError <: Exception
    func_name::Symbol
    type_missing::Type
end

# 2. Overload showerror for clean REPL output
function Base.showerror(io::IO, e::MethodNotImplementedError)
    print(io, "MethodNotImplementedError: Solver requires that `", e.func_name, 
              "` be implemented for `", e.type_missing, "`.")
end

"""
    calc_eqm(state, policy, obj; mode=SolveMode())

Compute equilibrium conditions at `state` given `policy`. Interface stub —
model packages must add methods dispatching on `ModelSolution{<:MyModel}`.

Expected return: a NamedTuple with fields `vals` (value outputs), `fx` (residuals),
`trans` (next-period transient state), and `addvars` (extra simulation variables).
"""
function calc_eqm(state, policy, mdl::AbstractModel, np::NextPeriod, mode::EvaluationMode) 
    # Not Implemented Error 
    throw(MethodNotImplementedError(:calc_eqm, typeof(mdl)))
end

function compute_next_period_values(next, mdl::AbstractModel)::NextPeriod
    throw(MethodNotImplementedError(:compute_next_period_values, typeof(mdl)))
end

function Exogenv(mdl::AbstractModel)::Exogenv
    throw(MethodNotImplementedError(:Exogenv, typeof(mdl)))
end

function define_endog_grids(mdl::AbstractModel, ::AbstractSteadyState)
    throw(MethodNotImplementedError(:define_endog_grids, typeof(mdl)))
end

# This is only necessary if steadystate is not implemented.
# In that case, the fallback steadystate will use this to construct next period vals.
# Users can also use it in their own calc_eqm to create the IterateMode() branch output.
function vals_from_policy(state, policy, mdl::AbstractModel, args...)
    throw(MethodNotImplementedError(:vals_from_policy, typeof(mdl)))
end

# This is only necessary if steady state is not implemented. Must return a tuple (policy_guess, endog_guess, exogenv)
# of NamedTuples keyed by the policy variable names and endogenous state variable names, respectively.
# Exogenv should contain the 1-state MC in deterministic steady state.
# It can also be nothing, in which case the default Exogenv(mdl) will be used to construct a constant MC at the ergodic mean.
function steadystate_guess(mdl::AbstractModel)
    throw(MethodNotImplementedError(:steadystate_guess, typeof(mdl)))
end

