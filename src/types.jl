"Abstract base for multidimensional state space grids."
abstract type AbstractStateSpaceGrid end

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
struct StateSpaceGrid{N,T,S} <: AbstractStateSpaceGrid
    data::NamedTuple{N,T}

    function StateSpaceGrid(data::NamedTuple{N,T}) where {N,T}
        S = singleton_keys(typeof(data))
        new{N,T,S}(data)
    end
end

"Return the names of scalar-valued (singleton) fields in `ssg`. Singletons are excluded from interpolation but appear in state NamedTuples."
singletons(::StateSpaceGrid{N,T,S}) where {N,T,S} = S

StateSpaceGrid(; kwargs...) = StateSpaceGrid(NamedTuple{keys(kwargs)}(values(kwargs)))

"Return the NamedTuple of univariate grid vectors keyed by state variable name."
@inline unigrids(ssg::AbstractStateSpaceGrid) = ssg.data

"Return the total number of grid points (product of all dimension lengths)."
@inline npt(ssg::AbstractStateSpaceGrid) = prod(dims(ssg))

@inline ndims(ssg::AbstractStateSpaceGrid) = length(keys(unigrids(ssg)))

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
    Y = [MVector(values(x)) for x in returnvals]

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

function (fct::Interpolant)(statevars::NamedTuple)
    statevars = Base.structdiff(statevars, NamedTuple{singletons(statespacegrid(fct))})
    return _evaluate(fct, statevars)
end

function (fct::Interpolant)(pts::Union{AbstractArray,Base.Generator})
    return StructArray(_evaluate(fct,x) for x in pts)
end

function _evaluate(fct, x)
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

"Wrap values `x` into a NamedTuple keyed by the output names of `fct`."
@inline name_values(x, fct::AbstractInterpolant) = NamedTuple{varnames(fct)}(x)

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
end

"Return the number of variables in each exogenous state vector."
@inline nex(exogenv::AbstractExogenv) = length(exogenv.pts[1])

"Return the number of discrete exogenous states."
@inline exnpt(exogenv::AbstractExogenv) = size(exogenv.pts,1)

@inline points(exogenv::AbstractExogenv) = exogenv.pts
@inline transition_matrix(exogenv::AbstractExogenv) = exogenv.mtrans

"Return the steady-state index of the exogenous environment."
@inline ssidx(exogenv::AbstractExogenv) = exogenv.ssidx

"Return a StructArray of exogenous points given an AbstractArray of exogenous state integers."
function points(exogenv::AbstractExogenv, states::AbstractArray{<:Integer})
    return StructArray(exogenv.pts[states])
end

"Abstract base for user-supplied economic model logic. Implement `calc_eqm` for each subtype."
abstract type AbstractModel end

"""
    calc_eqm(state, policy, obj; mode=SolveMode())

Compute equilibrium conditions at `state` given `policy`. Interface stub —
model packages must add methods dispatching on `ModelSolution{<:MyModel}`.

Expected return: a NamedTuple with fields `vals` (value outputs), `fx` (residuals),
`trans` (next-period transient state), and `addvars` (extra simulation variables).
"""
function calc_eqm end  # interface stub; model packages add methods dispatching on ModelSolution{M}

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
            VF<:AbstractInterpolant} <: AbstractModelSolution
    model::M
    exogenv::E
    ssg::S
    Pfct::PF
    Vfct::VF
end

"Abstract dispatch tag for evaluation context."
abstract type EvaluationMode end

"Evaluation context: Solving for equilibrium policy functions."
struct SolveMode <: EvaluationMode end

"Evaluation context: Iterating (updating) the value function."
struct IterateMode <: EvaluationMode end

"Evaluation context: Simulating the model."
struct SimulateMode <: EvaluationMode end

function StateSpaceGrid(exogenv::AbstractExogenv, endog_grids::NamedTuple)
    return StateSpaceGrid(; exst = eachindex(exogenv.pts), endog_grids...)
end
