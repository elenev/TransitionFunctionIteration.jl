### Transformations Into the Real Line
"""
    Bounded

Abstract type for variable bounds. Used to specify the transformation of policy function variables 
into an unconstrained space for the solver. Custom concrete types can be defined. They must
implement  `to_constrained` and `to_unconstrained` for the transformations.
"""
abstract type Bounded end
Base.broadcastable(b::Bounded) = Ref(b)

"""
    Unbounded <: Bounded
No bounds; identity transformation.
"""
struct Unbounded <: Bounded end
@inline Base.:(==)(::Unbounded, ::Unbounded) = ( true )

"""
    LowerBounded{L<:Real} <: Bounded
Bounded below by `lower`. Transformed to unconstrained space by `to_constrained(x, b) = exp(x) + b.lower`.
"""
struct LowerBounded{L<:Real} <: Bounded
    lower::L
end

"""
    PositiveBounded() := LowerBounded(0.0)

Bounded below by `0`. Transformed to unconstrained space by `to_constrained(x, b) = exp(x)`.
"""
PositiveBounded() = LowerBounded(0.0)

@inline Base.:(==)(b1::LowerBounded, b2::LowerBounded) = (b1.lower == b2.lower)

"""
    UpperBounded{U<:Real} <: Bounded
Bounded above by `upper`. Transformed to unconstrained space by `to_constrained(x, b) = -exp(x) + b.upper`.
"""
struct UpperBounded{U<:Real} <: Bounded
    upper::U
end
"""
    NegativeBounded() := UpperBounded(0.0)

Bounded above by `0`. Transformed to unconstrained space by `to_constrained(x, b) = -exp(x)`.
"""
NegativeBounded() = UpperBounded(0.0)
@inline Base.:(==)(b1::UpperBounded, b2::UpperBounded) = (b1.upper == b2.upper)

"""
    IntervalBounded{L<:Real = 0.0, U<:Real = 1.0} <: Bounded

Bounded between `lower` and `upper`. Transformed to unconstrained space by `to_constrained(x, b) = (b.upper - b.lower) / (1 + exp(-x)) + b.lower`.
"""
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


"""
    NextPeriod(; expectations, continuations, prnext=SVector{1}(1.0))

Container for next-period information returned by `compute_next_period_values`.

- `expectations`: probability-weighted averages (e.g., `E[β u′(C′)]`); scalar for a given grid point
- `continuations`: per-exogenous-state raw values (e.g., capital payoff each state); stacked by `compute_expectations` into an exnpt-length vector
- `prnext`: transition probability vector (default `[1.0]` for a single state, overridden by `Exogenv`)
"""
@kwdef struct NextPeriod{E,C,V<:AbstractVector}
    expectations::E
    continuations::C
    prnext::V  = SVector{1}(1.0)
end

_mean(f, x, w) = mapreduce((wi, xi) -> wi * f(xi), +, w, x)
@inline Statistics.mean(x, np::NextPeriod) = dot(np.prnext, x)
@inline Statistics.mean(f::Function, x, np::NextPeriod) = _mean(f, x, np.prnext)


"""
    ModelSolution(; model, exogenv, ssg, Pfct, Vfct)

Assembled solution container for an economic model.

- `model`: user-supplied `AbstractModel` instance (economics logic)
- `exogenv`: `Exogenv` for the Markov exogenous process
- `ssg`: combined `StateSpaceGrid` (exogenous + endogenous axes)
- `Pfct`: policy function `Interpolant`
- `Vfct`: value function `Interpolant`
- `Tfct`: transition function `Interpolant`
- `Ffct`: forecasting function `Interpolant`, or `nothing` (default) if the model has none.
  Forecasting values are current-period choices that next period's expectations depend on
  beyond the next-period state (e.g. the portfolio entering a next-period consumption stage).
  Unlike `Vfct`, `Ffct` is evaluated at the current state; its values are attached to every
  next-period point by `forecast_next`. See `forecasts_from_policy`.
"""
@kwdef struct ModelSolution{M<:AbstractModel,
            E<:AbstractExogenv,
            S<:AbstractStateSpaceGrid,
            PF<:AbstractInterpolant,
            VF<:AbstractInterpolant,
            TF<:AbstractInterpolant,
            PB<:NamedTuple,
            TB<:NamedTuple,
            FF<:Union{Nothing,AbstractInterpolant}} <: AbstractModelSolution
    model::M
    exogenv::E
    ssg::S
    Pfct::PF
    Vfct::VF
    Tfct::TF
    Pfct_bounds::PB   = NamedTuple{varnames(Pfct)}(ntuple(_ -> Unbounded(), nof(Pfct)))
    Tfct_bounds::TB   = NamedTuple{varnames(Tfct)}(ntuple(_ -> Unbounded(), nof(Tfct)))
    Ffct::FF          = nothing
end

_isapprox_or_nothing(::Nothing, ::Nothing; kwargs...) = true
_isapprox_or_nothing(a, b; kwargs...) = !isnothing(a) && !isnothing(b) && Base.isapprox(a, b; kwargs...)

function Base.:(==)(obj1::ModelSolution, obj2::ModelSolution)
    return (
        (obj1.model == obj2.model) &&
        (obj1.exogenv == obj2.exogenv) &&
        (obj1.ssg == obj2.ssg) &&
        (obj1.Pfct == obj2.Pfct) &&
        (obj1.Vfct == obj2.Vfct) &&
        (obj1.Tfct == obj2.Tfct) &&
        (obj1.Pfct_bounds == obj2.Pfct_bounds) &&
        (obj1.Tfct_bounds == obj2.Tfct_bounds) &&
        (obj1.Ffct == obj2.Ffct)
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
        (obj1.Tfct_bounds == obj2.Tfct_bounds) &&
        _isapprox_or_nothing(obj1.Ffct, obj2.Ffct; kwargs...)
    )
end



# Steady state solver modes

"""
    AnalyticalSteadyStateMode <: SteadyStateMode

Dispatch tag: steady state is computed by an analytic formula in `steadystate(model, ::AnalyticalSteadyStateMode)`.
To use this, implement a method for `steadystate` dispatching on your `AbstractModel` subtype and `AnalyticalSteadyStateMode()`.
"""
struct AnalyticalSteadyStateMode <: SteadyStateMode end

_exogenv_from_point(pt) = Exogenv(pts = StaticStructArray(SVector(pt)), mtrans = [1.;;], ssidx=1)


"""
    FixedPointSteadyStateMode{P<:NamedTuple, T<:NamedTuple, E<:AbstractExogenv, G<:NamedTuple} <: SteadyStateMode

Dispatch tag: steady state is found by iterating on the dynamic model equations and (constant) law of motion 
for endogenous variables with a starting guess for policy and transition functions."

- `policies`: `NamedTuple` of policy function values at the steady state (e.g., `(C = 0.8,)`)
- `transitions`: `NamedTuple` of transition function values at the steady state (e.g., `(K = 3.2,)`)
- `exogenv`: `Exogenv` representing the exogenous environment at the steady state. Can also be a 
  `NamedTuple` of exogenous variable values (e.g., `(A = 1.0,)`), which will be converted to a degenerate `Exogenv` with a single point.
- `eqmvalguesses` (keyword, default empty): guesses for Vfct variables that `vals_from_policy`
  cannot compute because they depend on expectations (e.g. a value function). They become extra
  unknowns; their residuals are the guess minus the value `calc_eqm` returns in `vals`.
"""
struct FixedPointSteadyStateMode{P<:NamedTuple, T<:NamedTuple, E<:AbstractExogenv, G<:NamedTuple} <: SteadyStateMode
    policies::P
    transitions::T
    exogenv::E
    eqmvalguesses::G

    function FixedPointSteadyStateMode(policies::P,
        transitions::T, exogenv::E; eqmvalguesses::G = (;)) where {P<:NamedTuple, T<:NamedTuple, E<:AbstractExogenv, G<:NamedTuple}
        N = exnpt(exogenv)
        if N > 1
            exog_mean = mean(exogenv)
            exogenv = _exogenv_from_point(exog_mean)
        end
        return new{P,T,typeof(exogenv),G}(policies, transitions, exogenv, eqmvalguesses)
    end
end

function FixedPointSteadyStateMode(policies::NamedTuple, transitions::NamedTuple, exogenv::NamedTuple; eqmvalguesses = (;))
    return FixedPointSteadyStateMode(policies, transitions, _exogenv_from_point(exogenv); eqmvalguesses)
end

"""
    SystemSteadyStateMode{S<:NamedTuple, B<:NamedTuple} <: SteadyStateMode

Dispatch tag: steady state is found by solving a steadystate specific nonlinear system (e.g., system of equilibrium conditions).
Must implement `steadystate_system(guess, model, mode)` which returns the system residuals given a guess for the steady state values.
"""
struct SystemSteadyStateMode{S<:NamedTuple, B<:NamedTuple} <: SteadyStateMode
    ssguess::S
    bounds::B
end
function SystemSteadyStateMode(ssguess::NamedTuple)
    return SystemSteadyStateMode(ssguess, map(_ -> Unbounded(), ssguess))
end

"""
    SteadyState(; model, exog, endog, values, guess)

Steady-state description returned by `steadystate`. Holds model-specific values at the non-stochastic
steady state and serves as the starting point for `main_setup`.

- `model`: the `AbstractModel` instance
- `exog`: exogenous state at the steady state (e.g., `(A=1.0,)`)
- `endog`: endogenous state at the steady state (e.g., `(K=3.2,)`)
- `values`: additional reported values (e.g., `(C=0.8,)`)
- `guess`: initial policy guess for the solver
"""
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


"Evaluation context: Solving for equilibrium policy functions."
struct SolveMode <: EvaluationMode end

"Evaluation context: Iterating (updating) the value function."
struct IterateMode <: EvaluationMode end

"Evaluation context: Simulating the model."
struct SimulateMode <: EvaluationMode end

"Solution method: Transition Function Iteration (TFI) of Elenev, Landvoigt, and Van Nieuwerburgh (2021)."
struct TFIMethod <: SolutionMethod end

"Solution method: Simultaneous Transition and Policy Function Iteration (STPFI) of Cao, Luo, and Nie."
struct STPFIMethod <: SolutionMethod end


