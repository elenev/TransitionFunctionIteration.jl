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
    Pfct_bounds::PB   = NamedTuple{varnames(Pfct)}(ntuple(_ -> Unbounded(), nof(Pfct)))
    Tfct_bounds::TB   = NamedTuple{varnames(Tfct)}(ntuple(_ -> Unbounded(), nof(Tfct)))
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



# Steady state solver modes

struct AnalyticalSteadyStateMode <: SteadyStateMode end

_exogenv_from_point(pt) = Exogenv(pts = StaticStructArray(SVector(pt)), mtrans = [1.;;], ssidx=1)
struct FixedPointSteadyStateMode{P<:NamedTuple, T<:NamedTuple, E<:AbstractExogenv} <: SteadyStateMode
    policies::P
    transitions::T
    exogenv::E

    function FixedPointSteadyStateMode(policies::P,
        transitions::T, exogenv::E) where {P<:NamedTuple, T<:NamedTuple, E<:AbstractExogenv}
        N = exnpt(exogenv)
        if N > 1
            exog_mean = mean(exogenv)
            exogenv = _exogenv_from_point(exog_mean)
        end
        return new{P,T,typeof(exogenv)}(policies, transitions, exogenv)
    end
end

function FixedPointSteadyStateMode(policies::NamedTuple, transitions::NamedTuple, exogenv::NamedTuple)
    return FixedPointSteadyStateMode(policies, transitions, _exogenv_from_point(exogenv))
end

struct SystemSteadyStateMode{S<:NamedTuple, B<:NamedTuple} <: SteadyStateMode
    ssguess::S
    bounds::B
end
function SystemSteadyStateMode(ssguess::NamedTuple)
    return SystemSteadyStateMode(ssguess, map(_ -> Unbounded(), ssguess))
end

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

"Solution method: Transition Function Iteration (TFI) of Elenev and Landvoigt."
struct TFIMethod <: SolutionMethod end

"Solution method: Simultaneous Transition and Policy Function Iteration (STPFI) of Cao, Luo, and Nie."
struct STPFIMethod <: SolutionMethod end


# Interface stubs requiring NextPeriod, Exogenv, or SteadyStateMode types

"""
    calc_eqm(state, policy, mdl, np, mode)

Compute equilibrium conditions at `state` given `policy`. Interface stub --
model packages must add methods dispatching on their `AbstractModel` subtype.

Expected return: a NamedTuple with fields `fx` (residuals), `trans` (next-period
transient state), and `vals` (value outputs and extra simulation variables in SimulateMode).
"""
function calc_eqm(state, policy, mdl::AbstractModel, np::NextPeriod, mode::EvaluationMode)
    throw(MethodNotImplementedError(:calc_eqm, typeof(mdl)))
end

function compute_next_period_values(next, mdl::AbstractModel, mode::EvaluationMode)::NextPeriod
    throw(MethodNotImplementedError(:compute_next_period_values, typeof(mdl)))
end

function Exogenv(mdl::AbstractModel)
    throw(MethodNotImplementedError(:Exogenv, typeof(mdl)))
end

function steadystate(mdl::AbstractModel, ::AnalyticalSteadyStateMode)
    throw(MethodNotImplementedError(:steadystate, typeof(mdl), " for AnalyticalSteadyStateMode"))
end
