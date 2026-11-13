## Interface stubs: model packages must implement these for their AbstractModel subtype, or the default method will throw a MethodNotImplementedError.

"""
    define_endog_grids(model, stv) -> NamedTuple

Return a `NamedTuple` of endogenous grid vectors (e.g., `(K = range(0.5, 10.0, 50),)`).
Called by `main_setup` to build the `StateSpaceGrid`. Interface stub -- implement for your
`AbstractModel` subtype when the default steady-state-proportional grid is not appropriate.
"""
function define_endog_grids(mdl::AbstractModel, ::AbstractSteadyState)
    throw(MethodNotImplementedError(:define_endog_grids, typeof(mdl)))
end

"""
    steadystate(model, mode) -> (SteadyState, Bool)

Compute the steady state for `model` using the given `mode` (e.g., `AnalyticalSteadyStateMode()`).
Interface stub -- model packages can add methods dispatching on their subtype to overload default
steady-state solving behaviors.

Returns `(stv, converged)` where `stv::SteadyState` and `converged::Bool` indicates success.
"""
function steadystate(mdl::AbstractModel, ::AnalyticalSteadyStateMode)
    throw(MethodNotImplementedError(:steadystate, typeof(mdl), " for AnalyticalSteadyStateMode"))
end


"""
    steadystate_system(model, mode) -> (SteadyState, Bool)

Provides an alternative to calc_eqm for the model equations that characterize the steady state. This is useful when the \
steady state can be characterized by a smaller system of equations.

Expected return: depends on the evaluation mode:
- in `SolveMode`: a vector `resid` of residuals for the steady state system of equations, given a guess for the steady state values.
- in `IterateMode`: a tuple of `(resid, stv)` where `resid` is the same as in `SolveMode` and `stv` is a `SteadyState` object 
containing the steady state values implied by the guess.
"""
function steadystate_system(::NamedTuple, mdl::AbstractModel, mode::EvaluationMode)
   throw(MethodNotImplementedError(:steadystate_system, typeof(mdl)))
end


"""
    calc_eqm(state, policy, mdl, np, mode)

Compute equilibrium conditions at `state` given `policy`. Interface stub --
model packages must add methods dispatching on their `AbstractModel` subtype.

Expected return: depends on the evaluation mode:
- in `SolveMode`: a NamedTuple with fields `fx` (residuals)
- in `IterateMode`: a NamedTuple with fields `trans` (next-period endogenous state values for 
  each  realization of exogenous state), and `vals` (this-period values of Vfct variables)
- in `SimulateMode`: a NamedTuple with fields `trans` (same as in IterateMode) and `simvars` 
  (value outputs and extra simulation variables in SimulateMode).
"""
function calc_eqm(state, policy, mdl::AbstractModel, np::NextPeriod, mode::EvaluationMode)
    throw(MethodNotImplementedError(:calc_eqm, typeof(mdl)))
end

"""
    compute_next_period_values(next, model, mode) -> NextPeriod

Compute next-period auxiliary quantities given next-period state `next`. Interface stub --
model packages must add methods dispatching on their `AbstractModel` subtype.

Return a `NextPeriod(; expectations, continuations)` where:
- `expectations`: probability-weighted scalars (averaged by `compute_expectations`)
- `continuations`: per-state raw values (stacked into an `SVector` across exogenous states)
"""
function compute_next_period_values(next, mdl::AbstractModel, mode::EvaluationMode)::NextPeriod
    throw(MethodNotImplementedError(:compute_next_period_values, typeof(mdl)))
end


"""
    forecasts_from_policy(state, policy, model) -> NamedTuple

Forecasting values: current-period quantities, other than the next-period state, that next
period's expectations depend on (e.g. a portfolio entering a next-period consumption stage).
`state` is the full current state (exogenous and endogenous variables), `policy` the named,
constrained policy. Optional: the default returns an empty NamedTuple (no forecasting
function).

If the model returns a non-empty NamedTuple, `main_setup` builds a forecasting function
`obj.Ffct` from it, and `forecast_next` attaches the forecasting values to every next-period
point, so `compute_next_period_values` sees them as fields of `next` (names must not clash
with state or Vfct names). Which values are attached:
- `STPFIMethod`, the fixed-point steady state, simulation and IRFs: the values implied by the
  policy being evaluated;
- `TFIMethod`: expectations are precomputed before the policy is solved, so they use the stored
  `obj.Ffct` from the previous iteration, which is then updated (and damped) like `Vfct`.
  The equilibrium conditions must then depend on the policy other than through the forecasts,
  or the system is singular in it; models whose choices act only through next period's
  expectations need `STPFIMethod`.
"""
forecasts_from_policy(state, policy, mdl::AbstractModel) = (;)


"""
    Exogenv(mdl::AbstractModel) -> Exogenv

Compute the exogenous environment for the given model. Interface stub --
model packages must add methods dispatching on their `AbstractModel` subtype.

Return a `Exogenv` instance representing the Markov process for the model's exogenous variables.
"""
function Exogenv(mdl::AbstractModel)
    throw(MethodNotImplementedError(:Exogenv, typeof(mdl)))
end