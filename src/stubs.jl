## Interface stubs: model packages must implement these for their AbstractModel subtype, or the default method will throw a MethodNotImplementedError.

"""
    define_endog_grids(model, stv) -> NamedTuple

Return one grid vector per endogenous state in a NamedTuple, for example
`(K = range(0.5, 10.0, 50),)`. Implement this for your model type; `main_setup` uses the vectors
to build its state grid.
"""
function define_endog_grids(mdl::AbstractModel, ::AbstractSteadyState)
    throw(MethodNotImplementedError(:define_endog_grids, typeof(mdl)))
end

"""
    steadystate(model, ::AnalyticalSteadyStateMode) -> (SteadyState, Bool)

Implement an analytical steady state for your model. Return `(stv, converged)`, where
`stv::SteadyState` contains its states, values and policy guess, and `converged::Bool` indicates
success. For a numerical steady state, use [`SystemSteadyStateMode`](@ref) or
[`FixedPointSteadyStateMode`](@ref).
"""
function steadystate(mdl::AbstractModel, ::AnalyticalSteadyStateMode)
    throw(MethodNotImplementedError(:steadystate, typeof(mdl), " for AnalyticalSteadyStateMode"))
end


"""
    steadystate_system(guess::NamedTuple, model, mode)

Evaluate a reduced system of steady-state equations for [`SystemSteadyStateMode`](@ref).
`guess` is a NamedTuple of the unknowns named in that mode's initial guess.

Return `(resid, nothing)` in `SolveMode` and `(resid, stv)` in `IterateMode`, where `stv` is the
complete [`SteadyState`](@ref) implied by the guess. Implement both modes for your model type.
"""
function steadystate_system(::NamedTuple, mdl::AbstractModel, mode::EvaluationMode)
   throw(MethodNotImplementedError(:steadystate_system, typeof(mdl)))
end


"""
    calc_eqm(state, policy, mdl, np, mode)

Evaluate the equilibrium conditions at `state` for a trial `policy`, using next-period quantities
in `np`. Implement this for your model type.

Return the quantities needed for the requested mode:
- in `SolveMode`: the residuals `fx` themselves, as an `SVector` or a `NamedTuple`
- in `IterateMode`: a NamedTuple `(; fx, trans, vals)` with the residuals, the next-period
  endogenous state `trans` (each entry a scalar, or an `SVector` with one element per exogenous
  state next period), and this period's values `vals` of the `Vfct` variables
- in `SimulateMode`: a NamedTuple `(; fx, trans, simvars)`, with `fx` and `trans` as above and
  `simvars` the variables to record in the simulation
"""
function calc_eqm(state, policy, mdl::AbstractModel, np::NextPeriod, mode::EvaluationMode)
    throw(MethodNotImplementedError(:calc_eqm, typeof(mdl)))
end

"""
    compute_next_period_values(next, model, mode) -> NextPeriod

Calculate the next-period quantities needed in today's equations for a single exogenous state.
Implement this for your model type.

It is called once for each exogenous state that can follow the current one. `next` holds that
state's exogenous variables, the next-period endogenous state, the `Vfct` values there and any
forecasting values. Return a `NextPeriod(; expectations, continuations)` of this state's values,
normally scalars:
- `expectations`: terms to take the expectation of; `compute_expectations` weights them by the
  transition probabilities and sums, so `calc_eqm` receives scalars in `np.expectations`
- `continuations`: other values needed by state; `compute_expectations` stacks them, so `calc_eqm`
  receives one entry per next-period exogenous state in `np.continuations`

Do not weight by the probabilities or return values for several exogenous states.
"""
function compute_next_period_values(next, mdl::AbstractModel, mode::EvaluationMode)::NextPeriod
    throw(MethodNotImplementedError(:compute_next_period_values, typeof(mdl)))
end


"""
    forecasts_from_policy(state, policy, model) -> NamedTuple

Forecasting values implied by `policy` at `state`: current-period quantities, beyond the
next-period state, that next period's expectations depend on (e.g. a portfolio entering a
next-period consumption stage). Overload when the model has a forecasting function `Ffct`; the
names must not clash with state or Vfct names. Default: no forecasting values.
"""
forecasts_from_policy(state, policy, mdl::AbstractModel) = (;)


"""
    Exogenv(mdl::AbstractModel) -> Exogenv

Return an `Exogenv` representing the Markov process for the model's exogenous variables.
Implement this for your model type; [`ExogenvAR1`](@ref) constructs a discretized AR(1) process.
"""
function Exogenv(mdl::AbstractModel)
    throw(MethodNotImplementedError(:Exogenv, typeof(mdl)))
end
