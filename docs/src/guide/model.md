# Implementing a Model

## Overview

To use this package, you need two things:

1. A struct `MyModel <: AbstractModel` holding your model's parameters.
2. Implementations of several interface functions dispatching on `MyModel`.

The easiest way to get started is:

```julia
generate_template("MyModel", "/path/to/model_file.jl")
```

This writes a template file with shells for all required and optional functions.

## Required interface

| Function | Signature | Purpose |
|---|---|---|
| `Exogenv` | `Exogenv(mdl::MyModel)` | Build the exogenous Markov process |
| `define_endog_grids` | `define_endog_grids(mdl, stv)` | Return grid vectors for endogenous states |
| `steadystate` | `steadystate(mdl, mode)` | Compute the deterministic steady state |
| `compute_next_period_values` | `compute_next_period_values(next, mdl, mode)` | Compute next-period expectations and continuations |
| `calc_eqm` | `calc_eqm(state, policy, mdl, np, mode)` | Evaluate equilibrium conditions |

### Exogenous environment

`Exogenv(mdl)` should return an `Exogenv` representing the model's Markov process.
For a single AR(1) process use [`ExogenvAR1`](@ref):

```julia
function Exogenv(mdl::MyModel)
    ExogenvAR1(ρ=mdl.ρ, σ=mdl.σ, n=mdl.N, transform=exp, jensen=true, shock_name=:A)
end
```

### Endogenous grids

`define_endog_grids(mdl, stv)` returns a `NamedTuple` of grid vectors -- one per endogenous
state variable. `stv` is the `SteadyState` computed earlier, so grids can be scaled
relative to steady-state values:

```julia
function define_endog_grids(mdl::MyModel, stv::AbstractSteadyState)
    K_grid = range(0.8, 1.2, length=40) * stv.endog.K
    return (; K = K_grid)
end
```

### `calc_eqm` return contract

The return value depends on `mode`:

| Mode | Return |
|---|---|
| `SolveMode` | `SVector` (or `NamedTuple`) of residuals `fx` |
| `IterateMode` | `(; fx, trans, vals)` |
| `SimulateMode` | `(; fx, trans, simvars)` |

- `trans`: `NamedTuple` of next-period endogenous state values. Each value can be a scalar
  (same across exogenous states) or an `SVector` (one entry per exogenous state). Scalars
  are automatically broadcast by `expand_transitions`.
- `vals`: `NamedTuple` of variables for `Vfct` -- the quantities whose next-period values
  enter expectations in equilibrium conditions.
- `simvars`: `NamedTuple` of all output variables to record in simulation.

Example dispatch pattern:

```julia
function calc_eqm(state, policy, mdl::MyModel, np::NextPeriod, mode::EvaluationMode)
    # ... compute Y, I, K′, MU, fx ...
    mode isa SolveMode    && return fx
    mode isa IterateMode  && return (; fx, trans=(; K=K′), vals=(; C))
    mode isa SimulateMode && return (; fx, trans=(; K=K′), simvars=(; Y, I, C))
end
```

### `compute_next_period_values` return contract

Returns a `NextPeriod(; expectations, continuations)` where:

- `expectations`: a `NamedTuple` of probability-weighted scalars (one per quantity), computed
  by `compute_expectations` as `E_t[f(x_{t+1})]`. These enter `calc_eqm` via `np.expectations`.
- `continuations`: a `NamedTuple` of per-state raw values (one entry for each exogenous state
  realization). These enter `calc_eqm` via `np.continuations` and are available in `SimulateMode`
  for the realized next-period state.

```julia
function compute_next_period_values(next, mdl::MyModel, ::EvaluationMode)
    (; A, K, C) = next
    (; α, β, δ, γ) = mdl
    MPK = α * A * K^(α-1)
    MU  = C^(-γ)
    return NextPeriod(;
        expectations  = (; Ek = β * MU * (MPK + 1 - δ)),
        continuations = (; MPK, MU),
    )
end
```

## Steady-state computation

Three modes are available, chosen by passing a `SteadyStateMode` to `main_setup`:

### `AnalyticalSteadyStateMode` (default)

Implement `steadystate(mdl, ::AnalyticalSteadyStateMode)` returning `(SteadyState(...), true)`:

```julia
function steadystate(mdl::MyModel, ::AnalyticalSteadyStateMode)
    # ... closed-form solution ...
    return SteadyState(; model=mdl, exog=(; A=1.0), endog=(; K), values=(; C), guess=(; C)), true
end
```

Fields of `SteadyState`:
- `exog`: exogenous state at the steady state
- `endog`: endogenous state at the steady state (used to center grids)
- `values`: quantities stored in `Vfct` at initialization
- `guess`: initial solver guess for policy variables

### `SystemSteadyStateMode`

Use when the steady state satisfies a reduced system (smaller than the full dynamic system).
Implement `steadystate_system(guess, mdl, mode)` and pass the mode to `main_setup`:

```julia
obj, stv = main_setup(mdl; mode=SystemSteadyStateMode((; K=3.0, C=0.8)))
```

### `FixedPointSteadyStateMode`

Use when you cannot solve for the steady state analytically or via a reduced system. The solver
fixes the law of motion to the initial guess and finds a fixed point of the dynamic equations:

```julia
obj, stv = main_setup(mdl; mode=FixedPointSteadyStateMode(
    (C=0.8,),          # initial policy guess
    (K=3.0,),          # initial transition guess
    (A=1.0,),          # exogenous state at the steady state
))
```

## Optional interface

These functions have sensible defaults but can be overloaded:

| Function | Purpose |
|---|---|
| `bounds(stv::SteadyState{MyModel})` | Variable bounds for solver transformation (default: all `Unbounded()`) |
| `augment_simulation(run, mdl)` | Append extra variables to a simulation run |
| `forecasts_from_policy(state, policy, mdl)` | Current-period quantities that next period's expectations depend on (a forecasting function `Ffct`) |
| `augment_errors(obj, state, out, trans)` | Add model-specific convergence diagnostics to the error output |
| `normalize_next_period_values(np, mdl, mode)` | Transform expectations before they enter `calc_eqm` |
| `normalize_fx_errors(mdl, fx_nt, simvars)` | Normalize residuals for comparability |
| `initial_guess_from_stv(stv, mdl)` | Adjust the steady-state values before initializing interpolants |
| `vals_from_policy(state, policy, mdl)` | Needed when `Vfct` variables are not a subset of policy variables |

### Variable bounds

Bounds are applied via a smooth transformation so the solver always operates in an unconstrained
space. Specify bounds by overloading `bounds(stv::SteadyState{MyModel})`:

```julia
function bounds(::SteadyState{MyModel})
    policy_bounds = (C = PositiveBounded(), q = IntervalBounded(0.0, 1.0))
    endog_bounds  = (K = PositiveBounded(),)
    return policy_bounds, endog_bounds
end
```

Available bound types: `Unbounded`, `LowerBounded`, `UpperBounded`, `IntervalBounded`
(documented in the API reference; none is exported, so `import TransitionFunctionIteration:
PositiveBounded, IntervalBounded` or qualify them). Convenience constructors:
`PositiveBounded() = LowerBounded(0.0)`, `NegativeBounded() = UpperBounded(0.0)`.

### Augmenting simulation output

`augment_simulation(run, mdl)` receives a vector of the T period-NamedTuples of one simulation run
and returns a `NamedTuple` of extra variables, one vector each. A vector of length T - d is a
variable that needs d periods of history (d = 0: a variable of each period; d = 1: a variable of
consecutive periods, such as a gross return). It is aligned with the last T - d periods and the
first d are `NaN`. `simulate` simulates `context_periods` periods before the output periods,
augments them with the rest and drops them, so no `NaN` remains if `context_periods` is at least
d. Declare the largest d in `dynamic_depth(mdl)`, which is the default of `context_periods` (and
of `simulate_irf`'s `hold_periods`). Return `nothing` (the default) for no augmentation:

```julia
function augment_simulation(run, mdl::MyModel)
    R = [run[t].Y / run[t-1].K for t in 2:length(run)]    # T-1 values
    return (; R)
end

dynamic_depth(::MyModel) = 1
```

### Forecasting functions and the choice of method

Decide by one question: *does a variable chosen at this grid point change what next period's expectations
are?* In an ordinary model the answer is no. Expectations at a grid point depend on the next-period
state, which the transition functions already give, so `TFIMethod` (the default) computes them once
before solving for the policy at that point. In some models the answer is yes: for example a portfolio
chosen today enters next period's consumption stage, so the expectations depend on the very choice being
solved for.

Such a model declares those choices as *forecasting values*:

```julia
forecasts_from_policy(state, policy, mdl::MyModel) = (; portfolio = policy.m)   # names must not clash with state or Vfct names
obj, stv = main_setup(mdl)
main_execute!(obj; solmethod = STPFIMethod())
```

- `forecasts_from_policy(state, policy, mdl)` returns a `NamedTuple` of current-period quantities beyond the
  next-period state. The solution then has a forecasting function `Ffct` over the state space, iterated
  together with `Pfct`, `Vfct` and `Tfct`; the forecasting values appear in `next` in
  `compute_next_period_values`.
- `STPFIMethod` (Simultaneous Transition and Policy Function Iteration) recomputes expectations from the
  solver's current guess instead of precomputing them. `TFIMethod` freezes them within a grid-point solve, so it
  misses any effect of the policy on next period's expectations: where the policy acts on the equilibrium
  only through those expectations, the choice has no effect in the solve.
- `Ffct` is saved with the solution and is part of its identity (`solution_id`). Switching a model
  between ordinary TFI and forecasting values changes the solution architecture; treat the two as
  different solutions.
- Convergence is reported separately for the value, transition and forecasting functions.

### Extra error diagnostics

`augment_errors(obj, state, out, trans)` is called for each simulated state alongside the
standard equation and transition errors. It returns a `(scalar_nt, vector_nt)` pair: scalars
are merged into the equation errors, and per-next-exogenous-state vectors are merged into the
transition errors. The default returns `(;), (;)` and costs nothing. A returned name that
matches an existing equation error replaces it.

```julia
function augment_errors(obj, state, out, trans)
    C_gap = out.simvars.C - obj.Vfct(state).C   # recomputed vs. interpolated
    return (; C_gap), (;)
end
```

## Re-solving after a parameter change

`update_model` creates a new model with updated parameters, and `main_setup` accepts a
previous solution as a warm start:

```julia
mdl2   = update_model(mdl; σ=0.03)
obj2, stv2 = main_setup(mdl2, obj)   # warm-start from previous solution
main_execute!(obj2)
```
