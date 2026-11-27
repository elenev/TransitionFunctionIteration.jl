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
Import the interface functions you extend (`import TransitionFunctionIteration: ...`), as in the
[Quick Start](@ref). On the solve path, keep numeric types generic so ForwardDiff dual numbers
work; use static vectors and matrices for small collections, and avoid shared mutable state
because grid points and simulation runs are evaluated on multiple threads.

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

The first grid axis is always `:exst`, the integer exogenous state index. Do not include it in
`define_endog_grids`. Endogenous grids may be unevenly spaced; interpolation is linear with flat
extrapolation. An out-of-range exogenous index throws an error. Check simulated bounds and errors
before relying on results outside the endogenous grid.

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
    if mode isa SolveMode
        return fx
    elseif mode isa IterateMode
        return (; fx, trans=(; K=K′), vals=(; C))
    elseif mode isa SimulateMode
        return (; fx, trans=(; K=K′), simvars=(; Y, I, C))
    end
end
```

### `compute_next_period_values` return contract

`compute_next_period_values(next, mdl, mode)` describes one next-period exogenous state. The package
calls it once for each exogenous state that can follow the current one, so it never deals with
probabilities or with more than one state at a time.

`next` is that state: its exogenous variables (`A`), the next-period endogenous state from the
transition functions (`K`), the value functions there (`C`), and the forecasting values if the model
has any. The function returns `NextPeriod(; expectations, continuations)`, two `NamedTuple`s of
this state's values, normally scalars:

- `expectations`: the terms whose conditional expectation enters the equilibrium conditions. The
  package multiplies each by the probability of the state and sums over states, so `calc_eqm`
  receives `np.expectations.Ek`, a scalar equal to the expectation.
- `continuations`: other next-period values the equilibrium conditions need state by state. The
  package collects them across states, so `calc_eqm` receives `np.continuations.MU`, a vector with
  one entry per next-period exogenous state; `np.prnext` holds the probabilities of those states, and
  `np.continuations` also holds the variables of `next` by state.

So do not weight by probabilities and do not return a vector over exogenous states: both are done
for you.

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

Import `steadystate_system` from TransitionFunctionIteration. In `SolveMode`, return
`(residuals, nothing)` for the NamedTuple `guess`. In `IterateMode`, return `(residuals, stv)`, where `stv` is the
complete `SteadyState` implied by that guess. Both modes must accept dual numbers.

### `FixedPointSteadyStateMode`

Use when you cannot solve for the steady state analytically or via a reduced system. The solver
solves the dynamic equations jointly for policies and endogenous states, imposing that transitions
return the same state:

```julia
obj, stv = main_setup(mdl; mode=FixedPointSteadyStateMode(
    (C=0.8,),          # initial policy guess
    (K=3.0,),          # initial transition guess
    (A=1.0,),          # exogenous state at the steady state
))
```

By default `vals_from_policy` supplies the policy itself as the initial `Vfct` values. Override it
for other values that can be computed from state and policy. Values that depend on expectations
can instead be extra unknowns: pass `eqmvalguesses = (; V = 1.0)` to the mode, using names returned
in `calc_eqm`'s `vals`. Their equations equate the guesses with those returned values. To bound the
fixed-point unknowns, overload `bounds(mdl, policies, transitions)` or, with equilibrium-value
guesses, `bounds(mdl, policies, transitions, eqmvalguesses)`, returning bounds for the policy,
transition and equilibrium-value guesses respectively.

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
| `vals_from_policy(state, policy, mdl)` | `Vfct` values as functions of the policy, for the fixed-point steady state (default: the policy itself) |

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

### Forecasting functions

Next period's expectations sometimes depend on current-period quantities besides next period's
state. For example, the stochastic discount factor `β (C′/C)^(-γ)` needs today's consumption `C`
as well as next period's `C′`. Declare such quantities with `forecasts_from_policy`, which returns a
`NamedTuple` computed from the state and the policy; each name then appears in `next` in
`compute_next_period_values`:

```julia
import TransitionFunctionIteration: forecasts_from_policy

forecasts_from_policy(state, policy, mdl::MyModel) = (; Cnow = policy.C)   # names must differ from state and Vfct names

function compute_next_period_values(next, mdl::MyModel, ::EvaluationMode)
    (; A, K, C, Cnow) = next
    # ... β * (C / Cnow)^(-γ) * (MPK + 1 - δ) ...
end
```

The solution then has a forecasting function `Ffct`, which stores these values at every grid point.
After each iteration it is updated from the newly solved policy, with the same damping as the value
and transition functions, and the iteration converges only when it has stopped changing too
(`distF` in the result of `main_execute!`). Both solution methods support forecasting functions;
they differ in which forecasting values a grid-point solve uses, as described next.

### Choosing a solution method

Each iteration solves the equilibrium conditions at every grid point, then updates the value,
transition and forecasting functions. The two methods differ in what the solver varies at a grid
point:

- `TFIMethod` (the default) computes the expectations once, before solving, from the next-period
  states given by the current transition functions and from the stored forecasting values. The
  solver then varies only the policy. Its effect on next period's state and on the forecasting
  values enters the expectations at the next iteration, when the updated functions are used.
- `STPFIMethod` (Cao, Luo, and Nie) makes the next-period states unknowns of the solve alongside
  the policy, with the transitions computed by `calc_eqm` as extra equations, and recomputes the
  forecasting values from the solver's current policy. The expectations are recomputed at every
  evaluation, so the policy's effect on them is part of the solve. The system at each point is
  larger: the policy variables plus one next-period value of each endogenous state per exogenous
  state.

```julia
main_execute!(obj; solmethod = STPFIMethod())
```

At a converged solution both methods satisfy the same equilibrium conditions, and TFI is cheaper
at each point. TFI fails when a policy variable affects the residuals of `calc_eqm` only through
next period's state or the forecasting values. Within the solve those are fixed, so the residuals
do not depend on that variable, the system is singular, and the solver fails at every point: the
variable stays at its initial guess. In the example above, if `calc_eqm` uses the expectation as
it is, `1 - E[β (C′/Cnow)^(-γ) (MPK′ + 1 - δ)]`, the residual depends on today's `C` only through
`Cnow`; TFI cannot solve it and `STPFIMethod` can. Writing the residual so that `C` also enters
directly, `1 - E[...] (Cnow/C)^(-γ)`, gives the same equation at the solution, where `Cnow = C`, and
TFI solves it. When TFI fails in this way, `main_execute!` stops after the first
iteration with the message "Solver failed at every grid point." and `status = :all_failed`.

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
