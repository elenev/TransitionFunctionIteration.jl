# Quick Start

Solve a neoclassical growth model, check its accuracy, and use the solution to study consumption,
capital and output. The model has capital `K` and productivity `A` as states. The household chooses
consumption `C`; output is `A * K^α`, and unconsumed output adds to next period's capital after
depreciation. Productivity follows a discretized AR(1) process, and the consumption policy
satisfies the Euler equation.

Install `TransitionFunctionIteration` and `StaticArrays` in your Julia project, then run the blocks
below in order, or run the complete model in `examples/ngm.jl`. StaticArrays supplies the small
fixed-size vectors used in the model equations. The later guides use the same model.

## Step 1: Define the model

A model is a struct that subtypes `AbstractModel` and holds all parameters:

```@example qs
using TransitionFunctionIteration, StaticArrays

import TransitionFunctionIteration: Exogenv, calc_eqm, define_endog_grids,
                                     steadystate, compute_next_period_values, bounds,
                                     PositiveBounded

@kwdef struct NeoclassicalGrowthModel <: AbstractModel
    α::Float64 = 0.33   # capital share
    β::Float64 = 0.96   # discount factor
    δ::Float64 = 0.10   # depreciation
    γ::Float64 = 2.0    # CRRA risk aversion
    σ::Float64 = 0.02   # shock volatility
    ρ::Float64 = 0.9    # shock persistence
    N_A::Int   = 5      # exogenous grid points
end
nothing # hide
```

## Step 2: Implement the model interface

Define productivity's process, the steady state, the capital grid and the Euler equation. The
equation uses expected discounted marginal utility from the next period; TFI computes that
expectation from the values returned by `compute_next_period_values`.
[Implementing a Model](@ref) describes the functions in more detail.

```@example qs
# Exogenous process: discretized AR(1) log productivity
function Exogenv(mdl::NeoclassicalGrowthModel)
    ExogenvAR1(ρ=mdl.ρ, σ=mdl.σ, n=mdl.N_A, transform=exp, jensen=true, shock_name=:A)
end

# Deterministic steady state, in closed form
function steadystate(model::NeoclassicalGrowthModel, ::AnalyticalSteadyStateMode)
    (; α, β, δ) = model
    K = (α * β / (1 - β * (1 - δ)))^(1/(1-α))
    C = K^α - δ * K
    return SteadyState(; model, exog=(; A=1.0), endog=(; K), values=(; C), guess=(; C)), true
end

# Grid for the endogenous state, capital
function define_endog_grids(model::NeoclassicalGrowthModel, stv::AbstractSteadyState)
    K_grid = range(0.8, 1.2, length=20) * stv.endog.K
    return (; K = K_grid)
end

# Consumption and capital must be positive
function bounds(::SteadyState{NeoclassicalGrowthModel})
    return (C = PositiveBounded(),), (K = PositiveBounded(),)
end

# Quantities in one next-period exogenous state
function compute_next_period_values(next, mdl::NeoclassicalGrowthModel, ::EvaluationMode)
    (; A, K, C) = next
    (; α, β, δ, γ) = mdl
    MPK     = α * A * K^(α-1)
    MU      = C^(-γ)
    Kpayoff = MPK + (1 - δ)
    return NextPeriod(; expectations=(; Ek = β * MU * Kpayoff), continuations=(; Kpayoff, MPK, MU))
end

# Equilibrium conditions at one point of the state space
function calc_eqm(state, policy, mdl::NeoclassicalGrowthModel, np::NextPeriod, mode::EvaluationMode)
    (; A, K) = state;  (; C) = policy;  (; α, δ, γ) = mdl
    Y  = A * K^α
    I  = Y - C
    K′ = (1 - δ) * K + I
    MU = C^(-γ)
    fx = SVector(1 - np.expectations.Ek / MU)     # Euler equation residual
    if mode isa SolveMode
        return fx
    elseif mode isa IterateMode
        return (; fx, trans=(; K=K′), vals=(; C))
    elseif mode isa SimulateMode
        return (; fx, trans=(; K=K′), simvars=(; Y, I, C))
    end
end
nothing # hide
```

## Step 3: Solve

```@example qs
ngm = NeoclassicalGrowthModel()
obj, stv = main_setup(ngm)
res = main_execute!(obj; maxiter=1000, logging=COMPLETE)
(res.status, res.iter, res.distV, res.distT)
```

`main_setup` returns the initial solution `obj` and steady state `stv`. It starts policies and
values from the steady state and transitions from unchanged endogenous states. `main_execute!`
updates `obj` until the functions change by less than the requested tolerance, or it stops for
one of the reasons below:

| `res.status` | Meaning |
|:--|:--|
| `:converged` | Changes in the value, transition and any forecasting functions are below `10^tol`. |
| `:maxiter` | The iteration limit was reached before convergence. |
| `:all_failed` | The solver failed at every grid point; see [Choosing a solution method](@ref). |

The default `tol = -5` means `1e-5`. The displayed distances `distV` and `distT` are in `log10`,
so a value below `-5` meets that tolerance. `damp` is the weight on the previous iterate, 0.5 by
default: 0 gives a full update, and 1 freezes the functions. `logging = COMPLETE` prints only the final
message; the default prints every iteration.

## Step 4: Check accuracy

Simulate the model to check how well its interpolated solution satisfies the equations at the
states the economy visits:

```@example qs
simseries, errs = simulate(obj, stv)
errstats = calc_error_stats(errs)
```

There is one row per equilibrium condition (`fx1`, the Euler equation) and one per transition
(`K`), and the columns are the mean and percentiles of the absolute error, with `p100` the
largest. Here the errors in the table above are of the order of `1e-5`, which suggests an
accurate solution. Large errors
can call for a finer grid where the simulation goes. `simulate` also reports each state's range
against its grid; widen a grid if the simulation leaves it.

## Step 5: Moments

```@example qs
stats = describe(simseries).moments.uncond
stats((:C, :K, :Y), (:mean, :std, :ac))
```

`stats` holds one row per variable and one column per statistic: `:mean`, `:std`, `:minimum`,
`:maximum`, `:ac` (first-order autocorrelation), the percentiles `:p5` to `:p95`, and the
correlation with each endogenous state (`:cor_K_lag0`). Read single values by name:

```@example qs
stats(:C, :mean), stats(:K, :std)
```

## Step 6: Impulse responses

```@example qs
irf = simulate_irf(obj, simseries)
propertynames(irf.means)
```

The default experiment starts at the mean of the simulated states and has one shock for each
exogenous state other than the starting one (the steady-state index, 3 here). Under `shock_5`, the
exogenous state moves to 5, the highest productivity, in period 2. `uncond` draws the exogenous
state every period as usual and is the no-shock baseline. The response to a shock is the shock path
compared with `uncond`:

```@example qs
irf.means.shock_5.C .- irf.means.uncond.C
```

Name the shocks yourself by giving each a forced path of exogenous state indices:

```@example qs
irf_named = simulate_irf(obj, simseries; shocks = (up = 5, down = 1))
propertynames(irf_named.means)
```

## Step 7: All steps at once

`solve_and_simulate` runs steps 3 to 6 and returns their results:

```@example qs
result = solve_and_simulate(ngm; execute_kwargs = (; logging = COMPLETE))
keys(result)
```

`result.stats` is the `NamedTuple` of moments by sample that `describe` returns, so
`result.stats.uncond` is the `stats` of step 5.

## Next steps

- [Simulation and Analysis](@ref): simulation options, conditional moments, impulse responses and
  plots.
- [Moments and Comparing Economies](@ref): tables of moments across economies.
- [Tables and Provenance](@ref): passing results to any Tables.jl-compatible writer, and saving them.
