# Quick Start

This page solves the neoclassical growth model (NGM), checks the accuracy of the solution,
simulates it and computes moments and impulse responses. The same model is in `examples/ngm.jl`,
and the later pages of this guide continue with it. The NGM is simple enough that its global
solution is close to the log-linear one, which makes it easy to check.

Install `TransitionFunctionIteration` and `StaticArrays` in your Julia project, then run the blocks
in order. StaticArrays supplies the small fixed-size vectors used in the model equations.

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

These functions describe the economics. [Implementing a Model](@ref) explains each one.

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

`main_setup` computes the steady state, builds the grid, initializes policies and values from it,
and initializes transitions to leave endogenous states unchanged.
`main_execute!` iterates until the largest change between successive iterates of the value and
transition functions (and forecasting functions, when present) is below `10^tol`; the default
`tol = -5` means `1e-5`. It stops in one of
three states, returned as `status`: `:converged`; `:maxiter`, when `maxiter` iterations ran without
convergence; or `:all_failed`, when the solver failed at every grid point (see
[Choosing a solution method](@ref) for a common cause). The output above is the status, the number
of iterations and the final distances of the value and transition functions in `log10`.
`damp` is the weight on the previous iterate (0 is a full update, 1 freezes the functions).
`logging = COMPLETE` prints only the final message; the default prints every iteration.

## Step 4: Check accuracy

Convergence of the iteration does not mean the solution is accurate. Simulate the model and look at
the equation errors along the simulated path:

```@example qs
simseries, errs = simulate(obj, stv)
errstats = calc_error_stats(errs)
```

There is one row per equilibrium condition (`fx1`, the Euler equation) and one per transition
(`K`), and the columns are the mean and percentiles of the absolute error, with `p100` the
largest. Errors of the order of `1e-5` suggest the solution is accurate; if they are large, use a finer
grid where the simulation goes. `simulate` also reports the simulated range of each state against
its grid.

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
