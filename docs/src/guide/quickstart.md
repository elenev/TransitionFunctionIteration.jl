# Quick Start

This page walks through the complete TFI workflow using the Neoclassical Growth Model (NGM),
the canonical reference example in `examples/ngm.jl`. The NGM is simple enough that its
global solution closely matches the log-linear approximation, which makes it easy to verify.

## Step 0: Define the model struct

A model is a plain Julia struct that inherits from `AbstractModel` and holds all parameters:

```julia
using TransitionFunctionIteration, StaticArrays

import TransitionFunctionIteration: Exogenv, calc_eqm, define_endog_grids,
                                     steadystate, compute_next_period_values, bounds

@kwdef struct NeoclassicalGrowthModel <: AbstractModel
    α::Float64 = 0.33   # capital share
    β::Float64 = 0.96   # discount factor
    δ::Float64 = 0.10   # depreciation
    γ::Float64 = 2.0    # CRRA risk aversion
    σ::Float64 = 0.02   # shock volatility
    ρ::Float64 = 0.9    # shock persistence
    N_A::Int   = 5      # exogenous grid points
end
```

## Step 1: Implement the model interface

Four functions are required. See [Implementing a Model](@ref) for the full interface.

```julia
# Exogenous process: discretized AR(1) log-productivity
function Exogenv(mdl::NeoclassicalGrowthModel)
    ExogenvAR1(ρ=mdl.ρ, σ=mdl.σ, n=mdl.N_A, transform=exp, jensen=true, shock_name=:A)
end

# Analytical steady state
function steadystate(model::NeoclassicalGrowthModel, ::AnalyticalSteadyStateMode)
    (; α, β, δ) = model
    K = (α * β / (1 - β * (1 - δ)))^(1/(1-α))
    C = K^α - δ * K
    return SteadyState(; model, exog=(; A=1.0), endog=(; K), values=(; C), guess=(; C)), true
end

# Endogenous grid: capital
function define_endog_grids(model::NeoclassicalGrowthModel, stv)
    K_grid = range(0.8, 1.2, length=20) * stv.endog.K
    return (; K = K_grid)
end

# Bounds: consumption and capital must be positive
function bounds(::SteadyState{NeoclassicalGrowthModel})
    return (C = PositiveBounded(),), (K = PositiveBounded(),)
end

# Next-period values: marginal utility and capital return for Euler equation
function compute_next_period_values(next, mdl::NeoclassicalGrowthModel, ::EvaluationMode)
    (; A, K, C) = next
    (; α, β, δ, γ) = mdl
    MPK    = α * A * K^(α-1)
    MU     = C^(-γ)
    Kpayoff = MPK + (1 - δ)
    return NextPeriod(; expectations=(; Ek = β * MU * Kpayoff), continuations=(; Kpayoff, MPK, MU))
end

# Equilibrium conditions
function calc_eqm(state, policy, mdl::NeoclassicalGrowthModel, np::NextPeriod, mode::EvaluationMode)
    (; A, K) = state;  (; C) = policy;  (; α, δ, γ) = mdl
    Y  = A * K^α
    I  = Y - C
    K′ = (1 - δ) * K + I
    MU = C^(-γ)
    fx = SVector(1 - np.expectations.Ek / MU)
    mode isa SolveMode    && return fx
    mode isa IterateMode  && return (; fx, trans=(; K=K′), vals=(; C))
    mode isa SimulateMode && return (; fx, trans=(; K=K′), simvars=(; Y, I, C))
end
```

## Step 2: Solve

```julia
ngm = NeoclassicalGrowthModel()
obj, stv = main_setup(ngm)
main_execute!(obj; maxiter=1000)
```

`main_setup` computes the steady state, builds the state space grid, and initializes the
interpolants. `main_execute!` iterates to convergence; default tolerance is `1e-5` on the
sup-norm distance between successive policy function iterates.

## Step 3: Simulate and summarize

```julia
simseries, errs = simulate(obj, stv)
stats    = describe(simseries).moments.uncond
errstats = calc_error_stats(errs)
```

`simseries` is a `SimulationPanel`; `stats` is a `MomentArray` of unconditional moments
(mean, standard deviation, autocorrelation, etc.).

```julia
# Access specific moments
stats(:C, :mean)          # mean consumption
stats(:K, :std)           # std of capital
```

## Step 4: Convenience wrapper

For a clean one-liner when you do not need to inspect intermediate objects:

```julia
result = solve_and_simulate(ngm)
# result.obj       -- ModelSolution
# result.stv       -- SteadyState
# result.simseries -- SimulationPanel
# result.stats     -- MomentArray (unconditional)
# result.errstats  -- error statistics
# result.irf       -- IRFResult
```

## Step 5: Impulse responses

```julia
irf = simulate_irf(obj, simseries)
# irf.means.A.C  -- mean C response to an A shock
# irf.means.A.K  -- mean K response to an A shock
```

See [Simulation and Analysis](@ref) for more on plotting, conditional moments, and
cross-economy comparison.
