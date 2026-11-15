# Quick Start

This page walks through the complete TFI workflow using the Neoclassical Growth Model (NGM),
the canonical reference example in `examples/ngm.jl`. The NGM is simple enough that its
global solution closely matches the log-linear approximation, which makes it easy to verify.

## Step 0: Define the model struct

A model is a plain Julia struct that inherits from `AbstractModel` and holds all parameters:

```julia
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
function define_endog_grids(model::NeoclassicalGrowthModel, stv::AbstractSteadyState)
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
interpolants. `main_execute!` iterates until the largest distance between successive iterates of the
value and transition functions (and the forecasting function, if the model has one) is below `tol`,
which is a `log10` distance: the default `tol = -5` means `1e-5`. `damp` is the weight on the
previous iterate (0 is a full update, 1 freezes the functions), and `logging = COMPLETE` suppresses the
iteration logs.

## Step 3: Simulate and summarize

```julia
simseries, errs = simulate(obj, stv)
stats    = describe(simseries).moments.uncond
errstats = calc_error_stats(errs)
```

`simseries` is a `SimulationPanel`; `stats` is a `MomentArray` of unconditional moments with
one row per variable and one column per statistic (`:mean`, `:std`, `:minimum`, `:maximum`, `:ac`
for the autocorrelation, the percentiles `:p5` to `:p95`, and correlations with the endogenous
states, e.g. `:cor_K_lag0`). `errstats` has one row per equilibrium equation and transition.

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
# result.stats     -- describe(simseries).moments: a NamedTuple of MomentArrays by sample (uncond, ...)
# result.errstats  -- error statistics (MomentArray)
# result.irf       -- IRFResult
```

## Step 5: Impulse responses

```julia
irf = simulate_irf(obj, simseries)
propertynames(irf.means)   # (:uncond, :shock_1, :shock_2, :shock_4, :shock_5)
irf.means.shock_5.C        # mean path of C when the exogenous state jumps to state 5 in period 2
irf.means.uncond.C         # mean path of C without a forced shock

# name the shocks yourself: a forced path of exogenous state indices
irf = simulate_irf(obj, simseries; shocks = (up = 5, down = 1, uncond = Int[]))
irf.means.up.K
```

By default there is one shock per exogenous state other than the starting one (`shock_<index>`),
plus the baseline `uncond`, which draws every period as usual. A response is a shock path
compared with `uncond`.

## Step 6: Save the results

The result of `solve_and_simulate` has the shape that OutputVault.jl saves, with the model's parameters
added (see its documentation for the configuration):

```julia
using OutputVault, Parquet2
set_config!(Config(output_root = "/path/to/outputs", table_format = OutputVault.Formats.ParquetFormat()))
save_result(RunHandle("baseline"; experdef = "demo"), merge(result, (; params = collect_parameters(result.obj))))
```

See [Simulation and Analysis](@ref) for more on plotting, conditional moments, and
cross-economy comparison.
