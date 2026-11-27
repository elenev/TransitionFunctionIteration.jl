# TransitionFunctionIteration.jl

[![Build Status](https://github.com/elenev/TransitionFunctionIteration.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/elenev/TransitionFunctionIteration.jl/actions/workflows/CI.yml)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://www.vadimelenev.com/TransitionFunctionIteration.jl/stable/)
[![CompatHelper](https://github.com/elenev/TransitionFunctionIteration.jl/actions/workflows/CompatHelper.yml/badge.svg)](https://github.com/elenev/TransitionFunctionIteration.jl/actions/workflows/CompatHelper.yml)
[![TagBot](https://github.com/elenev/TransitionFunctionIteration.jl/actions/workflows/TagBot.yml/badge.svg)](https://github.com/elenev/TransitionFunctionIteration.jl/actions/workflows/TagBot.yml)

TransitionFunctionIteration.jl solves discrete-time dynamic equilibrium models over a grid of
states. You provide the model's equations, parameters and grids; the package solves for policies
and equilibrium quantities, simulates the economy, and computes impulse responses and moments.

It supports Transition Function Iteration (TFI, Elenev, Landvoigt, and Van Nieuwerburgh, 2021) and
Simultaneous Transition and Policy Function Iteration (STPFI, Cao, Luo, and Nie, 2023). Results can
be plotted, compared across economies, or passed directly to Tables.jl-compatible file writers.

## Installation

```julia-repl
julia> using Pkg
julia> Pkg.add("TransitionFunctionIteration")
julia> Pkg.add(["StaticArrays", "StructArrays"])   # used by the example models and the template
```

Julia 1.12 or later.

## Getting started

To run the bundled growth model, including a solution, simulation and error checks:

```julia
using TransitionFunctionIteration
include(joinpath(pkgdir(TransitionFunctionIteration), "examples", "ngm.jl"))
```

To implement your own model, define a subtype of `AbstractModel` holding its parameters and
provide the functions below. Generate a file with their signatures to fill in:

```julia-repl
julia> generate_template("MyModel", "/path/to/model_file.jl")
```

| Function | Provides |
|:--|:--|
| `Exogenv(mdl)` | the discretized Markov process of the exogenous states |
| steady state: one of three options | the deterministic steady state, in closed form, from a reduced system, or as a fixed point of the dynamic equations (see [steady-state computation](https://www.vadimelenev.com/TransitionFunctionIteration.jl/stable/guide/model/#Steady-state-computation)) |
| `define_endog_grids(mdl, stv)` | grid points for the endogenous states |
| `compute_next_period_values(next, mdl, mode)` | expectations and continuations at each next-period state |
| `calc_eqm(state, policy, mdl, np, mode)` | the equilibrium conditions, and (by `mode`) the new `Vfct` and `Tfct` values or simulation variables |

Once the model is implemented, create `mdl = MyModel()` and solve and analyze it:

```julia
obj, stv = main_setup(mdl)                    # steady state, grids, initial functions
main_execute!(obj; maxiter = 1000)            # iterate (tol = -5 means a distance of 1e-5); .status says how it ended

simseries, errs = simulate(obj, stv)          # stochastic simulation
stats    = describe(simseries).moments.uncond # MomentArray of moments
errstats = calc_error_stats(errs)             # equation errors on the simulated states: judge accuracy by these
irf      = simulate_irf(obj, simseries)       # impulse responses
```

The [quick start](https://www.vadimelenev.com/TransitionFunctionIteration.jl/stable/guide/quickstart/)
walks through the growth model from its equations to simulation and analysis. Use the simulated
equation and transition errors to assess accuracy, and check that simulated states stay inside
the grid.

## The solution

A `ModelSolution` holds functions of the state, linearly interpolated on a rectangular grid whose
nodes may be unevenly spaced:

- `Pfct` (policy functions): the unknowns of the equilibrium conditions at each grid point.
- `Vfct` (value functions): quantities whose next-period values enter today's equilibrium
  conditions.
- `Tfct` (transition functions): the next-period value of each endogenous state, one value per
  realization of the next-period exogenous state.
- `Ffct` (forecasting functions, optional): current-period quantities that next period's
  expectations depend on besides the next-period state.

The first state axis is the index of the exogenous Markov chain. You choose the endogenous axes
and their grid points.

The default `TFIMethod()` solves for policies while holding expectations fixed within each
grid-point solve. Use `STPFIMethod()` when a policy affects the equations only through next-period
states or forecasting values. See
[choosing a solution method](https://www.vadimelenev.com/TransitionFunctionIteration.jl/stable/guide/model/#Choosing-a-solution-method)
for the distinction and an example.

## Examples

The `examples` folder contains two complete models:

1. `ngm.jl`: the neoclassical growth model, the minimal model. It is easy to solve, and its log-linear
   approximation is close to the global solution. Start here to learn the model interface.
2. `intermedap.jl`: the advanced example, an intermediary asset pricing model with two agents,
   occasionally binding constraints, trend growth and Epstein-Zin preferences (`intermedap.pdf` has the
   model). It shows a steady state from a smaller system, bounds, normalization of expectations and
   errors, augmented simulation variables, and `solve_and_simulate` with options for each stage,
   conditional moments and IRF paths. From a checkout of the repository, run it with
   `julia --project=examples examples/intermedap.jl`.

## Saving results

Use [OutputVault.jl](https://www.vadimelenev.com/OutputVault.jl/stable/) to keep a solved model together with
its simulations, moments and IRFs, record the code that produced them, and compare saved runs.
Install OutputVault and Parquet2 alongside TFI, then follow the
[saving and comparing guide](https://www.vadimelenev.com/TransitionFunctionIteration.jl/stable/guide/outputvault/).

## Development

From a checkout, run the tests and build the documentation with the commands below. Both need
OutputVault checked out at `../OutputVault`.

```sh
JULIA_NUM_THREADS=2 julia --project -e 'using Pkg; Pkg.test()'
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```
