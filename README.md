# TransitionFunctionIteration.jl

[![Build Status](https://github.com/elenev/TransitionFunctionIteration.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/elenev/TransitionFunctionIteration.jl/actions/workflows/CI.yml)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://elenev.github.io/TransitionFunctionIteration.jl/dev/)
[![CompatHelper](https://github.com/elenev/TransitionFunctionIteration.jl/actions/workflows/CompatHelper.yml/badge.svg)](https://github.com/elenev/TransitionFunctionIteration.jl/actions/workflows/CompatHelper.yml)
[![TagBot](https://github.com/elenev/TransitionFunctionIteration.jl/actions/workflows/TagBot.yml/badge.svg)](https://github.com/elenev/TransitionFunctionIteration.jl/actions/workflows/TagBot.yml)

A package for globally solving macroeconomic models with Transition Function Iteration (Elenev,
Landvoigt, and Van Nieuwerburgh, 2021) or, optionally, Simultaneous Transition and Policy Function
Iteration (STPFI, Cao, Luo, and Nie, 2023).

A solved model is a set of functions of the state variables, linearly interpolated on a (potentially
irregular) rectangular grid. The first state is always the index of the Markov chain that drives the
exogenous variables.

- **`Pfct`**: policies, the solution of the system of equilibrium conditions at each grid point.
- **`Vfct`**: equilibrium quantities and prices whose next-period values enter today's conditions.
- **`Tfct`**: next-period values of the endogenous states, for each realization of the exogenous state.
- **`Ffct`** (optional): forecasting values, current-period quantities that next-period expectations
  depend on.

The solver starts from the deterministic steady state and iterates to convergence. The package also
simulates, computes impulse responses, summarizes and compares results, and plots them. Results are
Tables.jl tables with provenance in their metadata.

## Installation

```julia-repl
julia> using Pkg
julia> Pkg.add("TransitionFunctionIteration")
julia> Pkg.add("StaticArrays")   # used by the example models
```

Julia 1.12 or later.

## Getting started

To run the bundled growth model, including a solution, simulation and error checks:

```julia
using TransitionFunctionIteration
include(joinpath(pkgdir(TransitionFunctionIteration), "examples", "ngm.jl"))
```

A model is a struct that inherits from `AbstractModel` and holds its parameters, plus five
functions with its economics. The quickest start is a template with a shell for each:

```julia-repl
julia> generate_template("MyModel", "/path/to/model_file.jl")
```

| Function | Provides |
|:--|:--|
| `Exogenv(mdl)` | the discretized Markov process of the exogenous states |
| `steadystate(mdl, ::AnalyticalSteadyStateMode)` | the deterministic steady state (or use `SystemSteadyStateMode` or `FixedPointSteadyStateMode`) |
| `define_endog_grids(mdl, stv)` | grid points for the endogenous states |
| `compute_next_period_values(next, mdl, mode)` | expectations and continuations at each next-period state |
| `calc_eqm(state, policy, mdl, np, mode)` | the equilibrium conditions, and (by `mode`) the new `Vfct` and `Tfct` values or simulation variables |

With `mdl = MyModel()`:

```julia
obj, stv = main_setup(mdl)                    # steady state, grids, initial functions
main_execute!(obj; maxiter = 1000)            # iterate (tol = -5 means a distance of 1e-5); .status says how it ended

simseries, errs = simulate(obj, stv)          # stochastic simulation
stats    = describe(simseries).moments.uncond # MomentArray of moments
errstats = calc_error_stats(errs)             # equation errors on the simulated states: judge accuracy by these
irf      = simulate_irf(obj, simseries)       # impulse responses
```

**Choosing the method.** `TFIMethod()` (the default) computes expectations at each grid point once,
from the stored transition and forecasting functions, and then solves for the policy.
`STPFIMethod()` solves for the policy and the next-period states together, recomputing expectations
at every step. Both support forecasting functions and have the same fixed point. TFI cannot solve for a
policy variable that enters the equilibrium conditions solely through next period's state or forecasting
values; `main_execute!` then stops with `status = :all_failed`, and such a model needs
`main_execute!(obj; solmethod = STPFIMethod())`.

The [documentation](https://elenev.github.io/TransitionFunctionIteration.jl/stable/) has a quick start
that builds and solves the neoclassical growth model step by step, the full model interface
(optional hooks, bounds, steady-state modes), simulation and accuracy, moments and comparing
economies, and tables and provenance.

## Examples

The `examples` folder holds two models, and is the reference for how a model is written.

1. `ngm.jl`: the neoclassical growth model, the minimal model. It is easy to solve, and its log-linear
   approximation is close to the global solution; the test suite uses the log-linear solution as a check.
   Start here.
2. `intermedap.jl`: the advanced example, an intermediary asset pricing model with two agents,
   occasionally binding constraints, trend growth and Epstein-Zin preferences (`intermedap.pdf` has the
   model). It shows a steady state from a smaller system, bounds, normalization of expectations and
   errors, augmented simulation variables, and `solve_and_simulate` with options for each stage,
   conditional moments and IRF paths. It runs as shipped in about 20 seconds with one thread
   (`julia --project=. examples/intermedap.jl`).

## See also

- [OutputVault.jl](https://github.com/elenev/OutputVault.jl) saves solutions, simulations, moments and IRFs of research runs with the git state that produced them, and queries them across runs; TransitionFunctionIteration's support for it loads with `using OutputVault`.

The [saving and comparing guide](https://elenev.github.io/TransitionFunctionIteration.jl/stable/guide/outputvault/)
walks through saving a run, reloading its solution and comparing saved economies. It requires
OutputVault 0.4 and Parquet2; the integration lives in TFI's extension.

## Development

Tests and documentation examples exercise the OutputVault integration. Check out OutputVault
beside this repository at `../OutputVault`, then run:

```sh
JULIA_NUM_THREADS=2 julia --project -e 'using Pkg; Pkg.test()'
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```
