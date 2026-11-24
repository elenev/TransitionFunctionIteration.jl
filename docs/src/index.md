# TransitionFunctionIteration.jl

`TransitionFunctionIteration.jl` solves discrete-time dynamic equilibrium models globally, using
the Transition Function Iteration (TFI) algorithm of Elenev, Landvoigt, and Van Nieuwerburgh (2021)
or the Simultaneous Transition and Policy Function Iteration (STPFI) algorithm of Cao, Luo, and Nie
(2023). [Choosing a solution method](@ref) explains when each applies; TFI is the default.

## What you provide and what you get

You write the economics of your model: its parameters, the exogenous Markov process, the
equilibrium conditions at a point of the state space, the next-period quantities that enter
expectations, a steady state, and the grids. The package does the rest: it builds the state space,
iterates the solution to convergence, simulates, computes impulse responses and moments, compares
economies, and plots.

A solved model is represented by functions of the state, linearly interpolated on a rectangular
(possibly irregular) grid:

- `Pfct` (policy functions): the unknowns of the equilibrium conditions at each grid point.
- `Vfct` (value functions): quantities whose next-period values enter today's expectations.
- `Tfct` (transition functions): the next-period value of each endogenous state variable, one
  value per realization of the exogenous state.
- `Ffct` (forecasting functions, optional): current-period quantities that next period's
  expectations depend on besides the next-period state.

The first state variable is always the index of the exogenous state in its Markov chain. The
solver starts all functions at the deterministic steady state and iterates with damped updates.

Every result (simulations, moments, impulse responses, error statistics) is a
[Tables.jl](https://tables.juliadata.org) table that records where it came from, so any
Tables.jl-compatible writer can save it; see [Tables and Provenance](@ref).

## Installation

```julia
using Pkg
Pkg.add("TransitionFunctionIteration")
```

The package requires Julia 1.12 or later.

## Getting started

The [Quick Start](@ref) solves, checks, simulates and analyzes the neoclassical growth model.
[Implementing a Model](@ref) describes the interface your model provides,
[Simulation and Analysis](@ref) the simulation, accuracy checks, impulse responses and plots, and
[Moments and Comparing Economies](@ref) tables of moments across economies.

## References

- Elenev, V., Landvoigt, T., and Van Nieuwerburgh, S. (2021). "A Macroeconomic Model with
  Financially Constrained Producers and Intermediaries." *Econometrica*, 89(3), 1361–1418.
- Cao, D., Luo, W., and Nie, G. (2023). "Global DSGE Models." Working paper.
