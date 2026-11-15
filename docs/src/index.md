# TransitionFunctionIteration.jl

`TransitionFunctionIteration.jl` is a Julia package for globally solving discrete-time dynamic
equilibrium models using the Transition Function Iteration (TFI) algorithm of Elenev, Landvoigt,
and Van Nieuwerburgh (2021), or optionally the Simultaneous Transition and Policy Function
Iteration (STPFI) algorithm of Cao, Luo, and Nie (2023).

## Choosing a solution method

`TFIMethod` computes expectations at each grid point before solving for the policy, so it misses
any effect of the policy on next period's expectations. Models where a policy acts only through
those expectations (e.g. a portfolio entering a next-period consumption stage) need
`STPFIMethod`, which recomputes expectations from the solver's current guess.

## What the package does

A solved model is represented by three sets of functions of the state variables, linearly
interpolated on a (possibly irregular) rectangular grid:

- **`Pfct`** (policy functions): solutions to the nonlinear system of equations defining
  equilibrium at each grid point.
- **`Vfct`** (value functions): equilibrium quantities and prices whose expected next-period
  values enter today's equilibrium conditions.
- **`Tfct`** (transition functions): the next-period value of each endogenous state variable,
  one value per realization of the exogenous state.

The first state variable is always an integer index into the Markov chain governing the
exogenous state dynamics. The solver initializes all three sets of functions at their
deterministic steady-state values and iterates to convergence using damped updates.

The package also provides utilities to simulate stochastic paths, compute conditional and
unconditional moments, simulate impulse responses, and plot results.

## Installation

```julia
using Pkg
Pkg.add("TransitionFunctionIteration")
```

The package requires Julia 1.12 or later.

## Getting started

See the [Quick Start](@ref) guide for a minimal working example using the Neoclassical
Growth Model, and [Implementing a Model](@ref) for a full description of the interface a
model package must provide.

## References

- Elenev, V., Landvoigt, T., and Van Nieuwerburgh, S. (2021). "A Macroeconomic Model with
  Financially Constrained Producers and Intermediaries." *Econometrica*, 89(3), 1361–1418.
- Cao, D., Luo, W., and Nie, G. (2023). "Global DSGE Models." Working paper.
