# TransitionFunctionIteration.jl

TransitionFunctionIteration.jl solves discrete-time dynamic equilibrium models over a grid of
states. Use it to compute policy functions, simulate a model, measure solution errors, and compare
moments and impulse responses across economies.

You provide a model's parameters, equilibrium equations, exogenous Markov process, steady state
and endogenous grids. The [Quick Start](@ref) shows how to define these for a neoclassical growth
model, solve it and analyze the results.

## Installation

```julia
using Pkg
Pkg.add("TransitionFunctionIteration")
```

Requires Julia 1.12 or later. Install StaticArrays and StructArrays as well to run the model
examples and the file that `generate_template` writes.

## How a solution is represented

A solved model is represented by functions of the state, linearly interpolated on a rectangular
grid whose nodes may be unevenly spaced:

- `Pfct` (policy functions): the unknowns of the equilibrium conditions at each grid point.
- `Vfct` (value functions): quantities whose next-period values enter today's equilibrium
  conditions.
- `Tfct` (transition functions): the next-period value of each endogenous state, one value per
  realization of the next-period exogenous state.
- `Ffct` (forecasting functions, optional): current-period quantities that next period's
  expectations depend on besides the next-period state.

The first state variable is the index of the exogenous Markov chain. You choose the grid for each
endogenous state. The solver iterates on the functions until their changes fall below the requested
tolerance. Check simulated equation errors and state bounds before using a solution for analysis;
[Checking solution accuracy](@ref) explains how.

The default method is Transition Function Iteration (TFI). Simultaneous Transition and Policy
Function Iteration (STPFI) solves for policies and next-period states jointly;
[Choosing a solution method](@ref) explains when to use it. Both support the same model interface.

## Getting started

- [Quick Start](@ref): implement and solve the neoclassical growth model.
- [Implementing a Model](@ref): write your own equations, grids and model-specific calculations.
- [Simulation and Analysis](@ref): simulate a solution, check errors and plot impulse responses.
- [Moments and Comparing Economies](@ref): build tables of results across economies.
- [Tables and Provenance](@ref): read result tables and save them with their metadata.
- [Saving and Comparing Results with OutputVault](@ref): save solutions and compare runs on disk.

## References

- Elenev, V., Landvoigt, T., and Van Nieuwerburgh, S. (2021). "A Macroeconomic Model with
  Financially Constrained Producers and Intermediaries." *Econometrica*, 89(3), 1361–1418.
- Cao, D., Luo, W., and Nie, G. (2023). "Global DSGE Models." Working paper.
