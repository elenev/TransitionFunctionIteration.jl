# TransitionFunctionIteration.jl

Generic Julia infrastructure for solving discrete-time dynamic equilibrium models globally via transition function iteration.

Provides:
- **State space / interpolation**: `StateSpaceGrid`, `Interpolant`, `Exogenv` — multidimensional grids with Cartesian products and gridded linear interpolation
- **Solver**: `main_execute!` — fixed-point iteration over policy and value functions using trust-region nonlinear solvers with multithreading and damping
- **Simulation**: `simulate`, `simulate_markov_chain` — Markov chain simulation with equation-error diagnostics
- **Logging / plotting**: `ExecutionLogger`, `Plots.plot` extension for `Interpolant`

Model-specific economics are supplied by the user: define a type `M <: AbstractEconomics` and extend `TransitionFunctionIteration.calc_eqm` dispatching on `GenericModel{M}`.
