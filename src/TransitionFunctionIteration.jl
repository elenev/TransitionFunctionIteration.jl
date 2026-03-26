module TransitionFunctionIteration

using Interpolations, NonlinearSolve, ForwardDiff, StaticArrays, StructArrays
using RecipesBase
using Logging, LinearAlgebra, Random, Statistics, StatsBase, PrettyTables

# Types
export AbstractModel, AbstractModelSolution, AbstractExogenv, AbstractStateSpaceGrid, AbstractInterpolant, AbstractSteadyState
export StateSpaceGrid, Interpolant, Exogenv, ModelSolution, SteadyState, NextPeriod
export EvaluationMode, SolveMode, IterateMode, SimulateMode
export SteadyStateMode, AnalyticalSteadyStateMode, SystemSteadyStateMode, FixedPointSteadyStateMode
export points, dims, npt, nof, elsize, varnames, unigrids, singletons
export coefficients, knots, statespacegrid, name_values

# Interface stub: model packages extend this
export calc_eqm, compute_next_period_values, steadystate, define_endog_grids

# Logging
export ExecutionLogger
export COMPLETE, ITERATION, ITERATIONFAILED, FAILEDPOINT
#export @complete, @iteration, @iterationfailed, @failedpoint

# Exogenous environment tools
export ExogenvAR1, ergodic_distro, ssidx, exnpt, nex, transition_matrix, conditional_distro

# Solver
export main_setup, main_execute!, make_point, compute_expectations, TFIMethod, STPFIMethod, fisher_burmeister

# Simulation
export simulate, simulate_markov_chain, sample_shocks, calc_error_stats, check_bounds, describe, linearize, E
export SimulationPanel, SimulationErrorPanel, endog_states, states

# MomentArray
export MomentArray, row_names, col_names

include("types.jl")
include("exog.jl")
include("utils.jl")
include("logging.jl")
include("plotting.jl")

include("setup.jl")
include("execution.jl")
include("simulation.jl")
include("economics.jl")

end
