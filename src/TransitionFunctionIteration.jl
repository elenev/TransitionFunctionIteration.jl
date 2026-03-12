module TransitionFunctionIteration

using Interpolations, NonlinearSolve, ForwardDiff, StaticArrays, StructArrays
using RecipesBase
using Logging, LinearAlgebra, Random, Statistics

# Types
export AbstractModel, AbstractModelSolution, AbstractExogenv, AbstractStateSpaceGrid, AbstractInterpolant
export StateSpaceGrid, Interpolant, Exogenv, ModelSolution
export EvaluationMode, SolveMode, IterateMode, SimulateMode
export points, dims, npt, nof, varnames, unigrids, singletons, ssidx, exnpt, nex
export coefficients, knots, statespacegrid, name_values
# Interface stub: model packages extend this
export calc_eqm

# Logging
export ExecutionLogger
export COMPLETE, ITERATION, ITERATIONFAILED, FAILEDPOINT
export @complete, @iteration, @iterationfailed, @failedpoint

# Solver
export main_execute!, forecast_next, forecast_at_point, make_point

# Simulation
export simulate, simulate_markov_chain, sample_shocks, calc_error_stats, check_bounds

include("types.jl")
include("logging.jl")
include("plotting.jl")
include("execution.jl")
include("simulation.jl")

end
