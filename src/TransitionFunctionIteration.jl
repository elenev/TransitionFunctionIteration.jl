module TransitionFunctionIteration

using Interpolations, NonlinearSolve, ForwardDiff, StaticArrays, StructArrays
using Plots
import Plots: plot
using Logging, LinearAlgebra, Random, Statistics

# Types
export AbstractEconomics, AbstractModel, AbstractExogenv, AbstractStateSpaceGrid, AbstractInterpolant
export StateSpaceGrid, Interpolant, Exogenv, GenericModel
export EvaluationMode, SolveMode, IterateMode, SimulationMode, SimulateMode
export points, dims, npt, nof, names, unigrids, singletons, ssidx, exnpt, nex
export coefficients, knots, statespacegrid, name_values
# Interface stub: model packages extend this
export calc_eqm

# Logging
export ExecutionLogger
export COMPLETE, ITERATION, ITERATIONFAILED, FAILEDPOINT
export @complete, @iteration, @iterationfailed, @failedpoint

# Solver
export main_execute!

# Simulation
export simulate, simulate_markov_chain, sample_shocks, calc_error_stats, check_bounds

include("types.jl")
include("logging.jl")
include("plotting.jl")
include("execution.jl")
include("simulation.jl")

end
