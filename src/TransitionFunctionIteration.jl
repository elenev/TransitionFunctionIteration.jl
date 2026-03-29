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
export simulate, simulate_markov_chain, sample_shocks, calc_error_stats, check_bounds, describe, linearize, E, augment_simulation
export ReportedMoment, compare_economies
export SimulationPanel, SimulationErrorPanel, endog_states, states

# IRFs
export IRFResult, Transformations, ReportedVariable, IRFPage, simulate_irf

# MomentArray
export MomentArray, row_names, col_names

# Main interface 
export solve_and_simulate

include("abstract_types.jl")
include("named_interpolations.jl")
include("solution_types.jl")
include("statistics_types.jl")
include("exog.jl")
include("utils.jl")
include("logging.jl")
include("plotting.jl")

include("setup.jl")
include("execution.jl")
include("simulation.jl")
include("economics.jl")

include("irf.jl")
include("irfplotting.jl")


function solve_and_simulate(model::AbstractModel;
    setup_kwargs=(;), execute_kwargs=(;), simulate_kwargs=(;), irf_kwargs=(;),
    run_simulate=true, run_irf=true)

    # Setup
    @info "Setting up the model..."
    obj, stv = main_setup(model; setup_kwargs...)

    # Execute
    @info "Solving..."
    main_execute!(obj; execute_kwargs...)

    if !run_simulate
        return obj, stv
    end

    # Simulate
    @info "Simulating to produce ergodic moments..."
    simseries, errs = simulate(obj, stv; simulate_kwargs...)

    stats, _ = describe(simseries).moments.uncond, nothing
    errstats = calc_error_stats(errs)

    if !run_irf
        return obj, stv, simseries, errs
    end

    # IRF
    @info "Running IRFs..."
    irf = simulate_irf(obj, simseries; irf_kwargs...)

    return (; obj, stv, simseries, errs, stats, errstats, irf)
end

end