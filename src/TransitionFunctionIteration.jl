module TransitionFunctionIteration

using Interpolations, NonlinearSolve, ForwardDiff, StaticArrays, StructArrays
using RecipesBase
using Logging, LinearAlgebra, Random, Statistics, StatsBase, PrettyTables
using Polyester: @batch
using Tables: Tables
using MappedArrays: mappedarray
import DataAPI, JSON, SHA

# Types
export AbstractModel, AbstractModelSolution, AbstractExogenv, AbstractStateSpaceGrid, AbstractInterpolant, AbstractSteadyState
export StateSpaceGrid, Interpolant, Exogenv, ModelSolution, SteadyState, NextPeriod
export EvaluationMode, SolveMode, IterateMode, SimulateMode
export SteadyStateMode, AnalyticalSteadyStateMode, SystemSteadyStateMode, FixedPointSteadyStateMode
export points, dims, npt, nof, elsize, varnames, unigrids, singletons
export coefficients, knots, statespacegrid, name_values, named_coefficients

# Interface stub: model packages extend this
export calc_eqm, compute_next_period_values, steadystate, define_endog_grids, forecasts_from_policy

# Logging
export ExecutionLogger, TerminalLogger, CSVLogger
export COMPLETE, ITERATION, ITERATIONFAILED, FAILEDPOINT

# Exogenous environment tools
export ExogenvAR1, ergodic_distro, ssidx, exnpt, nex, transition_matrix, conditional_distro

# Solver
export main_setup, main_execute!, make_point, compute_expectations, TFIMethod, STPFIMethod, fisher_burmeister

# Simulation
export simulate, simulate_markov_chain, sample_shocks, calc_error_stats, check_bounds, describe, linearize, augment_simulation, augment_errors
export AbstractMoment, SimulationMoment, ImpulseMoment, EconomyResults, MomentTable, getmoment, moment_keys, compare_economies
export SimulationPanel, SimulationErrorPanel, endog_states, states

# Tables and provenance
export Provenance, SimulationSettings, IRFSettings, solution_id, irf_paths

# IRFs
export IRFResult, Transformations, ReportedVariable, IRFPage, simulate_irf

# MomentArray
export MomentArray, @MomentArray, row_names, col_names, with_rowkey

# Main interface
export solve_and_simulate, generate_template, update_model

include("abstract_types.jl")
include("named_interpolations.jl")
include("solution_types.jl")

include("transformations.jl")
include("statistics_types.jl")

include("stubs.jl")

include("exog.jl")
include("utils.jl")
include("logging.jl")
include("plotting.jl")

include("setup.jl")
include("execution.jl")
include("simulation.jl")
include("economics.jl")
include("template.jl")

include("irf.jl")
include("moments.jl")
include("irfplotting.jl")
include("provenance.jl")
include("tables.jl")


"""
    solve_and_simulate(model; setup_kwargs=(;), execute_kwargs=(;), simulate_kwargs=(;),
                       irf_kwargs=(;), run_simulate=true, run_irf=true)

Convenience wrapper: calls `main_setup`, `main_execute!`, `simulate`, and `simulate_irf`
in sequence. Returns a NamedTuple `(; obj, stv, simseries, errs, stats, errstats, irf)`.

Set `run_simulate=false` to stop after solving; `run_irf=false` to skip IRF computation.
"""
function solve_and_simulate(model::AbstractModel;
    setup_kwargs=(;), execute_kwargs=(;), simulate_kwargs=(;), irf_kwargs=(;),
    describe_kwargs=(;),
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

    stats = describe(simseries; describe_kwargs...).moments
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