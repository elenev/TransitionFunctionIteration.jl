module TransitionFunctionIteration

using Interpolations, NonlinearSolve, ForwardDiff, StaticArrays, StructArrays
using RecipesBase
using Logging, LinearAlgebra, Random, Statistics, StatsBase, PrettyTables
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
export simulate, simulate_markov_chain, sample_shocks, calc_error_stats, check_bounds, describe, linearize, augment_simulation, dynamic_depth, augment_errors
export AbstractMoment, SimulationMoment, ImpulseMoment, ParameterMoment, CalculatedMoment, @CalculatedMoment, EconomyResult, getmoment, compare_economies, collect_parameters
export SimulationPanel, SimulationErrorPanel, endog_states, states

# Tables and provenance
export Provenance, SimulationSettings, IRFSettings, solution_id, irf_paths

# IRFs
export IRFResult, Transformations, ReportedVariable, IRFPage, simulate_irf

# MomentArray
export MomentArray, @MomentArray, row_names, col_names, row_key

# Main interface
export solve_and_simulate, generate_template, update_model

# Saved results (methods with OutputVault)
export save_result, load_result

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
include("portable.jl")


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

## Saved results

# Saving and loading results. The methods come with the package extension
# TransitionFunctionIterationOutputVaultExt, loaded with `using OutputVault`.

"""
    save_result(rf, result; start = :ss, dirty_repo_check = :warn) -> String

Save the output of [`solve_and_simulate`](@ref) (or a NamedTuple with some of its properties) in
the OutputVault run `rf`. With `obj`, the run starts over: `save_solution(rf, (; obj, stv))`, then
the parameters at `params/params` ([`collect_parameters`](@ref)). Then, for each property present:
`simseries` and `errs` at `simulation/simseries` and `simulation/errs`, `errstats` at
`errstats/errstats`, `stats` at `stats/<sample>`, the IRF `irf` at `irf/<start>/stats` and
`irf/<start>/paths`, and the interpolants of `recursions` (e.g. `claims`, `welfare`) at
`recursions/<name>`. Without `obj` the run must exist, and only the other properties are saved.
Returns the run's folder. If the new solution cannot be moved into place (on Windows, open files can
block it), `save_solution` warns and nothing else is saved; the staging folder holding the solution
is returned and the old run is left as it was. Needs `using OutputVault`.
"""
function save_result end

"""
    load_result(rf, item::Symbol; model = nothing, native = false)
    load_result(rf; model = nothing, native = false) -> NamedTuple

Load a saved item of the OutputVault run `rf` as TFI's types: `:solution` (`(; obj, stv)`),
`:params` and `:errstats` (`MomentArray`s), `:simulation` (`(; simseries, errs)`), `:stats` (moments by
sample), `:irf` (a NamedTuple of `IRFResult`s by starting point) or `:recursions` (a NamedTuple of
`Interpolant`s). Without `item`, every item the run holds, in one NamedTuple with the properties
of `solve_and_simulate`'s output.

The portable solution files store the model's parameters ([`to_portable`](@ref)), not the model: `model` is a model (its
parameters are compared with the stored ones, and differences are a warning) or a model type (the
model is rebuilt with [`from_portable`](@ref)). `native = true` reads the solution's native
copy instead, which needs no model but the code that saved it. Needs `using OutputVault`.
"""
function load_result end

"""
    migrate_irf!(rf)
    migrate_irf!(pattern; move = false)

Move the IRF files of runs saved with OutputVault 0.3 (`irf/irf` and `irf/irf_paths`) to the
current paths `irf/ss/stats` and `irf/ss/paths`, so that queries find them. For a `RunPattern`,
every run found is listed with what would move, and the files move only with `move = true`.
Needs `using OutputVault`. A transition aid for runs saved with OutputVault 0.3: not exported,
and to be removed in a later version.
"""
function migrate_irf! end

end