# API Reference

## Core workflow

```@docs
solve_and_simulate
main_setup
main_execute!
```

## Model interface (stubs)

Model packages implement the required economics by dispatching on their `AbstractModel` subtype.
Required stubs throw `MethodNotImplementedError`; augmentation hooks have defaults that add nothing.

```@docs
calc_eqm
compute_next_period_values
steadystate
TransitionFunctionIteration.steadystate_system
define_endog_grids
augment_simulation
dynamic_depth
augment_errors
```

## Solution types

```@docs
ModelSolution
SteadyState
NextPeriod
```

### Evaluation modes

```@docs
SolveMode
IterateMode
SimulateMode
```

### Steady-state modes

```@docs
AnalyticalSteadyStateMode
FixedPointSteadyStateMode
SystemSteadyStateMode
```

### Solution methods

```@docs
TFIMethod
STPFIMethod
```

### Variable bounds

```@docs
TransitionFunctionIteration.Bounded
TransitionFunctionIteration.Unbounded
TransitionFunctionIteration.LowerBounded
TransitionFunctionIteration.UpperBounded
TransitionFunctionIteration.IntervalBounded
```

## State space and interpolation

```@docs
StateSpaceGrid
Interpolant
points
```

## Exogenous environment

```@docs
ExogenvAR1
ergodic_distro
```

## Simulation

```@docs
simulate_markov_chain
sample_shocks
compute_expectations
make_point
```

## Impulse responses

```@docs
simulate_irf
IRFResult
ReportedVariable
IRFPage
```

## Moments and comparison

```@docs
MomentArray
row_names
col_names
row_key
stack_moments
AbstractMoment
TransitionFunctionIteration.transformation
TransitionFunctionIteration.compare_transformation
TransitionFunctionIteration.Transformations.requires_basis
SimulationMoment
ImpulseMoment
ParameterMoment
CalculatedMoment
@CalculatedMoment
TransitionFunctionIteration.lifted
EconomyResult
collect_parameters
getmoment
compare_economies
describe(::SimulationPanel)
calc_error_stats
```

## Saved results (with OutputVault)

```@docs
save_result
load_result
TransitionFunctionIteration.migrate_irf!
TransitionFunctionIteration.to_portable
TransitionFunctionIteration.from_portable
```

## Simulation results and provenance

```@docs
SimulationPanel
SimulationErrorPanel
irf_paths
solution_id
Provenance
SimulationSettings
IRFSettings
```

## Logging

```@docs
ExecutionLogger
TerminalLogger
CSVLogger
COMPLETE
ITERATION
ITERATIONFAILED
FAILEDPOINT
```

## Economics utilities

```@docs
fisher_burmeister
```

## Template

```@docs
generate_template
```

## All other documented symbols

```@autodocs
Modules = [TransitionFunctionIteration]
Filter = t -> begin
    already = [
        TransitionFunctionIteration.solve_and_simulate,
        TransitionFunctionIteration.main_setup,
        TransitionFunctionIteration.main_execute!,
        TransitionFunctionIteration.calc_eqm,
        TransitionFunctionIteration.compute_next_period_values,
        TransitionFunctionIteration.steadystate,
        TransitionFunctionIteration.steadystate_system,
        TransitionFunctionIteration.define_endog_grids,
        TransitionFunctionIteration.augment_simulation,
        TransitionFunctionIteration.augment_errors,
        TransitionFunctionIteration.ModelSolution,
        TransitionFunctionIteration.SteadyState,
        TransitionFunctionIteration.NextPeriod,
        TransitionFunctionIteration.SolveMode,
        TransitionFunctionIteration.IterateMode,
        TransitionFunctionIteration.SimulateMode,
        TransitionFunctionIteration.AnalyticalSteadyStateMode,
        TransitionFunctionIteration.FixedPointSteadyStateMode,
        TransitionFunctionIteration.SystemSteadyStateMode,
        TransitionFunctionIteration.TFIMethod,
        TransitionFunctionIteration.STPFIMethod,
        TransitionFunctionIteration.Bounded,
        TransitionFunctionIteration.Unbounded,
        TransitionFunctionIteration.LowerBounded,
        TransitionFunctionIteration.UpperBounded,
        TransitionFunctionIteration.IntervalBounded,
        TransitionFunctionIteration.StateSpaceGrid,
        TransitionFunctionIteration.Interpolant,
        TransitionFunctionIteration.points,
        TransitionFunctionIteration.ExogenvAR1,
        TransitionFunctionIteration.ergodic_distro,
        TransitionFunctionIteration.simulate_markov_chain,
        TransitionFunctionIteration.sample_shocks,
        TransitionFunctionIteration.compute_expectations,
        TransitionFunctionIteration.make_point,
        TransitionFunctionIteration.simulate_irf,
        TransitionFunctionIteration.IRFResult,
        TransitionFunctionIteration.ReportedVariable,
        TransitionFunctionIteration.IRFPage,
        TransitionFunctionIteration.MomentArray,
        TransitionFunctionIteration.row_names,
        TransitionFunctionIteration.col_names,
        TransitionFunctionIteration.row_key,
        TransitionFunctionIteration.stack_moments,
        TransitionFunctionIteration.AbstractMoment,
        TransitionFunctionIteration.SimulationMoment,
        TransitionFunctionIteration.ImpulseMoment,
        TransitionFunctionIteration.ParameterMoment,
        TransitionFunctionIteration.CalculatedMoment,
        TransitionFunctionIteration.var"@CalculatedMoment",
        TransitionFunctionIteration.lifted,
        TransitionFunctionIteration.transformation,
        TransitionFunctionIteration.compare_transformation,
        TransitionFunctionIteration.Transformations.requires_basis,
        TransitionFunctionIteration.EconomyResult,
        TransitionFunctionIteration.collect_parameters,
        TransitionFunctionIteration.getmoment,
        TransitionFunctionIteration.compare_economies,
        TransitionFunctionIteration.save_result,
        TransitionFunctionIteration.load_result,
        TransitionFunctionIteration.migrate_irf!,
        TransitionFunctionIteration.to_portable,
        TransitionFunctionIteration.from_portable,
        TransitionFunctionIteration.describe,
        TransitionFunctionIteration.calc_error_stats,
        TransitionFunctionIteration.SimulationPanel,
        TransitionFunctionIteration.SimulationErrorPanel,
        TransitionFunctionIteration.irf_paths,
        TransitionFunctionIteration.solution_id,
        TransitionFunctionIteration.Provenance,
        TransitionFunctionIteration.SimulationSettings,
        TransitionFunctionIteration.IRFSettings,
        TransitionFunctionIteration.ExecutionLogger,
        TransitionFunctionIteration.TerminalLogger,
        TransitionFunctionIteration.CSVLogger,
        TransitionFunctionIteration.COMPLETE,
        TransitionFunctionIteration.ITERATION,
        TransitionFunctionIteration.ITERATIONFAILED,
        TransitionFunctionIteration.FAILEDPOINT,
        TransitionFunctionIteration.fisher_burmeister,
        TransitionFunctionIteration.generate_template,
    ]
    !(t in already)
end
```
