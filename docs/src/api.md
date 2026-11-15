# API Reference

## Core workflow

```@docs
solve_and_simulate
main_setup
main_execute!
```

## Model interface (stubs)

These functions have default implementations that throw `MethodNotImplementedError`.
Model packages must add methods dispatching on their `AbstractModel` subtype.

```@docs
calc_eqm
compute_next_period_values
steadystate
TransitionFunctionIteration.steadystate_system
define_endog_grids
augment_simulation
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
AbstractMoment
SimulationMoment
ImpulseMoment
EconomyResults
MomentTable
getmoment
moment_keys
compare_economies
describe
calc_error_stats
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
        TransitionFunctionIteration.AbstractMoment,
        TransitionFunctionIteration.SimulationMoment,
        TransitionFunctionIteration.ImpulseMoment,
        TransitionFunctionIteration.EconomyResults,
        TransitionFunctionIteration.MomentTable,
        TransitionFunctionIteration.getmoment,
        TransitionFunctionIteration.moment_keys,
        TransitionFunctionIteration.compare_economies,
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
