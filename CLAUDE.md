# TransitionFunctionIteration.jl

Generic Julia infrastructure for solving discrete-time dynamic equilibrium models globally via
Transition Function Iteration (TFI). Developed by Vadim Elenev. Julia 1.11+.

---

## What the package does

Implements a time iteration algorithm over policy, value, and transition functions on a
discretized state space. The user supplies model-specific economics (equilibrium conditions,
exogenous process, steady state); the package handles everything else: grid construction,
interpolation, the nonlinear solver loop, simulation, IRFs, and moment reporting.

Two solution methods:
- `TFIMethod` (default): pre-computes expectations from the transition function before solving policy.
- `STPFIMethod`: solves policy and transition functions simultaneously (Cao, Luo, Nie).

Which method a model needs: `TFIMethod` computes expectations before solving for the policy,
so it misses any effect of the policy on next period's expectations. Models where a policy acts
only through those expectations need `STPFIMethod`.

---

## Source layout

```
src/
  abstract_types.jl       -- Abstract type hierarchy + MethodNotImplementedError
  named_interpolations.jl -- StateSpaceGrid, Interpolant, Exogenv
  solution_types.jl       -- Bounded/bounds types, NextPeriod, ModelSolution, SteadyState,
                             EvaluationMode subtypes, SolutionMethod subtypes
  stubs.jl                -- Interface stubs (throw MethodNotImplementedError by default)
  transformations.jl      -- Transformations module (Level, Pct, DiffPct, AnnPct, etc.)
  statistics_types.jl     -- MomentArray, SimulationPanel, SimulationErrorPanel
  moments.jl              -- AbstractMoment (SimulationMoment, ImpulseMoment, ParameterMoment, CalculatedMoment, @CalculatedMoment), EconomyResult, collect_parameters, momentvalues, getmoment
  setup.jl                -- main_setup, steadystate dispatches, bounds fallback
  execution.jl            -- main_execute!, perform_iteration!, solve_at_point, EquilibriumSystem
  simulation.jl           -- simulate, simulate_markov_chain, describe, linearize, compare_economies
  economics.jl            -- fisher_burmeister, crra/logutility/powerutility, ces, cobb_douglas,
                             epstein_zin_ce, return_on_ex_price, growth_rate
  exog.jl                 -- ExogenvAR1
  irf.jl                  -- simulate_irf, IRFResult, IRFPage
  irfplotting.jl          -- Plots recipes for IRFs
  logging.jl              -- ExecutionLogger, TerminalLogger, CSVLogger + log levels
  plotting.jl             -- Plots recipe for Interpolant
  tables.jl               -- Tables.jl interfaces and DataAPI metadata of MomentArray, panels, IRFs; restore constructors
  provenance.jl           -- Provenance, SimulationSettings/IRFSettings, solution_id, tfi_version
  template.jl             -- generate_template (model skeleton)
  utils.jl                -- merge_structarrays, StaticStructArray helpers

examples/
  ngm.jl                  -- Neoclassical Growth Model (canonical reference implementation)
  intermedap.jl           -- Intermediary AP model (more complex, in-progress on rev202603)

test/
  runtests.jl             -- Test suite entry point
```

---

## Type hierarchy

```
AbstractModel               -- user subtype; carries all model parameters
AbstractModelSolution
  ModelSolution{M,E,S,PF,VF,TF,PB,TB,FF}
                            -- assembled solution: .model, .exogenv, .ssg, .Pfct, .Vfct, .Tfct,
                               .Pfct_bounds, .Tfct_bounds, .Ffct (forecasting function or nothing)
AbstractSteadyState
  SteadyState               -- .model, .guess, .endog, .exog, .values
AbstractStateSpaceGrid
  StateSpaceGrid            -- named Cartesian product grid; :exst is the exog index axis
AbstractInterpolant
  Interpolant               -- gridded linear interp over StateSpaceGrid, multiple outputs
AbstractExogenv
  Exogenv                   -- discrete Markov chain: .pts (StructArray), .mtrans, .ssidx

EvaluationMode
  SolveMode     -- compute residuals (for NonlinearSolve)
  IterateMode   -- compute transitions + value outputs (for updating interpolants)
  SimulateMode  -- compute transitions + all simulation variables

SteadyStateMode
  AnalyticalSteadyStateMode   -- user provides closed-form formula
  FixedPointSteadyStateMode   -- solver finds fixed point of dynamic equations
  SystemSteadyStateMode       -- user provides a reduced steady-state system

Bounded (for policy/transition variable constraints)
  Unbounded
  LowerBounded{L}   -- PositiveBounded() = LowerBounded(0.0)
  UpperBounded{U}   -- NegativeBounded() = UpperBounded(0.0)
  IntervalBounded{L,U}  -- IntervalBounded() = [0,1]
```

---

## Interface: what a model package must implement

Implement methods dispatching on `MyModel <: AbstractModel`:

| Function | Required? | Notes |
|---|---|---|
| `Exogenv(mdl::MyModel)` | Yes | Return an `Exogenv` (use `ExogenvAR1` for AR(1)) |
| `calc_eqm(state, policy, mdl, np, mode)` | Yes | Core equilibrium conditions; dispatch on `mode` |
| `compute_next_period_values(next, mdl, mode)` | Yes | Return `NextPeriod(; expectations, continuations)` |
| `steadystate(mdl, ::AnalyticalSteadyStateMode)` | Yes (or other mode) | Return `(SteadyState(...), true)` |
| `define_endog_grids(mdl, stv)` | Yes | Return `NamedTuple` of grid vectors |
| `bounds(stv::SteadyState{MyModel})` | No | Return `(Pfct_bounds, Tfct_bounds)` NTs; default is all `Unbounded()` |
| `augment_simulation(run, mdl)` | No | Return `(static_nt, dynamic_nt)` for extra sim variables |
| `vals_from_policy(state, policy, mdl)` | No | Vfct vars that are functions of the policy but not policies (FixedPointSteadyState); Vfct vars that depend on expectations are `eqmvalguesses` of `FixedPointSteadyStateMode` instead, with residual guess minus `calc_eqm`'s `vals` |
| `forecasts_from_policy(state, policy, mdl)` | No | Forecasting values (current choices next period's expectations depend on); builds `obj.Ffct`, attached to every `next` point by `forecast_next` |

### `calc_eqm` return contract

- `SolveMode`: return `SVector` (or NamedTuple) of residuals `fx`
- `IterateMode`: return `(; fx, trans, vals)` where `trans` is next-period endog state NamedTuple and `vals` is Vfct NamedTuple
- `SimulateMode`: return `(; fx, trans, simvars)` where `simvars` is a NamedTuple of all output variables

Transition values in `trans` can be scalars (deterministic) or `SVector`s (state-contingent).
`expand_transitions` promotes scalars to SVectors automatically.

### `compute_next_period_values` return contract

```julia
return NextPeriod(; expectations, continuations)
# expectations: NamedTuple of probability-weighted scalars (e.g., E[β u′(C′) R′])
# continuations: NamedTuple of per-state values (stacked across exog states by compute_expectations)
```

---

## Main workflow

```julia
obj, stv = main_setup(model)            # builds ModelSolution with initial guess
main_execute!(obj; maxiter=200, tol=-5) # iterates to convergence (modifies obj in-place)
simseries, errs = simulate(obj, stv)    # stochastic simulation -> SimulationPanel
stats = describe(simseries).moments.uncond
errstats = calc_error_stats(errs)
irf = simulate_irf(obj, simseries)
```

Or use the convenience wrapper:
```julia
result = solve_and_simulate(model)
# result = (; obj, stv, simseries, errs, stats, errstats, irf)
```

---

## Testing

```
test/runtests.jl   -- loads all test_*.jl files
```

Run with `julia --project -e "using Pkg; Pkg.test()"` or `]test` in the package REPL.
Each test file covers one subsystem (statespace, exogenv, bounded, economics, interpolant,
execution, simulation, momentarray, logging, plotting, irfplotting, linearize, ngm).

---

## Coding conventions

- State space: always includes `:exst` (integer exog state index) as first axis. The `StateSpaceGrid`
  is constructed as `StateSpaceGrid(exogenv, endog_grids)`.
- Solver operates in unconstrained space; bounds are applied via `to_constrained`/`to_unconstrained`
  at the boundary. `Pfct` coefficients are stored unconstrained.
- `Interpolant` evaluates by NamedTuple lookup (`:exst` uses `NoInterp()`, others use `Gridded(Linear())`).
  Extrapolation: `Flat()` for endogenous, `Throw()` for exst.
- `StructArray` and `StaticStructArray` (SVector-backed) used throughout for performance.
- `StaticArrays` (`SVector`, `SMatrix`) used wherever size is compile-time known; avoids allocation
  inside the hot solver loop.
- `threaded_foreach` (`Threads.@threads`) for multithreading over grid points and simulation runs; it
  limits BLAS to one thread while threaded, to avoid oversubscription.
- Coefficients arrays are indexed as `coefs[cartesian_index]` -> `SVector` of nof values.

---

## Claude Conventions

### Testing
- Unless explicitly told otherwise, run tests at the end of your work.
- When asked to design unit tests, if tests fail, don't just update the test until it passes. Consider if the failed test identifies a bug in the code and propose a fix.
- Tests should ensure type stability and performance (no degradation in allocations) in addition to correct output.

### Math
- Standalone mathematical derivations files are usually in Typst. If updating those, make sure to use Typst syntax, not latex. Watch out for curly braces and backslashes. Those are Latex syntax, not Typst.
- When 


### Memory
- At the end of every session, write a summary of completed/discussed work into `.claude/logs/` with the date and time in the filename.
- At the start of every session, review the logs to guide your code review.


---

## Interface Notes

### `MomentArray` access

```julia
getmoment(ma, :varname, :stat)              # scalar
getmoment(ma, (:v1, :v2), :mean)            # NamedTuple
getmoment(ma, :varname, (:mean, :std))      # NamedTuple
getmoment(ma, Val((:v1,:v2)), Val((:mean,))) # MomentArray sub-view (compile-time)
ma(:varname, :stat)                         # callable shorthand
```


Row keys: `MomentArray{ROWS,COLS,T,M,ROWKEY}`. ROWKEY names the column of row names when the array is
a table (`__row` by default; creators set `variable`, `equation`, `moment`, `parameter`); it is in the
type so that `Tables.columns` is inferred. Arithmetic/vcat/hcat/sub-selection keep it.

### `describe`, moments and `compare_economies`

```julia
# Conditional moments:
out = describe(simseries; filters = (uncond = Returns(true), crisis = x -> x.K < threshold))
out.moments.crisis  # MomentArray for crisis filter (out.stderrs: standard errors)

# Moments are lazy specs; momentvalues(source, moment) -> (value, basis) is the one reading function and
# getmoment(source, moment) applies the transformation. Sources: describe output (SimulationMoment), IRFResult
# (ImpulseMoment), EconomyResult(; stats, irf, params), NamedTuple of economies, or a project's own source type
# (OutputVault: a Dict per economy of saved rows). CalculatedMoment reads from any source serving its inputs.
moment_defs = (
    MeanC = SimulationMoment(:C),
    StdK  = SimulationMoment(:K; stat=:std),
    irf_C = ImpulseMoment(:C, :shock_5, 2; transformation = Transformations.DiffLevel()),
    beta  = ParameterMoment(:β),
    mean_ratio = @CalculatedMoment mean(:C) / mean(:K),   # with mean(v) = SimulationMoment(v); see the moments guide
)
compare_economies(moment_defs, (; baseline = EconomyResult(; stats = stats_b, irf = irf_b, params = collect_parameters(obj_b)),
                                  alt      = EconomyResult(; stats = stats_a, irf = irf_a, params = collect_parameters(obj_a))))
```

### Logging

Log levels (in order of verbosity):
- `COMPLETE` -- convergence messages only
- `ITERATION` -- per-iteration distances (default: displayed)
- `ITERATIONFAILED` -- grid points that failed to converge
- `FAILEDPOINT` -- per-point solver failure details

Pass `logging=COMPLETE` to `main_execute!` to suppress iteration logs.
Pass an `ExecutionLogger` or `CSVLogger` for structured logging to file.


### Exogenous environment helpers

```julia
ExogenvAR1(; ρ, σ, n, transform=exp, jensen=true, shock_name=:A)
ergodic_distro(exogenv)
ssidx(exogenv)          # steady-state index
exnpt(exogenv)          # number of discrete states
conditional_distro(exogenv, state_idx)
```

---

## Active work (branch rev202603)

- `examples/intermedap.jl` is an intermediary asset pricing model under active development.
- `augment_simulation` and related code is being updated and tested.

---