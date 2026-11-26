# TransitionFunctionIteration.jl

Global solution of discrete-time dynamic equilibrium models by Transition Function Iteration:
time iteration over policy, value, transition and (optional) forecasting functions on a
discretized state space. The package supplies grids, interpolation, the solver loop, simulation,
IRFs, moments and result tables. A model package supplies the economics by dispatching on a
subtype of `AbstractModel`.

## Commands

```sh
julia --project -e 'using Pkg; Pkg.test()'    # tests (CI uses JULIA_NUM_THREADS=2)
julia --project=docs docs/make.jl             # documentation, including the executed examples
```

## Where things are

| Topic | Location |
|---|---|
| Module, exports, `solve_and_simulate` | `src/TransitionFunctionIteration.jl` |
| Model interface | `src/stubs.jl`, `docs/src/guide/model.md` |
| Setup and steady state | `src/setup.jl` |
| Solver | `src/execution.jl` |
| Simulation, IRFs, moments | `src/simulation.jl`, `src/irf.jl`, `src/moments.jl` |
| Reference model | `examples/ngm.jl` |
| Moments as a long table (`assemble_moments`, `MomentRows`) | `src/moment_tables.jl` |
| Saving with OutputVault: stubs, model metadata | `src/storage.jl`, `src/model_metadata.jl` |
| OutputVault extension: storage of TFI's types, portable files, `save_result`, queries | `ext/TransitionFunctionIterationOutputVaultExt.jl` |

Read the interface contracts in the model guide and the docstrings before changing them;
`generate_template` (`src/template.jl`) must stay consistent with them.

## Model interface

Required: `Exogenv`, `steadystate` (or `steadystate_system`, or a fixed-point mode),
`define_endog_grids`, `compute_next_period_values` and `calc_eqm`.

- `calc_eqm(state, policy, mdl, np, mode)` returns residuals `fx` in `SolveMode`,
  `(; fx, trans, vals)` in `IterateMode` and `(; fx, trans, simvars)` in `SimulateMode`.
  `trans` entries are scalars or `SVector`s with one element per exogenous state.
- `compute_next_period_values(next, mdl, mode)` describes a single next-period exogenous state
  and returns `NextPeriod(; expectations, continuations)`. The package applies the probabilities
  and stacks the states.

## Solution methods

`TFIMethod` (default) fixes expectations and forecasts during each grid-point solve and varies only
the policy. `STPFIMethod` also solves for next-period states and recomputes expectations within the
solve. Both have the same fixed point. TFI cannot solve for a policy variable that enters the
residuals only through next-period states or forecasts: every point fails and `main_execute!`
returns `status = :all_failed`. Use `STPFIMethod`, or write the residual so that the variable also
enters directly.

Judge accuracy by simulated equation and transition errors (`calc_error_stats`, `check_bounds`),
not by `status` alone.

## Invariants

- `:exst`, the integer exogenous state, is the first state axis. It is not interpolated and throws
  out of range; endogenous axes are linear with flat extrapolation.
- The solver works in unconstrained space and `Pfct` coefficients are stored unconstrained.
  Convert with `to_constrained` and `to_unconstrained` at the boundary.
- The solver differentiates with ForwardDiff, so code on the solve path must accept dual numbers.
- Code on the solve path must be type stable and avoid heap allocations (`SVector`, `SMatrix`,
  `StructArray`).
- Grid points and simulation runs are threaded (`threaded_foreach`), so the functions they call
  must be thread safe.

## OutputVault extension

- OutputVault does all file I/O. The extension only defines `to_storage`, `from_storage`, `write_portable_jld2`
  and `read_portable_jld2` for TFI's types, and conveniences built on OutputVault's functions.
- The portable `.jld2` layout of solutions, steady states and interpolants is documented at the top
  of the extension. Bump its schema version only for changes that break readers of earlier files.
- `assemble_moments` gives the same table for results in memory and for saved runs; a test checks it.
  Change both methods together.
- Runs saved with OutputVault 0.3 must keep loading with `load_result`.

## Code style

- Short, modular functions with brief section comments. Keep diffs to existing code minimal.
- Early `return`, `continue` and `break` go in explicit `if` blocks, never `cond || return`.
- Fix a bug where it arises rather than adding optional arguments or helpers around it.
- Tabular results are exposed through Tables.jl and DataAPI.jl only. Do not use or mention
  DataFrames.jl anywhere, including docs and examples.

## Testing

- Run the tests before finishing a change.
- Cover type stability and allocations as well as output (see `test/test_performance.jl`).
- A failing test may point to a bug: diagnose it before changing the test.
- Include order in `test/runtests.jl` matters: `test_ngm.jl` defines the solved model that several
  later files use.

## Documentation and releases

- Update docstrings and `docs/src` with any change in behavior. Write them for users, in terms of
  observable behavior.
- Derivations are in Typst (`examples/intermedap.typ`); use Typst syntax there, not LaTeX.
- Do not change the version in `Project.toml` unless asked.
