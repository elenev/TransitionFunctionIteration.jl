# Saving and Comparing Results with OutputVault

[OutputVault.jl](https://www.vadimelenev.com/OutputVault.jl/stable/) keeps each solved economy and its results
in a run folder, along with the git commit and any uncommitted changes. Use it to reload a solution
in a later session or compare moments across saved economies.

This page continues with `ngm` from the [Quick Start](@ref). Install OutputVault and Parquet2
in the same Julia project with `Pkg.add(url = "https://github.com/elenev/OutputVault.jl")` and
`Pkg.add("Parquet2")`.
OutputVault's own [quick start](https://www.vadimelenev.com/OutputVault.jl/stable/guide/quickstart/)
introduces runs and queries without TFI.

```@setup vault
using TransitionFunctionIteration
include(joinpath(pkgdir(TransitionFunctionIteration), "examples", "ngm.jl"))
```

## Setup

```@example vault
using TransitionFunctionIteration, OutputVault, Parquet2
const TFI = TransitionFunctionIteration

set_config!(Config(output_root = mktempdir(),            # use a permanent folder for real runs
                   table_format = OutputVault.Formats.ParquetFormat()))   # TFI's tables as Parquet
```

Use a permanent output folder for results you want to keep. A `RunHandle` names the variant,
experiment and branch of a run. By default, the branch and recorded git state come from the
active Julia project's repository. Pass `branch` explicitly on a detached HEAD or outside a
repository; it chooses the output folder without switching branches.

## Saving

```@example vault
result = solve_and_simulate(ngm; execute_kwargs = (; logging = COMPLETE))
rf = RunHandle("base"; branch = "docs", experdef = "v1")   # a fixed branch; usually omitted
save_result(rf, result; start = :ergodic)
list_run(rf)
```

`start = :ergodic` labels the saved IRF, which starts at the simulated mean. `start` is only the
name of the folder the IRF is saved under (`:ss` by default): the starting state is chosen when you
call `simulate_irf`, and `save_result` does not check the label against it.

[`save_result`](@ref) starts the run over with the solution and then saves each output:

| Item | Path | Stored as |
|:--|:--|:--|
| `obj`, `stv` | `solution/obj`, `solution/stv` | portable `.jld2` files, and both in `native.jld2` |
| parameters ([`collect_parameters`](@ref)) | `params/params` | a table, key `parameter`, column `value` |
| `simseries`, `errs` | `simulation/simseries`, `simulation/errs` | tables |
| `errstats` | `errstats/errstats` | a table, key `equation` |
| `stats` (moments by sample) | `stats/<sample>` | one table per sample, key `variable` |
| `irf` | `irf/<start>/stats`, `irf/<start>/paths` | two tables; `start = :ss` by default |
| `recursions` (interpolants) | `recursions/<name>` | portable `.jld2` files |

You can add results later. For example, compute an IRF from the deterministic steady state and
save it alongside the one from the simulated mean:

```@example vault
irf_ss = simulate_irf(result.obj, result.stv.endog)
save_result(rf, (; irf = irf_ss); start = :ss)
nothing # hide
```

Pass only the new result to keep the saved solution. Including `obj` starts the run over and
deletes its previous items. To add to an existing run:

- an IRF: `save_result(rf, (; irf); start = :label)`, as above;
- one sample of moments: `save_result_item(rf, "stats/low_K", sample_moments)`. Saving `stats`
  through `save_result` replaces the whole sample collection;
- interpolants: `save_result(rf, (; recursions = (; claims, welfare)))`;
- any other value: OutputVault's `save_result_item`, at any unreserved path.

With `run_simulate = false` or `run_irf = false`, `solve_and_simulate` returns an ordinary tuple.
Name its members before passing them to `save_result`, for example
`obj, stv = solve_and_simulate(ngm; run_simulate = false); save_result(rf, (; obj, stv))`.

## Loading

To reopen a run in another Julia session, configure the same output folder and create a
`RunHandle` with the same branch, experiment and variant. Load the model's definitions, then pass
a model instance to restore its solution:

```@example vault
loaded = load_result(rf; model = ngm)
solution_id(loaded.obj) == solution_id(result.obj)
```

```@example vault
keys(loaded.irf)                                 # both starting-point labels
```

Load an individual part when you only need the solution, statistics or one IRF experiment:

```julia
load_result(rf, :solution; model = ngm)          # (; obj, stv) from the portable files
load_result(rf, :solution; model = NeoclassicalGrowthModel)  # rebuild the model from parameters
load_result(rf, :solution; native = true)        # the native copy
load_result(rf, :stats)                         # moments by sample, in alphabetical order
load_result(rf, :irf).ergodic                    # one IRFResult
```

The full result contains the saved solution, simulation and summaries, including `params` when
present. Use `loaded.irf.ergodic` or `loaded.irf.ss` to choose an IRF experiment. The `irf` returned
by `solve_and_simulate` is a single `IRFResult`; the loader groups them by their saved labels.

Read custom items with OutputVault's `load_result_item`. You can also use it for a single solution
object: `load_result_item(rf, "solution/obj", ModelSolution; context = ngm)`. To save just one object
outside a run, use `write_object_file(path, obj)` and reload it with
`read_object_file(path, ModelSolution; context = ngm)`.

## Portable files

Solutions, steady states and interpolants are saved as numbers, strings and arrays, so their data
can be read without the original Julia type definitions, including with Python's h5py. To restore
a solution or steady state as a TFI object, supply a model instance or model type with `model`.

With a model instance, the loader compares its parameters with the saved ones and warns about
differences. With a model type, it rebuilds the model from the saved parameters. For a quick reload
with compatible Julia code, use `native = true` to read the original objects from `native.jld2`.

For custom storage, [`TransitionFunctionIteration.to_portable`](@ref) returns the data and
[`TransitionFunctionIteration.from_portable`](@ref) rebuilds an object from it. If your model has
derived fields, define `to_portable(model)` to select its parameters and
`from_portable(M, parameters)` to reconstruct it. Import these functions from
TransitionFunctionIteration when adding methods.

## Comparing economies

[`compare_economies`](@ref) takes saved runs as economies: a NamedTuple of `RunHandle`s, keyed by
the economies' names, or a `RunPattern`, whose runs are named by their variant. It takes the same
keywords as for results in memory and gives the same table for the same results:

```@example vault
ngm_hi = update_model(ngm; σ = 0.03)
hi = solve_and_simulate(ngm_hi; execute_kwargs = (; logging = COMPLETE))
rf_hi = RunHandle("high_vol"; branch = "docs", experdef = "v1")
save_result(rf_hi, hi; start = :ergodic)

using TransitionFunctionIteration.Transformations: RatioPct
moments = (mean_K = SimulationMoment(:K),
           resp_C = ImpulseMoment(:C, :shock_5, 3; transformation = RatioPct()),
           beta = ParameterMoment(:β))
compare_economies(moments, RunPattern(branch = "docs", experdef = "v1");
                  start = :ergodic, base = :base, parameters = :distinct)
```

For explicit economy names, pass `(; baseline = rf, high_vol = rf_hi)` instead of the pattern.
The comparison reads the required statistics, IRFs and parameters directly from the selected runs.
A number missing from one run produces a `missing` cell; a run with none of the requested numbers
is omitted with a warning. If a variant name occurs in more than one branch or experiment, narrow
the pattern or supply explicit economy names. Use `start` to select the saved IRF experiment
(`:ss` by default).

## Comparing IRFs from older runs

To compare IRFs saved with OutputVault 0.3, move their files to the current layout with
[`TransitionFunctionIteration.migrate_irf!`](@ref). Preview the moves first:

```julia
TFI.migrate_irf!(RunPattern(branch = "**"))             # lists what would move
TFI.migrate_irf!(RunPattern(branch = "**"); move = true) # moves
```

This moves `irf/irf` and `irf/irf_paths` to `irf/ss/stats` and `irf/ss/paths`, where comparison
queries can find them. To load an older IRF without moving files, use `load_result(rf, :irf).ss`.
