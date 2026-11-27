# Saving and Comparing Results with OutputVault

[OutputVault.jl](https://github.com/elenev/OutputVault.jl) saves the results of research runs, one
folder per run with the git state of the code that produced them, and queries them across runs.
TransitionFunctionIteration's support for it is a package extension, loaded with `using OutputVault`:
it says how TFI's types are stored, and adds the functions on this page.

This page continues with `ngm` from the [Quick Start](@ref). Install OutputVault 0.4 and Parquet2
in the model's project (`Pkg.add(url = "https://github.com/elenev/OutputVault.jl"); Pkg.add("Parquet2")`).
OutputVault's own [quick start](https://github.com/elenev/OutputVault.jl/blob/main/docs/src/guide/quickstart.md)
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

A project package typically calls `set_config!` in its `__init__`. A run is named by a
`RunHandle(variant; branch, experdef)`; the branch defaults to the git branch of the active project's
repository, whose commit and diffs the run records.
Pass `branch` explicitly on a detached HEAD or outside a git repository; it chooses a folder
without switching branches. Provenance still comes from the active project's repository.

## Saving

```@example vault
result = solve_and_simulate(ngm; execute_kwargs = (; logging = COMPLETE))
rf = RunHandle("base"; branch = "docs", experdef = "v1")
save_result(rf, result; start = :ergodic)
list_run(rf)
```

The default IRF experiment starts at the simulated mean. Here `start = :ergodic` labels that
experiment in storage; `start` never changes how an IRF was computed. The default storage label
`:ss` is only a name. To compute an IRF from the deterministic steady state, call
`simulate_irf(result.obj, result.stv.endog)` and save that IRF with `start = :ss`.

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

A NamedTuple without `obj` adds its items to the existing run and leaves the solution alone:

```@example vault
irf_ss = simulate_irf(result.obj, result.stv.endog)
save_result(rf, (; irf = irf_ss); start = :ss)
nothing # hide
```

Saving a new solution deletes all previous items of the run. Saving `stats` through `save_result`
replaces the whole sample collection; to add just one sample, use
`save_result_item(rf, "stats/low_K", sample_moments)`. Interpolants computed later can be added with
`save_result(rf, (; recursions = (; claims, welfare)))`. Other values can be saved at any unreserved
path with OutputVault's `save_result_item`.

With `run_simulate = false` or `run_irf = false`, `solve_and_simulate` returns an ordinary tuple.
Name its members before passing them to `save_result`, for example
`obj, stv = solve_and_simulate(ngm; run_simulate = false); save_result(rf, (; obj, stv))`.

## Portable files

Solutions, steady states and interpolants are stored as their portable data, [`TransitionFunctionIteration.to_portable`](@ref):
nested NamedTuples of numbers, strings and numeric arrays, without Julia objects, so they read
without the code that wrote them (h5py reads the files too). [`TransitionFunctionIteration.from_portable`](@ref) rebuilds them;
both work without OutputVault, for any storage of plain data. The
model is not stored as an object, only its parameters, `to_portable(model)`. When a portable
solution is loaded, the model is either supplied, and its parameters are compared with the stored
ones (differences are a warning), or rebuilt from the stored parameters with
`from_portable(M, parameters)` when a model type `M` is supplied. A model package whose model has
fields that are derived objects adds methods of both, imported from TransitionFunctionIteration
(they are not exported).

The native copy, `native.jld2`, holds the same solution as Julia objects: it reloads fastest, with
the same code.

## Loading

```@example vault
loaded = load_result(rf; model = ngm)
solution_id(loaded.obj) == solution_id(result.obj)
```

```@example vault
keys(loaded.irf)                                 # both starting-point labels
```

Load only what you need, without reading the solution when you want statistics or IRFs:

```julia
load_result(rf, :solution; model = ngm)          # (; obj, stv) from the portable files
load_result(rf, :solution; model = NeoclassicalGrowthModel)  # rebuild the model from parameters
load_result(rf, :solution; native = true)        # the native copy
load_result(rf, :stats)                         # moments by sample, in alphabetical order
load_result(rf, :irf).ergodic                    # one IRFResult
```

The full loader returns the recognized TFI items that are present, including `params`.
Its `irf` is a NamedTuple by starting-point label, while `solve_and_simulate` returns a single
`IRFResult`. Custom items saved at other paths are loaded separately through OutputVault.

[`load_result`](@ref) knows which TFI type each path holds. OutputVault's `load_result_item` reads a
single item, as data or as a type: `load_result_item(rf, "solution/obj", ModelSolution; context = ngm)`.
`write_object_file(path, obj)` and `read_object_file(path, ModelSolution; context = ngm)` do the same
for one file outside any run.

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
No simulation or solution objects are reconstructed to build the comparison table.

The numbers the moments need are read with one query per table (statistics, IRFs, parameters), which
opens only the files of the selected runs. A run that holds none of them is left out with a warning,
and a number that a run lacks is a `missing` cell. Variants that repeat across
branches or experdefs are an error: select runs with a narrower pattern. IRFs are read from the
starting point `start` (keyword, default `:ss`).

## Runs saved with OutputVault 0.3

Runs saved with OutputVault 0.3 load with `load_result`. Their IRF is `irf/irf` and `irf/irf_paths`,
which `load_result(rf, :irf)` returns as starting point `ss` but queries do not find.
[`TransitionFunctionIteration.migrate_irf!`](@ref), a transition aid to be removed in a later version, moves them to `irf/ss/stats` and `irf/ss/paths`:

```julia
TFI.migrate_irf!(RunPattern(branch = "**"))             # lists what would move
TFI.migrate_irf!(RunPattern(branch = "**"); move = true) # moves
```
