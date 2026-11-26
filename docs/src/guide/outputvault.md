# Saving and Comparing Results with OutputVault

[OutputVault.jl](https://github.com/elenev/OutputVault.jl) saves the results of research runs, one
folder per run with the git state of the code that produced them, and queries them across runs.
TransitionFunctionIteration's support for it is a package extension, loaded with `using OutputVault`:
it says how TFI's types are stored, and adds the functions on this page.

## Setup

```julia
using TransitionFunctionIteration, OutputVault, Parquet2
const TFI = TransitionFunctionIteration

set_config!(Config(output_root = "/path/to/outputs",
                   table_format = OutputVault.Formats.ParquetFormat()))   # TFI's tables as Parquet
```

A project package typically calls `set_config!` in its `__init__`. A run is named by a
`RunHandle(variant; branch, experdef)`; the branch defaults to the git branch of the active project's
repository, whose commit and diffs the run records.

## Saving

```julia
result = solve_and_simulate(mdl)
rf = RunHandle("base"; experdef = "v1")
save_result(rf, result)
```

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

```julia
save_result(rf, (; recursions = (; claims, welfare)))     # computed later from the solution
save_result(rf, (; irf = irf_boom); start = :boom)        # an IRF from another starting point
```

Any other value can be saved at a path of your choice with OutputVault's `save_result_item`.

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

```julia
load_result(rf, :solution; model = mdl)          # (; obj, stv) from the portable files
load_result(rf, :solution; model = MyModel)      # the model rebuilt from its parameters
load_result(rf, :solution; native = true)        # the native copy
load_result(rf, :stats)                          # moments by sample, in alphabetical order
load_result(rf, :irf)                            # (ss = IRFResult, ...), by starting point
load_result(rf; model = mdl)                     # everything the run holds, as solve_and_simulate returns it
```

[`load_result`](@ref) knows which TFI type each path holds. OutputVault's `load_result_item` reads a
single item, as data or as a type: `load_result_item(rf, "solution/obj", ModelSolution; context = mdl)`.
`write_object_file(path, obj)` and `read_object_file(path, ModelSolution; context = mdl)` do the same
for one file outside any run.

## Comparing economies

[`compare_economies`](@ref) takes saved runs as economies: a NamedTuple of `RunHandle`s, keyed by
the economies' names, or a `RunPattern`, whose runs are named by their variant. It takes the same
keywords as for results in memory and gives the same table for the same results:

```julia
moments = (mean_K = SimulationMoment(:K), resp = ImpulseMoment(:K, :pos, 3), beta = ParameterMoment(:β))
compare_economies(moments, RunPattern(experdef = "v1"); base = :base, parameters = :distinct)
compare_economies(moments, (; frm = RunHandle("frm"; experdef = "v1"), arm = RunHandle("arm"; experdef = "v1")))
```

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
