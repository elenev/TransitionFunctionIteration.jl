# Tables and Provenance

Pass simulations, moments, error statistics and impulse responses directly to a
[Tables.jl](https://tables.juliadata.org)-compatible writer or database. Each result exposes named
columns and carries metadata describing its type and origin. A file format that preserves that
metadata lets you reload the result as its original TFI type.

This page uses `obj`, `simseries`, `stats` and `errstats` from the neoclassical growth model in
the [Quick Start](@ref). You can also obtain them by running `examples/ngm.jl`.

```@setup tab
using TransitionFunctionIteration
include(joinpath(pkgdir(TransitionFunctionIteration), "examples", "ngm.jl"))
irf = simulate_irf(obj, simseries; path_vars = (:C, :K))
```

## Reading columns

```@example tab
using Tables
Tables.istable(simseries)
```

```@example tab
Tables.columnnames(Tables.columns(simseries))
```

| Type | Columns |
|:--|:--|
| `MomentArray` | the row names, then one column per statistic |
| `SimulationPanel`, `SimulationErrorPanel` | `__run`, `__period`, then one column per variable |
| `IRFResult` | `__shock`, `__stat`, `__period`, then one column per variable |
| `irf_paths(irf)` | `__shock`, `__run`, `__period`, then the variables in `path_vars` |

Names, such as the row names of a `MomentArray` or the shock names of an IRF, are written as strings,
since file formats such as Parquet have no symbol type, and are turned back into symbols when a table
is restored.

A `MomentArray` uses its first column for row names. The column is named for the kind of result;
for example, error statistics have an `equation` column:

```@example tab
Tables.columnnames(Tables.columns(errstats))
```

| Created by | Name of the row column |
|:--|:--|
| `describe`, `linearize` | `variable` |
| `calc_error_stats` | `equation` |
| `compare_economies` | `moment`, including rows added by `parameters` |
| `MomentArray(data, rows, cols)` | `__row`, or the `rowkey` keyword |

`row_key(m)` returns the name. No statistic column may share the row column's name.

## Metadata

Read a result's file metadata through [DataAPI.jl](https://github.com/JuliaData/DataAPI.jl):

```@example tab
import DataAPI
DataAPI.metadata(simseries)
```

| Key | Meaning |
|:--|:--|
| `tfi.kind` | the type of the table (`MomentArray`, `SimulationPanel`, ...) |
| `tfi.rowkey` | a `MomentArray`'s row column |
| `tfi.endog` | a `SimulationPanel`'s endogenous state names (a JSON array) |
| `tfi.provenance` | the provenance of a simulation or an IRF (JSON) |

## Provenance

Check `.provenance` to identify the solution and settings used for a simulation or IRF. The
record contains:

- the solution's [`solution_id`](@ref), computed from its grid and function values, which can be
  compared with the id of a solution reloaded from a file;
- the version of TransitionFunctionIteration.jl;
- the settings of the run: for a simulation the seed, burn-in and starting states, for an IRF the
  seed, number of paths and shocks.

```@example tab
simseries.provenance
```

## Saving and restoring results

Pass a result directly to a Tables.jl writer to save it. To restore its original type later, use
a format that keeps the metadata, such as Parquet. With
[Parquet2.jl](https://gitlab.com/ExpandingMan/Parquet2.jl), pass the metadata to the writer:

```@example tab
using Parquet2
dir = mktempdir()

Parquet2.writefile(joinpath(dir, "stats.parquet"), stats; metadata = DataAPI.metadata(stats))
Parquet2.writefile(joinpath(dir, "sim.parquet"), simseries; metadata = DataAPI.metadata(simseries))
Parquet2.writefile(joinpath(dir, "irf.parquet"), irf; metadata = DataAPI.metadata(irf))
Parquet2.writefile(joinpath(dir, "irf_paths.parquet"), irf_paths(irf); metadata = DataAPI.metadata(irf_paths(irf)))
nothing # hide
```

Each type has a constructor from any Tables.jl table:

```@example tab
stats2 = MomentArray(Parquet2.Dataset(joinpath(dir, "stats.parquet")))
stats2 == stats
```

```@example tab
sim2 = SimulationPanel(Parquet2.Dataset(joinpath(dir, "sim.parquet")))
sim2.provenance == simseries.provenance
```

```@example tab
irf2 = IRFResult(Parquet2.Dataset(joinpath(dir, "irf.parquet")),
                 Parquet2.Dataset(joinpath(dir, "irf_paths.parquet")))
irf2.means.shock_5.C == irf.means.shock_5.C
```

Constructors use the stored `tfi.kind` to check the result type. If the file has no metadata,
as with a plain CSV, supply the row column or endogenous state names explicitly:

```julia
MomentArray(table; rowkey = :variable)   # the name of the row column
SimulationPanel(table, (:K,))            # the endogenous states
```

`MomentArray{ROWS,COLS}(table)` also checks that the table has the row names `ROWS` and the
statistics `COLS`.

Solutions, steady states and interpolants are not tables. To save them in a file format of your
choice, or in run folders with OutputVault, see [Portable files](@ref).

For reloading a solution with the same versions of Julia and its packages, you can also use
Julia's standard library `Serialization`:

```@example tab
using Serialization
serialize(joinpath(dir, "obj.jls"), obj)
obj2 = deserialize(joinpath(dir, "obj.jls"))
solution_id(obj2) == solution_id(obj)
```

For saving a solution together with its simulations, moments and IRFs, see
[Saving and Comparing Results with OutputVault](@ref).
