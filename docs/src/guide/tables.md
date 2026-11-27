# Tables and Provenance

Every result of the package is a [Tables.jl](https://tables.juliadata.org) table: moments and error
statistics (`MomentArray`), simulations (`SimulationPanel`, `SimulationErrorPanel`) and impulse
responses (`IRFResult`). Any package that accepts a Tables.jl table, such as a CSV or Parquet writer
or a database, accepts them as they are, without a copy of the data. Each table also carries
metadata ([DataAPI.jl](https://github.com/JuliaData/DataAPI.jl)) that records what it is and where
it came from, so that it can be restored as the same type from a file that keeps the metadata.

This page uses the solved neoclassical growth model of `examples/ngm.jl`, which the hidden first
block runs.

```@setup tab
using TransitionFunctionIteration
include(joinpath(pkgdir(TransitionFunctionIteration), "examples", "ngm.jl"))
irf = simulate_irf(obj, simseries; path_vars = (:C, :K))
```

## The tables

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

The first column of a `MomentArray` holds its row names, and its name says what the rows are:

```@example tab
Tables.columnnames(Tables.columns(errstats))
```

| Created by | Name of the row column |
|:--|:--|
| `describe`, `linearize` | `variable` |
| `calc_error_stats` | `equation` |
| `compare_economies` | `moment`, including rows added by `parameters` |
| `MomentArray(data, rows, cols)` | `__row`, or the `rowkey` keyword |

`row_key(m)` returns the name. A statistic cannot have the same name.

## Metadata

The metadata are strings, read with `DataAPI.metadata`:

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

A `SimulationPanel`, `SimulationErrorPanel` and `IRFResult` record where they came from:

- the id of the solution they were computed from, [`solution_id`](@ref). It is a hash of the grid and
  of the values of the solution's functions, so two solutions with the same values have the same id,
  including a solution restored from a file;
- the version of TransitionFunctionIteration.jl;
- the settings of the run: for a simulation the seed, burn-in and starting states, for an IRF the
  seed, number of paths and shocks.

```@example tab
simseries.provenance
```

## Saving and restoring results

The package does not write files: pass the tables to any Tables.jl writer. To restore a table as its
type later, use a format that keeps the metadata, such as Parquet. With
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

A constructor checks `tfi.kind` when the table has it, so a table of one kind cannot be restored as
another. A table without metadata, such as one read from a CSV file, needs what the metadata would
have said:

```julia
MomentArray(table; rowkey = :variable)   # the name of the row column
SimulationPanel(table, (:K,))            # the endogenous states
```

`MomentArray{ROWS,COLS}(table)` also checks that the table has the row names `ROWS` and the
statistics `COLS`.

A solution (`ModelSolution`) is not a table. For a portable representation of solutions, steady
states and interpolants, use `TransitionFunctionIteration.to_portable(x)` and rebuild with
`TransitionFunctionIteration.from_portable(T, data; model)`. These functions return and consume
plain data without doing file I/O. [Saving and Comparing Results with OutputVault](@ref) shows how
to save them as portable files and reload a model from its parameters.

Julia's standard library `Serialization` also saves and
restores it, in a file that is read back with the same versions of Julia and the packages:

```@example tab
using Serialization
serialize(joinpath(dir, "obj.jls"), obj)
obj2 = deserialize(joinpath(dir, "obj.jls"))
solution_id(obj2) == solution_id(obj)
```

[OutputVault.jl](https://github.com/elenev/OutputVault.jl) is a separate, optional package that
saves solutions and results of research runs in an organized layout of plain files, using the tables
and metadata described here ([Saving and Comparing Results with OutputVault](@ref)).
