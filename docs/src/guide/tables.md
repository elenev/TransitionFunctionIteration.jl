# Tables and Provenance

Every result type of the package is a [Tables.jl](https://tables.juliadata.org) column table, so
it can go straight into a `DataFrame`, a CSV or Parquet writer, or DuckDB, without copying its
data: the columns are views of the stored arrays. Each table also carries table-level metadata
(DataAPI.jl), which records what it is and where it came from, so that it can be restored as the
same type from a file that keeps metadata, such as Parquet.

```julia
using DataFrames
DataFrame(stats)       # a MomentArray
DataFrame(simseries)   # a SimulationPanel
DataFrame(irf)         # an IRFResult
```

## The tables

| Type | Columns |
|:--|:--|
| `MomentArray` | the row names, in the column named by the array's row key, then one column per statistic |
| `SimulationPanel`, `SimulationErrorPanel` | `__run`, `__period`, then one column per variable |
| `IRFResult` | `__shock`, `__stat`, `__period`, then one column per variable |
| `irf_paths(irf)` | `__shock`, `__run`, `__period`, then the saved variables |

Names are written as strings, since file formats such as Parquet have no symbol type, and turned
back into symbols when a table is read.

### Row keys

The row names of a `MomentArray` are written in a column of their own. Its name, the *row key*, is
part of the array's type (`MomentArray{ROWS,COLS,T,M,ROWKEY}`), and whoever creates the array
chooses it to say what the rows are:

| Created by | Row key |
|:--|:--|
| `describe` (moments and standard errors), `linearize` | `variable` |
| `calc_error_stats` | `equation` |
| `compare_economies` | `moment` |
| the parameter rows of `compare_economies` | `parameter` |
| `MomentArray(data, rows, cols)` | `__row` |

```julia
m = MomentArray(data, (:a, :b), (:mean, :std); rowkey = :variable)
row_key(m)                                   # :variable
Tables.columnnames(Tables.columns(m))        # (:variable, :mean, :std)
```

Arithmetic, `vcat`, `hcat` and sub-selection keep the key; `vcat` of arrays with different keys
gives `__row`. A column name cannot equal the row key. Because the key is in the type, writing the
table is type-stable and allocates nothing beyond the column tuple.

## Metadata

The table-level metadata are strings, ready for any writer, and DataAPI's `metadata` reads them:

| Key | Meaning |
|:--|:--|
| `tfi.kind` | the type that wrote the table (`MomentArray`, `SimulationPanel`, ...) |
| `tfi.rowkey` | a `MomentArray`'s row key |
| `tfi.endog` | a `SimulationPanel`'s endogenous state names (a JSON array) |
| `tfi.provenance` | the provenance of a simulation or an IRF (JSON) |

## Restoring from a table

Each type has a constructor from any Tables.jl source:

```julia
stats = MomentArray(table)                       # names and row key from the table and its metadata
stats = MomentArray(table; rowkey = :variable)   # a table without metadata: say where the names are
simseries = SimulationPanel(table)               # needs "tfi.endog", or: SimulationPanel(table, (:K,))
irf = IRFResult(stats_table, paths_table)
```

The restore checks `tfi.kind` when the source has it, so a table of one kind cannot be read as
another. For a `MomentArray` you can fix more of the type to get a restore whose result type is
known at compile time:

```julia
MomentArray{ROWS,COLS}(table)                    # the table's names are checked against ROWS and COLS
MomentArray{ROWS,COLS,Float64,SMatrix{3,2,Float64,6},:variable}(table)   # inferred, no allocation
```

A static matrix type is filled straight from the table's columns.

## Provenance

A `SimulationPanel`, `SimulationErrorPanel` and `IRFResult` carry a `Provenance`:

- the id of the solution they were computed from, [`solution_id`](@ref), a SHA-1 hash of the grid
  and of the values of the policy, value, transition and forecasting functions, so that solutions with
  identical values have identical ids, including one restored from a file;
- the TFI version;
- the settings of the run (`SimulationSettings`: seed, burn-in, alignment and starting states;
  `IRFSettings`: seed, number of draws, shocks).

It is published as JSON under `tfi.provenance`, which DuckDB reads with `json_extract`.

## Saving to disk

The package defines the tables and their metadata and does not write files. The companion package
OutputVault.jl saves solutions, simulations, moments and IRFs in a fixed layout, with the git
commit that produced them, in plain formats that other languages read.
