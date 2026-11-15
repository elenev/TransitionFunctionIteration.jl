# Moments and Comparing Economies

A *moment* is a lazy specification of one number per economy: it says which statistic of which
variable, or which point of which impulse response, and how to transform it. Nothing is
computed when a moment is defined, so the same list of moments can be read from economies held in
memory or from saved results, and tabulated with [`compare_economies`](@ref).

## Defining moments

```julia
using TransitionFunctionIteration
using TransitionFunctionIteration.Transformations: DiffLevel, Pct

moments = (
    mean_C = SimulationMoment(:C),                          # mean of C in the unconditional sample
    std_K  = SimulationMoment(:K; stat = :std),
    ac_Y   = SimulationMoment(:Y; stat = :ac, transformation = Pct()),
    mean_K_crisis = SimulationMoment(:K; sample = :crisis),
    irf_C  = ImpulseMoment(:C, :shock_5, 2; transformation = DiffLevel()),
)
```

- [`SimulationMoment`](@ref)`(var; sample = :uncond, stat = :mean, transformation, compare_transformation)`
  is a statistic of a simulated variable in one sample of `describe`.
- [`ImpulseMoment`](@ref)`(var, shock, period; stat = :mean, transformation, compare_transformation)`
  is a point on an impulse response: `stat` of `var` under `shock`, `period` periods into the IRF
  (`period` is the index into the IRF lines, so 1 is the starting period). Its transformation is
  applied with the `uncond` path at the same period as the basis, so `DiffLevel()` is the response
  and `RatioPct()` the percentage response. The IRF is assumed to start at the stochastic steady
  state, where `uncond` is the no-shock baseline.

Both are subtypes of [`AbstractMoment`](@ref). A new kind of moment is a new subtype with a
[`moment_keys`](@ref) and a `momentvalues` method.

### Transformations

`transformation` is applied when the moment is reported, and `compare_transformation` when a
non-base economy is compared with the base one (below). They come from the `Transformations`
module and are called as `t(x; basis = ...)`:

| Family | Types | Result |
|:--|:--|:--|
| levels | `Level`, `AnnLevel`, `AnnInvLevel`, `AnnQuadLevel` | `x`, times 4, divided by 4, times 2 |
| percent | `Pct`, `AnnPct`, `AnnInvPct`, `AnnQuadPct` | the same, times 100 |
| differences | `DiffLevel`, `DiffAnnLevel`, ..., `DiffPct`, `DiffAnnPct`, ... | `x - basis` with the scaling of the prefix |
| ratios | `RatioPct`, `RatioAnnPct` | `(x / basis - 1) * 100` |
| from the start | `DiffInitLevel`, `DiffInitPct` | difference from the first element |

## Reading a moment

[`getmoment`](@ref)`(source, moment)` reads a moment from a source of results, and which property
it reads follows from the pair (source, moment type):

| Source | Moments it can serve |
|:--|:--|
| `describe(...).moments` (a `NamedTuple` of `MomentArray`s by sample) | `SimulationMoment` |
| an `IRFResult` | `ImpulseMoment` |
| [`EconomyResults`](@ref)`(stats, irf)` | both |
| a [`MomentTable`](@ref) (below), for one economy | both |

```julia
simseries, errs = simulate(obj, stv)
irf = simulate_irf(obj, simseries)
results = EconomyResults(describe(simseries).moments, irf)

getmoment(results, moments.std_K)
getmoment(results, moments.irf_C)
```

A missing sample, variable, statistic or shock is a `KeyError` that names it, a period outside the
IRF is a `BoundsError`, and an `ImpulseMoment` on results without an IRF is an `ArgumentError`.

## Comparing economies

[`compare_economies`](@ref) tabulates moments (rows) across economies (columns) in a `MomentArray`
whose row key is `moment`:

```julia
economies = (baseline = EconomyResults(stats_base, irf_base),
             alt      = EconomyResults(stats_alt,  irf_alt))
compare_economies(moments, economies)
compare_economies(moments, economies; base = :alt)
```

- In the base economy (the first by default) a moment is its `transformation` of the value.
- In every other economy it is the `compare_transformation` of the value with the base economy's
  value as the basis. With `compare_transformation = DiffPct()` the table shows the difference from
  the base in percent; with the default, `Level()`, it shows the value itself.
- A bare variable name in place of a moment means `SimulationMoment(name)`: `(C = :C, K = :K)`.
- With `models = (; baseline = mdl_base, alt = mdl_alt)`, the parameters that differ between the
  models are added as rows, so that the table shows what was varied.

The economies may also be `describe` outputs (a `NamedTuple` of `MomentArray`s by sample) when
all moments are `SimulationMoment`s.

## From a table of statistics

Results saved to disk are read back as a long table, one row per statistic, for example by
OutputVault's DuckDB query. A [`MomentTable`](@ref) wraps any such Tables.jl source with the
columns `economy`, `kind`, `group`, `variable`, `stat`, `period` and `value`:

| `kind` | `group` | `period` |
|:--|:--|:--|
| `simulation` | the sample (`uncond`, ...) | 0 |
| `impulse` | the shock (`up`, ..., and the baseline `uncond`) | the IRF period |

[`moment_keys`](@ref)`(moment)` lists the rows a moment needs, so a query can fetch only those.
With the table, the same `getmoment` and `compare_economies` apply:

```julia
table = MomentTable(rows)                    # rows: any Tables.jl source
getmoment(table, moments.std_K, :baseline)   # one economy, as for any other source
compare_economies(moments, (:baseline, :alt), table)   # a MomentArray as above
```

Because both routes go through `momentvalues`, an in-memory comparison and one from saved
results give the same numbers.
