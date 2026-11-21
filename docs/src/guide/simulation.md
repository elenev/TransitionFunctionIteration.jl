# Simulation and Analysis

## Stochastic simulation

After solving with `main_execute!`, simulate a stochastic panel with:

```julia
simseries, errs = simulate(obj, stv)
```

`simseries` is a `SimulationPanel` -- a `StructArray` of per-period NamedTuples covering all
state and simulation variables defined in `calc_eqm` under `SimulateMode`. `errs` is a
`SimulationErrorPanel` containing the equilibrium residuals at each simulated point, useful
for assessing solution accuracy.

Key keyword arguments for `simulate`:

| Argument | Default | Description |
|---|---|---|
| `sz` | `10_000 => 12` | Periods per run `=>` number of runs |
| `burn` | `500` | Burn-in periods discarded |
| `seed` | `100` | RNG seed for reproducibility |
| `exog_start` | steady-state index | Exogenous state in the first period |
| `endog_start` | `stv.endog` | Endogenous state in the first period |
| `context_periods` | `dynamic_depth(model)` | Periods simulated after the burn-in and dropped after augmentation, so that variables defined from consecutive periods have no `NaN` (see `augment_simulation`) |
| `report_bounds` | `true` | Report the ranges of the simulated states against the grid |

A simulation records where it came from (the solution's id, the TFI version and the settings
above) in its `provenance`; see [Tables and Provenance](@ref).

Runs are simulated in parallel across threads, so a model's `calc_eqm`, `augment_simulation`
and `augment_errors` must be thread-safe (no shared mutable state or global RNG use).

## Computing moments

`describe` summarizes a `SimulationPanel` and returns `(; moments, stderrs)`, each a `NamedTuple`
of `MomentArray`s by sample:

```julia
out = describe(simseries)
out.moments.uncond   # MomentArray of unconditional moments
```

The statistics are `:mean`, `:std`, `:minimum`, `:maximum`, `:ac` (first-order autocorrelation), the
percentiles `:p5`, `:p25`, `:p50`, `:p75`, `:p95`, and `:cor_<state>_lag0` for each endogenous state
(`correlation_vars` chooses the states and lags). Pass `additional_stats` to add statistics, either
functions of one variable's series or `MultivariateStatFun`s that also read another variable.

### Conditional moments

Pass `filters` as a `NamedTuple` of predicates, one per regime:

```julia
out = describe(simseries;
    filters = (
        uncond = Returns(true),
        crisis = x -> x.K < 0.9 * stv.endog.K,
    )
)
out.moments.uncond   # unconditional MomentArray
out.moments.crisis   # crisis-regime MomentArray
```

Pass `dropnan=true` to exclude every observation containing a NaN from all filters. A filter
that selects no observations yields NaN moments rather than an error.

### `MomentArray` access

`MomentArray` rows are variable names and columns are statistics. The names are part of the type,
so lookup by name is resolved at compile time. Several access patterns are available:

```julia
stats = out.moments.uncond
stats(:C, :mean)                      # scalar
stats((:C, :K), :mean)                # NamedTuple of means
stats(:C, (:mean, :std))              # NamedTuple of stats for C
stats((:C, :K), (:mean, :std))        # MomentArray of the selected rows and columns
stats(Val((:C, :K)), Val((:mean,)))   # the same with the names fixed at compile time
```

A `MomentArray` also has a name for its column of row names when it is written as a table (its
*row key*): `variable` for the statistics of `describe`, `equation` for `calc_error_stats`, `moment`
for `compare_economies`. See [Tables and Provenance](@ref).

## Comparing economies

`compare_economies` builds a table of moments across economies. A moment is a lazy specification,
`SimulationMoment` for a statistic of a simulated variable and `ImpulseMoment` for a point on an
impulse response; the economies are `describe` outputs, or `EconomyResult`s that also hold an IRF and the parameters.

```julia
moment_defs = (
    MeanC = SimulationMoment(:C),
    StdK  = SimulationMoment(:K; stat=:std),
    MeanY = SimulationMoment(:Y),
)

desc_base = describe(simseries_base)
desc_alt  = describe(simseries_alt)

compare_economies(moment_defs, (; baseline=desc_base.moments, alt=desc_alt.moments))
```

See [Moments and Comparing Economies](@ref) for impulse moments, transformations, the base
economy, and building the same table from saved results.

## Checking solution accuracy

`calc_error_stats` summarizes the `SimulationErrorPanel` returned by `simulate`:

```julia
errstats = calc_error_stats(errs)
```

The resulting table reports max and mean absolute residuals for each equilibrium equation,
evaluated at the simulated states: one row per equation (and per transition), with the mean and
the percentiles `:p50`, `:p75`, `:p90`, `:p95`, `:p99` and the maximum `:p100` of the absolute
residual. Large errors (above ~1e-4) indicate the solver may not have converged sufficiently or
the grid is too coarse near the ergodic distribution. Judge accuracy by these errors, not by the
convergence flag of the solver.

## Impulse response functions

`simulate_irf` computes impulse responses by comparing simulation paths with and without
a shock, holding the initial endogenous state fixed:

```julia
irf = simulate_irf(obj, simseries)
```

Key keyword arguments:

| Argument | Default | Description |
|---|---|---|
| `T` | `20` | IRF horizon (periods) |
| `N` | `5_000` | Number of draws for averaging |
| `seed` | `100` | RNG seed |
| `shocks` | one per other exogenous state, plus `uncond` | `NamedTuple` of forced exogenous paths: a state index or a vector of indices, e.g. `(up = 5, up2 = [5, 5])`; `uncond = Int[]` forces nothing |
| `exog_start` | steady-state index | Exogenous state before the shock |
| `hold_periods` | `dynamic_depth(model)` | Periods the exogenous state is held at `exog_start` before the IRF; they fill the variables of `augment_simulation` in period 1 and are not output |
| `path_vars` | `()` | Variables whose individual draws are kept (`irf_paths`) |
| `other_stats` | `(;)` | Statistics besides the mean and standard deviation |

`endogstart`, the second argument, is a `NamedTuple` of starting values, or a `SimulationPanel`, in
which case the ergodic mean is used. `irf` is an `IRFResult`. Paths are indexed by shock:

```julia
propertynames(irf.means)   # the shock names, e.g. (:uncond, :shock_1, ...)
irf.means.shock_5.C        # mean C path after the exogenous state moves to state 5 (length T)
irf.stdevs.shock_5.K       # standard deviation of K across draws
```

The first element of a path is the starting period, before the shock. The response to a shock
is the shock path compared with the `uncond` path; `ImpulseMoment` and the plot recipes do this.

## Plotting

The package provides Plots.jl recipes for the main output types. Load a backend
(e.g., `using Plots`) before calling `plot`.

### Simulation panel scatter matrix

```julia
plot(simseries)
```

Produces a P×P matrix of the endogenous state variables. Diagonal panels are histograms;
off-diagonal panels scatter the column variable (x) against the row variable (y).
`legend` is off by default and can be overridden.

Keyword arguments:

| keyword      | default                | meaning |
|:-------------|:-----------------------|:--------|
| `varnames`   | all endogenous states  | variables to plot, in the order of rows and columns |
| `exst`       | all                    | keep only simulation rows whose exogenous state equals this value |
| `maxpoints`  | `nothing` (all points) | randomly subsample scatter panels to at most this many points; histograms always use the full data |
| `stategrid`  | `nothing`              | an `AbstractStateSpaceGrid`; draws thin grey lines at the grid nodes of every plotted endogenous state, within the data range plus the nearest node on each side |

```julia
plot(simseries; exst=2, maxpoints=5_000)
```

For long simulations set `maxpoints`: interactive backends such as `plotly()` serialize every
plotted point and can be very slow on the full data.

### Simulation error grid

```julia
plot(errs, simseries)
```

Rows correspond to exogenous states, columns to error variables. Shows the distribution
of residuals across the ergodic simulation. With one endogenous state on the x-axis, each
panel scatters the error against that state, with y-limits shared across panels. With two,
each panel scatters the two states against each other, colored by the absolute error, with
shared color limits and one colorbar in an extra column. Titles appear only in the top row,
x-axis labels and ticks only in the bottom row, and y-axis labels and ticks only in the
leftmost column. `legend` is off by default.

Keyword arguments:

| keyword       | default                         | meaning |
|:--------------|:--------------------------------|:--------|
| `errnames`    | all fields of `errs`            | tuple of error names to plot (the columns) |
| `exogstates`  | all values in `simseries.exst`  | exogenous state values to plot (the rows) |
| `endogstates` | the first two endogenous states | one or two endogenous state names for the axes |
| `stategrid`   | `nothing`                       | an `AbstractStateSpaceGrid`; draws thin grey lines at its nodes for the selected states, within the data range |
| `maxpoints`   | `nothing` (all points)          | randomly subsample each panel to at most this many points; axis and color limits still use the full data, and one sample is shared across the error columns of a row |

```julia
plot(errs, simseries; endogstates=[:K], maxpoints=2_000)
```

### Policy/value function slice

```julia
plot(obj.Pfct, :C, (; exst=3))
```

Plots one output variable against the remaining free state dimensions, holding all others
fixed. The third argument pins values for all but one or two states; those remaining
states are swept over their full grid. Produces a line plot for one free dimension and a
filled-contour plot for two. For a vector-valued output, pass `(:C, i)` instead of `:C` to
plot element `i`. The title is the variable name, and `legend` is off by default.

Keyword argument:

| keyword     | default    | meaning |
|:------------|:-----------|:--------|
| `transform` | `identity` | function applied elementwise to the plotted values, e.g. `transform=log` |

### IRF plots

Define the variables and transformations to display using `ReportedVariable` and group them
into pages with `IRFPage`:

```julia
using TransitionFunctionIteration.Transformations

page = IRFPage([
    ReportedVariable(:C; label="Consumption", transformation=Pct()),
    ReportedVariable(:K; label="Capital",     transformation=Pct()),
    ReportedVariable(:Y; label="Output",      transformation=Pct()),
])

plot(irf, page, :A)   # IRF to the A shock
```

`Transformations` provides: `Level()`, `Pct()`, `DiffLevel()`, `AnnPct()`, and others.

To plot one variable without building a page, pass its name or a `ReportedVariable`:

```julia
plot(irf, :C, :A)                                          # level of C, labeled "C"
plot(irf, ReportedVariable(:C; transformation=Pct()), :A)  # with a transformation and label
```

A bare name uses `ReportedVariable(:C)`, so the transformation is `Level()` and the label is
the variable name. The keyword arguments below apply to a single variable and to a page.

| keyword     | default                  | meaning |
|:------------|:-------------------------|:--------|
| `pers`      | `0:T-1`                  | periods to plot, counting the starting period as 0 |
| `basis`     | first key of `irf.means` | key of `irf.means` (e.g. `:uncond`) whose path of the same variable is passed as `basis` to the transformation; used by the `Diff*` and `Ratio*` transformations and ignored by the others |
| `init`      | `nothing`                | name of a variable in the shock's mean path; its first value is passed to the transformation as `init`, which is used by `DiffInitLevel` and `DiffInitPct`. When `nothing`, the transformation uses the first value of the plotted variable |
| `trend`     | `nothing`                | name of a variable in the shock's mean path, passed to the transformation as `trend`; none of the built-in transformations use it |
| `sd`        | `false`                  | `true` or a multiplier `m`: shade a ribbon of `m` standard deviations (`true` means 1) around the mean, using `irf.stdevs` |
| `fan`       | `false`                  | `true` or an integer number of bands: plot a fan chart of the simulated paths in `irf.paths` instead of the `sd` ribbon, with the mean path as a line |
| `n_bands`   | `10`                     | number of bands when `fan=true` |
| `tail_pct`  | `1.0`                    | lower tail percentile where the outermost band begins; bands are symmetric, so the default spans the 1st to 99th percentile |
| `fan_alpha` | `0.8`                    | maximum cumulative opacity at the center of the fan |
| `legendon`  | `false`                  | page only: which subplot shows the legend (`false` none, `true` the first, an integer that subplot) |

The page layout follows the shape of the `IRFPage` matrix, with one subplot per
`ReportedVariable`. Each panel is titled with the variable's `label` and its y-axis is
labeled by the transformation.

```julia
plot(irf, page, :A; pers=0:20, sd=2)
plot(irf, page, :A; fan=true, tail_pct=5)
```
