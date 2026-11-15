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
| `N_runs` | `1` | Number of simulation runs (paths) |
| `N_periods` | `1000` | Length of each run |
| `burn` | `200` | Burn-in periods discarded |
| `seed` | `100` | RNG seed for reproducibility |

Runs are simulated in parallel across threads, so a model's `calc_eqm`, `augment_simulation`
and `augment_errors` must be thread-safe (no shared mutable state or global RNG use).

## Computing moments

`describe` summarizes a `SimulationPanel` and returns a `NamedTuple` of `MomentArray`s:

```julia
out = describe(simseries)
out.moments.uncond   # MomentArray of unconditional moments
```

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

`MomentArray` rows are variable names and columns are statistics (`:mean`, `:std`, `:auto`,
etc.). Several access patterns are available:

```julia
stats = out.moments.uncond
stats(:C, :mean)                      # scalar
stats((:C, :K), :mean)                # NamedTuple of means
stats(:C, (:mean, :std))              # NamedTuple of stats for C
stats[Val((:C,:K)), Val((:mean,))]    # MomentArray sub-view (compile-time dimensions)
```

## Comparing economies

`compare_economies` produces a formatted table of selected moments across multiple solved
models. Define moments with `SimulationMoment` (a statistic of a simulated variable) or `ImpulseMoment` (a point on an impulse response), then pass a `NamedTuple` of economies: `describe` outputs, or `EconomyResults(describe(...).moments, irf)` when impulse moments are included:

```julia
moment_defs = (
    MeanC = SimulationMoment(:C),
    StdK  = SimulationMoment(:K; stat=:std),
    MeanY = SimulationMoment(:Y),
)

desc_base = describe(simseries_base)
desc_alt  = describe(simseries_alt)

compare_economies(moment_defs, (; baseline=desc_base, alt=desc_alt))
```

## Checking solution accuracy

`calc_error_stats` summarizes the `SimulationErrorPanel` returned by `simulate`:

```julia
errstats = calc_error_stats(errs)
```

The resulting table reports max and mean absolute residuals for each equilibrium equation,
evaluated at the simulated states. Large errors (above ~1e-4) indicate the solver may not
have converged sufficiently or the grid is too coarse near the ergodic distribution.

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
| `shocks` | `(;)` | `NamedTuple` of shock sizes (default: one-std deviations) |

`irf` is an `IRFResult`. Access paths as:

```julia
irf.means.A.C    # mean C response path to an A shock (length T)
irf.stdevs.A.K   # std of K response path to an A shock
```

## Plotting

The package provides Plots.jl recipes for the main output types. Load a backend
(e.g., `using Plots`) before calling `plot`.

### Simulation panel scatter matrix

```julia
plot(simseries)
```

Produces a P×P scatter/histogram matrix of the endogenous state variables.

### Simulation error grid

```julia
plot(errs, simseries)
```

Rows correspond to exogenous states, columns to error variables. Shows the distribution
of residuals across the ergodic simulation.

### Policy/value function slice

```julia
plot(obj.Pfct, :C, (; exst=3))
```

Plots one output variable against the remaining free state dimensions, holding all others
fixed. Produces a line plot for one free dimension and a filled-contour plot for two.

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
