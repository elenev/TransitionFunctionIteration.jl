# Moments and Comparing Economies

Use `compare_economies` to report means, volatilities, impulse responses and parameters side by
side. Each row specifies a number to read from the economies' results and how to compare it with
a baseline.

The example compares two versions of the neoclassical growth model with different productivity
volatility. Start with `obj`, `stv` and `simseries` from the [Quick Start](@ref), or obtain them by
running `examples/ngm.jl`.

```@setup mom
using TransitionFunctionIteration
include(joinpath(pkgdir(TransitionFunctionIteration), "examples", "ngm.jl"))
```

## Two economies

Increase productivity volatility from `0.02` to `0.03`, using the baseline solution as the starting
guess for the new solve:

```@example mom
ngm_hi = update_model(ngm; σ = 0.03)
obj_hi, stv_hi = main_setup(ngm_hi, obj)      # warm start from the baseline solution
main_execute!(obj_hi; logging = COMPLETE)
sim_hi, errs_hi = simulate(obj_hi, stv_hi)
calc_error_stats(errs_hi)
```

For each economy, collect its moments (with a sample of low-capital periods), its impulse responses
and its parameters in an [`EconomyResult`](@ref):

```@example mom
samples = (uncond = Returns(true), low_K = x -> x.K < 0.97 * stv.endog.K)

economies = (
    baseline = EconomyResult(; stats  = describe(simseries; filters = samples).moments,
                               irf    = simulate_irf(obj, simseries),
                               params = collect_parameters(obj)),
    high_vol = EconomyResult(; stats  = describe(sim_hi; filters = samples).moments,
                               irf    = simulate_irf(obj_hi, sim_hi),
                               params = collect_parameters(obj_hi)),
)
nothing # hide
```

## A table of moments

```@example mom
using TransitionFunctionIteration.Transformations: RatioPct, DiffPct

moments = (
    mean_C      = SimulationMoment(:C),
    std_C       = SimulationMoment(:C; stat = :std, compare_transformation = RatioPct()),
    mean_K_lowK = SimulationMoment(:K; sample = :low_K),
    beta        = ParameterMoment(:β),
)
compare_economies(moments, economies; parameters = :distinct)
```

Each entry of `moments` specifies a number to read from each economy's results. The rows of the
table are the moments and its columns the economies:

- `mean_C` is the mean of `C` in the unconditional sample.
- `std_C` is the standard deviation of `C`. In the baseline economy it is the value itself; in
  `high_vol` it is reported relative to the baseline, as `100 (value / baseline - 1)`, so consumption
  is that many percent more volatile.
- `mean_K_lowK` is the mean of `K` in the `low_K` sample defined above.
- `beta` is the parameter `β`.
- `parameters = :distinct` adds a row above the moments for each parameter that differs between
  the economies, here `σ`, so the table shows what was varied. The parameters are those of each
  `EconomyResult`; `parameters = (:β, :σ)` shows the ones named instead.

A bare variable name in place of a moment means `SimulationMoment(name)`, so `(C = :C, K = :K)` is a
valid collection of moments.

## Kinds of moments

- [`SimulationMoment`](@ref)`(var; sample = :uncond, stat = :mean)` is a statistic of a simulated
  variable in one sample of `describe` (see [Computing moments](@ref) for the statistics).
- [`ImpulseMoment`](@ref)`(var, shock, period; stat = :mean)` is a point of an impulse response:
  `stat` of `var` under `shock` in `period`, where period 1 is the starting period and the forced
  exogenous path starts in period 2. (Plots count the starting period as 0, so `period = 2` is the
  point plotted at 1.)
- [`ParameterMoment`](@ref)`(name)` is a model parameter, as [`collect_parameters`](@ref) returns it.
- [`CalculatedMoment`](@ref) is a function of other moments (below).

Each takes `transformation` to choose the units of the reported value and `compare_transformation`
to compare that value with another economy; the next section explains both.

## Transformations

Choose a transformation for the units you want to report. Differences and ratios use a reference
value called the *basis*:

| Family | Types | Result |
|:--|:--|:--|
| levels | `Level`, `AnnLevel`, `AnnInvLevel`, `AnnQuadLevel` | `x`, `4x`, `x / 4`, `2x` |
| percent | `Pct`, `AnnPct`, `AnnInvPct`, `AnnQuadPct` | the same, times 100 |
| differences | `DiffLevel`, `DiffAnnLevel`, `DiffAnnInvLevel`, `DiffAnnQuadLevel` | `x - basis`, with the scaling of the level type |
| | `DiffPct`, `DiffAnnPct`, `DiffAnnInvPct`, `DiffAnnQuadPct` | the same, times 100 |
| ratios | `RatioPct`, `RatioAnnPct` | `100 (x / basis - 1)` |
| from the start | `DiffInitLevel`, `DiffInitPct` | `x - init`, `100 (x - init) / init`, with `init` the first value |

The `Ann` prefixes annualize quarterly values: a flow is multiplied by 4, a standard deviation by 2
(the square root of 4), and `AnnInv` divides by 4. The types are in the `Transformations`
submodule, which exports no names: import the ones you use, as above.

`DiffPct` and `RatioPct` answer different questions. Comparing 2.2 with a basis of 2.0, `DiffPct`
gives `100 (2.2 - 2.0) = 20` and `RatioPct` gives `100 (2.2 / 2.0 - 1) = 10`. For a rate, such as an
interest rate of 0.022 against 0.020, `DiffPct` is the difference in percentage points (0.2) and
`RatioPct` the relative change in percent (10). For a level, such as consumption, the relative
change `RatioPct` is usually what you want.

A moment can first transform a value within an economy and then compare it with the baseline:

- `transformation` (default `Level()`) gives the reported value within each economy. For a
  `SimulationMoment` or a `ParameterMoment` the basis is the default of the transformation, so only
  the level and percent types are meaningful. For an `ImpulseMoment` the basis is the path named by
  `basis` (`:uncond` by default) at the same period, so `DiffLevel()` gives the response and
  `RatioPct()` the percentage response. For example,
  `ImpulseMoment(:C, :up, 2; transformation = RatioPct(), basis = :control)` reports the
  percentage response relative to the `:control` path. Level and percent transformations need only
  the requested shock path; the others raise `KeyError` if the basis path is missing. IRF plots
  default to the first path in `irf.means` instead, so pass the same `basis` to both to compare
  them.
- `compare_transformation` (default `Level()`) is applied in every economy other than the base, to
  its reported value with the base economy's reported value as the basis. With the default
  `Level()`, every economy shows its reported value.

The base economy is the first by default; choose another with `base`:

```@example mom
compare_economies(moments, economies; base = :high_vol)
```

Now `high_vol` shows the values themselves, and `std_C` in the baseline column is relative to
`high_vol`: consumption is that many percent less volatile in the baseline.

The `compare_transform` keyword of `compare_economies` switches all comparisons on or off; it is
not the same as a moment's `compare_transformation`, which chooses how one moment is compared. With
`compare_transform = false` every economy shows the moment's own value. `stack_moments` puts
the table in long form, columns `economy`, `moment` and `value` with one row per cell, to save or to
plot:

```@example mom
stack_moments(compare_economies(moments, economies; compare_transform = false))
```

## Impulse responses across economies

An `ImpulseMoment` with a `transformation` reports the response, here the percentage response of
`C` relative to the unconditional path, in every economy:

```@example mom
using TransitionFunctionIteration.Transformations: DiffLevel

response_C = ImpulseMoment(:C, :shock_5, 2; transformation = RatioPct())
response_diff = ImpulseMoment(:C, :shock_5, 2; transformation = RatioPct(), compare_transformation = DiffLevel())
compare_economies((; response_C, response_diff), economies)
```

A `compare_transformation` compares the reported values: `response_diff` is the baseline response
in the base column and, in the `high_vol` column, the difference of the `high_vol` response from the
baseline one, in percentage points.

## Calculated moments

[`CalculatedMoment`](@ref)`(f, inputs)` is a function of other moments. `inputs` is a `NamedTuple` of
moments, and `f` receives a `NamedTuple` with the same names and the inputs' values in one economy:

```@example mom
C_over_Y = CalculatedMoment(nt -> nt.C / nt.Y, (C = SimulationMoment(:C), Y = SimulationMoment(:Y)))
compare_economies((; C_over_Y), economies)
```

Each input is read with its own `transformation`, and its `compare_transformation` is ignored. The
calculated moment's own `transformation` and `compare_transformation` apply to the result, as for any
moment. Two moments built from equal anonymous functions are not equal (`==`); use a named function
if you need to compare or store them.

[`@CalculatedMoment`](@ref) writes the same thing as an expression. Every function call whose
arguments include a moment becomes a calculated moment, so any function of numbers works, operators
included:

```@example mom
avg(v) = SimulationMoment(v)
sd(v)  = SimulationMoment(v; stat = :std)

C_over_Y_expr = @CalculatedMoment avg(:C) / avg(:Y)
rel_vol       = @CalculatedMoment sd(:C) / sd(:Y)
log_C         = @CalculatedMoment log(avg(:C))
compare_economies((; C_over_Y_expr, rel_vol, log_C), economies)
```

With a collection of moments as the first argument, a name that is not a variable in scope is read
from the collection:

```@example mom
ratio = @CalculatedMoment moments std_C / mean_C
compare_economies((; ratio), economies)
```

Only function calls become calculated moments: `if`, `&&`, broadcasting and anonymous functions are
left as they are, and a moment passed as a keyword argument is not converted.

## Reading one moment

[`getmoment`](@ref)`(source, moment)` reads one moment, transformed, from one economy:

```@example mom
getmoment(economies.baseline, response_C)
```

Pass an `EconomyResult`, `describe(...).moments` for simulation moments, or an `IRFResult` for
impulse moments. You can also use these directly as the economies passed to `compare_economies`.
`getmoment(economies, moment, economy, base)` reads one cell of a table.

A missing sample, variable, statistic, shock or parameter is a `KeyError` that names it:

```@example mom
try
    getmoment(economies.baseline, SimulationMoment(:C; sample = :crisis))
catch err
    err
end
```

An out-of-range IRF period raises `BoundsError`. Requesting a part an `EconomyResult` lacks raises
`KeyError`; an unsupported source raises `ArgumentError`. In a comparison table, a missing input
or out-of-range period gives a `missing` cell so other economies can still be reported.

## Advanced: reading moments from your own results

To compare results stored in your own type, define `extract_moment_values(source, moment)` for
each kind of moment you need. Return `(value, basis)`, where `value` is the untransformed number
and `basis` is its reference value, or `nothing` to use the transformation's default.
Calculated moments can then combine those inputs as usual:

```@example mom
struct MyEconomy
    values::Dict{Symbol, Float64}
end

TransitionFunctionIteration.extract_moment_values(e::MyEconomy, m::SimulationMoment) =
    (e.values[Symbol(m.stat, "_", m.var)], nothing)

mine = (a = MyEconomy(Dict(:mean_C => 1.0, :std_C => 0.02)),
        b = MyEconomy(Dict(:mean_C => 1.1, :std_C => 0.03)))
compare_economies((mean_C = SimulationMoment(:C), cv_C = @CalculatedMoment sd(:C) / avg(:C)), mine)
```

For custom transformations, [`Transformations.requires_basis`](@ref) controls whether reading an
impulse moment needs a separate reference path. To define a new kind of moment, subtype
[`AbstractMoment`](@ref) and implement `extract_moment_values` for its sources. Supply fields or
methods named `transformation` and `compare_transformation` to choose how it is reported.

To build the same comparison from saved runs, see
[Saving and Comparing Results with OutputVault](@ref).
