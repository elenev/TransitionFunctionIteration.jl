# Moments and Comparing Economies

A typical table of results reports the same numbers (a mean, a volatility, a point of an impulse
response, a parameter) for several economies. This page builds such a table for two versions of the
neoclassical growth model of the [Quick Start](@ref), then explains its parts. The baseline economy
is the solution `obj` of `examples/ngm.jl`, which the hidden first block of this page runs.

```@setup mom
using TransitionFunctionIteration
include(joinpath(pkgdir(TransitionFunctionIteration), "examples", "ngm.jl"))
```

## Two economies

The second economy has more volatile productivity. Its solution starts from the baseline one:

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
compare_economies(moments, economies; models = (baseline = ngm, high_vol = ngm_hi))
```

Each entry of `moments` names one number per economy. Defining it computes nothing: it is a
recipe that says where to find the number, and `compare_economies` reads it from each economy. The
rows of the table are the moments and its columns the economies:

- `mean_C` is the mean of `C` in the unconditional sample.
- `std_C` is the standard deviation of `C`. In the baseline economy it is the value itself; in
  `high_vol` it is reported relative to the baseline, as `100 (value / baseline - 1)`, so consumption
  is that many percent more volatile.
- `mean_K_lowK` is the mean of `K` in the `low_K` sample defined above.
- `beta` is the parameter `β`.
- `models` adds a row for each parameter that differs between the models, here `σ`, so the table
  shows what was varied.

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

Each takes two more keywords, `transformation` and `compare_transformation`, described next.
An `ImpulseMoment` also accepts `basis = :uncond`, the name of its reference shock path.
For example, `ImpulseMoment(:C, :up, 2; transformation = RatioPct(), basis = :control)`
reports the percentage response relative to `:control`. The reference path is read only when
[`Transformations.requires_basis`](@ref) is true: levels and percentages need only the requested
shock, while differences and ratios need the reference path too. A missing required path is a
`KeyError`; path ordering does not affect the choice. Unlike moments, plots default to the first path.

## Transformations

A transformation turns a value `x` into the number reported, using a reference value, its *basis*:

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

A moment has two transformations, applied in different places:

- `transformation` (default `Level()`) is applied to the value within one economy, and only in the
  base economy of a table. For a `SimulationMoment` or a `ParameterMoment` the basis is the default
  of the transformation, so only the level and percent types are meaningful. For an
  `ImpulseMoment` the basis is the path named by `basis` (`:uncond` by default) at the same period, so `DiffLevel()` gives the
  response and `RatioPct()` the percentage response.
- `compare_transformation` (default `Level()`) is applied in every other economy, to that economy's
  value with the base economy's value as the basis. The `transformation` is not applied there.

The base economy is the first by default; choose another with `base`:

```@example mom
compare_economies(moments, economies; base = :high_vol)
```

Now `high_vol` shows the values themselves, and `std_C` in the baseline column is relative to
`high_vol`: consumption is that many percent less volatile in the baseline.

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

The *source* is the object the moment is read from. It can be an `EconomyResult`, whose parts are
each optional; the `moments` of a `describe` output, for `SimulationMoment`s; or an `IRFResult`, for
`ImpulseMoment`s. These can also be the economies of `compare_economies`.
`getmoment(economies, moment, economy, base)` reads one cell of a table.

A missing sample, variable, statistic, shock or parameter is a `KeyError` that names it:

```@example mom
try
    getmoment(economies.baseline, SimulationMoment(:C; sample = :crisis))
catch err
    err
end
```

A period outside the impulse response is a `BoundsError`. A moment that needs a part the
`EconomyResult` lacks, such as an `ImpulseMoment` without an IRF, is an `ArgumentError`, as is a
moment read from a source that cannot serve it.

## Advanced: reading moments from your own results

Results kept in your own format can be a source too. `getmoment` and `compare_economies` read every
moment through one function, `momentvalues(source, moment)`, which returns the pair
`(value, basis)`: the value before its transformation, and the basis of the transformation
(`nothing` for the default basis, which every moment but an `ImpulseMoment` uses). A source needs a
method for each kind of moment it serves; calculated moments then work with no further code, since
they read their inputs from the same source:

```@example mom
struct MyEconomy
    values::Dict{Symbol, Float64}
end

TransitionFunctionIteration.momentvalues(e::MyEconomy, m::SimulationMoment) =
    (e.values[Symbol(m.stat, "_", m.var)], nothing)

mine = (a = MyEconomy(Dict(:mean_C => 1.0, :std_C => 0.02)),
        b = MyEconomy(Dict(:mean_C => 1.1, :std_C => 0.03)))
compare_economies((mean_C = SimulationMoment(:C), cv_C = @CalculatedMoment sd(:C) / avg(:C)), mine)
```

A moment type with no method for a source is an `ArgumentError`. A new kind of moment is a subtype of
[`AbstractMoment`](@ref) with a `momentvalues` method for each source that can serve it, and methods
of `TransitionFunctionIteration.transformation` and `compare_transformation` if it does not have
fields of those names.

Saved runs are read through such a source, built from [`assemble_moments`](@ref): with OutputVault loaded,
`compare_economies` takes runs as economies ([Saving and Comparing Results with OutputVault](@ref)).
