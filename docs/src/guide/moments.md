# Moments and Comparing Economies

A *moment* is a lazy specification of one number per economy: it says which statistic of which
variable, which point of which impulse response, or which parameter, and how to transform it.
Nothing is computed when a moment is defined, so the same list of moments can be read from
economies held in memory or from saved results, and tabulated with [`compare_economies`](@ref).

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
    beta   = ParameterMoment(:β),
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
- [`ParameterMoment`](@ref)`(name; transformation, compare_transformation)` is a model parameter, as
  [`collect_parameters`](@ref) returns it.
- [`CalculatedMoment`](@ref)`(f, inputs; transformation, compare_transformation)` is a function of
  other moments. `inputs` is a `NamedTuple` of moments (which may be calculated moments too), and
  `f` receives a `NamedTuple` with the same keys and the inputs' values in one economy:

```julia
irf_ratio = CalculatedMoment(nt -> nt.irf_C / nt.irf_Y,
                             (irf_C = ImpulseMoment(:C, :shock_5, 2; transformation = DiffLevel()),
                              irf_Y = ImpulseMoment(:Y, :shock_5, 2; transformation = DiffLevel())))
```

  Each input is read with its own `transformation` (its `compare_transformation` is ignored); the
  calculated moment's own `transformation` and `compare_transformation` apply to the result, as
  for any moment, and its transformation gets no basis. It reads from any source that serves all
  its inputs. A lambda is convenient, but two specs with equal lambdas are not `==`; use a named
  function if the spec must be stored.

All are subtypes of [`AbstractMoment`](@ref). A new kind of moment is a new subtype with a
`momentvalues` method per source that can serve it, and the accessors
`TransitionFunctionIteration.transformation` and `compare_transformation` (the fields of those
names by default).

### Calculated moments as expressions

[`@CalculatedMoment`](@ref) writes a calculated moment as an expression of moments. Every function
call in it is made through [`lifted`](@ref TransitionFunctionIteration.lifted): the plain call if no
argument is a moment, otherwise a calculated moment of the moment arguments. So any function of
numbers works, operators included, with no list of supported functions:

```julia
mean(v) = SimulationMoment(v)
sd(v)   = SimulationMoment(v; stat = :std)
cor_r(v) = SimulationMoment(v; stat = :cor_r_lag0)

ols_r(v) = @CalculatedMoment cor_r(v) * sd(v) / sd(:r)        # sd(:r) is a plain call: a SimulationMoment
total    = @CalculatedMoment sqrt(sd(:C)^2 + 2 * mean(:K))
mean_Y_floored = @CalculatedMoment max(mean(:Y), 1.5)         # any function, constants carried along (a number: mean Y, at least 1.5)
```

With a collection of named moments as the first argument, a bare name that is not a variable in
scope is read from it:

```julia
@CalculatedMoment moments mean_C / std_K
```

Only calls are lifted: `if`, `&&`, broadcasting and anonymous functions are left as they are, and a
moment passed as a keyword argument is not lifted.

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
| [`EconomyResult`](@ref)`(; stats, irf, params)` | all three; each part is optional |

A calculated moment reads from any source that serves its inputs. Saved results are read through
OutputVault, which provides a source for the runs of one economy.

```julia
simseries, errs = simulate(obj, stv)
irf = simulate_irf(obj, simseries)
results = EconomyResult(; stats = describe(simseries).moments, irf, params = collect_parameters(obj))

getmoment(results, moments.std_K)
getmoment(results, moments.irf_C)
getmoment(results, moments.beta)
```

A missing sample, variable, statistic, shock or parameter is a `KeyError` that names it, a period
outside the IRF is a `BoundsError`, and a moment that needs a part the `EconomyResult` lacks (an
`ImpulseMoment` without an IRF, say) is an `ArgumentError`, as is reading a moment from a source that
cannot serve its type.

## Writing a source

A source is any object that can serve moments. [`getmoment`](@ref) and [`compare_economies`](@ref)
read through one function, `momentvalues(source, moment)`, which returns the pair `(value, basis)`:
the moment's value before its transformation, and the basis for the transformation (`nothing` for its
default, which is what every moment but an `ImpulseMoment` has). A source needs one method per kind of
moment it can serve, and a calculated moment then works on it with no further code, since it reads
its inputs from the same source:

```julia
struct MyEconomy
    values::Dict{Symbol, Float64}
end

TransitionFunctionIteration.momentvalues(e::MyEconomy, m::SimulationMoment) =
    (e.values[Symbol(m.stat, "_", m.var)], nothing)     # a missing key is a KeyError that names it

getmoment(MyEconomy(Dict(:mean_C => 1.0, :std_K => 0.2)), SimulationMoment(:K; stat = :std))
```

A moment type with no method for a source is an `ArgumentError`. OutputVault's saved results are a
source in this sense: one `Dict` per economy, read through `momentvalues`. Both its row layout and
its DuckDB queries are described in OutputVault's guide on querying with DuckDB.

## Comparing economies

[`compare_economies`](@ref) tabulates moments (rows) across economies (columns) in a `MomentArray`
whose row key is `moment`:

```julia
economies = (baseline = EconomyResult(; stats = stats_base, irf = irf_base),
             alt      = EconomyResult(; stats = stats_alt,  irf = irf_alt))
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
all moments are `SimulationMoment`s. `getmoment(economies, moment, economy, base)` reads one cell.
Any source with the `momentvalues` methods works as an economy, so results saved to disk give the
same numbers as the same economies in memory.
