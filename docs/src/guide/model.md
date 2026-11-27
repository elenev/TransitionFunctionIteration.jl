# Implementing a Model

Define your model's parameters in a subtype of `AbstractModel`, then write functions for its
equilibrium equations, exogenous process, steady state and grids. Start with a template:

```julia
generate_template("MyModel", "/path/to/model_file.jl")
```

Fill in the required functions and any optional calculations your model needs. Import each
function you extend with `import TransitionFunctionIteration: ...`, as in the [Quick Start](@ref).

## Required interface

| Function | Signature | Purpose |
|---|---|---|
| `Exogenv` | `Exogenv(mdl::MyModel)` | Build the exogenous Markov process |
| `define_endog_grids` | `define_endog_grids(mdl, stv)` | Return grid vectors for endogenous states |
| steady state: one of three options | `steadystate(mdl, mode)`, `steadystate_system(guess, mdl, mode)`, or none with `FixedPointSteadyStateMode` | Compute the deterministic steady state; see [Steady-state computation](@ref) |
| `compute_next_period_values` | `compute_next_period_values(next, mdl, mode)` | Compute next-period expectations and continuations |
| `calc_eqm` | `calc_eqm(state, policy, mdl, np, mode)` | Evaluate equilibrium conditions |

### Exogenous environment

`Exogenv(mdl)` should return an `Exogenv` representing the model's Markov process.
For a single AR(1) process use [`ExogenvAR1`](@ref):

```julia
function Exogenv(mdl::MyModel)
    ExogenvAR1(ρ=mdl.ρ, σ=mdl.σ, n=mdl.N, transform=exp, jensen=true, shock_name=:A)
end
```

### Endogenous grids

Return one grid vector per endogenous state in a NamedTuple. The steady state `stv` is available
to scale the grid; for example, place capital between 80% and 120% of its steady-state value:

```julia
function define_endog_grids(mdl::MyModel, stv::AbstractSteadyState)
    K_grid = range(0.8, 1.2, length=40) * stv.endog.K
    return (; K = K_grid)
end
```

The first grid axis is always `:exst`, the integer exogenous state index. Do not include it in
`define_endog_grids`. Endogenous grids may be unevenly spaced; interpolation is linear with flat
extrapolation. An out-of-range exogenous index throws an error. Check simulated bounds and errors
before relying on results outside the endogenous grid.

### `calc_eqm` return contract

`calc_eqm` evaluates the equilibrium equations at a given state and policy, using the next-period
quantities in `np`. For the growth model, the state is productivity and capital, the policy is
consumption, and the equation is the Euler condition. The requested `mode` determines which
results to return:

| Mode | Return |
|---|---|
| `SolveMode` | `SVector` (or `NamedTuple`) of residuals `fx` |
| `IterateMode` | `(; fx, trans, vals)` |
| `SimulateMode` | `(; fx, trans, simvars)` |

- `trans`: `NamedTuple` of next-period endogenous state values. Each value is a scalar, which
  applies to every next-period exogenous state, or an `SVector` with one entry per exogenous state.
- `vals`: `NamedTuple` of variables for `Vfct` -- the quantities whose next-period values
  enter expectations in equilibrium conditions.
- `simvars`: `NamedTuple` of all output variables to record in simulation.

Compute the common quantities once, then return the results needed for solving, updating the
functions or recording a simulation:

```julia
function calc_eqm(state, policy, mdl::MyModel, np::NextPeriod, mode::EvaluationMode)
    # ... compute Y, I, K′, MU, fx ...
    if mode isa SolveMode
        return fx
    elseif mode isa IterateMode
        return (; fx, trans=(; K=K′), vals=(; C))
    elseif mode isa SimulateMode
        return (; fx, trans=(; K=K′), simvars=(; Y, I, C))
    end
end
```

### `compute_next_period_values` return contract

Use `compute_next_period_values(next, mdl, mode)` to calculate the terms needed in today's
equations from one next-period state. In the growth model, `next` contains productivity `A`,
capital `K` from the transition functions, and consumption `C` from the value functions.
Models with forecasting functions also receive their current-period forecasting values.

Return `NextPeriod(; expectations, continuations)`, with NamedTuples for two kinds of quantities:

- `expectations`: terms whose conditional expectation enters the equilibrium conditions. In
  `calc_eqm`, `np.expectations.Ek` is a scalar: the probability-weighted sum over next-period
  states.
- `continuations`: other next-period values that the equilibrium conditions need state by state.
  In `calc_eqm`, `np.continuations.MU` is a vector with one entry per next-period exogenous state.

Return values for this single state, normally scalars; the package calls the function for each
next-period exogenous state and combines the results. `np` then holds `expectations`,
`continuations` (which also contains the variables of `next`, by state) and `prnext`, the
probabilities of the next-period exogenous states.

```julia
function compute_next_period_values(next, mdl::MyModel, ::EvaluationMode)
    (; A, K, C) = next
    (; α, β, δ, γ) = mdl
    MPK = α * A * K^(α-1)
    MU  = C^(-γ)
    return NextPeriod(;
        expectations  = (; Ek = β * MU * (MPK + 1 - δ)),
        continuations = (; MPK, MU),
    )
end
```

## Steady-state computation

Choose how to find the deterministic steady state with the `mode` keyword of `main_setup`.
Use an analytical solution, solve a reduced system of equations, or solve the dynamic equations
with the additional condition that endogenous states remain unchanged.

### `AnalyticalSteadyStateMode` (default)

Implement `steadystate(mdl, ::AnalyticalSteadyStateMode)` returning `(SteadyState(...), true)`:

```julia
function steadystate(mdl::MyModel, ::AnalyticalSteadyStateMode)
    # ... closed-form solution ...
    return SteadyState(; model=mdl, exog=(; A=1.0), endog=(; K), values=(; C), guess=(; C)), true
end
```

The returned `SteadyState` supplies the starting values for the solution:

- `exog`: exogenous state at the steady state
- `endog`: endogenous state at the steady state (used to center grids)
- `values`: quantities stored in `Vfct` at initialization
- `guess`: initial solver guess for policy variables

### `SystemSteadyStateMode`

Use when the steady state satisfies a reduced system (smaller than the full dynamic system).
Implement `steadystate_system(guess, mdl, mode)` and pass the mode to `main_setup`:

```julia
obj, stv = main_setup(mdl; mode=SystemSteadyStateMode((; K=3.0, C=0.8)))
```

Import `steadystate_system` from TransitionFunctionIteration. Its `guess` argument is a NamedTuple
with the names passed to `SystemSteadyStateMode`. Return `(residuals, nothing)` in `SolveMode` and
`(residuals, stv)` in `IterateMode`, where `stv` is the complete steady state implied by the guess.
Both modes must accept dual numbers.

### `FixedPointSteadyStateMode`

Use when you cannot solve for the steady state analytically or via a reduced system. The solver
solves the dynamic equations jointly for policies and endogenous states, imposing that transitions
return the same state:

```julia
obj, stv = main_setup(mdl; mode=FixedPointSteadyStateMode(
    (C=0.8,),          # initial policy guess
    (K=3.0,),          # initial transition guess
    (A=1.0,),          # exogenous state at the steady state
))
```

Three options control the `Vfct` values and bounds of the fixed-point problem:

- `vals_from_policy(state, policy, mdl)` supplies the initial `Vfct` values; by default, the policy
  itself. Override it for other values that can be computed from the state and policy.
- Values that depend on expectations can instead be extra unknowns: pass
  `eqmvalguesses = (; V = 1.0)` to the mode, using names returned in `calc_eqm`'s `vals`. Their
  equations equate the guesses with those returned values.
- To bound the unknowns, overload `bounds(mdl, policies, transitions)`, or with equilibrium-value
  guesses `bounds(mdl, policies, transitions, eqmvalguesses)`. Return bounds for the policy,
  transition and equilibrium-value guesses, in that order.

## Optional interface

Add these methods when your model needs bounds, extra simulation variables, forecasting values
or custom error measures:

| Function | Purpose |
|---|---|
| `bounds(stv::SteadyState{MyModel})` | Variable bounds for solver transformation (default: all `Unbounded()`) |
| `augment_simulation(run, mdl)` | Append extra variables to a simulation run |
| `forecasts_from_policy(state, policy, mdl)` | Current-period quantities that next period's expectations depend on (a forecasting function `Ffct`) |
| `augment_errors(obj, state, out, trans)` | Add model-specific convergence diagnostics to the error output |
| `normalize_next_period_values(np, mdl, mode)` | Transform expectations before they enter `calc_eqm` |
| `normalize_fx_errors(mdl, fx_nt, simvars)` | Normalize residuals for comparability |
| `initial_guess_from_stv(stv, mdl)` | Adjust the steady-state values before initializing interpolants |
| `vals_from_policy(state, policy, mdl)` | `Vfct` values as functions of the policy, for the fixed-point steady state (default: the policy itself) |

### Variable bounds

Specify admissible values for policies and endogenous states with `bounds`. Here consumption and
capital must be positive, and `q` must lie between 0 and 1:

```julia
function bounds(::SteadyState{MyModel})
    policy_bounds = (C = PositiveBounded(), q = IntervalBounded(0.0, 1.0))
    endog_bounds  = (K = PositiveBounded(),)
    return policy_bounds, endog_bounds
end
```

Available bound types: `Unbounded`, `LowerBounded`, `UpperBounded`, `IntervalBounded`
(documented in the API reference; none is exported, so `import TransitionFunctionIteration:
PositiveBounded, IntervalBounded` or qualify them). Convenience constructors:
`PositiveBounded() = LowerBounded(0.0)`, `NegativeBounded() = UpperBounded(0.0)`.

Bounds are enforced through smooth transformations. The solver and stored `Pfct` coefficients use
unconstrained values; the model equations receive policies in their original units.

### Augmenting simulation output

Use `augment_simulation` for variables computed from several periods of a simulated path.
It receives a vector of one run's period-by-period NamedTuples and returns a NamedTuple of new
series. Here `R` divides current output by the previous period's capital:

```julia
function augment_simulation(run, mdl::MyModel)
    R = [run[t].Y / run[t-1].K for t in 2:length(run)]    # T-1 values
    return (; R)
end

dynamic_depth(::MyModel) = 1
```

`R` has one fewer observation than `run` because it needs a lag. Setting `dynamic_depth` to 1
asks `simulate` to include one extra period before the reported sample, calculate `R`, and then
drop that context period. It also sets the default `hold_periods` for `simulate_irf`.

More generally, a series of length `T - d` is aligned with the last `T - d` periods of a run of
length `T`. Its first `d` entries are padded with `NaN` before context periods are dropped.
Set `dynamic_depth` to the longest history any added variable needs, or choose `context_periods`
explicitly in `simulate`. Return `nothing` when there are no variables to add.

### Forecasting functions

Next period's expectations sometimes depend on current-period quantities besides next period's
state. For example, the stochastic discount factor `β (C′/C)^(-γ)` needs today's consumption `C`
as well as next period's `C′`. Declare such quantities with `forecasts_from_policy`, which returns a
`NamedTuple` computed from the state and the policy; each name then appears in `next` in
`compute_next_period_values`:

```julia
import TransitionFunctionIteration: forecasts_from_policy

forecasts_from_policy(state, policy, mdl::MyModel) = (; Cnow = policy.C)   # names must differ from state and Vfct names

function compute_next_period_values(next, mdl::MyModel, ::EvaluationMode)
    (; A, K, C, Cnow) = next
    # ... β * (C / Cnow)^(-γ) * (MPK + 1 - δ) ...
end
```

TFI stores these quantities in `Ffct` and updates them from the solved policy with the same damping
as the value and transition functions. Their convergence distance is `distF`; the iteration
continues until it, too, meets the tolerance.

### Choosing a solution method

Start with `TFIMethod()`, the default, which solves a smaller nonlinear system at each grid point.
Use `STPFIMethod()` when a policy variable affects the equilibrium equations only through
next-period states or forecasting values.

- `TFIMethod` computes expectations from the current transition and forecasting functions before
  solving for the policy. Those expectations stay fixed during the grid-point solve.
- `STPFIMethod` solves jointly for the policy and next-period endogenous states, recomputing
  forecasting values and expectations at each evaluation. There is one extra unknown for each
  endogenous state and possible next-period exogenous state.

```julia
main_execute!(obj; solmethod = STPFIMethod())
```

Both methods have the same fixed point. During a TFI grid-point solve, however, a policy that
appears only inside the fixed expectations cannot change the residuals. This makes the system
singular and can cause every grid point to fail, returning `status = :all_failed`.

For the consumption forecast above, the residual
`1 - E[β (C′/Cnow)^(-γ) (MPK′ + 1 - δ)]` depends on today's consumption only through `Cnow`.
Use STPFI to solve it in this form. Alternatively, write `1 - E[...] (Cnow/C)^(-γ)` so that the
current policy `C` enters the residual directly. At the fixed point `Cnow = C`, and the equation
is unchanged; TFI can solve this form.

### Extra error diagnostics

Add model-specific accuracy checks with `augment_errors(obj, state, out, trans)`. For example,
compare recomputed consumption with its interpolated value:

```julia
function augment_errors(obj, state, out, trans)
    C_gap = out.simvars.C - obj.Vfct(state).C   # recomputed vs. interpolated
    return (; C_gap), (;)
end
```

Return a pair of NamedTuples: scalar diagnostics are added to equation errors, and vectors with
one entry per next-period exogenous state are added to transition errors. A name matching an
existing equation error replaces it. The default returns two empty NamedTuples.

## Numerical requirements

TFI differentiates the model equations with ForwardDiff. Let intermediate calculations accept
dual numbers rather than requiring `Float64`. Use `SVector` and `SMatrix` for small collections in
the equations to keep grid-point solves efficient.

Grid points and simulation runs are evaluated on multiple threads. Model functions must be thread
safe: avoid shared mutable state and global random-number generators.

## Re-solving after a parameter change

`update_model` creates a new model with updated parameters, and `main_setup` accepts a
previous solution as a warm start. `update_model` calls the keyword constructor `MyModel(; ...)`
with the current field values and the new ones, so define the struct with `@kwdef`, as in the
[Quick Start](@ref), or write a keyword constructor:

```julia
mdl2   = update_model(mdl; σ=0.03)
obj2, stv2 = main_setup(mdl2, obj)   # warm-start from previous solution
main_execute!(obj2)
```
