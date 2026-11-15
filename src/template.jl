"""
    generate_template(typename, filename)

Write a Julia source file to `filename` containing shells for all interface functions
that a model package must (or may) implement to use TransitionFunctionIteration.jl.

`typename` is the name of the model struct (e.g., `"NeoclassicalGrowthModel"`).
Required functions (`Exogenv`, `define_endog_grids`, `steadystate`, `compute_next_period_values`,
`calc_eqm`) are uncommented; optional overloads (bounds, `vals_from_policy`, a steady-state system,
normalization, `augment_simulation`, `forecasts_from_policy`) are wrapped in block comments.
"""
function generate_template(typename::Union{String,Symbol}, filename::String)
    T = string(typename)

    content = """
# $T model
using TransitionFunctionIteration, StaticArrays, StructArrays

# Import all functions that will be overloaded for $T
import TransitionFunctionIteration: Exogenv, calc_eqm, define_endog_grids,
    steadystate, compute_next_period_values

# ── Model struct ─────────────────────────────────────────────────────────────

@kwdef struct $T <: AbstractModel
    # Add model parameters here
end

# ── Required interface functions ──────────────────────────────────────────────

# Construct the discretized exogenous Markov process for the model.
function Exogenv(mdl::$T)
    return ExogenvAR1(; ρ=mdl.ρ, σ=mdl.σ, n=mdl.N_Z, transform=exp, jensen=true, shock_name=:Z)
end

# Return a NamedTuple of grid vectors for each endogenous state variable.
function define_endog_grids(mdl::$T, stv::AbstractSteadyState)
    return (;)
end

# Compute the analytical steady state; return (SteadyState, converged::Bool).
# Replace with steadystate_system + SystemSteadyStateMode for a numerical approach.
function steadystate(mdl::$T, ::AnalyticalSteadyStateMode)
    exog   = (;)   # NamedTuple of exogenous steady-state values
    endog  = (;)   # NamedTuple of endogenous state values
    guess  = (;)   # NamedTuple of policies and prices that are solutions of the system of eqm conditions (Pfct)
    values = (;)   # NamedTuple of variables, whose next period values affect eqm conditions (Vfct)
    return SteadyState(; model=mdl, exog, endog, guess, values), true
end

# Compute next-period auxiliary quantities for each exogenous state realization.
# expectations: probability-weighted scalars (averaged by compute_expectations).
# continuations: per-state raw values (stacked into an SVector across exog states).
function compute_next_period_values(next, mdl::$T, ::EvaluationMode)
    expectations  = (;) # NamedTuple of variables, whose expectations enter equilibrium conditions
    continuations = (;) # NamedTuple of variables, whose per-state values are used to update transitions or in simulations
    return NextPeriod(; expectations, continuations)
end

# Compute equilibrium conditions at a grid point.
# fx names become the error variable names in simulation output -- use descriptive names.
function calc_eqm(state, policy, mdl::$T, np::NextPeriod, mode::EvaluationMode=SolveMode())
    # Unpack state and policy
    # (; var1, var2) = state
    # (; pol1, pol2) = policy

    # --- Equilibrium conditions ---
    # fx is a NamedTuple; field names label the residuals in error output
    fx = (;)

    if mode isa SolveMode
        return fx
    else
        trans   = (;)   # next-period endogenous state (scalars or vectors across exog states)
        vals    = (;)   # next-period policy/price variables (for Vfct update)

        if mode isa IterateMode
            return (; fx, trans, vals)
        elseif mode isa SimulateMode
            simvars = (;)   # all variables to record in simulation output
            return (; fx, trans, simvars)
        end
    end
end

# ── Optional overloads (uncomment and implement as needed) ────────────────────

#=
# Constrained domains for policy and endogenous-state variables (AnalyticalSteadyStateMode).
# Unbounded() is the default for all variables; override for positive/bounded quantities.
import TransitionFunctionIteration: bounds
function bounds(::SteadyState{$T})
    policy_bounds = (;)
    endog_bounds  = (;)
    return policy_bounds, endog_bounds
end
=#

#=
# Same as above but for FixedPointSteadyStateMode, which passes explicit policy/transition NTs.
import TransitionFunctionIteration: bounds
function bounds(::$T, policies::NamedTuple, transitions::NamedTuple)
    policy_bounds     = (;)
    transition_bounds = (;)
    return policy_bounds, transition_bounds
end
=#

#=
# Map policy variables to value-function variables when they differ.
# Only needed when Vfct tracks variables not directly in the policy vector.
import TransitionFunctionIteration: vals_from_policy
function vals_from_policy(state, policy, mdl::$T)
    return (;)
end
=#

#=
# Alternative steady-state system for SystemSteadyStateMode.
# Useful when the SS can be characterized by a smaller system than the full dynamic model.
# SolveMode: return (fx, nothing). IterateMode: return (fx, SteadyState(...)).
import TransitionFunctionIteration: steadystate_system
function steadystate_system(g::NamedTuple, mdl::$T, mode::EvaluationMode)
    fx = (;)
    if mode isa SolveMode
        return fx, nothing
    elseif mode isa IterateMode
        exog   = (;)
        guess  = (;)
        endog  = (;)
        values = vals_from_policy(nothing, guess, mdl)
        stv    = SteadyState(mdl, guess, endog, exog, values)
        return fx, stv
    end
end
=#

#=
# Post-expectation normalization of next-period values (e.g., Epstein-Zin CE adjustment).
# Called automatically after compute_expectations; default is identity.
import TransitionFunctionIteration: normalize_next_period_values
function normalize_next_period_values(np, mdl::$T, mode::EvaluationMode)
    return np
end
=#

#=
# Normalize residuals for simulation error output (e.g., divide by a level variable).
# Default returns fx_nt unchanged.
import TransitionFunctionIteration: normalize_fx_errors
function normalize_fx_errors(::$T, fx_nt, simvars)
    return fx_nt
end
=#

#=
# Append derived variables to the simulation panel (e.g., growth rates, excess returns).
# Return a NamedTuple of extra vectors (length T: a variable of each period; length T-1: of
# consecutive periods); each becomes a column in the simulation output.
import TransitionFunctionIteration: augment_simulation
function augment_simulation(sim::StructVector, mdl::$T)
    return (;)
end
=#

#=
# Adjust the steady state before using it as the initial interpolant guess.
# Useful when the deterministic SS differs meaningfully from the stochastic one.
import TransitionFunctionIteration: initial_guess_from_stv
function initial_guess_from_stv(stv::SteadyState, mdl::$T)
    return stv
end
=#

#=
# Forecasting values: current-period quantities, beyond the next-period state, that next period's
# expectations depend on (e.g. a portfolio entering next period's consumption stage). They are
# stored in a forecasting function Ffct and appear in `next` in compute_next_period_values.
# A model with forecasting values is solved with STPFIMethod:
#     main_execute!(obj; solmethod = STPFIMethod())
import TransitionFunctionIteration: forecasts_from_policy
function forecasts_from_policy(state, policy, mdl::$T)
    return (;)   # names must not clash with state or Vfct names
end
=#

# ── Usage ─────────────────────────────────────────────────────────────────────

# mdl = $T()
# obj, stv = main_setup(mdl)
# main_execute!(obj)        # tol = -5 means a distance of 1e-5; STPFIMethod for forecasting values
# simseries, errs = simulate(obj, stv)
# stats = describe(simseries).moments.uncond
"""

    write(filename, content)
    println("Template written to $filename")
    return nothing
end
