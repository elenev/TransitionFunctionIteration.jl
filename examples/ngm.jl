# Neoclassical Growth Model 
using TransitionFunctionIteration, StaticArrays

# Explicitly import functions that must be overloaded
import TransitionFunctionIteration: Exogenv, calc_eqm, define_endog_grids, steadystate, compute_next_period_values, vals_from_policy, bounds

# Import some other functions for convenience
import TransitionFunctionIteration: crra, cobb_douglas, cobb_douglas_prime


# For debugging:
# import TransitionFunctionIteration: StructArrays, ForwardDiff, NonlinearSolve

# Define model 
# can be a more nested, complicated struct if needed
# but if it's more complicated it's probably worth it to overload update_model()
@kwdef struct NeoclassicalGrowthModel <: AbstractModel
    α::Float64 = 0.33   # capital share (implicit labor endowment of 1)
    β::Float64 = 0.96   # discount factor
    δ::Float64 = 0.10   # depreciation rate
    γ::Float64 = 2.0    # risk aversion (CRRA coefficient)
    σ::Float64 = 0.02   # productivity shock volatility
    ρ::Float64 = 0.9    # productivity shock persistence
    N_A::Int   = 5      # number of points in discretized productivity grid
end

## Model 

## Model specific functions

# [Optional] Only needed if steadystate() is not overloaded.
vals_from_policy(_, policy, ::NeoclassicalGrowthModel) = (C = policy.C,)

steadystate_guess(::NeoclassicalGrowthModel) = (; C = 1.0, K = 3.0)

# Analytic steady state solution for initial guess
function steadystate(model::NeoclassicalGrowthModel)
    (; α, β, δ) = model
    A = 1.0
    K = (α * β * A / (1 - β * (1 - δ)))^(1/(1-α))
    C = A * K^α - δ * K

    exog = (; A,)
    values = (; C,)
    guess = (; C,)
    endog = (; K)

    return SteadyState(; model, exog, endog, values, guess)
end

# If solver is guessing quantities whose domain is constrained, it's best to specify
# the constrained domain. Then the solver will operate in an unconstrained space and transform guesses and interpolant evaluations to the constrained space as needed. 
function bounds(::SteadyState{NeoclassicalGrowthModel})
    policy_bounds = (C = TFI.PositiveBounded(),)
    endog_bounds = (K = TFI.PositiveBounded(),)
    return policy_bounds, endog_bounds
end
    

# Grids for endogenous variables
function define_endog_grids(model::NeoclassicalGrowthModel, ::AbstractSteadyState)
    # Maximum attainable K is given by the solution to the equation 
    # A * K^α - δ * K = 0, which defines the point where all output is invested and there is no consumption. 
    # We can use this to set a reasonable range for the capital grid around the steady state value.
    # Set A = 1. Then,
    # (; α, δ) = model    
    # Kmax = (δ)^(1/(α-1))
    K_grid = range(0.8, 1.2, length=20) * steadystate(model).endog.K
    return (; K = K_grid)
end

# Create exogenous environment
function Exogenv(mdl::NeoclassicalGrowthModel)
    # Discretize AR(1) process for log productivity
    return ExogenvAR1(ρ=mdl.ρ, σ=mdl.σ, n=mdl.N_A, transform=exp, jensen=true, shock_name=:A)
end

# Compute next period values
function compute_next_period_values(next, mdl::NeoclassicalGrowthModel)
    (; A, K, C) = next
    (; α, β, δ, γ) = mdl
    MPK = A * cobb_douglas_prime(K, 1.0, α, 1)  # Marginal Product of Capital at next period's state
    MU = crra.u′(C, γ)  # Marginal Utility at next period's state
    Kpayoff = MPK + (1 - δ)  # Next period's payoff from capital
    expectations = (; Ek = β * MU * Kpayoff,)
    continuations = (; Kpayoff, MPK, MU)
    return NextPeriod(; expectations, continuations)
end

# Compute equilibrium
function calc_eqm(state, policy, mdl::NeoclassicalGrowthModel, np::NextPeriod, mode::EvaluationMode=SolveMode())
    (; A, K) = state
    (; C) = policy
    (; α, δ, γ) = mdl

    # Budget constraint 
    Y = A * cobb_douglas(K, 1.0, α)  # Output
    I = Y - C  # Investment
    K′ =  (1 - δ) * K  + I # Next period's capital from current state and policy
    
    # SDF
    MU = crra.u′(C, γ)
    
    # Euler equation residual: 1 = β E[u′(C′) / u′(C) * (MPK′ + 1 - δ)] 
    #fx = (euler = 1 - np.expectations.Ek / MU,)
    fx = SVector(1 - np.expectations.Ek / MU,)

    if mode isa SolveMode
        return fx
    else
        trans = (K = K′,)
        vals = (; C,)

        if mode isa IterateMode
            return (; fx, trans, vals)
        elseif mode isa SimulateMode
            addvars = (; Y, I, C)
            return (; fx, trans, addvars)
        end

    end
end

# Build initial solution object
ngm = NeoclassicalGrowthModel()
obj, stv = main_setup(ngm)
main_execute!(obj, maxiter=1000)
simseries, errs = simulate(obj, stv)
stats, _ = describe(simseries)
errstats = calc_error_stats(errs)
