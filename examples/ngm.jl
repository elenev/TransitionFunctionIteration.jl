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
    policy_bounds = (C = TransitionFunctionIteration.PositiveBounded(),)
    endog_bounds = (K = TransitionFunctionIteration.PositiveBounded(),)
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

# Compare with linearized solution
using LinearAlgebra

"""
Solves the linear expectations model:

    y' = W y + M z + Π η
    z' = Φ z + ϵ

where y = [x; c] is the stacked vector of state variables and controls,
z is the vector of exogenous states, and ϵ is the vector of shocks. 
Π = [0; I] and η = c - E[c] captures the forecast errors, which are linear in ϵ.
The solution is of the form:

    y = Ψ_x x + Ψ_z z

where Ψ_x and Ψ_z are the policy function matrices for the endogenous 
and exogenous states, respectively.

Uses Klein (2000) method to solve for the policy function matrices Ψ_x and Ψ_z.
"""
function klein(W::AbstractMatrix, M::AbstractVecOrMat, Φ::Union{Real,AbstractMatrix}, nx::Integer)
    # 1. Real Schur decomposition of the transition matrix
    F_schur = schur(W)
    
    # 2. Order the Schur decomposition (stable roots top-left)
    select = abs.(F_schur.values) .< 1.0
    F_ordered = ordschur(F_schur, select)
    
    T = F_ordered.T
    Z = F_ordered.Z  # Orthogonal matrix (W = Z * T * Z')
    
    # 3. Partition the orthogonal matrix Z
    Z11 = Z[1:nx, 1:nx]
    Z12 = Z[1:nx, nx+1:end]
    Z21 = Z[nx+1:end, 1:nx]
    Z22 = Z[nx+1:end, nx+1:end]
    
    # 4. Solve for Psi_x (Policy function on endogenous states)
    # The / operator in Julia efficiently solves X * Z11 = Z21
    Psi_x = Z21 / Z11  
    
    # 5. Transform the shock mapping matrix M into the Schur basis
    # M_rot = Z' * M, and extract the explosive block (bottom rows)
    M_rot = Z' * M
    M2 = M_rot[nx+1:end, :]
    
    # 6. Solve for the explosive block mapping F using the Sylvester equation
    # Equation: -T22 * F + F * Phi = M2
    T22 = T[nx+1:end, nx+1:end]
    if length(T22) == 1 && length(Φ) == 1
        # If T22 and Φ are scalars, solve directly
        T22 = T22[1]
        F = -M2 / (T22 - Φ)
    else
        # Otherwise, solve the Sylvester equation
        F = sylvester(-T22, Φ, -M2)
    end
    
    # 7. Recover Psi_z (Policy function on exogenous states)
    Psi_z = (Z22 - Psi_x * Z12) * F
    
    return Psi_x, Psi_z
end

function klein(; XX, XC, CX, CC, XZ, CZ, Φ)
    W = vcat( hcat(XX, XC), hcat(CX, CC) )
    M = vcat( XZ, CZ )
    nx = size(XX, 1)
    # if Φ isa Number
    #     Φ = Φ * I(1)
    # end
    return klein(W, M, Φ, nx)
end

function solve_linearized(model::NeoclassicalGrowthModel)
    # Linearize the model around the steady state
    # For the NGM, the state variables are K and A, and the control variable is C.
    # The exogenous state is A, which follows an AR(1) process.
    
    # Steady state values
    stv = steadystate(model)
    Kbar = stv.endog.K
    Cbar = stv.values.C
    Abar = 1.0

    # Define some notation for the linearized system
    (; α, β, γ, ρ) = model
    η = Abar * Kbar^(α-1)  # Marginal product of capital at steady state
    κ = α * β * η / γ
    
    # Define coefficients
    KK = 1/β;          KC = -Cbar/Kbar;              KA = η;
    CK = κ*(α-1)/β;    CC = 1 - κ*(α-1)*Cbar/Kbar;   CA = κ*((α-1)*η + ρ);

    # Solve for policy function matrices
    Ck, Ca = klein(XX=KK, XC=KC, CX=CK, CC=CC, XZ=KA, CZ=CA, Φ=ρ)

    Ck = length(Ck) == 1 ? Ck[1] : Ck
    Ca = length(Ca) == 1 ? Ca[1] : Ca

    Kk = KK + KC * Ck
    Ka = KA + KC * Ca

    # log C - log Cbar = Ck * (log K - log Kbar) + Ca * (log A - log Abar)
    # log K - log Kbar = Kk * (log K - log Kbar) + Ka * (log A - log Abar)

    vcoefs = TransitionFunctionIteration.MomentArray([0.0 Ca Ck], (:C,), (:constant, :A, :K))
    tcoefs = TransitionFunctionIteration.MomentArray([0.0 Ka Kk 0.0], (:K′,), (:constant, :A, :K, :A′))

    return vcoefs, tcoefs #, pt -> (C = Cbar * (pt.K/Kbar)^Ck * (pt.A/Abar)^Ca, 
                          #        K = Kbar * (pt.K/Kbar)^Kk * (pt.A/Abar)^Ka)

end

Vlin, Tlin = solve_linearized(ngm)
Vsim, Tsim = linearize(simseries, stv)

# Exclude constant -- even for a linear model, the stochastic steady state will be a bit different than deterministic
# so we don't expect the constant to be zero
isapprox(Vlin[:,2:end], Vsim[:,2:end], atol=1e-3)
isapprox(Tlin[:,2:end], Tsim[:,2:end], atol=1e-3)
