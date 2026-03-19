# Neoclassical Growth Model 
using TransitionFunctionIteration, StaticArrays

# Explicitly import functions that must be overloaded
import TransitionFunctionIteration: Exogenv, calc_eqm, define_endog_grids, steadystate, compute_next_period_values, vals_from_policy, bounds

# Import some other functions for convenience
import TransitionFunctionIteration: crra, cobb_douglas, cobb_douglas_prime

@kwdef struct IntermediaryAssetPricingModel <: AbstractModel
    α::Float64 = 0.33   # capital share (implicit labor endowment of 1)
    β::Float64 = 0.96   # discount factor
    δ::Float64 = 0.10   # depreciation rate
    γ::Float64 = 2.0    # risk aversion (CRRA coefficient)
    σ::Float64 = 0.02   # productivity shock volatility
    ρ::Float64 = 0.9    # productivity shock persistence
    π::Float64 = 0.07   # payout target
    ξ::Float64 = 0.5    # equity issuance cost
    ϕ::Float64 = 0.5    # capital adjustment cost
    θ::Float64 = 0.8    # leverage constraint
    N_A::Int   = 5      # number of points in discretized productivity grid
end

# Equilibrium conditions

# Household
# Stochastic Discount Factor: Mₜ₊₁ ≡ β * (Cₜ₊₁/Cₜ)^(-γ)
# Euler equation for bonds: 1 =  Eₜ[Mₜ₊₁ * (1 + rₜ)]
# Budget constraint: Cₜ + qₜ Bₜ₊₁ = wₜ + Bₜ + Dₜ

# Intermediary
# Stochastic Discount Factor: MIₜ₊₁ ≡ 
# Investment optimality: qₜ = ϕ * (Iₜ/Kₜ - δ) + 1
# Euler equation for bonds: pₜ = λₜ + Eₜ[MIₜ₊₁]
# Euler equation for capital: qₜ = λₜ θ + Eₜ[MIₜ₊₁ * (α * Aₜ₊₁ * Kₜ₊₁^(α-1) + (1 - δ) * qₜ₊₁)]
# Capital accumulation: Kₜ₊₁ = (1 - δ) Kₜ + Iₜ
# Net worth: Nₜ =  Aₜ * Kₜ^α + (1 - δ) Kₜ - Bₜ
# Objective: Vₜ = π Nₜ - Xₜ - ξ/2 * Xₜ² + E[MIₜ₊₁ * Vₜ₊₁] where Xₜ is equity issuance and Vₜ is the value function
# Investment adjustment cost: Φₜ = ϕ/2 * (Iₜ/Kₜ - δ)² * K
# Budget constraint: (1-π) Nₜ + Xₜ + pₜ Bₜ₊₁ = qₜ Kₜ₊₁
