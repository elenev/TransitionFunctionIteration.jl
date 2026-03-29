# Neoclassical Growth Model 
using TransitionFunctionIteration, StaticArrays, StructArrays

# Explicitly import functions that must be overloaded
import TransitionFunctionIteration: Exogenv, calc_eqm, define_endog_grids, steadystate_system, 
    compute_next_period_values, vals_from_policy, bounds, 
    normalize_next_period_values, normalize_fx_errors

# Import some other functions for convenience
import TransitionFunctionIteration: crra, cobb_douglas, cobb_douglas_prime, ces, epstein_zin_ce,
    PositiveBounded, Unbounded

@kwdef struct IntermediaryAssetPricingModel <: AbstractModel
    α::Float64 = 0.33   # capital share (implicit labor endowment of 1)
    β::Float64 = 0.99   # discount factor
    δ::Float64 = 0.10   # depreciation rate
    γ::Float64 = 2.0    # risk aversion 
    ψ::Float64 = 0.5    # elasticity of intertemporal substitution
    gbar::Float64 = 0.02   # aggregate productivity growth rate
    σ::Float64 = 0.02   # productivity shock volatility
    ρ::Float64 = 0.65   # productivity shock persistence
    η::Float64 = 0.07   # payout target
    ξ::Float64 = 0.5    # equity issuance cost
    ϕ::Float64 = 2.5    # capital adjustment cost
    θ::Float64 = 0.9    # leverage constraint
    N_G::Int   = 5      # number of points in discretized productivity grid
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

function Exogenv(model::IntermediaryAssetPricingModel)
    (; ρ, σ, gbar) = model
    return ExogenvAR1(; σ=σ, ρ=ρ, μ=gbar, jensen=true, n=model.N_G, transform=exp, shock_name=:G)
end

function bounds(::IntermediaryAssetPricingModel, policies::NamedTuple, transitions::NamedTuple) 
    policy_bounds = (C = PositiveBounded(), q = PositiveBounded(), x = Unbounded(), μ = PositiveBounded(), p = PositiveBounded(), λ = Unbounded())
    transition_bounds = (K = PositiveBounded(), L = PositiveBounded())
    return policy_bounds, transition_bounds
end

function steadystate_system(g::NamedTuple, mdl::IntermediaryAssetPricingModel, mode::EvaluationMode)
    (; α, β, δ, gbar, ψ, ξ, η, θ) = mdl
    (; x) = g 

    G = exp(gbar)
    p = β * exp(-gbar / ψ)

    # Step 1: 𝓋 = v^F/μ as functions of x
    𝓋 = (η - x)*(1 - ξ*x) + (1 - η - ξ/2 * x^2)

    # Step 2: normalized leverage multiplier from bond Euler
    λ = p*(1 - 𝓋)

    # Step 3: gross return on capital from capital Euler
    RK = (1 - λ*θ) / (p * 𝓋)

    # Step 4: capital-to-net-worth ratio from net worth evolution (n=1)
    k = G / (RK - θ/p)

    # Residual: budget constraint (1-θ)k = 1 - η + x - ξ/2 x²
    resid = (1 - θ)*k - (1 - η + x - ξ/2 * x^2)

    # Remaining steady-state values given x
    q = 1.0
    rK = RK - (1 - δ)
    K = ( rK / α )^(1/(α-1))
    K′ = G*K
    N = K′ / k

    w = (1-α) * K^α
    C = (η - x)*N + (1/p - G)*θ*K + w

    L = N / (q*K)

    fx = (resid, )
    if mode isa SolveMode
        return fx, nothing
    elseif mode isa IterateMode
        policy_sol = (; C, x, q, p, λ) 
        transition_sol = (; K, L)
        exog = (; G)
        vals = vals_from_policy(nothing, policy_sol, mdl)
        stv = SteadyState(mdl, policy_sol, transition_sol, exog, vals)
        return fx, stv
    end
end

function vals_from_policy(state, policy, model::IntermediaryAssetPricingModel)
    (; C, q, x) = policy
    (; ξ, η, β, ψ, gbar) = model
    μ = 1 / (1 - ξ*x)
    Q = η - x + μ * (1-η - ξ/2 * x^2)
    
    if ψ ≈ 1.0
        statadj = exp( gbar * β/(1-β) )
    else
        SDF_hat = β * exp(gbar*(1-1/ψ))
        statadj = ( (1-β)/(1-SDF_hat) )^(1/(1-1/ψ))
    end

    V = C * statadj

    return (; C, q, Q, V)
end

function define_endog_grids(model::IntermediaryAssetPricingModel, stv::SteadyState)
    K_grid = range(0.5, 2.5, length=15) * stv.endog.K
    L_grid = range(0.10, 0.125, length=9)
    return (; K=K_grid, L=L_grid)
end

function calc_production(vars, model::IntermediaryAssetPricingModel)
    (; K, q) = vars 
    (; δ, gbar, α, ϕ) = model
    δ_tilde = δ + exp(gbar) - 1
    #ℓ = G / K

    Y = cobb_douglas(K, 1.0, α)
    MPK = cobb_douglas_prime(K, α, Y)
    
    #  q = 1 +  ϕ (i - δ)
    i = δ_tilde + (q - 1) / ϕ
    rK = MPK + (i*(q-1) - ϕ/2 * (i - δ_tilde)^2)
    RK = rK + 1 - δ
    return (; RK, Y, MPK, i, rK, δ_tilde)
end

function compute_next_period_values(next, model::IntermediaryAssetPricingModel, mode::EvaluationMode)
    (; γ, β, ψ) = model
    (; C, G, Q, K, V) = next
    # Household

    # The stochastic part of Epstein Zin SDF
    MU = crra.u′(C, 1/ψ)
    Λ = β * G^(1-γ) * MU * V^(1/ψ - γ)   
    BHpayoff = 1/G
    SHpayoff = Q
    
    # Firm
    (; q, K, Q) = next
    (; RK) = calc_production((K=K, q=q), model)

    KIpayoff = RK / G
    BIpayoff = 1  / G

    #

    # Expectations
    ΛF = Λ * Q
    expectations = (
                    KI = ΛF * KIpayoff, 
                    BI = ΛF * BIpayoff, 
                    BH = Λ  * BHpayoff,
                    SH = Λ  * SHpayoff,

                    CE = (G*V)^(1-γ)
                    )

    # Q and C added to continuations only in SimulateMode -- not needed in solve/iterate hot loops
    if mode isa SimulateMode
        continuations = (; RK, G, q, V, Q, C)
    else
        continuations = (; RK, G, q, V)
    end

    return NextPeriod(; expectations, continuations)
end

function normalize_next_period_values(np, mdl::IntermediaryAssetPricingModel, mode::EvaluationMode)
    (; γ, ψ) = mdl
    #@debug "np.continuations.V: $(np.continuations.V)"
    #@debug "exp(g) np.continuations.V: $( exp(mdl.gbar) * np.continuations.V )"
    CE = np.expectations.CE^(1/(1-γ))
    #@debug "CE: $CE"
    CEterm = CE^(γ - 1/ψ)
    #@debug "CEterm: $CEterm"

    expectations = (
                    KI = np.expectations.KI * CEterm, 
                    BI = np.expectations.BI * CEterm, 
                    BH = np.expectations.BH * CEterm,
                    SH = np.expectations.SH * CEterm,

                    CE = CE
                    )
    return NextPeriod(expectations = expectations, continuations = np.continuations,
            prnext = np.prnext)
end

function calc_eqm(state, policy, model::IntermediaryAssetPricingModel, np::NextPeriod, mode::EvaluationMode=SolveMode())
    (; K, L) = state
    (; C, x, q, p, λ) = policy
    (; α, δ, η, ξ, gbar, ϕ, ψ, β, θ) = model

    (; RK, Y, i) = calc_production((K=K, q=q), model)

    d = η - x
    μ = 1 / (1 - ξ*x)

    Q = d + μ * (1 - η - ξ/2 * x^2)

    N = L*q*K
    B = RK*K - N

    D = d * N

    w = cobb_douglas_prime(K, 1.0, α, Y, 2)
    BH = ( B + D + w - C ) / p
    
    K′ = (1 - δ + i) * K
    BI = (q*K′ - ((1-η) + x - ξ/2 * x^2) * N) / p

    MU = crra.u′(C, 1/ψ)
    
    Λ = MU 
    ΛF = Λ * μ

    fx = (EE_BH = p - np.expectations.BH / Λ,
          EE_KI = q - np.expectations.KI / ΛF   - λ*θ,
          EE_BI = p - np.expectations.BI / ΛF   - λ,
          MC_B  = BH - BI,
          LC    = fisher_burmeister(λ, θ*q*K′ - p*BI))

    if mode isa SolveMode
        return fx
    else
        N′ = @. K′*np.continuations.RK - BI
        trans = (K = K′ ./ np.continuations.G, 
                 L = N′ ./ (q.*K′))
        V = ces(C, np.expectations.CE, 1-β, ψ)
        vals = (; C, q, Q, V)
        
        if mode isa IterateMode
            return (; fx, trans, vals)
        elseif mode isa SimulateMode
            # np.continuations.Q: equity multiplier from next period (stacked per exog state)
            E_Q = mean(np.continuations.Q, np)

            r = 1/p - 1
            E_rS  = E_Q / (Q - d) - 1
            E_xrS = E_rS - r
            simvars = (; Y, RK, i, I = i*K, N, BH, BI, D, w, d, x,
                bind = λ > 0, r, E_rS, E_xrS, E_Q)
            simvals = merge(vals, simvars)

            return (; fx, trans, vals = simvals)
        end
    end

end

function normalize_fx_errors(::IntermediaryAssetPricingModel, fx_nt, addvars)
    return (; EE_BH = fx_nt.EE_BH / addvars.p, 
              EE_KI = fx_nt.EE_KI / addvars.q, 
              EE_BI = fx_nt.EE_BI / addvars.p, 
              MC_B  = fx_nt.MC_B  / addvars.BH, 
              LC    = fx_nt.LC    / addvars.BH)
end

mdl = IntermediaryAssetPricingModel()

# Configuration
xguess = 0.018
setup_opts = (mode = SystemSteadyStateMode((x=xguess,)), check_steadystate=true)
exec_opts  = (; maxiter=1000, logging=COMPLETE)
irf_opts   = (; path_vars=(:C, :K),)
out = solve_and_simulate(mdl; setup_kwargs=setup_opts, execute_kwargs=exec_opts, irf_kwargs=irf_opts)


# Steady state fixed point 
# _, stvguess = steadystate_system((x=xguess,), mdl, IterateMode())
# configfp = FixedPointSteadyStateMode(stvguess.guess, stvguess.endog, stvguess.exog)
# obj, stv = main_setup(mdl, mode=configfp)