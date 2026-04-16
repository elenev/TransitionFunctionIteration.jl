# Neoclassical Growth Model 
using TransitionFunctionIteration, StaticArrays, StructArrays

# Explicitly import functions that must be overloaded
import TransitionFunctionIteration: Exogenv, calc_eqm, define_endog_grids, steadystate_system, 
    initial_guess_from_stv, compute_next_period_values, vals_from_policy, bounds, 
    normalize_next_period_values, normalize_fx_errors, augment_simulation

# Import some other functions for convenience
import TransitionFunctionIteration: crra, cobb_douglas, cobb_douglas_prime, ces, epstein_zin_ce,
    PositiveBounded, Unbounded, growth_rate, return_on_cum_price, return_on_ex_price

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

# function initial_guess_from_stv(stv::SteadyState, model::IntermediaryAssetPricingModel)
#     (; β, γ, ψ, gbar, σ, ρ) = model
#     μ = gbar
#     ρ = 0.5
#     σ = σ * sqrt(1 - ρ^2)
#     vadj = TransitionFunctionIteration.epstein_zin_jensen_correction(; β, γ, ψ, μ, σ, ρ)
#     V = stv.values.V * exp(vadj)
#     values = merge(stv.values, (; V))
#     return SteadyState(; model=stv.model, guess=stv.guess, 
#         endog=stv.endog, exog=stv.exog, values)
# end

function define_endog_grids(model::IntermediaryAssetPricingModel, stv::SteadyState)
    K_grid = range(1.8, 9.0, 20) #range(0.5, 2.5, length=15) * stv.endog.K
    L_grid = range(0.09, 0.12, length=15)
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
                    CE = (G*V)^(1-γ),
                    q  = q,
                    )

    continuations = (; RK, G, q, V)

    # Additional next-period quantities for simulation output.
    if mode isa SimulateMode
        QG = Q * G
        expectations  = merge(expectations, (; QG))
        continuations = merge(continuations, (; C))
    end

    return NextPeriod(; expectations, continuations)
end

function normalize_next_period_values(np, mdl::IntermediaryAssetPricingModel, mode::EvaluationMode)
    (; γ, ψ) = mdl
    CE     = np.expectations.CE^(1/(1-γ))
    CEterm = CE^(γ - 1/ψ)

    base_exp = (
        KI = np.expectations.KI * CEterm,
        BI = np.expectations.BI * CEterm,
        BH = np.expectations.BH * CEterm,
        SH = np.expectations.SH * CEterm,
        CE = CE,
        q  = np.expectations.q,
    )
    # QG is the physical-measure E[Q*G]; it is not SDF-weighted so CEterm does not apply
    expectations = mode isa SimulateMode ?
        merge(base_exp, (; QG = np.expectations.QG, )) : base_exp

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
          LC    = fisher_burmeister(λ, θ*q*K′ - p*BI)) #np.expectations.

    if mode isa SolveMode
        return fx
    else
        N′ = @. K′*np.continuations.RK - BI
        trans = (K = K′ ./ np.continuations.G, 
                 L = N′ ./ (np.continuations.q.*K′))
        V = ces(C, np.expectations.CE, 1-β, ψ)
        vals = (; C, q, Q, V)
        
        if mode isa IterateMode
            return (; fx, trans, vals)
        elseif mode isa SimulateMode
            r     = 1/p - 1
            E_rS  = np.expectations.QG / (Q - d) - 1
            E_xrS = E_rS - r
            E_rK   = RK - 1
            simvars = (; Y, RK, i, I = i*K, N, BH, BI, D, w, d, x,
                bind = λ > 0, r, E_rS, E_xrS)
            simvars = merge(vals, simvars)

            return (; fx, trans, simvars)
        end
    end

end

function normalize_fx_errors(::IntermediaryAssetPricingModel, fx_nt, simvars)
    return (; EE_BH = fx_nt.EE_BH / simvars.p, 
              EE_KI = fx_nt.EE_KI / simvars.q, 
              EE_BI = fx_nt.EE_BI / simvars.p, 
              MC_B  = fx_nt.MC_B  / simvars.BH, 
              LC    = fx_nt.LC    / simvars.BH)
end


function augment_simulation(sim::StructVector, mdl::IntermediaryAssetPricingModel)
    growth_rate_vars = (:Y, :K, :C, :D, :I)
    growth_rates = TransitionFunctionIteration.apply_function_to_variables((v,s) -> growth_rate(v,s.G), 
        growth_rate_vars, sim, v -> Symbol("$(v)_gr"))

    rS = return_on_cum_price(sim.Q, sim.d, sim.G)
    r = @view sim.r[2:end]
    xrS = rS .- r

    (; MPK) = StructArray(calc_production.(sim, Ref(mdl)))
    rK = return_on_ex_price(sim.K, MPK, sim.G)
    xrK = rK .- r

    return (; growth_rates..., rS, xrS, rK, xrK)
end

mdl = IntermediaryAssetPricingModel()

# Configuration
xguess = 0.018
setup_opts = (mode = SystemSteadyStateMode((x=xguess,)), check_steadystate=true)
exec_opts  = (; maxiter=1000, logging=COMPLETE)
stat_opts  = (; filters = (; neggrowth = x -> x.Y_gr < 0))
irf_opts   = (; path_vars=(:C, :K),)
out = solve_and_simulate(mdl; setup_kwargs=setup_opts, execute_kwargs=exec_opts, irf_kwargs=irf_opts, describe_kwargs=stat_opts)

# mdl2 = IntermediaryAssetPricingModel(ϕ=5.0)
# out2 = solve_and_simulate(mdl2; setup_kwargs=setup_opts, execute_kwargs=exec_opts, irf_kwargs=irf_opts, describe_kwargs=stat_opts)

# Steady state fixed point 
# _, stvguess = steadystate_system((x=xguess,), mdl, IterateMode())
# configfp = FixedPointSteadyStateMode(stvguess.guess, stvguess.endog, stvguess.exog)
# obj, stv = main_setup(mdl, mode=configfp)