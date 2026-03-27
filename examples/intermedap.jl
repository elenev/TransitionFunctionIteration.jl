# Neoclassical Growth Model 
using TransitionFunctionIteration, StaticArrays, StructArrays

# Explicitly import functions that must be overloaded
import TransitionFunctionIteration: Exogenv, calc_eqm, define_endog_grids, steadystate_system, 
    compute_next_period_values, vals_from_policy, bounds

# Import some other functions for convenience
import TransitionFunctionIteration: crra, cobb_douglas, cobb_douglas_prime, PositiveBounded, Unbounded

@kwdef struct IntermediaryAssetPricingModel <: AbstractModel
    α::Float64 = 0.33   # capital share (implicit labor endowment of 1)
    β::Float64 = 0.99   # discount factor
    δ::Float64 = 0.10   # depreciation rate
    γ::Float64 = 2.0    # risk aversion (CRRA coefficient)
    gbar::Float64 = 0.02   # aggregate productivity growth rate
    σ::Float64 = 0.02   # productivity shock volatility
    ρ::Float64 = 0.65    # productivity shock persistence
    η::Float64 = 0.07   # payout target
    ξ::Float64 = 0.5    # equity issuance cost
    ϕ::Float64 = 0.5    # capital adjustment cost
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
    (; α, β, δ, gbar, γ, ξ, η, θ) = mdl
    (; x) = g 

    G = exp(gbar)
    p = β * exp(-gbar * γ)

    # Step 1: μ, v^F/μ as functions of x
    μ = 1 / (1 - ξ*x)
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
    (; ξ, η) = model
    μ = 1 / (1 - ξ*x)
    Q = η - x + μ * (1-η - ξ/2 * x^2)
    return (; C, q, Q)
end

function define_endog_grids(model::IntermediaryAssetPricingModel, stv::SteadyState)
    K_grid = range(0.7, 1.8, length=11) * stv.endog.K
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

function compute_next_period_values(next, model::IntermediaryAssetPricingModel)
    # Household
    (; γ, β) = model
    (; C, G, Q, K) = next
    MU = crra.u′(C, γ)    
    Λ = β * G^(1-γ) * MU
    BHpayoff = 1/G
    SHpayoff = Q
    
    # Firm
    (; q, K, Q) = next
    (; RK) = calc_production((K=K, q=q), model)

    KIpayoff = RK / G
    BIpayoff = 1 / G

    # Expectations
    ΛF = Λ * Q
    expectations = (KI = ΛF * KIpayoff, 
                    BI = ΛF * BIpayoff, 
                    BH = Λ  * BHpayoff,
                    SH = Λ  * SHpayoff)

    continuations = (RK = RK, G = G, q=q)

    #if mode isa SimulateMode
        addvars = (; RK, Q, C, q)
    #else
    #    addvars = nothing
    #end
    
    return NextPeriod(; expectations, continuations, addvars)
end

function calc_eqm(state, policy, model::IntermediaryAssetPricingModel, np::NextPeriod, mode::EvaluationMode=SolveMode())
    (; K, L) = state
    (; C, x, q, p, λ) = policy
    (; α, δ, η, ξ, gbar, ϕ, γ, θ) = model

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

    MU = crra.u′(C, γ)
    
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
        # @debug "K" K
        # @debug "K′" K′
        # @debug "trans.K" trans.K
        vals = (; C, q, Q)
        
        if mode isa IterateMode
            return (; fx, trans, vals)
        elseif mode isa SimulateMode
            evars = map(x -> E(x,np), components(np.addvars))
            enames = map(x -> Symbol("E_$x"), keys(components(np.addvars)))
            evars = NamedTuple{enames}(values(evars))


            addvars = (; Y, RK, i, I = i*K, N, BH, BI, D, w, d, x, Q, q,
                bind = λ > 0,
                r = 1/p - 1, E_rS = evars.E_Q / (Q - d) - 1)
            addvars = merge(addvars, evars)
            
            return (; fx, trans, addvars)
        end
    end

end

mdl = IntermediaryAssetPricingModel()

# Steady state system 
xguess = 0.018
configsys = SystemSteadyStateMode((x=0.018,))
obj, stv = main_setup(mdl;mode=configsys, check_steadystate=true)

main_execute!(obj; maxiter=1000)
simseries, errs = simseries, errs = simulate(obj, stv)
stats, _ = describe(simseries)
errstats = calc_error_stats(errs)

# Steady state fixed point 
# _, stvguess = steadystate_system((x=xguess,), mdl, IterateMode())
# configfp = FixedPointSteadyStateMode(stvguess.guess, stvguess.endog, stvguess.exog)
# obj, stv = main_setup(mdl, mode=configfp)

using Plots
function calc_at_guess(mdl, x=0.0056)
    resid, g,s,e=steadystate_guess(mdl,x)
    np = compute_next_period_values((q = 1.0, K=s.K, vals_from_policy(nothing,g,mdl).Q, C = g.C, G = e.G), mdl)
    out = calc_eqm((K = s.K, L=s.L, G = 1.02, exst=1), (C=g.C, x=g.x, μ=g.μ, q=1.0, p=g.p, λ=g.λ), mdl, np, IterateMode())
    stateerrs = (K = s.K - out.trans.K, L = s.L - out.trans.L)
    return merge(out.fx, stateerrs)
end


function plot_errs(errs::SimulationErrorPanel, sim::SimulationPanel;
                    errnames=propertynames(errs),
                    exogstates=sort(unique(sim.exst)),
                    endogstates=endog_states(sim)[1:2])

    N = length(exogstates)
    M = length(errnames)
    n_endog = length(endogstates)

    if exogstates == :all
        _get_exst = _ -> trues(length(sim.exst))
    else
        _get_exst = exst -> sim.exst .== exst
    end

    # Global error limits; for 2D, color by absolute error magnitude.
    # Use extrema(f, v) to find each panel's range in one pass without allocating.
    _efn = n_endog == 2 ? abs : identity
    err_lims = mapreduce(
        t -> extrema(_efn, getproperty(errs.data[_get_exst(t[1])], t[2])),
        (a, b) -> (min(a[1], b[1]), max(a[2], b[2])),
        Iterators.product(exogstates, errnames)
    )
    function _plot_err_1D(mask, errname, statename; kw...)
        scatter(getproperty(sim.data[mask], statename), getproperty(errs.data[mask], errname);
            ylims=err_lims, markerstrokewidth=0, markerstrokecolor=:transparent, kw...)
    end

    function _plot_err_2D(mask, errname, statename1, statename2; kw...)
        abs_errs = abs.(getproperty(errs.data[mask], errname))
        scatter(getproperty(sim.data[mask], statename1), getproperty(sim.data[mask], statename2);
            zcolor=abs_errs, colormap=:jet, clims=err_lims, colorbar=false,
            markerstrokewidth=0, markerstrokecolor=:transparent, kw...)
    end

    if n_endog == 1
        _plot = (mask, errname; kw...) -> _plot_err_1D(mask, errname, endogstates[1]; kw...)
    elseif n_endog == 2
        _plot = (mask, errname; kw...) -> _plot_err_2D(mask, errname, endogstates[1], endogstates[2]; kw...)
    else
        error("Plotting only implemented for 1 or 2 endog states")
    end

    # Build as a 1D vector in row-major order so Plots.jl places them correctly:
    # row r = exogstates[r], col c = errnames[c]; title on top row, xlabel on bottom, ylabel on left
    plots = [_plot(_get_exst(exogstates[r]), errnames[c];
                   title      = r == 1 ? string(errnames[c]) : "",
                   xlabel     = r == N ? string(endogstates[1]) : "",
                   ylabel     = c == 1 && n_endog == 2 ? string(endogstates[2]) : "",
                   xformatter = r == N ? :auto : _ -> "",
                   yformatter = c == 1 ? :auto : _ -> "",
                   colorbar   = false,
                   legend     = false)
             for r in 1:N for c in 1:M]

    min_width = 200
    min_height = 200
    total_size = (M * min_width, N * min_height)

    if n_endog == 2
        # Invisible scatter panel: no visible body, colorbar strip only
        cbar_plot = scatter([0.5], [0.5]; zcolor=[mean(err_lims)], clims=err_lims,
                            color=:jet, ms=0, markerstrokewidth=0, colorbar=true,
                            framestyle=:none, xlims=(0,1), ylims=(0,1), label="")
        l = @layout [grid(N, M) a{0.12w}]
        plot(plots..., cbar_plot, layout=l, size=(total_size[1] + 80, total_size[2]))
    else
        plot(plots..., layout=(N, M), size=total_size)
    end
end