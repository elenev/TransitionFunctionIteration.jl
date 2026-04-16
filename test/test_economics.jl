@testset "Utility functions — CRRA and log" begin

    # γ=1 limit: crra dispatches to log utility
    @test TFI.crra.u(1.0, 1.0)  ≈ log(1.0)
    @test TFI.crra.u(ℯ, 1.0)   ≈ 1.0
    @test TFI.crra.u(2.0, 1.0)  ≈ log(2.0)

    # γ≠1: standard CRRA formula (c^(1-γ)-1)/(1-γ)
    @test TFI.crra.u(1.0, 2.0)  ≈ (1.0^(1-2) - 1) / (1-2)   # = 0.0
    @test TFI.crra.u(2.0, 2.0)  ≈ (2.0^(-1)  - 1) / (-1)     # = 0.5

    # marginal utility: c^(-γ) for all γ
    @test TFI.crra.u′(2.0, 1.0) ≈ 2.0^(-1.0)
    @test TFI.crra.u′(2.0, 2.0) ≈ 2.0^(-2.0)
    @test TFI.crra.u′(3.0, 0.5) ≈ 3.0^(-0.5)

    # logutility constants
    @test TFI.logutility.u(ℯ)   ≈ 1.0
    @test TFI.logutility.u(1.0) ≈ 0.0
    @test TFI.logutility.u′(2.0) ≈ 0.5
    @test TFI.logutility.u′(1.0) ≈ 1.0

    # powerutility at γ=1 matches logutility
    @test TFI.powerutility.u(2.0, 1.0)  ≈ TFI.logutility.u(2.0)
    @test TFI.powerutility.u′(2.0, 1.0) ≈ TFI.logutility.u′(2.0)

    # crra_strict skips the γ=1 check; powerutility_strict likewise
    @test TFI.crra_strict.u(2.0, 2.0)         ≈ TFI.crra.u(2.0, 2.0)
    @test TFI.powerutility_strict.u(2.0, 2.0) ≈ TFI.powerutility.u(2.0, 2.0)

end

@testset "CES — two-input scalar" begin

    # Unit inputs with equal weights always return 1
    @test TFI.ces(1.0, 1.0, 0.5, 2.0) ≈ 1.0
    @test TFI.ces(1.0, 1.0, 0.5, 0.5) ≈ 1.0

    # ε=1 limit: CES collapses to Cobb-Douglas
    for (x1, x2, w) in [(2.0, 3.0, 0.4), (1.5, 2.5, 0.7)]
        @test TFI.ces(x1, x2, w, 1.0) ≈ TFI.cobb_douglas(x1, x2, w)
    end

    # Exact formula check: ces_strict(2,3,0.5,2) = (0.5*2^0.5 + 0.5*3^0.5)^2
    a = (2.0 - 1) / 2.0   # a = 0.5
    expected = (0.5 * 2.0^a + 0.5 * 3.0^a)^(1/a)
    @test TFI.ces_strict(2.0, 3.0, 0.5, 2.0) ≈ expected

    # Symmetry: swapping inputs with swapped weights returns same value
    @test TFI.ces(2.0, 3.0, 0.4, 2.0) ≈ TFI.ces(3.0, 2.0, 0.6, 2.0)

end

@testset "CES — vector inputs" begin

    x  = [2.0, 3.0]
    w  = [0.5, 0.5]
    xs = SVector(2.0, 3.0)
    ws = SVector(0.5, 0.5)

    # ε=1: vector CES equals Cobb-Douglas
    @test TFI.ces(x,  w,  1.0) ≈ TFI.cobb_douglas(x,  w)
    @test TFI.ces(xs, ws, 1.0) ≈ TFI.cobb_douglas(xs, ws)

    # AbstractVector and SVector give same result
    @test TFI.ces(xs, ws, 2.0) ≈ TFI.ces(x, w, 2.0)

    # Unit inputs → 1 regardless of ε
    @test TFI.ces([1.0, 1.0], [0.5, 0.5], 3.0) ≈ 1.0

end

@testset "CES — marginal products" begin

    x = [2.0, 3.0]; w = [0.4, 0.6]; ε = 1.5
    δ = 1e-6

    for i in 1:2
        # Finite-difference derivative vs analytical
        x_hi = copy(x); x_hi[i] += δ
        fd_prime = (TFI.ces(x_hi, w, ε) - TFI.ces(x, w, ε)) / δ
        @test TFI.ces_prime(x, w, ε, i) ≈ fd_prime  rtol=1e-5
    end

    # Version with Y pre-computed gives same result
    Y = TFI.ces(x, w, ε)
    @test TFI.ces_prime(x, w, ε, Y, 1) ≈ TFI.ces_prime(x, w, ε, 1)

end

@testset "Cobb-Douglas — scalar, SVector, AbstractVector" begin

    # Two-input scalar
    @test TFI.cobb_douglas(2.0, 3.0, 0.5) ≈ sqrt(6.0)
    @test TFI.cobb_douglas(4.0, 1.0, 0.5) ≈ 2.0

    # SVector and AbstractVector should agree
    xs = SVector(2.0, 3.0); ws = SVector(0.5, 0.5)
    xv = [2.0, 3.0];        wv = [0.5, 0.5]
    @test TFI.cobb_douglas(xs, ws) ≈ sqrt(6.0)
    @test TFI.cobb_douglas(xv, wv) ≈ sqrt(6.0)
    @test TFI.cobb_douglas(xs, ws) ≈ TFI.cobb_douglas(xv, wv)

    # Exponents summing to 1 (standard CRS): x₁^w₁ * x₂^w₂
    @test TFI.cobb_douglas(SVector(3.0, 2.0), SVector(0.3, 0.7)) ≈ 3.0^0.3 * 2.0^0.7

end

@testset "Cobb-Douglas — marginal products" begin

    x1, x2, w = 2.0, 3.0, 0.5
    Y = TFI.cobb_douglas(x1, x2, w)
    δ = 1e-7

    # Analytical MP₁ = w₁ * Y / x₁
    @test TFI.cobb_douglas_prime(x1, x2, w, 1) ≈ w * Y / x1   rtol=1e-10
    @test TFI.cobb_douglas_prime(x1, x2, w, 2) ≈ (1-w) * Y / x2  rtol=1e-10

    # Version with Y pre-computed agrees
    @test TFI.cobb_douglas_prime(x1, x2, w, Y, 1) ≈ TFI.cobb_douglas_prime(x1, x2, w, 1)

    # Finite-difference consistency for vector version
    xv = [2.0, 3.0, 1.5]; wv = [0.3, 0.5, 0.2]
    for i in 1:3
        xv_hi = copy(xv); xv_hi[i] += δ
        fd = (TFI.cobb_douglas(xv_hi, wv) - TFI.cobb_douglas(xv, wv)) / δ
        @test TFI.cobb_douglas_prime(xv, wv, i) ≈ fd  rtol=1e-5
    end

    # Invalid index → ArgumentError
    @test_throws ArgumentError TFI.cobb_douglas_prime(x1, x2, w, 3)

end

@testset "Fisher-Burmeister complementarity function" begin

    # f(a,b) = a + b - sqrt(a²+b²)
    @test TFI.fisher_burmeister(3.0, 4.0) ≈ 3.0 + 4.0 - 5.0  # known 3-4-5 triangle
    @test TFI.fisher_burmeister(0.0, 0.0) ≈ 0.0
    @test TFI.fisher_burmeister(1.0, 0.0) ≈ 1.0 + 0.0 - 1.0   # = 0
    @test TFI.fisher_burmeister(0.0, 1.0) ≈ 0.0 + 1.0 - 1.0   # = 0
    # Complementarity: if a≥0, b≥0, a*b=0, result=0
    @test abs(TFI.fisher_burmeister(2.0, 0.0)) < 1e-14
    @test abs(TFI.fisher_burmeister(0.0, 5.0)) < 1e-14

end

@testset "Epstein-Zin certainty equivalent" begin

    x = [1.0, 2.0, 3.0, 4.0]
    p = fill(0.25, 4)  # uniform weights

    # Scalar passthrough: any single value returns itself unchanged
    @test TFI.epstein_zin_ce(2.5, p, 2.0) == 2.5
    @test TFI.epstein_zin_ce(1.0, p, 1.0) == 1.0

    # γ≈1: E-Z CE reduces to geometric mean exp(E[log x])
    geom_mean = exp(sum(p .* log.(x)))
    @test TFI.epstein_zin_ce(x, p, 1.0) ≈ geom_mean

    # γ≠1: power mean — verify via the E-function overload
    # E[f,x] overload with explicit function: E(f,x) = sum(p .* f.(x))
    E_fn(f, x) = sum(p .* f.(x))
    γ = 2.0
    expected = E_fn(xi -> xi^(1-γ), x)^(1/(1-γ))
    @test TFI.epstein_zin_ce(x, E_fn, γ) ≈ expected

    # γ≈1 via E-function overload
    geom_E = exp(E_fn(log, x))
    @test TFI.epstein_zin_ce(x, E_fn, 1.0) ≈ geom_E

    # Both overloads agree for γ≠1
    @test TFI.epstein_zin_ce(x, p, γ) ≈ TFI.epstein_zin_ce(x, E_fn, γ)

end

@testset "Epstein-Zin Jensen correction" begin

    # γ=1: no uncertainty correction regardless of other parameters
    @test TFI.epstein_zin_jensen_correction(γ=1.0, σ=0.05) == 0.0
    @test TFI.epstein_zin_jensen_correction(γ=1.0, ψ=2.0, σ=0.05) == 0.0

    # ψ≈1: multiplier preferences -- analytical formula, result is negative (uncertainty reduces welfare)
    corr_mult = TFI.epstein_zin_jensen_correction(β=0.99, γ=2.0, ψ=1.0, σ=0.02, ρ=0.9)
    @test corr_mult isa Float64
    @test corr_mult != 0.0
    # Correction is negative: uncertainty lowers expected utility relative to deterministic SS
    @test corr_mult < 0.0

    # General EZ (ψ≠1): nonlinear solve -- result is also non-zero and finite
    corr_ez = TFI.epstein_zin_jensen_correction(β=0.99, γ=2.0, ψ=2.0, σ=0.02, ρ=0.9)
    @test isfinite(corr_ez)
    @test corr_ez != 0.0

    # σ=0: no uncertainty → correction should be zero (or very close) for both cases
    corr_zero = TFI.epstein_zin_jensen_correction(β=0.99, γ=2.0, ψ=1.0, σ=0.0)
    @test corr_zero ≈ 0.0  atol=1e-10

end

@testset "Asset returns — vector overloads" begin

    # return_on_ex_price: R_t = (p_{t+1} + d_{t+1}) / p_t * G_{t+1} - 1
    p  = [1.0, 1.1, 1.2]
    d  = [0.0, 0.05, 0.1]
    G  = [1.0, 1.0, 1.0]

    rets_ex = TFI.return_on_ex_price(p, d, G)
    @test length(rets_ex) == 2
    @test rets_ex[1] ≈ TFI.return_on_ex_price(p[1], p[2], d[2], G[2])
    @test rets_ex[2] ≈ TFI.return_on_ex_price(p[2], p[3], d[3], G[3])

    # return_on_cum_price: R_t = p_{t+1} / (p_t - d_t) * G_{t+1} - 1
    rets_cum = TFI.return_on_cum_price(p, d, G)
    @test length(rets_cum) == 2
    @test rets_cum[1] ≈ TFI.return_on_cum_price(p[1], p[2], d[1], G[2])
    @test rets_cum[2] ≈ TFI.return_on_cum_price(p[2], p[3], d[2], G[3])

    # Default d=0 and G=1 vector overloads
    rets_ex_def = TFI.return_on_ex_price(p)
    @test rets_ex_def[1] ≈ p[2]/p[1] - 1

    rets_cum_def = TFI.return_on_cum_price(p)
    @test rets_cum_def[1] ≈ p[2]/p[1] - 1

end

@testset "Growth rate — scalar and vector" begin

    # Scalar: growth_rate = return_on_ex_price with d=0
    @test TFI.growth_rate(1.0, 1.1) ≈ 0.1
    @test TFI.growth_rate(2.0, 3.0) ≈ 0.5
    @test TFI.growth_rate(1.0, 1.0, 1.0) ≈ 0.0

    # Vector overload: consecutive growth rates
    x = [1.0, 1.1, 1.21]
    g = TFI.growth_rate(x)
    @test length(g) == 2
    @test g[1] ≈ 0.1  rtol=1e-10
    @test g[2] ≈ 0.1  rtol=1e-10

    # With G (renormalization) vector
    G = [1.0, 2.0, 4.0]
    g2 = TFI.growth_rate(x, G)
    @test g2[1] ≈ TFI.growth_rate(x[1], x[2], G[2])

end

@testset "CES marginal product — single-input overload" begin

    # ces_prime(xi, w, ε, Y) = w * (Y/xi)^(1/ε)
    xi, w, ε, Y = 2.0, 0.5, 1.5, 4.0
    expected = w * (Y / xi)^(1 / ε)
    @test TFI.ces_prime(xi, w, ε, Y) ≈ expected

    # Consistent with vector version
    x = [2.0, 3.0]; wv = [0.5, 0.5]
    Y_vec = TFI.ces(x, wv, ε)
    @test TFI.ces_prime(x[1], wv[1], ε, Y_vec) ≈ TFI.ces_prime(x, wv, ε, 1)

end

@testset "Economics — ForwardDiff compatibility" begin

    # Verify that ForwardDiff.gradient/derivative runs without error.
    # This tests that all functions are analytically differentiable and
    # that no non-differentiable branches are hit.

    # Cobb-Douglas gradient w.r.t. inputs
    g_cd = ForwardDiff.gradient(x -> TFI.cobb_douglas(x[1], x[2], 0.5), [2.0, 3.0])
    @test g_cd ≈ [TFI.cobb_douglas_prime(2.0, 3.0, 0.5, 1),
                  TFI.cobb_douglas_prime(2.0, 3.0, 0.5, 2)]  rtol=1e-8

    # CES (scalar) derivative w.r.t. x₁
    d_ces = ForwardDiff.derivative(x -> TFI.ces(x, 3.0, 0.5, 1.5), 2.0)
    @test d_ces ≈ TFI.ces_prime([2.0, 3.0], [0.5, 0.5], 1.5, 1)  rtol=1e-8

    # CRRA marginal utility derivative matches -γ c^(-γ-1)
    d_mu = ForwardDiff.derivative(c -> TFI.crra.u′(c, 2.0), 2.0)
    @test d_mu ≈ -2.0 * 2.0^(-3.0)  rtol=1e-8

end
