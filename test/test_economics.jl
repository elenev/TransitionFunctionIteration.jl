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
