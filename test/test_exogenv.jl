@testset "Exogenv — ExogenvAR1 constructor" begin

    # basic construction (iid, n=3)
    @test exnpt(e3) == 3
    @test nex(e3)   == 1
    @test ssidx(e3) == 2                          # middle point
    @test propertynames(points(e3)) == (:z,)      # default shock name
    @test size(e3.mtrans) == (3, 3)
    @test all(≈(1.0), sum(e3.mtrans, dims=2))     # rows sum to 1

    # custom shock name
    ea = ExogenvAR1(n=3, shock_name=:a)
    @test propertynames(points(ea)) == (:a,)

    # ssidx for larger n
    e5 = ExogenvAR1(n=5)
    @test ssidx(e5) == 3

    # shock_volatility=true produces wider state support than =false
    e_sv_false = ExogenvAR1(ρ=0.9, σ=1.0, n=5, shock_volatility=false)
    e_sv_true  = ExogenvAR1(ρ=0.9, σ=1.0, n=5, shock_volatility=true)
    @test maximum(abs, collect(points(e_sv_true).z)) >
          maximum(abs, collect(points(e_sv_false).z))

    # jensen=true shifts μ so ergodic mean ≈ exp(μ)=1
    e_log = ExogenvAR1(n=5, σ=0.2, μ=0.0, jensen=true)
    @test all(>(0), points(e_log).z)
    @test mean(e_log).z ≈ exp(0.0)  atol=0.05

end

@testset "Exogenv — product-space blending constructor" begin

    e1 = ExogenvAR1(n=3, shock_name=:z1)
    e2 = ExogenvAR1(n=5, shock_name=:z2)
    eb = Exogenv(e1, e2)

    @test exnpt(eb) == 15
    @test size(eb.mtrans) == (15, 15)
    @test all(≈(1.0), sum(eb.mtrans, dims=2))

    @test :z1 in propertynames(points(eb))
    @test :z2 in propertynames(points(eb))

    # ssidx: linear index of (ssidx(e1), ssidx(e2)) in vec(product(1:3, 1:5))
    combos       = vec(collect(Iterators.product(1:3, 1:5)))
    expected_idx = findfirst(==((ssidx(e1), ssidx(e2))), combos)
    @test ssidx(eb) == expected_idx

    # transition matrix is Kronecker product
    @test eb.mtrans ≈ kron(e2.mtrans, e1.mtrans)

    # duplicate field names → error
    @test_throws ErrorException Exogenv(e1, e1)

end

@testset "Exogenv — ergodic statistics" begin

    # iid symmetric chain: ergodic distribution should be uniform
    e = ExogenvAR1(ρ=0.0, σ=0.5, μ=0.0, n=3)

    pvec = TransitionFunctionIteration.ergodic_distro(e)
    @test length(pvec) == 3
    @test sum(pvec) ≈ 1.0
    @test all(>=(0.0), pvec)

    @test mean(e).z  ≈ 0.0  atol=1e-10
    @test var(e).z   >= 0.0
    @test std(e).z   ≈ sqrt(var(e).z)

    idx_mode, pt_mode = TransitionFunctionIteration.mode(e)
    @test 1 <= idx_mode <= exnpt(e)
    @test pt_mode == points(e)[idx_mode]

    idx_med, pt_med = median(e)
    @test idx_med == 2
    @test pt_med == points(e)[2]

    # even n → warns and returns (nothing, nothing)
    e_even = ExogenvAR1(n=4)
    @test_warn r"Even number" begin
        idx2, pt2 = median(e_even)
        @test isnothing(idx2)
        @test isnothing(pt2)
    end

end

@testset "Exogenv — equality and accessors" begin

    e5 = ExogenvAR1(n=5)

    @test e3 == e3
    @test e3 != e5
    @test isapprox(e3, e3)
    @test !isapprox(e3, e5)

    # transition_matrix returns the raw mtrans field
    @test transition_matrix(e3) === e3.mtrans

    # conditional_distro: each row should sum to 1
    for i in 1:exnpt(e3)
        @test sum(conditional_distro(e3, i)) ≈ 1.0
    end

end

@testset "Exogenv — ExogenvAR1 with one state" begin
    e1 = ExogenvAR1(ρ=0.724, σ=0.0136, μ=0.0335, n=1, shock_volatility=true, shock_name=:r)
    @test exnpt(e1) == 1
    @test ssidx(e1) == 1
    @test points(e1).r[1] == 0.0335
    @test e1.mtrans == [1.0;;]
    @test conditional_distro(e1, 1) == [1.0]
    # no Jensen correction: the single point keeps the transformed mean exp(μ)
    e1log = ExogenvAR1(σ=0.2, μ=0.1, n=1, jensen=true)
    @test points(e1log).z[1] ≈ exp(0.1)
end
