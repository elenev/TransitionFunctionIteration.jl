@testset "Bounded — construction and equality" begin

    @test TFI.Unbounded()           == TFI.Unbounded()
    @test TFI.LowerBounded(0.0)     == TFI.LowerBounded(0.0)
    @test TFI.LowerBounded(0.0)     != TFI.LowerBounded(1.0)
    @test TFI.UpperBounded(1.0)     == TFI.UpperBounded(1.0)
    @test TFI.UpperBounded(1.0)     != TFI.UpperBounded(2.0)
    @test TFI.IntervalBounded(0.0, 1.0) == TFI.IntervalBounded(0.0, 1.0)
    @test TFI.IntervalBounded(0.0, 1.0) != TFI.IntervalBounded(0.0, 2.0)

    # convenience constructors
    @test TFI.PositiveBounded() == TFI.LowerBounded(0.0)
    @test TFI.NegativeBounded() == TFI.UpperBounded(0.0)
    @test TFI.IntervalBounded() == TFI.IntervalBounded(0.0, 1.0)

end

@testset "Bounded — to_constrained / to_unconstrained round-trips" begin

    # For each bound type, verify to_constrained(to_unconstrained(y, b), b) ≈ y
    # over a set of representative values within the valid range.

    for y in [2.5, 0.0, -1.0]
        b = TFI.Unbounded()
        @test TFI.to_constrained(TFI.to_unconstrained(y, b), b) ≈ y
    end

    for y in [1.5, 0.6, 10.0]     # y > lower = 0.5
        b = TFI.LowerBounded(0.5)
        @test TFI.to_constrained(TFI.to_unconstrained(y, b), b) ≈ y
    end

    for y in [1.0, 2.5, -5.0]     # y < upper = 3.0  (avoid 0.0: isapprox(≈1e-16, 0) needs atol)
        b = TFI.UpperBounded(3.0)
        @test TFI.to_constrained(TFI.to_unconstrained(y, b), b) ≈ y
    end

    for y in [0.1, 0.4, 0.9]      # y ∈ (0, 1)
        b = TFI.IntervalBounded(0.0, 1.0)
        @test TFI.to_constrained(TFI.to_unconstrained(y, b), b) ≈ y
    end

end

@testset "Bounded — correct images at x=0" begin

    # to_constrained(0, b) has an analytically known value for each bound type

    # Unbounded: identity → 0
    @test TFI.to_constrained(0.0, TFI.Unbounded()) == 0.0

    # LowerBounded(lower): exp(0) + lower = 1 + lower
    @test TFI.to_constrained(0.0, TFI.LowerBounded(0.0)) ≈ 1.0
    @test TFI.to_constrained(0.0, TFI.LowerBounded(2.0)) ≈ 3.0

    # UpperBounded(upper): -exp(0) + upper = upper - 1
    @test TFI.to_constrained(0.0, TFI.UpperBounded(1.0)) ≈ 0.0
    @test TFI.to_constrained(0.0, TFI.UpperBounded(5.0)) ≈ 4.0

    # IntervalBounded(lo, hi): logistic(0) maps 0 to midpoint (lo+hi)/2
    @test TFI.to_constrained(0.0, TFI.IntervalBounded(0.0, 1.0)) ≈ 0.5
    @test TFI.to_constrained(0.0, TFI.IntervalBounded(2.0, 4.0)) ≈ 3.0

end

@testset "Bounded — monotonicity" begin

    # to_constrained should be strictly monotone in x for bounded types
    xs = [-2.0, -1.0, 0.0, 1.0, 2.0]

    # LowerBounded: increasing in x
    b_lo = TFI.LowerBounded(0.0)
    ys_lo = TFI.to_constrained.(xs, b_lo)
    @test issorted(ys_lo)
    @test all(>=(0.0), ys_lo)     # output always above lower bound

    # UpperBounded: decreasing in x (exp flipped)
    b_up = TFI.UpperBounded(5.0)
    ys_up = TFI.to_constrained.(xs, b_up)
    @test issorted(reverse(ys_up))
    @test all(<=(5.0), ys_up)     # output always below upper bound

    # IntervalBounded: increasing in x, output ∈ (lo, hi)
    b_iv = TFI.IntervalBounded(0.0, 1.0)
    ys_iv = TFI.to_constrained.(xs, b_iv)
    @test issorted(ys_iv)
    @test all(x -> 0.0 < x < 1.0, ys_iv)

end
