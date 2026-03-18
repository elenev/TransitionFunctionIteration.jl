@testset "Interpolant — constructors and accessors" begin

    @test varnames(fct1) == (:y, :z)
    @test nof(fct1)      == 2
    @test npt(fct1)      == 5
    @test statespacegrid(fct1) === ssg1

    # coefficients are reshaped to match grid dims; each element holds nof values
    @test size(coefficients(fct1)) == (5,)
    @test size(coefficients(fct2)) == (3, 3)

    # elsize: shape of one coefficient element (nof-length vector for multi-output)
    @test elsize(fct1) == (2,)
    @test elsize(fct2) == (1,)

    # knots: tuple of univariate grid vectors
    @test first(knots(fct1)) ≈ xgrid

    # name_values: wraps a raw value into a NamedTuple keyed by fct's output names
    coef1 = first(coefficients(fct1))   # SVector or similar with nof elements
    nt = name_values(coef1, fct1)
    @test nt isa NamedTuple
    @test keys(nt) == (:y, :z)

    # NamedTuple constructor: fills all grid points with constant values
    fct_const = Interpolant(ssg1, (y=0.0, z=0.0))
    @test varnames(fct_const) == (:y, :z)
    @test nof(fct_const)      == 2
    @test fct_const((; x=3.0)).y == 0.0

end

@testset "Interpolant — evaluation" begin

    # Exact evaluation at a grid point
    @test fct1((; x=3.0)).y ≈ 6.0
    @test fct1((; x=3.0)).z ≈ 9.0

    # Off-grid: Gridded(Linear()) → y(2.5) = (4+6)/2 = 5
    @test fct1((; x=2.5)).y ≈ 5.0

    # Keyword syntax
    @test fct1(; x=3.0).y ≈ 6.0

    # Extrapolation: Flat() clamps to boundary → y(10.0) = y(5.0) = 10.0
    @test fct1((; x=10.0)).y ≈ 10.0

    # Multi-point: generator input → StructArray output
    results = fct1(points(ssg1))
    @test length(results) == 5
    @test results[3].y ≈ 6.0

    # 2D exact at a grid point
    @test fct2((; x=2.0, y=20.0)).f ≈ 22.0

    # 2D interpolated: bilinear → exact for affine functions
    @test fct2((; x=1.5, y=15.0)).f ≈ 16.5

    # Exogenous state (NoInterp): integer dispatch, no cross-state interpolation
    @test fct_ex((; exst=1, x=2.0)).y ≈ 2.0    # 1 * 2 = 2
    @test fct_ex((; exst=2, x=2.5)).y ≈ 5.0    # 2 * 2.5 = 5

    # Singleton fields are stripped before evaluation
    ssg_s   = StateSpaceGrid(; x=collect(1.0:5.0), s=5.0)
    fvals_s = StructArray([(y=2.0*xi,) for xi in 1.0:5.0])
    fct_s   = Interpolant(ssg_s, fvals_s)
    @test fct_s((; x=3.0, s=5.0)).y ≈ 6.0

end

@testset "Interpolant — equality and isapprox" begin

    @test fct1 == fct1
    @test fct2 == fct2
    @test fct1 != fct2         # different grids and output names

    @test isapprox(fct1, fct1)
    @test !isapprox(fct1, fct2)

    # Perturbed coefficients: still approx-equal within loose tolerance,
    # but not equal and not approx-equal at default tolerance
    fvals1_perturbed = StructArray([(y=2.0*xi + 1e-8, z=xi^2) for xi in xgrid])
    fct1_p = Interpolant(ssg1, fvals1_perturbed)
    @test fct1_p != fct1
    @test isapprox(fct1_p, fct1; atol=1e-6)
    @test !isapprox(fct1_p, fct1; atol=1e-10)

end
