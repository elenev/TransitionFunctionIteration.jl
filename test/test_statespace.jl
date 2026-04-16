@testset "StateSpaceGrid" begin

    # ── accessors ─────────────────────────────────────────────────────────────

    @test varnames(ssg1) == (:x,)
    @test dims(ssg1)     == (x=5,)
    @test npt(ssg1)      == 5
    @test unigrids(ssg1).x == xgrid

    @test varnames(ssg2) == (:x, :y)
    @test dims(ssg2)     == (x=3, y=3)
    @test npt(ssg2)      == 9

    # ── points() ──────────────────────────────────────────────────────────────
    # Cartesian product: x varies fastest (first arg to Iterators.product)

    pts2 = collect(points(ssg2))
    @test length(pts2) == 9
    @test pts2[begin] == (; x=1.0, y=10.0)
    @test pts2[2]     == (; x=2.0, y=10.0)
    @test pts2[end]   == (; x=3.0, y=30.0)

    # ── singletons ────────────────────────────────────────────────────────────
    # Scalar-valued fields are metadata, not iterable axes

    @test singletons(ssg1) == ()
    ssg_s = StateSpaceGrid(; x=collect(1.0:3.0), s=5.0)
    @test singletons(ssg_s) == (:s,)

    # ── equality and isapprox ─────────────────────────────────────────────────

    @test ssg1 == ssg1
    @test ssg1 != ssg2
    @test isapprox(ssg1, ssg1)
    @test !isapprox(ssg1, ssg2)

    # ── StateSpaceGrid(exogenv, endog_nt) constructor ─────────────────────────
    # Combines exst integer range with endogenous axes

    exog_pts = StructArray([(z=v,) for v in [0.9, 1.0, 1.1]])
    mtrans   = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    exv      = Exogenv(pts=exog_pts, mtrans=mtrans, ssidx=2)
    ssg_from_exog = StateSpaceGrid(exv, (K=collect(1.0:4.0),))
    @test varnames(ssg_from_exog) == (:exst, :K)
    @test dims(ssg_from_exog)     == (exst=3, K=4)

end
