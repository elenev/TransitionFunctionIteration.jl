using Test
using TransitionFunctionIteration
using StructArrays, StaticArrays, RecipesBase, Statistics

# ─── shared fixtures ──────────────────────────────────────────────────────────

# 1D: y(x) = 2x, z(x) = x²
xgrid = collect(1.0:5.0)
ssg1  = StateSpaceGrid(; x=xgrid)
fvals1 = StructArray([(y=2.0*xi, z=xi^2) for xi in xgrid])
fct1  = Interpolant(ssg1, fvals1)

# 2D: f(x,y) = x + y  (bilinear, so linear interpolation is exact)
xgrid2 = collect(1.0:3.0)
ygrid2 = collect(10.0:10.0:30.0)
ssg2   = StateSpaceGrid(; x=xgrid2, y=ygrid2)
fvals2 = StructArray([(f=xi+yj,) for xi in xgrid2, yj in ygrid2])
fct2   = Interpolant(ssg2, fvals2)

# With exogenous state (NoInterp integer indexing) + endogenous (linear): f(exst,x) = exst*x
ssg_ex  = StateSpaceGrid(; exst=1:2, x=collect(1.0:3.0))
fvals_ex = StructArray([(y=Float64(ei)*xi,) for ei in 1:2, xi in 1.0:3.0])
fct_ex  = Interpolant(ssg_ex, fvals_ex)

# ─── tests ────────────────────────────────────────────────────────────────────

@testset "TransitionFunctionIteration" begin

@testset "StateSpaceGrid" begin
    # Accessors on 1D grid
    @test varnames(ssg1) == (:x,)
    @test dims(ssg1) == (x=5,)    # dims() returns map(length, unigrids(ssg)), i.e. a NamedTuple
    @test npt(ssg1) == 5
    @test unigrids(ssg1).x == xgrid

    # Accessors on 2D grid
    @test varnames(ssg2) == (:x, :y)
    @test dims(ssg2) == (x=3, y=3)
    @test npt(ssg2) == 9

    # points(): Cartesian product via Iterators.product(xgrid, ygrid)
    # x varies fastest (first argument to product), so second element has x=xgrid[2]
    pts2 = collect(points(ssg2))
    @test length(pts2) == 9
    @test pts2[begin] == (; x=1.0, y=10.0)   # first x, first y
    @test pts2[2]     == (; x=2.0, y=10.0)   # second x, still first y
    @test pts2[end]   == (; x=3.0, y=30.0)   # last x, last y

    # singletons(): fields whose type <: Number are singletons (scalar metadata, not iterable axes)
    @test singletons(ssg1) == ()
    ssg_s = StateSpaceGrid(; x=collect(1.0:3.0), s=5.0)
    @test singletons(ssg_s) == (:s,)
end

@testset "Interpolant — constructors and accessors" begin
    @test varnames(fct1) == (:y, :z)
    @test nof(fct1) == 2
    @test npt(fct1) == 5
    @test statespacegrid(fct1) === ssg1

    # coefficients are reshaped to match grid dims
    @test size(coefficients(fct1)) == (5,)
    @test size(coefficients(fct2)) == (3, 3)

    # knots(): first element is the first grid axis
    @test first(knots(fct1)) ≈ xgrid

    # NamedTuple constructor: fills all grid points with the given constant values
    fct_const = Interpolant(ssg1, (y=0.0, z=0.0))
    @test varnames(fct_const) == (:y, :z)
    @test nof(fct_const) == 2
    @test fct_const((; x=3.0)).y == 0.0
end

@testset "Interpolant — evaluation" begin
    # Exact evaluation at a grid point
    @test fct1((; x=3.0)).y ≈ 6.0
    @test fct1((; x=3.0)).z ≈ 9.0

    # Off-grid: Gridded(Linear()) interpolates between adjacent points
    # y(2.5) = (y(2) + y(3)) / 2 = (4 + 6) / 2 = 5
    @test fct1((; x=2.5)).y ≈ 5.0

    # Keyword syntax delegates to the NamedTuple dispatch
    @test fct1(; x=3.0).y ≈ 6.0

    # Extrapolation: Flat() clamps to the nearest boundary value
    # y(10.0) is clamped to y(5.0) = 10.0
    @test fct1((; x=10.0)).y ≈ 10.0

    # Multi-point: evaluating over a generator returns a StructArray
    results = fct1(points(ssg1))
    @test length(results) == 5
    @test results[3].y ≈ 6.0

    # 2D exact at a grid point
    @test fct2((; x=2.0, y=20.0)).f ≈ 22.0

    # 2D interpolated: f is bilinear so linear interpolation is exact
    @test fct2((; x=1.5, y=15.0)).f ≈ 16.5

    # Exogenous state (NoInterp): integer indexing, no interpolation across exst
    @test fct_ex((; exst=1, x=2.0)).y ≈ 2.0   # 1 * 2 = 2
    @test fct_ex((; exst=2, x=2.5)).y ≈ 5.0   # 2 * 2.5 = 5, interpolated between x=2 and x=3

    # Singleton fields: stripped from the NamedTuple before evaluation, so they are ignored
    ssg_s  = StateSpaceGrid(; x=collect(1.0:5.0), s=5.0)
    fvals_s = StructArray([(y=2.0*xi,) for xi in 1.0:5.0])
    fct_s  = Interpolant(ssg_s, fvals_s)
    @test fct_s((; x=3.0, s=5.0)).y ≈ 6.0   # s is stripped; y(3.0) = 6
end

@testset "Utility functions" begin
    # make_point: merges exogenous state index with transient NamedTuple
    @test make_point(2, (; k=1.5))        == (; exst=2, k=1.5)
    @test make_point(1, (; a=0.5, b=0.3)) == (; exst=1, a=0.5, b=0.3)

    # fisher_burmeister: a + b - sqrt(a² + b²)
    # Pythagorean triple 3-4-5:  3 + 4 - 5 = 2
    @test TransitionFunctionIteration.fisher_burmeister(3.0, 4.0) ≈ 2.0
    @test TransitionFunctionIteration.fisher_burmeister(0.0, 0.0) == 0.0
end

@testset "sample_shocks" begin
    shocks = sample_shocks(5=>10, 42, 0)

    # Shape: first(sz) = N_runs, last(sz)+burn = N_periods
    @test size(shocks) == (5, 10)

    # Values are uniform [0, 1)
    @test all(0 .<= shocks .< 1)

    # Deterministic: same seed always produces the same matrix
    @test sample_shocks(5=>10, 42, 0) == shocks

    # Different seed → different draws
    @test sample_shocks(5=>10, 99, 0) != shocks

    # Burn: adds extra columns (trimmed by simulate, but stored here)
    @test size(sample_shocks(5=>10, 42, 3)) == (5, 13)
end

@testset "simulate_markov_chain" begin
    # 3-state chain with an identity transition matrix: each state is absorbing
    exog_pts = StructArray([(z=v,) for v in [0.9, 1.0, 1.1]])
    mtrans   = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    exogenv  = Exogenv(pts=exog_pts, mtrans=mtrans, ssidx=2)
    shocks   = fill(0.5, 3, 5)   # 3 runs, 5 periods

    # Shape: N_runs × (N_periods + 1) because t=0 is included
    states_sc = simulate_markov_chain(exogenv, 1, shocks)
    @test size(states_sc) == (3, 6)

    # Scalar start: all 3 runs start at absorbing state 1 and stay there
    @test all(states_sc .== 1)

    # Vector start: each run stays at its own absorbing starting state
    states_v = simulate_markov_chain(exogenv, [1, 2, 3], shocks)
    @test states_v[:, 1]    == [1, 2, 3]     # initial conditions
    @test all(states_v[1,:] .== 1)            # run 1 stays at 1
    @test all(states_v[2,:] .== 2)            # run 2 stays at 2
    @test all(states_v[3,:] .== 3)            # run 3 stays at 3
    @test all(1 .<= states_v .<= 3)           # all states within valid range
end

@testset "Plotting recipe" begin
    # RecipesBase.apply_recipe is the low-level function @recipe compiles to.
    # It runs the recipe body and returns a Vector{RecipeData}, each element having:
    #   .args            — the data tuple (e.g. (xvec, yvec) for 1D, (x,y,z) for 2D)
    #   .plotattributes  — attributes set by --> (suggest, no-overwrite) or := (force)
    # This lets us verify recipe outputs without loading a Plots backend.
    # Custom kwargs like `transform` are passed via the incoming dict and consumed
    # by pop!(plotattributes, :transform, identity) inside the recipe.

    # 1D: all dims free, identity transform (default)
    rd = first(RecipesBase.apply_recipe(Dict{Symbol,Any}(), fct1, :y, (;)))
    @test rd.args[1] ≈ xgrid               # x-axis is the grid
    @test rd.args[2] ≈ 2 .* xgrid          # fvec: y(x) = 2x at each grid point
    @test rd.plotattributes[:title]  == "y" # --> title from varname
    @test rd.plotattributes[:xlabel] == "x" # --> xlabel from the free dimension name
    @test get(rd.plotattributes, :legend, false) == false

    # transform kwarg: passed in the dict, consumed by pop!(), applied to fvec
    rd2 = first(RecipesBase.apply_recipe(Dict{Symbol,Any}(:transform => x -> -x), fct1, :y, (;)))
    @test rd2.args[2] ≈ -2 .* xgrid
    @test !haskey(rd2.plotattributes, :transform)   # consumed, not left as a stray attribute

    # 2D grid with one dim fixed: fixedstates=(;y=20.0) removes y from free dims → 1D plot
    # _make_point injects y=20.0 for each evaluation point
    rd3 = first(RecipesBase.apply_recipe(Dict{Symbol,Any}(), fct2, :f, (; y=20.0)))
    @test rd3.args[1] ≈ xgrid2
    @test rd3.args[2] ≈ xgrid2 .+ 20.0    # f(x, 20.0) = x + 20
    @test rd3.plotattributes[:xlabel] == "x"

    # 2D: both dims free → contourf recipe with (x1vec, x2vec, fmat)
    # _points(grid) uses Iterators.product(xgrid2, ygrid2): xgrid2 varies fastest,
    # so the comprehension produces an n1×n2 matrix: fmat[i,j] = f(xgrid2[i], ygrid2[j])
    rd4 = first(RecipesBase.apply_recipe(Dict{Symbol,Any}(), fct2, :f, (;)))
    @test rd4.args[1] ≈ xgrid2
    @test rd4.args[2] ≈ ygrid2
    @test size(rd4.args[3]) == (3, 3)                    # n1 × n2
    @test rd4.args[3][2, 2] ≈ 22.0                      # f(xgrid2[2]=2.0, ygrid2[2]=20.0) = 22
    @test rd4.plotattributes[:seriestype] == :contourf   # forced via :=
    @test rd4.plotattributes[:title]      == "f"
    @test rd4.plotattributes[:xlabel]     == "x"
    @test rd4.plotattributes[:ylabel]     == "y"

    # >2D: recipe hits the else branch and throws an error
    ssg3   = StateSpaceGrid(; x=collect(1.0:2.0), y=collect(1.0:2.0), z=collect(1.0:2.0))
    fvals3 = StructArray([(f=1.0,) for _ in 1:2, _ in 1:2, _ in 1:2])
    fct3   = Interpolant(ssg3, fvals3)
    @test_throws ErrorException RecipesBase.apply_recipe(Dict{Symbol,Any}(), fct3, :f, (;))
end

@testset "Exogenv — ExogenvAR1 constructor" begin
    # --- basic construction (iid, n=3) ---
    e3 = ExogenvAR1(n=3)
    @test exnpt(e3) == 3
    @test nex(e3) == 1
    @test ssidx(e3) == 2                        # middle point for n=3
    @test propertynames(points(e3)) == (:z,)    # default shock name
    @test size(e3.mtrans) == (3, 3)
    @test all(≈(1.0), sum(e3.mtrans, dims=2))  # rows sum to 1

    # --- custom shock_name ---
    ea = ExogenvAR1(n=3, shock_name=:a)
    @test propertynames(points(ea)) == (:a,)

    # --- ssidx for odd vs larger odd n ---
    e5 = ExogenvAR1(n=5)
    @test ssidx(e5) == 3

    # --- shock_volatility flag ---
    # With shock_volatility=false (default), σ is the *state* volatility → adjusted internally.
    # For ρ=0 the adjustment is trivial (σ_shock = σ), so test with ρ≠0.
    # The state std under shock_volatility=false should match the requested σ=1.
    e_sv_false = ExogenvAR1(ρ=0.9, σ=1.0, n=5, shock_volatility=false)
    e_sv_true  = ExogenvAR1(ρ=0.9, σ=1.0, n=5, shock_volatility=true)
    # shock_volatility=true uses σ as-is, so it produces wider states
    pts_false = collect(points(e_sv_false).z)
    pts_true  = collect(points(e_sv_true).z)
    @test maximum(abs, pts_true) > maximum(abs, pts_false)

    # --- jensen correction ---
    # Without jensen, Exogenv state values are in log space and expectation is exp(μ + σ²/2).
    # With jensen=true, the mean correction shifts μ so that the ergodic mean ≈ exp(μ).
    e_log  = ExogenvAR1(n=5, σ=0.2, μ=0.0, jensen=true)
    @test propertynames(points(e_log)) == (:z,)
    @test all(>(0), points(e_log).z)                       # exp transform → all positive
    @test mean(e_log).z ≈ exp(0.0)  atol=0.05             # ergodic mean ≈ exp(μ)=1
end

@testset "Exogenv — product-space blending constructor" begin
    e1 = ExogenvAR1(n=3, shock_name=:z1)
    e2 = ExogenvAR1(n=5, shock_name=:z2)
    eb = Exogenv(e1, e2)

    # State count and transition matrix size
    @test exnpt(eb) == 3 * 5
    @test size(eb.mtrans) == (15, 15)
    @test all(≈(1.0), sum(eb.mtrans, dims=2))

    # Both field names present
    @test :z1 in propertynames(points(eb))
    @test :z2 in propertynames(points(eb))

    # Steady-state index: linear position of (ssidx(e1), ssidx(e2)) = (2,3) in
    # vec(Iterators.product(1:3, 1:5)) where first index varies fastest
    combos = vec(collect(Iterators.product(1:3, 1:5)))
    expected_ssidx = findfirst(==((ssidx(e1), ssidx(e2))), combos)
    @test ssidx(eb) == expected_ssidx

    # Transition matrix is kron(P2, P1) — rows sum to 1 already checked above
    P1 = e1.mtrans
    P2 = e2.mtrans
    @test eb.mtrans ≈ kron(P2, P1)

    # Duplicate field names → error
    @test_throws ErrorException Exogenv(e1, e1)
end

@testset "Exogenv — ergodic statistics" begin
    # Use a 3-state iid chain (ρ=0, μ=0) as a simple fixture.
    # Its ergodic distribution is uniform (all rows of P are equal by symmetry of Rouwenhorst).
    e = ExogenvAR1(ρ=0.0, σ=0.5, μ=0.0, n=3)

    pvec = TransitionFunctionIteration.ergodic_distro(e)
    @test length(pvec) == 3
    @test sum(pvec) ≈ 1.0
    @test all(>=(0.0), pvec)

    # Symmetric chain → ergodic mean ≈ 0
    @test mean(e).z  ≈ 0.0  atol=1e-10

    # Variance is non-negative
    @test var(e).z >= 0.0

    # std = sqrt(var)
    @test std(e).z ≈ sqrt(var(e).z)

    # mode: returns a valid index and the corresponding state NamedTuple
    idx_mode, pt_mode = TransitionFunctionIteration.mode(e)
    @test 1 <= idx_mode <= exnpt(e)
    @test pt_mode == points(e)[idx_mode]

    # median: odd n → returns middle-index state
    idx_med, pt_med = median(e)
    @test idx_med == 2
    @test pt_med == points(e)[2]

    # median: even n → warns and returns (nothing, nothing)
    e2 = ExogenvAR1(n=4)
    @test_warn r"Even number" begin
        idx2, pt2 = median(e2)
        @test isnothing(idx2)
        @test isnothing(pt2)
    end
end

end  # @testset "TransitionFunctionIteration"
