@testset "IRF plotting" begin

    using RecipesBase
    using StructArrays
    using TransitionFunctionIteration: Transformations, IRFResult, ReportedVariable, IRFPage

    # RecipesBase.is_key_supported has no methods outside a plotting backend;
    # define a stub so apply_recipe can be called in isolation.
    RecipesBase.is_key_supported(::Symbol) = true

    T = 5
    # Unconditional baseline: C rises slowly
    uncond_C = [1.0, 1.05, 1.1, 1.15, 1.2]
    # Shock path: C rises faster
    shock_C  = [1.0, 1.15, 1.3, 1.4, 1.5]

    means = (
        uncond  = StructArray((C = uncond_C,)),
        shock_1 = StructArray((C = shock_C,)),
    )
    stdevs = (
        uncond  = StructArray((C = fill(0.05, T),)),
        shock_1 = StructArray((C = fill(0.05, T),)),
    )
    irf_mock = IRFResult(means, stdevs, (;), (uncond=(;), shock_1=(;)))

    @testset "Transformations module" begin
        x     = [2.0, 3.0, 4.0]
        basis = [1.0, 1.0, 1.0]

        @test Transformations.Level()(x)               == x
        @test Transformations.DiffLevel()(x; basis=basis) ≈ [1.0, 2.0, 3.0]
        @test Transformations.DiffPct()(x; basis=basis) ≈ [100.0, 200.0, 300.0]
        @test Transformations.RatioPct()(x; basis=basis) ≈ [100.0, 200.0, 300.0]
        @test Transformations.DiffInitLevel()(x)        ≈ [0.0, 1.0, 2.0]
        @test Transformations.DiffInitPct()(x)          ≈ [0.0, 50.0, 100.0]
    end

    @testset "IRF recipe: DiffLevel uses uncond as default basis" begin
        rv = ReportedVariable(:C; transformation=Transformations.DiffLevel())
        # Apply the recipe -- returns a list of plot specs
        specs = RecipesBase.apply_recipe(Dict{Symbol,Any}(), irf_mock, rv, :shock_1)
        @test length(specs) >= 1
        # First spec carries (x_data, y_data) as the positional series args
        y = specs[1].args[2]
        expected = shock_C .- uncond_C
        @test y ≈ expected
    end

    @testset "IRF recipe: Level transformation ignores basis" begin
        rv = ReportedVariable(:C; transformation=Transformations.Level())
        specs = RecipesBase.apply_recipe(Dict{Symbol,Any}(), irf_mock, rv, :shock_1)
        y = specs[1].args[2]
        @test y ≈ shock_C
    end

    @testset "IRF recipe: RatioPct uses uncond as default basis" begin
        rv = ReportedVariable(:C; transformation=Transformations.RatioPct())
        specs = RecipesBase.apply_recipe(Dict{Symbol,Any}(), irf_mock, rv, :shock_1)
        y = specs[1].args[2]
        expected = (shock_C ./ uncond_C .- 1) .* 100
        @test y ≈ expected
    end

    @testset "IRF recipe: explicit basis override works" begin
        rv = ReportedVariable(:C; transformation=Transformations.DiffLevel())
        # Use shock_1 as its own basis -> should give zeros
        specs = RecipesBase.apply_recipe(
            Dict{Symbol,Any}(:basis => :shock_1), irf_mock, rv, :shock_1)
        y = specs[1].args[2]
        @test all(y .≈ 0.0)
    end

    @testset "IRFPage dispatches over ReportedVariables" begin
        page = IRFPage([
            ReportedVariable(:C),
            ReportedVariable(:C; transformation=Transformations.DiffLevel()),
        ])
        @test size(page) == (1, 2)                    # a vector is one row of panels
        @test page[1].varname == :C
        @test page[2].transformation isa Transformations.DiffLevel
        @test size(IRFPage([:C :K; :C :K])) == (2, 2)  # a matrix keeps its shape; names become ReportedVariables
        @test IRFPage([:C, :K])[2] == ReportedVariable(:K)
    end

    @testset "IRFPage recipe — layout follows the page shape" begin
        for (page, layout) in ((IRFPage([:C, :C, :C]), (1, 3)), (IRFPage(reshape([:C, :C, :C, :C], 2, 2)), (2, 2)))
            attrs = Dict{Symbol,Any}()
            specs = RecipesBase.apply_recipe(attrs, irf_mock, page, :shock_1)
            @test attrs[:layout] == layout
            @test [s.plotattributes[:subplot] for s in specs] == 1:length(page)
        end
    end

    @testset "IRF recipe — pers parameter" begin
        rv = ReportedVariable(:C)
        specs = RecipesBase.apply_recipe(
            Dict{Symbol,Any}(:pers => 1:3), irf_mock, rv, :shock_1)
        @test length(specs) >= 1
        # x-axis should cover the requested period range
        @test specs[1].args[1] == 1:3
    end

    @testset "IRF recipe — sd ribbon" begin
        rv = ReportedVariable(:C)
        # sd=true adds a ribbon; recipe should complete without error
        specs = RecipesBase.apply_recipe(
            Dict{Symbol,Any}(:sd => true), irf_mock, rv, :shock_1)
        @test length(specs) >= 1
        # ribbon attribute should be set
        @test haskey(specs[1].plotattributes, :ribbon)
    end

    @testset "IRF recipe — fan chart" begin
        # irf.paths.shock_1.C is T × N (periods by paths), as simulate_irf stores it. All paths
        # share the starting period, so the fan has zero width there.
        N_paths = 20
        paths_C = shock_C .+ [zeros(1, N_paths); randn(Random.Xoshiro(1), T - 1, N_paths)]   # T × N
        irf_with_paths = IRFResult(
            irf_mock.means,
            irf_mock.stdevs,
            irf_mock.other_stats,
            (shock_1 = StructArray((C = paths_C,)),),
        )

        rv = ReportedVariable(:C)
        specs = RecipesBase.apply_recipe(
            Dict{Symbol,Any}(:fan => true, :n_bands => 3, :tail_pct => 5.0),
            irf_with_paths, rv, :shock_1)
        # fan=true produces n_bands series + 1 mean line
        @test length(specs) == 4

        # each band spans quantiles across paths, period by period
        for (k, p_lo) in enumerate((0.05, 0.05 + 0.15, 0.05 + 0.30))
            pers, y_lo = specs[k].args
            y_hi = specs[k].plotattributes[:fillrange]
            @test pers == 0:T-1
            @test y_lo ≈ [quantile(paths_C[t, :], p_lo) for t in 1:T]
            @test y_hi ≈ [quantile(paths_C[t, :], 1 - p_lo) for t in 1:T]
            @test y_lo[1] ≈ y_hi[1] ≈ shock_C[1]
        end
    end

    @testset "IRFPage recipe — smoke test" begin
        page = IRFPage([
            ReportedVariable(:C),
            ReportedVariable(:C; transformation=Transformations.DiffLevel()),
        ])
        specs = RecipesBase.apply_recipe(Dict{Symbol,Any}(), irf_mock, page, :shock_1)
        # One @series per variable in the page
        @test length(specs) >= length(page)
    end

end
