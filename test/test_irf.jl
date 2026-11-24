# IRF tests — depend on obj, stv, simseries from test_ngm.jl (included first)

@testset "simulate_irf — basic call" begin
    irf = simulate_irf(obj, simseries)

    @test irf isa IRFResult
    @test !isempty(TFI.shock_names(irf))

    # means and stdevs have one key per shock
    @test propertynames(irf.means) == TFI.shock_names(irf)
    @test propertynames(irf.stdevs) == TFI.shock_names(irf)

    # each shock's mean path is a non-empty StructArray
    first_shock = first(TFI.shock_names(irf))
    @test irf.means[first_shock] isa StructArray
    @test length(irf.means[first_shock]) > 0

    # no path_vars → paths is a NamedTuple of empty NamedTuples
    @test all(isempty, irf.paths)
end

@testset "simulate_irf — explicit shock and path_vars" begin
    irf2 = simulate_irf(obj, simseries;
        shocks=(pos_shock=5,), path_vars=(:K, :C), T=10, N=100)

    @test only(TFI.shock_names(irf2)) == :pos_shock
    @test irf2.paths.pos_shock isa StructArray
    @test :K in propertynames(irf2.paths.pos_shock)
    @test :C in propertynames(irf2.paths.pos_shock)
end

@testset "fan chart of simulate_irf paths: bands are quantiles across paths by period" begin
    # T ≠ N, so reading the stored T × N paths the wrong way round cannot pass
    # (RecipesBase.is_key_supported is stubbed in test_irfplotting.jl)
    irf_f = simulate_irf(obj, simseries; shocks = (up = 5,), path_vars = (:C,), T = 10, N = 100)
    paths = irf_f.paths.up.C
    @test size(paths) == (10, 100)
    specs = RecipesBase.apply_recipe(Dict{Symbol,Any}(:fan => 2, :tail_pct => 5.0),
                                     irf_f, TFI.ReportedVariable(:C), :up)
    for (k, p_lo) in enumerate((0.05, 0.275))
        pers, y_lo = specs[k].args
        @test pers == 0:9
        @test y_lo ≈ [quantile(paths[t, :], p_lo) for t in 1:10]
        @test specs[k].plotattributes[:fillrange] ≈ [quantile(paths[t, :], 1 - p_lo) for t in 1:10]
    end
    # every path starts at the same state, so the fan has zero width in the starting period
    @test specs[1].args[2][1] ≈ specs[1].plotattributes[:fillrange][1]
end

@testset "simulate_irf — NamedTuple endogstart" begin
    irf3 = simulate_irf(obj, stv.endog)
    irf_ref = simulate_irf(obj, simseries)
    @test TFI.shock_names(irf3) == TFI.shock_names(irf_ref)
end

@testset "IRFResult constructor" begin
    sa = StructArray([(K=Float64(t), C=Float64(t)) for t in 1:5])
    lines = (shock1=(; mean=sa, std=sa), shock2=(; mean=sa, std=sa))
    paths = (;)
    irf_direct = IRFResult(lines, paths)

    @test TFI.shock_names(irf_direct) == (:shock1, :shock2)
    @test irf_direct.means.shock1 === sa
    @test irf_direct.stdevs.shock1 === sa
    @test isempty(irf_direct.other_stats)
    @test irf_direct.paths === paths
end

@testset "simulate_irf — hold periods fill the dynamic variables of period 1" begin
    kw = (exog_start = 2, shocks = (uncond = Int[], up = 3), T = 6, N = 5)
    irf = simulate_irf(obj_as, (x = 1.0,); kw...)          # default: dynamic_depth = 2
    @test irf.provenance.settings.hold_periods == 2
    for sh in (:uncond, :up)
        @test length(irf.means[sh]) == 6                      # the held periods are not output
    end
    @test !any(isnan, irf.means.up.dexst) && !any(isnan, irf.means.up.ddexst)

    # period 1 is the starting state for every shock; the shock arrives in period 2
    @test irf.means.up.exst[1] == irf.means.uncond.exst[1] == 2
    @test irf.means.up.exst[2] == 3
    @test irf.means.up.dexst[1:2] == [0, 1] && irf.means.up.ddexst[1:2] == [0, 1]

    # without held periods the history of period 1 is missing
    irf0 = simulate_irf(obj_as, (x = 1.0,); kw..., hold_periods = 0)
    @test length(irf0.means.up) == 6
    @test isnan(irf0.means.up.dexst[1]) && all(isnan, irf0.means.up.ddexst[1:2])
    @test irf0.means.up.exst[1:2] == [2, 3]
end
