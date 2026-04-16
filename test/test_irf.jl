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
