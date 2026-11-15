using Tables, DataAPI

@testset "Moments — SimulationMoment and ImpulseMoment" begin
    T = TFI.Transformations

    # describe-style statistics and an IRF built from tables (no metadata needed)
    stats(scale) = (uncond = TFI.MomentArray([scale 2scale; 3scale 4scale], (:C, :K), (:mean, :std)),
                    bust   = TFI.MomentArray([10scale 20scale; 30scale 40scale], (:C, :K), (:mean, :std)))
    irf_table(scale) = (
        __shock  = repeat(["hike", "uncond"], inner = 3, outer = 2),
        __stat   = repeat(["mean", "std"], inner = 6),
        __period = repeat(1:3, 4),
        C = scale .* [1.0, 2.0, 3.0, 1.0, 1.5, 2.0, 0.1, 0.2, 0.3, 0.1, 0.1, 0.1])
    irf(scale) = IRFResult(irf_table(scale))

    @testset "construction" begin
        sm = SimulationMoment(:C)
        @test sm isa AbstractMoment && sm.sample == :uncond && sm.stat == :mean
        im = ImpulseMoment(:C, :hike, 2)
        @test im isa AbstractMoment && im.stat == :mean && im.period == 2
        @test im.transformation isa T.Level
        @test_throws UndefKeywordError ImpulseMoment(var = :C, shock = :hike)
    end

    @testset "getmoment dispatches on the moment" begin
        res = EconomyResults(stats(1.0), irf(1.0))
        @test getmoment(res, SimulationMoment(:K; stat = :std)) == 4.0
        @test getmoment(res, SimulationMoment(:C; sample = :bust)) == 10.0
        @test getmoment(res, ImpulseMoment(:C, :hike, 2)) == 2.0                       # Level: the path value
        @test getmoment(res, ImpulseMoment(:C, :hike, 3; transformation = T.DiffLevel())) == 1.0   # minus uncond
        @test getmoment(res, ImpulseMoment(:C, :hike, 3; transformation = T.RatioPct())) ≈ 50.0
        @test getmoment(res, ImpulseMoment(:C, :hike, 2; stat = :std)) == 0.2
        # describe output alone serves simulation moments
        @test getmoment(stats(1.0), SimulationMoment(:C)) == 1.0
        @test getmoment(irf(1.0), ImpulseMoment(:C, :hike, 1)) == 1.0
        # errors name what is missing
        @test_throws ArgumentError getmoment(EconomyResults(stats(1.0)), ImpulseMoment(:C, :hike, 2))
        @test_throws KeyError getmoment(res, SimulationMoment(:nope))
        @test_throws KeyError getmoment(res, ImpulseMoment(:C, :cut, 2))
        @test_throws BoundsError getmoment(res, ImpulseMoment(:C, :hike, 9))
    end

    moments = (mean_C = SimulationMoment(:C),
               std_K = SimulationMoment(:K; stat = :std, compare_transformation = T.DiffLevel()),
               irf_C = ImpulseMoment(:C, :hike, 3; transformation = T.DiffLevel()))

    @testset "economies tuple" begin
        economies = (base = EconomyResults(stats(1.0), irf(1.0)), big = EconomyResults(stats(2.0), irf(2.0)))
        @test getmoment(economies, moments.mean_C, :big, :base) == 2.0
        @test getmoment(economies, moments.std_K, :big, :base) == 4.0     # DiffLevel against base: 8 - 4
        @test getmoment(economies, moments.std_K, :base, 1) == 4.0        # the base uses its transformation
        @test getmoment(economies, moments.irf_C, :big) == 6.0        # non-base: compare_transformation (Level) of the path value
        @test getmoment(economies, moments.irf_C, :base) == 1.0       # base: its own DiffLevel against uncond
        tab = compare_economies(moments, economies)
        @test tab isa TFI.MomentArray && size(tab) == (3, 2)
        @test tab(:irf_C, :big) == 6.0
        # legacy form: describe outputs and bare variable names
        legacy = compare_economies((C = :C,), (base = stats(1.0), big = stats(2.0)))
        @test legacy(:C, :big) == 2.0
    end

    @testset "long table gives the same table" begin
        rows = NamedTuple[]
        for (econ, scale) in (("base", 1.0), ("big", 2.0)), m in values(moments)
            res = EconomyResults(stats(scale), irf(scale))
            vals = (TFI.momentvalues(res, m)[1], TFI.momentvalues(res, m)[2])
            for (k, v) in zip(moment_keys(m), vals)
                push!(rows, (; economy = econ, kind = String(k.kind), group = String(k.group),
                             variable = String(k.variable), stat = String(k.stat), period = k.period, value = v))
            end
        end
        table = Tables.columntable(unique(rows))
        economies = (base = EconomyResults(stats(1.0), irf(1.0)), big = EconomyResults(stats(2.0), irf(2.0)))
        from_table = compare_economies(moments, (:base, :big), table)
        @test from_table.data == compare_economies(moments, economies).data
        @test row_names(from_table) == keys(moments) && col_names(from_table) == (:base, :big)
        @test compare_economies(moments, [:big, :base], table; base = :big)(:std_K, :base) == -4.0
        @test getmoment(TFI.MomentTable(table), moments.irf_C, :big) == 6.0
        @test_throws KeyError getmoment(TFI.MomentTable(table), SimulationMoment(:K; stat = :p99), :big)
        @test_throws ArgumentError TFI.MomentTable((economy = ["a"],))
        @test_throws ArgumentError TFI.MomentTable(Tables.columntable([rows[1], rows[1]]))
    end
end

@testset "with_rowkey: a MomentArray table with its own row-name column" begin
    m = TFI.MomentArray([1.0 2.0; 3.0 4.0], (:a, :b), (:mean, :std))
    k = with_rowkey(m, :variable)
    @test Tables.columnnames(Tables.columns(k)) == (:variable, :mean, :std)
    @test Tables.getcolumn(Tables.columns(k), :variable) == ["a", "b"]
    @test DataAPI.metadata(k, "tfi.rowkey") == "variable"
    @test DataAPI.metadata(k, "tfi.kind") == "MomentArray"
    @test TFI.MomentArray(Tables.columntable(k); rowkey = :variable) == m
    @test Tables.columnnames(Tables.columns(m)) == (:__row, :mean, :std)   # the default is unchanged
end
