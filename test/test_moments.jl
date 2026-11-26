using Tables, DataAPI

struct CustomBasisTransformation <: Transformations.AbstractTransformation end
(::CustomBasisTransformation)(x, basis, trend, basistrend, init) = x - basis
struct CustomNoBasisTransformation <: Transformations.AbstractTransformation end
Transformations.requires_basis(::CustomNoBasisTransformation) = false
(::CustomNoBasisTransformation)(x, basis, trend, basistrend, init) = 2x

@testset "Impulse moments — optional and named basis" begin
    T = Transformations
    table = (__shock = repeat(["up", "control", "uncond"], 2), __stat = repeat(["mean", "std"], inner = 3),
             __period = ones(Int, 6), C = [6.0, 2.0, 3.0, 0.6, 0.2, 0.3])
    irf = IRFResult(table)
    single = IRFResult(map(v -> v[[1, 4]], table))
    level = ImpulseMoment(:C, :up, 1; basis = :absent, compare_transformation = T.DiffLevel())
    diff = ImpulseMoment(:C, :up, 1; basis = :control, transformation = T.DiffLevel())
    ratio = ImpulseMoment(:C, :up, 1; basis = :control, transformation = T.RatioPct())
    @test getmoment(single, level) == 6.0
    @test getmoment(single, ImpulseMoment(:C, :up, 1; transformation = CustomNoBasisTransformation())) == 12.0
    @test getmoment(irf, diff) == 4.0
    @test getmoment(irf, ImpulseMoment(:C, :up, 1; stat = :std, basis = :control, transformation = T.DiffLevel())) ≈ 0.4
    @test getmoment(irf, ratio) == 200.0
    @test getmoment(irf, ImpulseMoment(:C, :up, 1; transformation = T.DiffLevel())) == 3.0
    @test getmoment(irf, ImpulseMoment(:C, :up, 1; transformation = T.RatioPct())) == 100.0
    @test getmoment(irf, ImpulseMoment(:C, :up, 1; basis = :up, transformation = T.DiffLevel())) == 0.0
    @test_throws KeyError getmoment(single, diff)
    @test_throws KeyError getmoment(single, ImpulseMoment(:C, :up, 1; transformation = T.RatioPct()))
    @test getmoment(irf, ImpulseMoment(:C, :up, 1; basis = :control, transformation = CustomBasisTransformation())) == 4.0
    reordered = IRFResult(map(v -> v[[3, 2, 1, 6, 5, 4]], table))
    @test getmoment(reordered, ratio) == getmoment(irf, ratio)
    nested = CalculatedMoment(nt -> nt.x * 2, (x = CalculatedMoment(nt -> nt.r, (r = ratio,)),))
    @test getmoment(irf, nested) == 400.0
    economies = (a = single, b = single)
    @test compare_economies((; level), economies).data == [6.0 0.0]
end

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
    params(scale) = TFI.MomentArray(reshape([0.9scale, 2.0scale], :, 1), (:β, :σ), (:value,); rowkey = :parameter)
    result(scale) = EconomyResult(; stats = stats(scale), irf = irf(scale), params = params(scale))

    @testset "construction" begin
        sm = SimulationMoment(:C)
        @test sm isa AbstractMoment && sm.sample == :uncond && sm.stat == :mean
        im = ImpulseMoment(:C, :hike, 2)
        @test im isa AbstractMoment && im.stat == :mean && im.period == 2
        @test im.transformation isa T.Level
        @test_throws UndefKeywordError ImpulseMoment(var = :C, shock = :hike)
    end

    @testset "getmoment dispatches on the moment" begin
        res = result(1.0)
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
        @test_throws ArgumentError getmoment(EconomyResult(; stats = stats(1.0)), ImpulseMoment(:C, :hike, 2))
        @test_throws KeyError getmoment(res, SimulationMoment(:nope))
        @test_throws KeyError getmoment(res, ImpulseMoment(:C, :cut, 2))
        @test_throws BoundsError getmoment(res, ImpulseMoment(:C, :hike, 9))
    end

    moments = (mean_C = SimulationMoment(:C),
               std_K = SimulationMoment(:K; stat = :std, compare_transformation = T.DiffLevel()),
               irf_C = ImpulseMoment(:C, :hike, 3; transformation = T.DiffLevel()))

    @testset "economies tuple" begin
        economies = (base = result(1.0), big = result(2.0))
        @test getmoment(economies, moments.mean_C, :big, :base) == 2.0
        @test getmoment(economies, moments.std_K, :big, :base) == 4.0     # DiffLevel against base: 8 - 4
        @test getmoment(economies, moments.std_K, :base, 1) == 4.0        # the base uses its transformation
        @test getmoment(economies, moments.irf_C, :big) == 2.0        # non-base: Level comparison of its own response (6 - 4)
        @test getmoment(economies, moments.irf_C, :base) == 1.0       # base: its own DiffLevel against uncond
        tab = compare_economies(moments, economies)
        @test tab isa TFI.MomentArray && size(tab) == (3, 2)
        @test tab(:irf_C, :big) == 2.0
        # legacy form: describe outputs and bare variable names
        legacy = compare_economies((C = :C,), (base = stats(1.0), big = stats(2.0)))
        @test legacy(:C, :big) == 2.0
    end

    @testset "CalculatedMoment" begin
        allocs(f, args...) = (f(args...); @allocated f(args...))
        res1 = result(1.0)
        economies = (base = res1, big = result(2.0))
        irfC = ImpulseMoment(:C, :hike, 3; transformation = T.DiffLevel())       # 1 * scale
        stdK = SimulationMoment(:K; stat = :std)                                  # 4 * scale
        total = CalculatedMoment(nt -> nt.irf_C + nt.std_K, (irf_C = irfC, std_K = stdK))
        nested = CalculatedMoment(nt -> 2nt.total - nt.mean_C, (total = total, mean_C = SimulationMoment(:C)))

        @test total isa AbstractMoment && total.transformation isa T.Level && total.compare_transformation isa T.Level
        @test TFI.transformation(total) === total.transformation && TFI.compare_transformation(total) === total.compare_transformation
        @test_throws ArgumentError CalculatedMoment(identity, (a = 1.0,))

        @testset "values and nesting" begin
            @test getmoment(res1, total) == 5.0                              # inputs use their own transformations
            @test getmoment(res1, nested) == 9.0
            @test getmoment(res1, CalculatedMoment(nt -> nt.a, (a = irfC,); transformation = T.Pct())) == 100.0
            @test getmoment(economies, total, :big, :base) == 10.0           # Level compare: the calculated value
            cmp = CalculatedMoment(total.f, total.inputs; compare_transformation = T.DiffLevel())
            @test getmoment(economies, cmp, :big, :base) == 5.0              # 10 - 5; child compare_transformations play no role
            @test getmoment(economies, cmp, :base) == 5.0
            child_cmp = CalculatedMoment(nt -> nt.a, (a = SimulationMoment(:K; stat = :std, compare_transformation = T.DiffLevel()),))
            @test getmoment(economies, child_cmp, :big, :base) == 8.0        # the child's DiffLevel is ignored
            tab = compare_economies((total = total, nested = nested, std_K = stdK), economies)
            @test tab(:total, :big) == 10.0 && tab(:nested, :base) == 9.0 && tab(:nested, :big) == 18.0
            @test_throws ArgumentError getmoment(EconomyResult(; stats = stats(1.0)), total)
            @test_throws KeyError getmoment(res1, CalculatedMoment(nt -> nt.a, (a = SimulationMoment(:nope),)))
        end

        @testset "inference and allocations" begin
            @test (@inferred getmoment(res1, total)) == 5.0
            @test (@inferred getmoment(res1, nested)) == 9.0
            @test (@inferred getmoment(economies, total, :big, :base)) == 10.0
            @test allocs(getmoment, res1, nested) <= allocs(getmoment, res1, irfC) + allocs(getmoment, res1, stdK) + allocs(getmoment, res1, SimulationMoment(:C))
        end
    end
end

struct ToyParametersModel <: TFI.AbstractModel
    β::Float64
    flag::Bool
    label::String
end

@testset "Moments — parameters and calculated expressions" begin
    T = TFI.Transformations
    stats = (uncond = TFI.MomentArray([1.0 2.0; 3.0 4.0], (:C, :K), (:mean, :std)),)
    params = TFI.MomentArray(reshape([0.9, 2.0], :, 1), (:β, :σ), (:value,); rowkey = :parameter)
    res = EconomyResult(; stats, params)

    @testset "collect_parameters" begin
        p = collect_parameters(ToyParametersModel(0.9, true, "x"))
        @test row_names(p) == (:β, :flag) && col_names(p) == (:value,) && p.data == reshape([0.9, 1.0], :, 1)
        @test row_key(p) == :parameter
    end

    @testset "ParameterMoment" begin
        @test getmoment(res, ParameterMoment(:β)) == 0.9
        @test getmoment(res, ParameterMoment(:σ; transformation = T.Pct())) == 200.0
        @test getmoment(EconomyResult(; params = (β = 0.5,)), ParameterMoment(:β)) == 0.5
        @test getmoment((a = res, b = EconomyResult(; params = (β = 0.5, σ = 1.0))), ParameterMoment(:σ), :b, :a) == 1.0
        @test_throws KeyError getmoment(res, ParameterMoment(:nope))
        @test_throws ArgumentError getmoment(EconomyResult(; stats), ParameterMoment(:β))
        @test_throws ArgumentError getmoment(EconomyResult(; params), SimulationMoment(:C))
        @test_throws ArgumentError getmoment(stats, ParameterMoment(:β))     # a describe output cannot serve it
    end

    C = SimulationMoment(:C)                      # 1
    K = SimulationMoment(:K; stat = :std)         # 4

    @testset "@CalculatedMoment" begin
        total = @CalculatedMoment C + K
        @test total isa CalculatedMoment && getmoment(res, total) == 5.0
        @test getmoment(res, @CalculatedMoment sqrt(K^2 + 2C)) == sqrt(18.0)
        @test getmoment(res, @CalculatedMoment -C) == -1.0
        @test getmoment(res, @CalculatedMoment log(1 + C / K)) ≈ log(1.25)
        @test getmoment(res, @CalculatedMoment total * 2) == 10.0                # nests
        @test getmoment(res, @CalculatedMoment C * ParameterMoment(:σ)) == 2.0   # any kind of moment
        scaled(x; by) = x * by
        @test getmoment(res, @CalculatedMoment scaled(K; by = 3.0)) == 12.0      # keywords are constants
        @test getmoment(res, @CalculatedMoment max(C, K, 10.0)) == 10.0          # constants among moments
        @test getmoment(res, @CalculatedMoment (a = 2; C * a)) == 2.0            # blocks and assignments pass through
        # a call with no moment among its arguments is the plain call
        @test (@CalculatedMoment SimulationMoment(:C)) isa SimulationMoment
        @test (@CalculatedMoment 1 + 2) == 3
        ratio_of(v) = @CalculatedMoment SimulationMoment(v; stat = :std) / SimulationMoment(v)
        @test getmoment(res, ratio_of(:K)) == 4 / 3
        # equal expressions are ==, and reading is inferable
        @test (@CalculatedMoment C + K) == (@CalculatedMoment C + K)
        @test (@inferred getmoment(res, total)) == 5.0
        @test getmoment((a = res, b = EconomyResult(; stats = (uncond = TFI.MomentArray([2.0 4.0; 6.0 8.0], (:C, :K), (:mean, :std)),))),
                        total, :b, :a) == 10.0
    end

    @testset "@CalculatedMoment with named moments" begin
        named = (mean_C = C, std_K = K)
        @test getmoment(res, @CalculatedMoment named std_K / mean_C) == 4.0
        x = 10.0
        @test (@CalculatedMoment named x + 1) == 11.0                              # a variable in scope is a variable
        @test getmoment(res, @CalculatedMoment named std_K - x) == -6.0
        @test getmoment(res, @CalculatedMoment named named.std_K + mean_C) == 5.0   # field access: only the object is a name
        @test_throws ArgumentError (@CalculatedMoment named nope / mean_C)
    end

    @testset "fallback" begin
        @test_throws ArgumentError TFI.momentvalues(1.0, C)
        err = try; getmoment(1.0, ParameterMoment(:β)); catch e; e; end
        @test occursin("ParameterMoment", err.msg) && occursin("Float64", err.msg)
    end
end

@testset "MomentArray row key" begin
    allocs(f, x) = (f(x); @allocated f(x))
    data = [1.0 2.0; 3.0 4.0; 5.0 6.0]
    plain = TFI.MomentArray(data, (:a, :b, :c), (:mean, :sd))
    keyed = TFI.MomentArray(data, (:a, :b, :c), (:mean, :sd); rowkey = :variable)

    @testset "default and keyword" begin
        @test row_key(plain) == :__row && row_key(keyed) == :variable
        @test Tables.columnnames(Tables.columns(plain)) == (:__row, :mean, :sd)
        @test Tables.columnnames(Tables.columns(keyed)) == (:variable, :mean, :sd)
        @test Tables.getcolumn(Tables.columns(keyed), :variable) == ["a", "b", "c"]
        @test DataAPI.metadata(keyed, "tfi.rowkey") == "variable" && DataAPI.metadata(plain, "tfi.rowkey") == "__row"
        @test keyed == plain                                              # equal as matrices
        @test_throws ArgumentError TFI.MomentArray(data, (:a, :b, :c), (:mean, :sd); rowkey = :mean)   # the key is a column name
        @test TFI.MomentArray{(:a, :b, :c), (:mean, :sd)}(data; rowkey = :equation) isa TFI.MomentArray
        @test @MomentArray(data, (:a, :b, :c), (:mean, :sd)) isa TFI.MomentArray
    end

    @testset "the key survives operations" begin
        @test row_key(keyed .* 2.0) == :variable && row_key(2 .* keyed .+ keyed) == :variable
        @test row_key(vcat(keyed, keyed[1:1, :] |> m -> TFI.MomentArray(m, (:d,), (:mean, :sd); rowkey = :variable))) == :variable
        @test row_key(hcat(keyed, TFI.MomentArray(data[:, 1:1], (:a, :b, :c), (:p50,); rowkey = :variable))) == :variable
        @test row_key(getmoment(keyed, (:a, :b), (:mean,))) == :variable
        @test row_key(getmoment(keyed, Val((:a,)), Val((:sd,)))) == :variable
        # arrays of different kinds of rows give the default
        other = TFI.MomentArray(data, (:x, :y, :z), (:mean, :sd); rowkey = :moment)
        @test row_key(vcat(keyed, other)) == :__row
    end

    @testset "creators name their rows" begin
        p = TFI.SimulationPanel(StructArray([(C = 1.0 + 0.1i, K = 2.0 + 0.1j) for i in 1:3, j in 1:10]), (:K,))
        d = describe(p)
        @test row_key(d.moments.uncond) == :variable && row_key(d.stderrs.uncond) == :variable
        cmp = compare_economies((m = SimulationMoment(:C),), (a = d.moments, b = d.moments))
        @test row_key(cmp) == :moment
    end

    @testset "reading" begin
        # the key is restored from the table's metadata; a foreign table needs it by keyword
        @test row_key(TFI.MomentArray(keyed)) == :variable && TFI.MomentArray(keyed) == keyed
        nt = Tables.columntable(keyed)                                    # no metadata
        @test_throws ArgumentError TFI.MomentArray(nt)
        @test row_key(TFI.MomentArray(nt; rowkey = :variable)) == :variable
        @test row_key(TFI.MomentArray{(:a, :b, :c), (:mean, :sd)}(nt; rowkey = :variable)) == :variable
        @test row_key(TFI.MomentArray{(:a, :b, :c), (:mean, :sd), Float64}(nt; rowkey = :variable)) == :variable
        # in the type: inferred, and the wrong key is an error
        restore(t) = TFI.MomentArray{(:a, :b, :c), (:mean, :sd), Float64, SMatrix{3, 2, Float64}, :variable}(t)
        @test (@inferred restore(nt)) == keyed && row_key(restore(nt)) == :variable
        @test_throws ArgumentError TFI.MomentArray{(:a, :b, :c), (:mean, :sd), Float64, Matrix{Float64}, :equation}(nt)
    end

    @testset "writing allocates nothing" begin
        for m in (keyed, plain, TFI.MomentArray(SMatrix{3, 2}(data), (:a, :b, :c), (:mean, :sd); rowkey = :equation))
            @test allocs(Tables.columns, m) == 0
            @test allocs(Tables.columntable, m) == 0
            @test allocs(row_key, m) == 0
        end
        @test Base.return_types(Tables.columns, (typeof(keyed),))[1] |> isconcretetype
        k2 = TFI.MomentArray(data, (:a, :b, :c), (:mean, :sd); rowkey = :variable)
        build(d) = TFI.MomentArray(d, (:a, :b, :c), (:mean, :sd); rowkey = :variable)
        @test (@inferred build(data)) isa TFI.MomentArray{(:a, :b, :c), (:mean, :sd), Float64, Matrix{Float64}, :variable}
    end
end
