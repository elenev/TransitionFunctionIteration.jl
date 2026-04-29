@testset "MomentArray — construction" begin

    data = [Float64(i + 10j) for i in 1:3, j in 1:4]
    m = TFI.MomentArray(data, (:a, :b, :c), (:mean, :p50, :p90, :p99))

    @test size(m)    == (3, 4)
    @test eltype(m)  == Float64
    @test m.data === data

    # Row and column index maps
    @test row_names(m) == (:a, :b, :c)
    @test col_names(m) == (:mean, :p50, :p90, :p99)

end

@testset "MomentArray — AbstractArray interface" begin

    data = reshape(Float64.(1:12), 3, 4)
    m = TFI.MomentArray(data, (:a, :b, :c), (:s1, :s2, :s3, :s4))

    # Cartesian indexing
    @test m[1, 1] == data[1, 1]
    @test m[2, 3] == data[2, 3]

    # Linear indexing
    @test m[3]  == data[3]
    @test m[12] == data[12]

    # setindex! pass-through
    m[1, 1] = 99.0
    @test data[1, 1] == 99.0

end

@testset "MomentArray — getmoment dispatch" begin

    data = [Float64(i + 10j) for i in 1:3, j in 1:4]
    m = TFI.MomentArray(data, (:a, :b, :c), (:mean, :p50, :p90, :p99))

    # Scalar / Scalar → Number
    @test m(:a, :mean) == data[1, 1]
    @test m(:c, :p99) == data[3, 4]

    # Scalar / Collection → NamedTuple (all stats for one variable)
    row_a = m(:a, (:mean, :p50))
    @test row_a isa NamedTuple
    @test row_a.mean == data[1, 1]
    @test row_a.p50  == data[1, 2]

    # Scalar / Colon → NamedTuple of all stats
    row_b_all = m(:b, :)
    @test row_b_all isa NamedTuple
    @test length(row_b_all) == 4

    # Collection / Scalar → (returns values for multiple rows)
    col_mean = m((:a, :b), :mean)
    @test length(col_mean) == 2

    # Colon / Scalar → all rows for one stat
    all_mean = m(:, :mean)
    @test length(all_mean) == 3

    # Tuple / Tuple → sub-MomentArray
    sub = m((:a, :b), (:mean, :p50))
    @test sub isa TFI.MomentArray
    @test size(sub) == (2, 2)
    @test m(:a, :mean) == m(:a, :mean)
    @test m(:b, :p50)  == m(:b, :p50)

end

@testset "ReportedMoment defaults" begin
    rm = ReportedMoment(var=:C)
    @test rm.var    == :C
    @test rm.sample == :uncond
    @test rm.stat   == :mean
end

@testset "StatisticsSet — construction and callable" begin

    # Construct from NamedTuple of functions
    ss = TFI.StatisticsSet((mean = mean, std = std))
    @test TFI.varnames(ss) == (:mean, :std)
    @test length(ss) == 2

    # Callable: applies each function to args and returns a NamedTuple
    x = [1.0, 2.0, 3.0, 4.0, 5.0]
    result = ss(x)
    @test result isa NamedTuple
    @test result.mean ≈ mean(x)
    @test result.std  ≈ std(x)

    # Construct from Tuple form
    ss2 = TFI.StatisticsSet((:p50, :p90), (x -> median(x), x -> quantile(x, 0.9)))
    @test TFI.varnames(ss2) == (:p50, :p90)
    result2 = ss2(x)
    @test result2.p50 ≈ median(x)

    # Identity: StatisticsSet(ss) returns the same object
    @test TFI.StatisticsSet(ss) === ss

end

@testset "StatisticsSet — merge" begin

    ss1 = TFI.StatisticsSet((mean = mean,))
    ss2 = TFI.StatisticsSet((std = std,))
    ss3 = merge(ss1, ss2)

    @test TFI.varnames(ss3) == (:mean, :std)
    x = [1.0, 2.0, 3.0]
    r = ss3(x)
    @test r.mean ≈ mean(x)
    @test r.std  ≈ std(x)

end

@testset "SimulationPanel — property delegation and metadata" begin

    data = StructArray([(C = 1.0*i, K = 2.0*i) for _ in 1:4, i in 1:10])
    panel = SimulationPanel(data, (:K,))

    # propertynames delegates to the underlying StructArray
    @test :C in propertynames(panel)
    @test :K in propertynames(panel)

    # getproperty delegation
    @test panel.C === data.C
    @test panel.K === data.K

    # endog_states / states
    @test TFI.endog_states(panel) == (:K,)
    @test TFI.states(panel) == (:exst, :K)

    # size matches underlying data
    @test size(panel) == (4, 10)

end

@testset "SimulationErrorPanel — property delegation" begin

    data = StructArray([(euler = 1e-5*i,) for _ in 1:4, i in 1:10])
    ep = TFI.SimulationErrorPanel(data)

    @test :euler in propertynames(ep)
    @test ep.euler === data.euler
    @test size(ep) == (4, 10)

end

@testset "MomentArray — show" begin

    data = reshape(Float64.(1:6), 2, 3)
    m = TFI.MomentArray(data, (:a, :b), (:mean, :std, :p90))

    # Compact show (no MIME)
    io = IOBuffer()
    show(io, m)
    s = String(take!(io))
    @test occursin("MomentArray", s)

    # text/plain show
    io2 = IOBuffer()
    show(io2, MIME"text/plain"(), m)
    s2 = String(take!(io2))
    @test occursin("MomentArray", s2)
    @test occursin("mean", s2)

end

@testset "compare_economies" begin

    # Build two minimal describe() outputs using synthetic data
    N_runs, N_periods = 3, 15
    Random.seed!(1)

    function _make_panel(scale)
        data = StructArray([(C = scale * (1 + 0.1*randn()), K = scale * (2 + 0.1*randn()))
                             for _ in 1:N_runs, _ in 1:N_periods])
        SimulationPanel(data, (:K,))
    end

    p1 = _make_panel(1.0)
    p2 = _make_panel(2.0)

    d1 = describe(p1)
    d2 = describe(p2)

    moment_defs = (
        mean_C = ReportedMoment(var=:C),
        mean_K = ReportedMoment(var=:K, stat=:mean),
    )
    economies = (baseline = d1.moments, scaled = d2.moments)

    table = compare_economies(moment_defs, economies)

    @test table isa TFI.MomentArray
    @test size(table) == (2, 2)              # 2 moments × 2 economies
    @test row_names(table) == (:mean_C, :mean_K)
    @test col_names(table) == (:baseline, :scaled)

    # Scaled economy has ~2× larger moments
    @test table(:mean_C, :scaled) > table(:mean_C, :baseline)
    @test table(:mean_K, :scaled) > table(:mean_K, :baseline)

end

# Lightweight model struct for parameter comparison tests
struct _TwoParamModel
    alpha::Float64
    beta::Float64
    label::String   # non-numeric: must be ignored by _distinct_parameters
end

@testset "MomentArray — vcat" begin
    m1 = TFI.MomentArray([1.0 2.0; 3.0 4.0], (:a, :b), (:x, :y))
    m2 = TFI.MomentArray([5.0 6.0],           (:c,),    (:x, :y))

    m = vcat(m1, m2)

    @test m isa TFI.MomentArray
    @test size(m) == (3, 2)
    @test row_names(m) == (:a, :b, :c)
    @test col_names(m) == (:x, :y)
    @test m[1, :] == [1.0, 2.0]
    @test m[3, :] == [5.0, 6.0]

    # Three-way vcat
    m3  = TFI.MomentArray([7.0 8.0], (:d,), (:x, :y))
    m123 = vcat(m1, m2, m3)
    @test size(m123) == (4, 2)
    @test row_names(m123) == (:a, :b, :c, :d)
end

@testset "MomentArray — hcat" begin
    m1 = TFI.MomentArray([1.0 2.0; 3.0 4.0],        (:a, :b), (:x, :y))
    m2 = TFI.MomentArray(reshape([5.0, 6.0], 2, 1),  (:a, :b), (:z,))

    m = hcat(m1, m2)

    @test m isa TFI.MomentArray
    @test size(m) == (2, 3)
    @test row_names(m) == (:a, :b)
    @test col_names(m) == (:x, :y, :z)
    @test m[:, 1] == [1.0, 3.0]
    @test m[:, 3] == [5.0, 6.0]

    # Three-way hcat
    m3   = TFI.MomentArray(reshape([7.0, 8.0], 2, 1), (:a, :b), (:w,))
    m123 = hcat(m1, m2, m3)
    @test size(m123) == (2, 4)
    @test col_names(m123) == (:x, :y, :z, :w)
end

@testset "_distinct_parameters" begin
    mdl1 = _TwoParamModel(0.3, 0.99, "base")
    mdl2 = _TwoParamModel(0.4, 0.99, "alt")
    models = (baseline = mdl1, alt = mdl2)

    params = TFI._distinct_parameters(models)

    @test :alpha in params      # differs across models
    @test :beta  ∉ params       # identical across models
    @test :label ∉ params       # non-numeric, always excluded

    # All numeric parameters identical → empty
    mdl3 = _TwoParamModel(0.3, 0.99, "dup")
    same_models = (baseline = mdl1, dup = mdl3)
    @test isempty(TFI._distinct_parameters(same_models))
end

@testset "_compare_parameters" begin
    mdl1 = _TwoParamModel(0.3, 0.99, "base")
    mdl2 = _TwoParamModel(0.4, 0.98, "alt")
    models = (baseline = mdl1, alt = mdl2)

    # Explicit param list
    ptable = TFI._compare_parameters((:alpha, :beta), models)
    @test ptable isa TFI.MomentArray
    @test size(ptable) == (2, 2)
    @test row_names(ptable) == (:alpha, :beta)
    @test col_names(ptable) == (:baseline, :alt)
    @test ptable(:alpha, :baseline) ≈ 0.3
    @test ptable(:alpha, :alt)      ≈ 0.4
    @test ptable(:beta,  :baseline) ≈ 0.99
    @test ptable(:beta,  :alt)      ≈ 0.98

    # Auto-detect distinct params via single-arg overload (both differ here)
    ptable2 = TFI._compare_parameters(models)
    @test row_names(ptable2) == (:alpha, :beta)
    @test ptable2(:alpha, :alt) ≈ 0.4
end

@testset "compare_economies — with models" begin
    N_runs, N_periods = 3, 15
    Random.seed!(2)

    function _make_panel2(scale)
        data = StructArray([(C = scale*(1 + 0.1*randn()), K = scale*(2 + 0.1*randn()))
                             for _ in 1:N_runs, _ in 1:N_periods])
        SimulationPanel(data, (:K,))
    end

    d1 = describe(_make_panel2(1.0))
    d2 = describe(_make_panel2(2.0))
    economies = (baseline = d1.moments, scaled = d2.moments)

    moment_defs = (
        mean_C = ReportedMoment(var=:C),
        mean_K = ReportedMoment(var=:K, stat=:mean),
    )

    # alpha differs, beta is the same → only alpha appears as a parameter row
    models = (
        baseline = _TwoParamModel(0.3, 0.99, "base"),
        scaled   = _TwoParamModel(0.4, 0.99, "alt"),
    )

    table = compare_economies(moment_defs, economies; models=models)

    # 1 param row + 2 moment rows, 2 economy columns
    @test table isa TFI.MomentArray
    @test size(table) == (3, 2)
    @test col_names(table) == (:baseline, :scaled)

    all_rows = row_names(table)
    @test :alpha  in all_rows
    @test :mean_C in all_rows
    @test :mean_K in all_rows

    # Parameter rows are prepended before moment rows
    @test all_rows[1]   == :alpha
    @test all_rows[2:3] == (:mean_C, :mean_K)

    # Parameter values are correct
    @test table(:alpha, :baseline) ≈ 0.3
    @test table(:alpha, :scaled)   ≈ 0.4

    # Moment values are unchanged vs. no-models call
    plain = compare_economies(moment_defs, economies)
    @test table(:mean_C, :baseline) ≈ plain(:mean_C, :baseline)
    @test table(:mean_K, :scaled)   ≈ plain(:mean_K, :scaled)
end
