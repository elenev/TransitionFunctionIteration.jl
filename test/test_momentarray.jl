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
