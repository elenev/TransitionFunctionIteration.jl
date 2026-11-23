@testset "value() — passthrough and Dual extraction" begin
    @test TFI.value(1.5) === 1.5
    @test TFI.value(0.0) === 0.0
    d = ForwardDiff.Dual(3.0, 1.0)
    @test TFI.value(d) === 3.0
end

@testset "flatten_structarray" begin
    # Build a StructArray whose columns hold SVector{2} elements
    sa = StructArray((x = [SVector(1.0,2.0), SVector(3.0,4.0)],
                      y = [SVector(5.0,6.0), SVector(7.0,8.0)]))
    flat = TFI.flatten_structarray(sa)

    # After flattening, the SVector dimension is prepended: shape (2,) → (2, 2)
    @test size(flat.x) == (2, 2)
    @test size(flat.y) == (2, 2)

    # Values are preserved (reinterpret, not copied)
    @test flat.x[1,1] == 1.0
    @test flat.x[2,1] == 2.0
    @test flat.x[1,2] == 3.0
end

@testset "StaticStructArray — round-trip" begin
    sa = StructArray((a = [1.0, 2.0, 3.0], b = [4.0, 5.0, 6.0]))
    ssa = TFI.StaticStructArray(sa)

    @test ssa.a == SVector(1.0, 2.0, 3.0)
    @test ssa.b == SVector(4.0, 5.0, 6.0)
end

@testset "StaticStructArray — from StaticArray of NamedTuples" begin
    # @generated overload: takes a StaticArray{Sz, NamedTuple}
    sv = SVector((a=1.0, b=10.0), (a=2.0, b=20.0), (a=3.0, b=30.0))
    ssa = TFI.StaticStructArray(sv)

    @test ssa.a == SVector(1.0, 2.0, 3.0)
    @test ssa.b == SVector(10.0, 20.0, 30.0)
end

@testset "merge_structarrays" begin
    sa1 = StructArray((x = [1.0, 2.0, 3.0],))
    sa2 = StructArray((y = [4.0, 5.0, 6.0],))
    sa3 = StructArray((z = [7.0, 8.0, 9.0],))

    # Two arrays
    merged2 = TFI.merge_structarrays(sa1, sa2)
    @test :x in propertynames(merged2)
    @test :y in propertynames(merged2)
    @test merged2.x == sa1.x
    @test merged2.y == sa2.y
    @test length(merged2) == 3

    # Three arrays
    merged3 = TFI.merge_structarrays(sa1, sa2, sa3)
    @test :z in propertynames(merged3)
    @test merged3.z == sa3.z
end

@testset "fmap — tuple and NamedTuple" begin
    x = 3.0

    # Tuple of functions
    result = TFI.fmap((sin, cos, sqrt), x)
    @test result isa Tuple
    @test result[1] ≈ sin(x)
    @test result[2] ≈ cos(x)
    @test result[3] ≈ sqrt(x)

    # NamedTuple of functions
    nt = (s = sin, c = cos)
    result_nt = TFI.fmap(nt, x)
    @test result_nt isa NamedTuple
    @test result_nt.s ≈ sin(x)
    @test result_nt.c ≈ cos(x)

    # Multiple args
    result2 = TFI.fmap((+, *), 2.0, 3.0)
    @test result2[1] ≈ 5.0
    @test result2[2] ≈ 6.0
end

@testset "select — StructArray column subset" begin
    sa = StructArray((a = [1.0, 2.0], b = [3.0, 4.0], c = [5.0, 6.0]))

    # Select two out of three fields
    sub = TFI.select(sa, Val((:a, :c)))
    @test propertynames(sub) == (:a, :c)
    @test sub.a == sa.a
    @test sub.c == sa.c
    @test !(:b in propertynames(sub))

    # Single field
    sub1 = TFI.select(sa, Val((:b,)))
    @test propertynames(sub1) == (:b,)
    @test sub1.b == sa.b
end

@testset "unrolled_any — NamedTuple and Tuple" begin
    @test !TFI.unrolled_any(isnan, (;))                    # no fields: false
    @test !TFI.unrolled_any(isnan, (a = 1.0, b = 2.0))
    @test TFI.unrolled_any(isnan, (a = NaN, b = 2.0))      # NaN in first field
    @test TFI.unrolled_any(isnan, (a = 1.0, b = NaN))      # NaN in last field
    @test TFI.unrolled_any(isnan, (1.0, NaN, 3))           # plain Tuple, mixed types
    @test !TFI.unrolled_any(isnan, (1, 2.5, true))         # Int, Float64, Bool
    @test TFI.unrolled_any(>(2), (a = 1, b = 3))           # arbitrary predicate
end

@testset "unrolled_any — type stability on a wide NamedTuple" begin
    wide     = NamedTuple{ntuple(i -> Symbol("v", i), 12)}(ntuple(float, 12))
    wide_nan = merge(wide, (v12 = NaN,))
    @test (@inferred TFI.unrolled_any(isnan, wide)) === false
    @test (@inferred TFI.unrolled_any(isnan, wide_nan)) === true
end

@testset "with_blas_threads and threaded_foreach" begin
    original = TFI.BLAS.get_num_threads()
    TFI.BLAS.set_num_threads(2)
    blas() = TFI.BLAS.get_num_threads()
    try
        # limited inside, restored after, also when f throws
        @test TFI.with_blas_threads(blas, 1) == 1
        @test blas() == 2
        @test_throws ErrorException TFI.with_blas_threads(() -> error("boom"), 1)
        @test blas() == 2

        # every element is visited, in parallel with one BLAS thread or in order with BLAS untouched
        for use_threads in (true, false)
            seen, threads = zeros(Int, 8), zeros(Int, 8)
            TFI.threaded_foreach(1:8; use_threads) do i
                seen[i] += 1
                threads[i] = blas()
            end
            @test seen == ones(Int, 8)
            @test threads == fill(use_threads ? 1 : 2, 8)
            @test blas() == 2
        end
    finally
        TFI.BLAS.set_num_threads(original)
    end
end
