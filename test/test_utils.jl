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
