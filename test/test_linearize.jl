# Uses `obj`, `stv`, `simseries` from the NGM, which are in global scope after
# test_ngm.jl is included. test_ngm.jl must precede this file in runtests.jl.

@testset "linearize(ModelSolution, SteadyState)" begin
    vcoefs, tcoefs = linearize(obj, stv)

    @test vcoefs isa TFI.MomentArray
    @test tcoefs isa TFI.MomentArray

    # NGM: 1 value var (C), regressors = (constant, A, K)
    @test size(vcoefs) == (1, 3)
    @test row_names(vcoefs) == (:C,)
    @test col_names(vcoefs) == (:constant, :A, :K)

    # NGM: 1 transition var (K′), regressors = (constant, A, K, A′)
    @test size(tcoefs) == (1, 4)
    @test row_names(tcoefs) == (:K′,)
    @test col_names(tcoefs) == (:constant, :A, :K, :A′)
end

@testset "linearize(simseries, SteadyState)" begin
    vcoefs, tcoefs = linearize(simseries, stv)

    @test vcoefs isa TFI.MomentArray
    @test tcoefs isa TFI.MomentArray
    @test size(vcoefs) == (1, 3)
    @test row_names(vcoefs) == (:C,)
    @test col_names(vcoefs) == (:constant, :A, :K)
    @test size(tcoefs) == (1, 4)
    @test row_names(tcoefs) == (:K′,)
    @test col_names(tcoefs) == (:constant, :A, :K, :A′)
end

@testset "linearize — loglin=false" begin
    # Level deviations: same structure, different coefficient values
    vl, tl = linearize(obj, stv; loglin=false)
    @test vl isa TFI.MomentArray
    @test tl isa TFI.MomentArray
    @test size(vl) == (1, 3)
    @test size(tl) == (1, 4)

    vl2, tl2 = linearize(simseries, stv; loglin=false)
    @test vl2 isa TFI.MomentArray
    @test tl2 isa TFI.MomentArray
    @test size(vl2) == (1, 3)
    @test size(tl2) == (1, 4)
end
