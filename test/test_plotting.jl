@testset "Plotting recipe" begin

    # RecipesBase.apply_recipe runs the @recipe body and returns Vector{RecipeData}.
    # Each element has:
    #   .args            — the data tuple (e.g. (xvec, yvec) for 1D)
    #   .plotattributes  — attributes set by --> (suggest) or := (force)
    # This lets us verify recipe outputs without loading a Plots backend.

    # 1D: all dims free, identity transform (default)
    rd = first(RecipesBase.apply_recipe(Dict{Symbol,Any}(), fct1, :y, (;)))
    @test rd.args[1] ≈ xgrid
    @test rd.args[2] ≈ 2 .* xgrid
    @test rd.plotattributes[:title]  == "y"
    @test rd.plotattributes[:xlabel] == "x"
    @test get(rd.plotattributes, :legend, false) == false

    # transform kwarg: consumed by the recipe, applied to the output vector
    rd2 = first(RecipesBase.apply_recipe(Dict{Symbol,Any}(:transform => x -> -x), fct1, :y, (;)))
    @test rd2.args[2] ≈ -2 .* xgrid
    @test !haskey(rd2.plotattributes, :transform)   # must be consumed, not left behind

    # 2D with one dim fixed: fixedstates removes y from free dims → 1D plot
    rd3 = first(RecipesBase.apply_recipe(Dict{Symbol,Any}(), fct2, :f, (; y=20.0)))
    @test rd3.args[1] ≈ xgrid2
    @test rd3.args[2] ≈ xgrid2 .+ 20.0
    @test rd3.plotattributes[:xlabel] == "x"

    # 2D: both dims free → contourf with (x1vec, x2vec, fmat)
    rd4 = first(RecipesBase.apply_recipe(Dict{Symbol,Any}(), fct2, :f, (;)))
    @test rd4.args[1] ≈ xgrid2
    @test rd4.args[2] ≈ ygrid2
    @test size(rd4.args[3]) == (3, 3)
    @test rd4.args[3][2, 2] ≈ 22.0               # f(2.0, 20.0) = 22
    @test rd4.plotattributes[:seriestype] == :contourf
    @test rd4.plotattributes[:title]      == "f"
    @test rd4.plotattributes[:xlabel]     == "x"
    @test rd4.plotattributes[:ylabel]     == "y"

    # >2D: recipe should throw
    ssg3   = StateSpaceGrid(; x=collect(1.0:2.0), y=collect(1.0:2.0), z=collect(1.0:2.0))
    fvals3 = StructArray([(f=1.0,) for _ in 1:2, _ in 1:2, _ in 1:2])
    fct3   = Interpolant(ssg3, fvals3)
    @test_throws ErrorException RecipesBase.apply_recipe(Dict{Symbol,Any}(), fct3, :f, (;))

end
