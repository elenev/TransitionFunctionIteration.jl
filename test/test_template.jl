@testset "generate_template" begin
    path = joinpath(mktempdir(), "mymodel.jl")
    redirect_stdout(devnull) do
        generate_template("MyModel", path)
    end
    src = read(path, String)
    ex = Meta.parseall(src)
    @test !any(e -> e isa Expr && e.head in (:error, :incomplete), ex.args)       # valid Julia
    @test occursin("struct MyModel <: AbstractModel", src)
    # the five required functions are live code, the optional ones are in block comments
    for f in ("Exogenv(mdl::MyModel)", "define_endog_grids(mdl::MyModel", "steadystate(mdl::MyModel",
              "compute_next_period_values(next, mdl::MyModel", "calc_eqm(state, policy, mdl::MyModel")
        @test occursin("function " * f, split(src, "#=")[1])
    end
    @test occursin("forecasts_from_policy(state, policy, mdl::MyModel)", src)
    @test !occursin("forecasts_from_policy", split(src, "#=")[1])
end
