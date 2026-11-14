# Solution files (JLD2 extension) and model metadata. Uses ngm, obj from test_ngm.jl and
# ForecastNGM from test_forecasts.jl.
using JLD2

struct MetaTestModel <: AbstractModel
    a::Float64
    f::Function
    big::Vector{Float64}
    small::Vector{Float64}
    nt::NamedTuple
end

@testset "model_metadata" begin
    @test model_metadata(ngm) == (α = ngm.α, β = ngm.β, δ = ngm.δ, γ = ngm.γ, σ = ngm.σ, ρ = ngm.ρ, N_A = ngm.N_A)
    @test model_from_metadata(NeoclassicalGrowthModel, model_metadata(ngm)) == ngm
    fm = ForecastNGM(NeoclassicalGrowthModel())
    @test model_metadata(fm) == (base = model_metadata(fm.base),)
    @test TFI.flatten_metadata(model_metadata(fm), "model")["model.base.β"] == string(fm.base.β)
    # derived objects are left out
    m = MetaTestModel(1.0, sin, zeros(10_000), [1.0, 2.0], (x = :s, g = cos))
    @test model_metadata(m) == (a = 1.0, small = [1.0, 2.0], nt = (x = :s,))
end

@testset "save_solution / load_solution" begin
    path = joinpath(mktempdir(), "sol.jld2")
    @test save_solution(path, obj; provenance = Dict("git.commit" => "abc123")) == path

    same(a, b) = a.Pfct == b.Pfct && a.Vfct == b.Vfct && a.Tfct == b.Tfct &&
                 a.ssg == b.ssg && isapprox(a.exogenv, b.exogenv) &&
                 a.Pfct_bounds == b.Pfct_bounds && a.Tfct_bounds == b.Tfct_bounds &&
                 (isnothing(a.Ffct) ? isnothing(b.Ffct) : a.Ffct == b.Ffct)

    o1 = load_solution(path)                                  # exact model
    @test same(o1, obj) && o1.model == ngm
    @test solution_id(o1) == solution_id(obj)                 # the id survives the file
    o2 = @test_logs load_solution(path, ngm)                  # supplied model, matching metadata
    @test same(o2, obj)
    @test_logs (:warn, r"metadata differ") match_mode = :any load_solution(path, NeoclassicalGrowthModel(β = 0.95))
    o3 = load_solution(path, NeoclassicalGrowthModel)         # model_from_metadata
    @test same(o3, obj) && o3.model == ngm

    # the restored solution evaluates and simulates like the original
    pt = (exst = 2, K = 3.0)
    @test o1.Pfct(pt) == obj.Pfct(pt)
    @test simulate(o1, stv; sz = 100 => 2)[1].data == simulate(obj, stv; sz = 100 => 2)[1].data

    # forecasting functions and a one-state exogenous environment
    fm = ForecastNGM(NeoclassicalGrowthModel(N_A = 1))
    objf, _ = main_setup(fm)
    main_execute!(objf; maxiter = 300, tol = -7, logging = COMPLETE)
    pathf = save_solution(joinpath(mktempdir(), "solf.jld2"), objf)
    of = load_solution(pathf)
    @test same(of, objf) && varnames(of.Ffct) == (:Cnow,)
    @test solution_id(of) == solution_id(objf) != solution_id(obj)
    @test same(load_solution(pathf, fm), objf)

    # stored layout: plain arrays over the grid axes, Tfct with a next-state axis
    jldopen(path, "r") do io
        @test io["meta/schema_version"] == 1
        @test io["grid/__axes"] == ["exst", "K"]
        @test size(io["Pfct/C"]) == (ngm.N_A, 20)
        @test size(io["Tfct/K"]) == (ngm.N_A, 20, ngm.N_A)
        @test io["meta/model/β"] == ngm.β
        @test io["meta/solution_id"] == solution_id(obj)
        @test io["meta/provenance/git.commit"] == "abc123"
    end
end
