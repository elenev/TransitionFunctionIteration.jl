# TFI's results with OutputVault (the package extension): portable solution, steady-state and
# interpolant files, tables, save_result and load_result, IRFs saved with OutputVault 0.3,
# moments from saved runs. Uses ngm, obj, stv, simseries, errs, stats, errstats from test_ngm.jl
# and ForecastNGM from test_forecasts.jl.
using OutputVault, Parquet2, JLD2, DataAPI, Tables
const OVExt = Base.get_extension(TransitionFunctionIteration, :TransitionFunctionIterationOutputVaultExt)

struct MetaTestModel <: AbstractModel
    a::Float64
    f::Function
    big::Vector{Float64}
    small::Vector{Float64}
    nt::NamedTuple
end

struct ParametersModel <: AbstractModel
    a::Float64
    n::Int
    flag::Bool
    label::String
    grid::Vector{Float64}
end

ov_git(dir, args...) = read(`git -C $dir -c user.name=Test -c user.email=test@example.com $args`, String)

# f(root) with an OutputVault configuration over a temporary root, provenance from a throwaway
# repository made the active project
function with_vault(f)
    repo = mktempdir()
    root = mktempdir()
    ov_git(repo, "init", "-q", "-b", "main")
    write(joinpath(repo, "Project.toml"), "name = \"Throwaway\"\n")
    ov_git(repo, "add", "Project.toml")
    ov_git(repo, "commit", "-q", "-m", "init")
    old = Base.ACTIVE_PROJECT[]
    Base.ACTIVE_PROJECT[] = joinpath(repo, "Project.toml")
    try
        with_config(() -> f(root), Config(; output_root = root, table_format = OutputVault.Formats.ParquetFormat()))
    finally
        Base.ACTIVE_PROJECT[] = old
        GC.gc(); GC.gc()
        rm(root; recursive = true, force = true)
        rm(repo; recursive = true, force = true)
    end
end


@testset "OutputVault extension" begin
    @testset "model metadata" begin
        @test to_portable(ngm) == (α = ngm.α, β = ngm.β, δ = ngm.δ, γ = ngm.γ, σ = ngm.σ, ρ = ngm.ρ, N_A = ngm.N_A)
        @test from_portable(NeoclassicalGrowthModel, to_portable(ngm)) == ngm
        fm = ForecastNGM(NeoclassicalGrowthModel())
        @test to_portable(fm) == (base = to_portable(fm.base),)
        m = MetaTestModel(1.0, sin, zeros(10_000), [1.0, 2.0], (x = :s, g = cos))   # derived objects left out
        @test to_portable(m) == (a = 1.0, small = [1.0, 2.0], nt = (x = :s,))
    end

    @testset "default parameters are the scalar real and Boolean fields" begin
        p = collect_parameters(ParametersModel(0.5, 3, true, "x", [1.0, 2.0]))
        @test row_names(p) == (:a, :n, :flag) && col_names(p) == (:value,)
        @test vec(p.data) == [0.5, 3.0, 1.0]
    end

    @testset "portable files" begin
        dir = mktempdir()
        path = write_object_file(joinpath(dir, "sol.jld2"), obj)
        o1 = read_object_file(path, ModelSolution; context = NeoclassicalGrowthModel)   # from_portable(M, parameters)
        @test same_solution(o1, obj) && o1.model == ngm && solution_id(o1) == solution_id(obj)
        o2 = @test_logs read_object_file(path, ModelSolution; context = ngm)            # matching parameters
        @test same_solution(o2, obj)
        @test_logs (:warn, r"parameters differ") match_mode = :any read_object_file(path, ModelSolution; context = NeoclassicalGrowthModel(β = 0.95))
        @test_throws ArgumentError read_object_file(path, ModelSolution)
        pt = (exst = 2, K = 3.0)
        @test o1.Pfct(pt) == obj.Pfct(pt)
        @test simulate(o1, stv; sz = 100 => 2)[1].data == simulate(obj, stv; sz = 100 => 2)[1].data

        jldopen(path, "r") do io
            @test io["meta/schema_version"] == 2
            @test io["Pfct/__bounds/__names"] == ["C"] && io["Pfct/__bounds/C"] == [0.0, Inf]
            @test io["grid/__axes"] == ["exst", "K"]
            @test size(io["Pfct/C"]) == (ngm.N_A, 20)
            @test size(io["Tfct/K"]) == (ngm.N_A, 20, ngm.N_A)
            @test io["meta/model/β"] == ngm.β
            @test io["meta/solution_id"] == solution_id(obj)
            @test !haskey(io, "value")                                     # plain data, not the object
        end

        # forecasting functions and a one-state exogenous environment
        fm = ForecastNGM(NeoclassicalGrowthModel(N_A = 1))
        objf, _ = main_setup(fm)
        main_execute!(objf; maxiter = 300, tol = -7, logging = COMPLETE)
        of = read_object_file(write_object_file(joinpath(dir, "solf.jld2"), objf), ModelSolution; context = fm)
        @test same_solution(of, objf) && varnames(of.Ffct) == (:Cnow,)

        spath = write_object_file(joinpath(dir, "stv.jld2"), stv)
        @test read_object_file(spath, SteadyState; context = ngm) == stv
        @test read_object_file(spath, SteadyState; context = NeoclassicalGrowthModel) == stv

        fpath = write_object_file(joinpath(dir, "fct.jld2"), obj.Vfct)
        @test read_object_file(fpath, Interpolant) == obj.Vfct

    end

    irf = simulate_irf(obj, simseries; shocks = (pos = 2, uncond = Int64[]), path_vars = (:K,), T = 5, N = 20)
    samples = describe(simseries; filters = (a_all = Returns(true),)).moments

    @testset "save_result and load_result" begin
        with_vault() do root
            rf = RunHandle("ngm"; branch = "main", experdef = "base")
            result = (; obj, stv, stats = samples, errstats, simseries, errs, irf)
            dir = save_result(rf, result)
            for f in ("native.jld2", ".provenance/git_commit.txt", "solution/obj.jld2", "solution/stv.jld2",
                      "params/params.parquet", "simulation/simseries.parquet", "simulation/errs.parquet",
                      "errstats/errstats.parquet", "stats/uncond.parquet", "stats/a_all.parquet",
                      "irf/ss/stats.parquet", "irf/ss/paths.parquet")
                @test isfile(joinpath(dir, f))
            end
            ds = Parquet2.Dataset(joinpath(dir, "simulation", "simseries.parquet"))
            @test DataAPI.metadata(ds, "tfi.kind") == "SimulationPanel"
            @test first(Tables.columnnames(Parquet2.Dataset(joinpath(dir, "errstats", "errstats.parquet")))) == :equation
            @test first(Tables.columnnames(Parquet2.Dataset(joinpath(dir, "params", "params.parquet")))) == :parameter

            sol = load_result(rf, :solution; native = true)
            @test same_solution(sol.obj, obj) && sol.stv == stv
            for model in (NeoclassicalGrowthModel, ngm)
                solp = load_result(rf, :solution; model)
                @test same_solution(solp.obj, obj) && solp.obj.model == ngm && solp.stv == stv
            end
            @test_throws ArgumentError load_result(rf, :solution)
            @test load_result(rf, :errstats) == errstats
            st = load_result(rf, :stats)
            @test keys(st) == (:a_all, :uncond) && st.uncond == samples.uncond && st.a_all == samples.a_all
            pars = load_result(rf, :params)
            @test pars == collect_parameters(obj)
            sim = load_result(rf, :simulation)
            @test sim.simseries.data == simseries.data && sim.simseries.provenance == simseries.provenance
            @test endog_states(sim.simseries) == endog_states(simseries)
            @test sim.errs.data == errs.data && sim.errs.provenance == errs.provenance
            r = load_result(rf, :irf)
            @test keys(r) == (:ss,)
            @test r.ss.means == irf.means && r.ss.stdevs == irf.stdevs && r.ss.paths.pos == irf.paths.pos
            @test r.ss.provenance == irf.provenance
            @test_throws ArgumentError load_result(rf, :other)

            # recursions saved later add to the run; the solution is untouched
            before = read(joinpath(dir, "solution", "obj.jld2"))
            save_result(rf, (; recursions = (; claims = obj.Vfct)))
            save_result(rf, (; recursions = (; welfare = obj.Vfct)))
            @test read(joinpath(dir, "solution", "obj.jld2")) == before
            rec = load_result(rf, :recursions)
            @test keys(rec) == (:claims, :welfare) && rec.claims == obj.Vfct
            pts = collect(points(obj.ssg))
            @test rec.welfare(pts) == obj.Vfct(pts)

            all_ = load_result(rf; model = ngm)
            @test keys(all_) == (:obj, :stv, :params, :simseries, :errs, :errstats, :stats, :irf, :recursions)
            @test all_.errstats == errstats

            # another starting point is another folder
            save_result(rf, (; irf); start = :boom)
            @test keys(load_result(rf, :irf)) == (:boom, :ss)
        end
    end

    @testset "IRFs saved with OutputVault 0.3" begin
        with_vault() do root
            rf = RunHandle("old"; branch = "main", experdef = "base")
            save_result(rf, (; obj))
            # the 0.3 layout: two table files in irf/
            dir = mkpath(joinpath(OutputVault.run_path(rf), "irf"))
            OutputVault.Formats.write_table(OutputVault.Formats.ParquetFormat(), joinpath(dir, "irf.parquet"), irf)
            OutputVault.Formats.write_table(OutputVault.Formats.ParquetFormat(), joinpath(dir, "irf_paths.parquet"), irf_paths(irf))
            r = load_result(rf, :irf)
            @test keys(r) == (:ss,) && r.ss.means == irf.means && r.ss.paths.pos == irf.paths.pos
            shape = OVExt.SHAPES.irf.query
            @test isempty(query_results([rf], shape).value)          # queries do not see the old files
            @test (@test_logs TFI.migrate_irf!(RunPattern(); move = false)) === nothing
            @test isfile(joinpath(OutputVault.run_path(rf), "irf", "irf.parquet"))
            TFI.migrate_irf!(RunPattern(); move = true)
            @test isfile(joinpath(OutputVault.run_path(rf), "irf", "ss", "stats.parquet"))
            @test load_result(rf, :irf).ss.means == irf.means
            @test !isempty(query_results([rf], shape).value)
            @test TFI.migrate_irf!(rf) == []
        end
    end

    @testset "moments from saved runs" begin
        with_vault() do root
            hi = NeoclassicalGrowthModel(σ = 2 * ngm.σ)
            runs = (base = RunHandle("base"; branch = "main", experdef = "v1"), hi = RunHandle("hi"; branch = "main", experdef = "v1"))
            save_result(runs.base, (; obj, stats = samples, irf))
            obj_hi = ModelSolution(; (f => getfield(obj, f) for f in fieldnames(typeof(obj)) if f != :model)..., model = hi)
            samples_hi = map(m -> MomentArray{row_names(m), col_names(m)}(2 .* m.data; rowkey = row_key(m)), samples)
            save_result(runs.hi, (; obj = obj_hi, stats = samples_hi, irf))

            moments = (mean_K = SimulationMoment(:K),
                       std_C_all = SimulationMoment(:C; sample = :a_all, stat = :std),
                       C = :C,
                       resp = ImpulseMoment(:K, :pos, 3; transformation = TransitionFunctionIteration.Transformations.DiffLevel()),
                       sig = ParameterMoment(:σ),
                       ratio = @CalculatedMoment SimulationMoment(:K) / ParameterMoment(:σ))
            in_memory = (base = EconomyResult(; stats = samples, irf, params = collect_parameters(obj)),
                         hi = EconomyResult(; stats = samples_hi, irf, params = collect_parameters(obj_hi)))

            # saved runs and results in memory are read the same way, with every option
            pattern = RunPattern(; experdef = "v1")
            for opts in ((;), (; base = :hi), (; compare_transform = false), (; parameters = :distinct),
                         (; parameters = (:σ, :β)))
                expected = compare_economies(moments, in_memory; opts...)
                @test isequal(compare_economies(moments, runs; opts...), expected)
                @test isequal(compare_economies(moments, pattern; opts...), expected)
            end
            @test row_names(compare_economies(moments, runs; parameters = :distinct))[1] == :σ
            @test isequal(stack_moments(compare_economies(moments, runs; parameters = :distinct)),
                          stack_moments(compare_economies(moments, in_memory; parameters = :distinct)))

            # a number a run lacks is a missing cell, in memory and saved alike
            lacking = (; moments..., p99 = SimulationMoment(:K; stat = :p99), late = ImpulseMoment(:K, :pos, 10_000))
            expected = compare_economies(lacking, in_memory)
            @test ismissing(expected(:p99, :base)) && ismissing(expected(:late, :hi))
            @test isequal(compare_economies(lacking, runs), expected)

            # a run with none of the numbers is left out with a warning
            save_result(RunHandle("bare"; branch = "main", experdef = "v1"), (; obj))
            t = @test_logs (:warn, r"left out: bare") match_mode = :any compare_economies((; m = SimulationMoment(:K)), pattern)
            @test col_names(t) == (:base, :hi)
            @test_throws ArgumentError compare_economies((; m = SimulationMoment(:K)), (; bare = RunHandle("bare"; branch = "main", experdef = "v1")))
            # variants that repeat across experdefs
            save_result(RunHandle("base"; branch = "main", experdef = "v2"), (; obj, stats = samples))
            err = try compare_economies((; m = SimulationMoment(:K)), RunPattern(; experdef = "*", variant = "base")); nothing catch e; e end
            @test err isa ArgumentError && occursin("narrower pattern", err.msg)

            p = query_results(collect(values(runs)), OVExt.SHAPES.params.query; parameter = :σ)
            @test p.value == [ngm.σ, hi.σ]
        end
    end
end
