# Portable data (src/portable.jl): solutions, steady states and interpolants as plain nested
# NamedTuples and back, without any storage. Uses obj, stv and ngm from test_ngm.jl and ForecastNGM
# from test_forecasts.jl.

same_solution(a, b) = a.Pfct == b.Pfct && a.Vfct == b.Vfct && a.Tfct == b.Tfct &&
    a.ssg == b.ssg && isapprox(a.exogenv, b.exogenv) &&
    a.Pfct_bounds == b.Pfct_bounds && a.Tfct_bounds == b.Tfct_bounds &&
    (isnothing(a.Ffct) ? isnothing(b.Ffct) : a.Ffct == b.Ffct)

# whether every leaf of a nested NamedTuple is plain data
plain(x::NamedTuple) = all(plain, values(x))
plain(x) = x isa Union{Number, AbstractString} || x isa AbstractArray{<:Union{Number, AbstractString}}

@testset "portable data" begin
    @testset "solution" begin
        d = to_portable(obj)
        @test plain(d)
        @test d.meta.schema_version == 2 && d.meta.solution_id == solution_id(obj) && d.meta.model.β == ngm.β
        @test d.grid.__axes == ["exst", "K"] && size(d.Pfct.C) == (ngm.N_A, 20) && size(d.Tfct.K) == (ngm.N_A, 20, ngm.N_A)
        @test d.Pfct.__bounds.__names == ["C"] && d.Pfct.__bounds.C == [0.0, Inf]

        o1 = from_portable(ModelSolution, d; model = NeoclassicalGrowthModel)   # model_from_metadata
        @test same_solution(o1, obj) && o1.model == ngm && solution_id(o1) == solution_id(obj)
        o2 = @test_logs from_portable(ModelSolution, d; model = ngm)
        @test same_solution(o2, obj)
        @test_logs (:warn, r"metadata differ") match_mode = :any from_portable(ModelSolution, d; model = NeoclassicalGrowthModel(β = 0.95))
        @test_throws ArgumentError from_portable(ModelSolution, d)
        @test_throws ArgumentError from_portable(ModelSolution, merge(d, (meta = merge(d.meta, (schema_version = 1,)),)); model = ngm)
        @test same_solution(from_portable(ModelSolution, merge(d, (extra = [1.0],)); model = ngm), obj)   # unknown keys ignored
    end

    @testset "forecasting functions and one exogenous state" begin
        fm = ForecastNGM(NeoclassicalGrowthModel(N_A = 1))
        objf, _ = main_setup(fm)
        main_execute!(objf; maxiter = 300, tol = -7, logging = COMPLETE)
        of = from_portable(ModelSolution, to_portable(objf); model = fm)
        @test same_solution(of, objf) && varnames(of.Ffct) == (:Cnow,)
    end

    @testset "steady state and interpolant" begin
        @test plain(to_portable(stv)) && plain(to_portable(obj.Vfct))
        @test from_portable(SteadyState, to_portable(stv); model = ngm) == stv
        @test from_portable(SteadyState, to_portable(stv); model = NeoclassicalGrowthModel) == stv
        @test from_portable(Interpolant, to_portable(obj.Vfct)) == obj.Vfct
    end

    @testset "bounds" begin
        TFI = TransitionFunctionIteration
        for b in (TFI.Unbounded(), TFI.LowerBounded(-1.0), TFI.UpperBounded(2.0), TFI.IntervalBounded(0.0, 3.0))
            @test TFI.interval_bound(TFI.bound_interval(b)) == b
        end
    end
end
