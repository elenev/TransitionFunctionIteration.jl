@testset "normalize_next_period_values — default passthrough" begin

    # The default implementation returns the NextPeriod unchanged
    expectations = (e = 1.5,)
    continuations = StructArray([(c = 0.9,), (c = 1.0,), (c = 1.1,)])
    prnext = SVector(1/3, 1/3, 1/3)
    np = NextPeriod(; expectations, continuations, prnext)

    # Default passthrough for any model
    np2 = TFI.normalize_next_period_values(np, DummyModel(), SolveMode())

    @test np2.expectations === np.expectations
    @test np2.continuations === np.continuations

end

@testset "track_history — captures per-iteration snapshots" begin

    # track_history is a NamedTuple of (fct_sym, var_name, grid_idx) triples.
    # The config key becomes the field name in the returned history StructVector.
    obj_copy = deepcopy(obj)
    result = main_execute!(obj_copy; maxiter=3, logging=COMPLETE,
                           track_history=(v_trace = (:Vfct, :v, 1),))

    @test haskey(result, :history)
    hist = result.history
    @test !isnothing(hist)
    @test hist isa StructArrays.StructVector

    # One slot pre-allocated per iteration (length == maxiter)
    @test length(hist) == 3

    # The config key :v_trace becomes the field name
    @test hasproperty(hist, :v_trace)

    # Each snapshot is a scalar (value of :v at grid index 1 for that iteration)
    @test hist.v_trace[1] isa Number

end

@testset "main_execute! — return type contains exec fields" begin

    obj_copy = deepcopy(obj)
    result = main_execute!(obj_copy; maxiter=2, logging=COMPLETE)

    @test result isa NamedTuple
    @test haskey(result, :maxiter)
    @test haskey(result, :iter)
    @test haskey(result, :history)
    @test isnothing(result.history)   # no track_history set

end
