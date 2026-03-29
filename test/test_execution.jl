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

    # Use the shared DummyModel obj from runtests.jl
    obj_copy = deepcopy(obj)
    result = main_execute!(obj_copy; maxiter=3, logging=COMPLETE,
                           track_history=(:Vfct, (:v,)))

    @test haskey(result, :history)
    hist = result.history
    @test !isnothing(hist)
    @test hist isa StructArrays.StructVector

    # Should have one entry per iteration (at most maxiter)
    @test length(hist) <= 3
    @test length(hist) >= 1

    # Each entry has the requested field :v
    @test hasproperty(hist, :v)

    # Each snapshot is an array of grid values (one per grid point)
    snap = hist.v[1]
    @test snap isa AbstractVector
    @test length(snap) == TFI.npt(obj_copy.Vfct)

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
