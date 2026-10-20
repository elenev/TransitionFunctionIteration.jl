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

    # History is trimmed to however many iterations actually ran (may be less
    # than maxiter if convergence happens early).
    @test length(hist) == result.iter

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

@testset "main_execute! — damping modes" begin

    # Same first iteration: both modes blend against the initial guess
    o1, o2 = deepcopy(obj), deepcopy(obj)
    main_execute!(o1; maxiter=1, damp=0.5, damping=:damped_iterate, logging=COMPLETE)
    main_execute!(o2; maxiter=1, damp=0.5, damping=:previous_solve, logging=COMPLETE)
    @test TFI.coefficients(o1.Vfct) ≈ TFI.coefficients(o2.Vfct)

    # Both reach the same fixed point
    o1, o2 = deepcopy(obj), deepcopy(obj)
    r1 = main_execute!(o1; maxiter=2000, tol=-8, damp=0.5, damping=:damped_iterate, logging=COMPLETE)
    r2 = main_execute!(o2; maxiter=2000, tol=-8, damp=0.5, damping=:previous_solve, logging=COMPLETE)
    @test r1.damping === :damped_iterate && r2.damping === :previous_solve
    @test r1.distV < -8 && r2.distV < -8
    @test all(isapprox.(TFI.coefficients(o1.Vfct), TFI.coefficients(o2.Vfct); atol=1e-6))

    # A fresh call per iteration under :previous_solve equals :damped_iterate
    o1, o2 = deepcopy(obj), deepcopy(obj)
    main_execute!(o1; maxiter=5, damp=0.5, damping=:damped_iterate, logging=COMPLETE)
    for _ in 1:5
        main_execute!(o2; maxiter=1, damp=0.5, damping=:previous_solve, logging=COMPLETE)
    end
    @test TFI.coefficients(o1.Vfct) ≈ TFI.coefficients(o2.Vfct)

    @test_throws ArgumentError main_execute!(deepcopy(obj); maxiter=1, damping=:unknown, logging=COMPLETE)

end
