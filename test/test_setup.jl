# Tests for setup.jl utility functions and main_setup branching.
# NGM fixtures (ngm, obj, stv) are brought into scope by test_ngm.jl which is included first.

@testset "setup — explicit_parameters and update_model" begin

    # explicit_parameters returns fieldnames
    params = TFI.explicit_parameters(ngm)
    @test params == fieldnames(NeoclassicalGrowthModel)
    @test :α ∈ params
    @test :β ∈ params

    # update_model with NamedTuple
    ngm2 = TFI.update_model(ngm, (; α = 0.40))
    @test ngm2 isa NeoclassicalGrowthModel
    @test ngm2.α ≈ 0.40
    @test ngm2.β == ngm.β   # unchanged fields preserved

    # update_model with keyword arguments
    ngm3 = TFI.update_model(ngm; σ = 0.03)
    @test ngm3.σ ≈ 0.03
    @test ngm3.α == ngm.α

end

@testset "setup — bounds fallback" begin

    # DummyModel does not override bounds, so the default (all Unbounded) applies.
    # Build a minimal SteadyState for DummyModel.
    dummy_stv = SteadyState(
        model  = DummyModel(),
        guess  = (p = 1.0,),
        endog  = (x = 1.0,),
        exog   = (;),
        values = (;),
    )

    pb, tb = TFI.bounds(dummy_stv)
    @test pb.p == TFI.Unbounded()
    @test tb.x == TFI.Unbounded()

end

@testset "setup — vals_from_policy fallback" begin

    # Default: returns policy unchanged
    state  = (x = 1.0, exst = 1)
    policy = (p = 0.5,)
    result = TFI.vals_from_policy(state, policy, DummyModel())
    @test result == policy

end

@testset "setup — main_setup with check_steadystate=false" begin

    # Should construct obj without running the steady-state consistency check
    obj2, stv2 = main_setup(ngm; check_steadystate=false)
    @test obj2 isa ModelSolution
    @test stv2 isa SteadyState

end

@testset "setup — main_setup with initobj (warm-start)" begin

    # Re-initializing from an existing solution re-interpolates coefficients onto the grid.
    # The resulting solution should have the same model and approximate interpolant values.
    obj_warm, stv_warm = main_setup(ngm, obj)
    @test obj_warm isa ModelSolution
    @test obj_warm.model == ngm
    # Interpolant values at the grid should be close to the original
    @test isapprox(obj_warm, obj; rtol=1e-6)

end

@testset "Plotting recipes — smoke tests (NGM)" begin

    # smoke: just verify apply_recipe completes without error and returns non-empty RecipeData.

    # SimulationPanel solo recipe (NGM simseries has K as endog state)
    specs_panel = RecipesBase.apply_recipe(Dict{Symbol,Any}(), simseries)
    @test !isempty(specs_panel)

    # SimulationErrorPanel + SimulationPanel recipe.
    # Use a synthetic 1-endog-state panel so errs and sim are guaranteed the same size,
    # and to verify the plotting.jl:29 fix (min(2, n_endog) clamp) works correctly.
    nruns, nperiods = 3, 20
    mock_data = StructArray([
        (K = 3.0 + 0.1*i, exst = mod(i, 3) + 1)
        for i in 1:nruns, j in 1:nperiods
    ])
    mock_panel = SimulationPanel(mock_data, (:K,))
    mock_errdata = StructArray([(euler = 1e-6 * Float64(i),) for i in 1:nruns, j in 1:nperiods])
    mock_errs = TFI.SimulationErrorPanel(mock_errdata)

    specs_err = RecipesBase.apply_recipe(Dict{Symbol,Any}(), mock_errs, mock_panel)
    @test !isempty(specs_err)

end

@testset "setup — FixedPointSteadyStateMode" begin

    # Construct mode from steady-state guess values and exog NamedTuple
    fp_mode = FixedPointSteadyStateMode(stv.guess, stv.endog, stv.exog)
    @test fp_mode isa FixedPointSteadyStateMode

    # Run main_setup with FixedPointSteadyStateMode; steady state should be found
    obj_fp, stv_fp = main_setup(ngm; mode=fp_mode, check_steadystate=false)
    @test obj_fp isa ModelSolution
    @test stv_fp isa SteadyState

    # Fixed-point steady state should be close to the analytical one
    @test stv_fp.endog.K ≈ stv.endog.K  rtol=1e-3

end
