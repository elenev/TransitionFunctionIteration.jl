# Forecasting functions (Ffct, forecasts_from_policy).
#
# Test model: the neoclassical growth model of examples/ngm.jl, written so that next period's
# expectations need current consumption: compute_next_period_values returns the Euler kernel
#     β (C′/C)^(-γ) (MPK′ + 1 - δ),
# with C′ from Vfct at the next-period state and current C supplied as a forecasting value
# (Cnow). calc_eqm rescales the kernel by (Cnow/C)^(-γ) with its own policy C, so the Euler
# equation is the standard NGM's whatever Cnow is. Under TFIMethod, Cnow is the stored forecast
# (frozen during the solve) and the policy still enters through C; under STPFIMethod Cnow = C.
# Both solution methods must reproduce `obj` from test_ngm.jl. Uses ngm, obj from test_ngm.jl.

import TransitionFunctionIteration: Exogenv, calc_eqm, define_endog_grids, steadystate,
    compute_next_period_values, vals_from_policy, bounds, forecasts_from_policy

struct ForecastNGM <: AbstractModel
    base::NeoclassicalGrowthModel
end

forecasts_from_policy(state, policy, ::ForecastNGM) = (Cnow = policy.C,)
vals_from_policy(_, policy, ::ForecastNGM) = (C = policy.C,)
Exogenv(m::ForecastNGM) = Exogenv(m.base)
define_endog_grids(m::ForecastNGM, stv::AbstractSteadyState) = define_endog_grids(m.base, stv)
function steadystate(m::ForecastNGM, mode::AnalyticalSteadyStateMode)
    s, ok = steadystate(m.base, mode)
    SteadyState(; model = m, exog = s.exog, endog = s.endog, values = s.values, guess = s.guess), ok
end
bounds(::SteadyState{ForecastNGM}) = ((C = TransitionFunctionIteration.PositiveBounded(),),
                                      (K = TransitionFunctionIteration.PositiveBounded(),))

function compute_next_period_values(next, m::ForecastNGM, ::EvaluationMode)
    (; A, K, C, Cnow) = next
    (; α, β, δ, γ) = m.base
    Kpayoff = A * TransitionFunctionIteration.cobb_douglas_prime(K, 1.0, α, 1) + 1 - δ
    return NextPeriod(; expectations = (; euler_kernel = β * (C / Cnow)^(-γ) * Kpayoff,),
                        continuations = (; Kpayoff, Cnow))
end

function calc_eqm(state, policy, m::ForecastNGM, np::NextPeriod, mode::EvaluationMode = SolveMode())
    (; A, K) = state
    (; C) = policy
    (; α, δ) = m.base
    Y = A * TransitionFunctionIteration.cobb_douglas(K, 1.0, α)
    Cnow = np.continuations.Cnow[1]          # the same in every next-period state
    fx = SVector(1 - np.expectations.euler_kernel * (Cnow / C)^(-m.base.γ),)
    mode isa SolveMode && return fx
    trans = (K = (1 - δ) * K + Y - C,)
    mode isa IterateMode && return (; fx, trans, vals = (; C))
    return (; fx, trans, simvars = (; Y, C))
end

# A variant whose Euler residual uses the expectation as it is, so today's C enters it only through
# the forecast Cnow. Within a TFI solve the forecast is fixed, so the residual does not depend on the
# policy and the solver fails at every point.
struct OnlyForecastNGM <: AbstractModel
    base::ForecastNGM
end
for f in (:forecasts_from_policy, :vals_from_policy)
    @eval $f(state, policy, m::OnlyForecastNGM) = $f(state, policy, m.base)
end
Exogenv(m::OnlyForecastNGM) = Exogenv(m.base)
define_endog_grids(m::OnlyForecastNGM, stv::AbstractSteadyState) = define_endog_grids(m.base.base, stv)
function steadystate(m::OnlyForecastNGM, mode::AnalyticalSteadyStateMode)
    s, ok = steadystate(m.base.base, mode)
    SteadyState(; model = m, exog = s.exog, endog = s.endog, values = s.values, guess = s.guess), ok
end
bounds(::SteadyState{OnlyForecastNGM}) = ((C = TransitionFunctionIteration.PositiveBounded(),),
                                           (K = TransitionFunctionIteration.PositiveBounded(),))
compute_next_period_values(next, m::OnlyForecastNGM, mode::EvaluationMode) = compute_next_period_values(next, m.base, mode)
function calc_eqm(state, policy, m::OnlyForecastNGM, np::NextPeriod, mode::EvaluationMode = SolveMode())
    out = calc_eqm(state, policy, m.base, np, mode)
    fx = SVector(1 - np.expectations.euler_kernel,)
    mode isa SolveMode && return fx
    return merge(out, (; fx))
end

@testset "Forecasting functions" begin
    fm = ForecastNGM(NeoclassicalGrowthModel())

    @testset "setup builds Ffct; models without forecasts get nothing" begin
        objf, stvf = main_setup(fm)
        @test !isnothing(objf.Ffct)
        @test varnames(objf.Ffct) == (:Cnow,)
        @test all(c -> c[1] ≈ stvf.values.C, coefficients(objf.Ffct))
        @test isnothing(obj.Ffct)                     # NGM from test_ngm.jl
        @test forecasts_from_policy(nothing, (C = 1.0,), ngm) == (;)
    end

    @testset "fixed-point steady state sees the forecasts" begin
        stv_an, _ = steadystate(fm, AnalyticalSteadyStateMode())
        stv_fp, ok = steadystate(fm, FixedPointSteadyStateMode((C = 0.9 * stv_an.values.C,),
                                                               (K = 1.1 * stv_an.endog.K,), (A = 1.0,)))
        @test ok
        @test isapprox(stv_fp.endog.K, stv_an.endog.K; rtol = 1e-8)
        @test isapprox(stv_fp.guess.C, stv_an.values.C; rtol = 1e-8)
    end

    @testset "TFIMethod reproduces the NGM solution" begin
        objf, _ = main_setup(fm)
        res = main_execute!(objf; maxiter = 1000, tol = -7, logging = COMPLETE, solmethod = TFIMethod())
        @test res.status == :converged
        @test res.distF < -7
        @test isapprox(coefficients(objf.Pfct), coefficients(obj.Pfct); rtol = 1e-5)
        # the stored forecasts equal the solved policy at every point (up to the iteration tolerance)
        @test all(i -> isapprox(coefficients(objf.Ffct)[i][1],
                       TransitionFunctionIteration.to_constrained(coefficients(objf.Pfct)[i][1], objf.Pfct_bounds.C);
                       atol = 1e-6),
                  eachindex(coefficients(objf.Ffct)))
    end

    @testset "main_execute! exit states" begin
        # every point fails: exits after the first iteration with status :all_failed, logged at COMPLETE
        objo, _ = main_setup(OnlyForecastNGM(fm))
        buf = IOBuffer()
        res = main_execute!(objo; maxiter = 50, solmethod = TFIMethod(), logging = TerminalLogger(COMPLETE; io = buf))
        @test res.status == :all_failed
        @test res.iter == 1
        @test contains(String(take!(buf)), "Solver failed at every grid point.")
        # A point that fails keeps its own policy, also when its guesses come from its neighbors in
        # the second pass. Start from a policy that differs at every point, so a neighbor's guess differs.
        objn, _ = main_setup(OnlyForecastNGM(fm))
        pc = coefficients(objn.Pfct)
        foreach(k -> pc[k] = pc[k] .+ 0.01k, eachindex(pc))
        pfct_before = copy(pc)
        resn = main_execute!(objn; maxiter = 1, npasses = 2, solmethod = TFIMethod(), logging = COMPLETE)
        @test resn.status == :all_failed
        @test coefficients(objn.Pfct) == pfct_before
        # maximum iterations
        objm, _ = main_setup(fm)
        @test main_execute!(objm; maxiter = 2, logging = COMPLETE).status == :maxiter
    end

    @testset "STPFIMethod reproduces the NGM solution" begin
        objf, stvf = main_setup(fm)
        main_execute!(objf; maxiter = 200, logging = COMPLETE, solmethod = STPFIMethod())
        @test isapprox(coefficients(objf.Pfct), coefficients(obj.Pfct); rtol = 1e-3)
        simf, errf = simulate(objf, stvf)
        @test simf isa TransitionFunctionIteration.SimulationPanel
    end

    @testset "solve_at_point with forecasts: allocation-free, type-stable" begin
        objf, _ = main_setup(fm)
        i = 1
        state = objf.ssg[i]
        pol = coefficients(objf.Pfct)[i]
        val = coefficients(objf.Vfct)[i]
        trans = coefficients(objf.Tfct)[i]
        ff = coefficients(objf.Ffct)[i]
        old_tfi = TransitionFunctionIteration.FunctionValuesAtPoint(state, pol, val, trans, ff)
        old_stpfi = TransitionFunctionIteration.FunctionValuesAtPoint(state, (pol, trans), val, trans, ff)
        @test allocs_after_warmup(TransitionFunctionIteration.solve_at_point, old_tfi, objf, TFIMethod()) == 0
        @test allocs_after_warmup(TransitionFunctionIteration.solve_at_point, old_stpfi, objf, STPFIMethod()) == 0
        res, _ = @inferred(TransitionFunctionIteration.solve_at_point(old_tfi, objf, TFIMethod()))
        @test res.ffct isa SVector
    end
end
