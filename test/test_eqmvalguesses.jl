# Fixed-point steady state with a Vfct variable that depends on expectations.
#
# Test model: the NGM of examples/ngm.jl with the value function V = u(C) + β E[V′] added to
# Vfct. vals_from_policy can return C but not V, so V is an equilibrium-value guess. In the
# steady state V = u(C)/(1 - β). Uses NeoclassicalGrowthModel from test_ngm.jl.

import TransitionFunctionIteration: Exogenv, calc_eqm, define_endog_grids, steadystate,
    compute_next_period_values, vals_from_policy, bounds

struct ValueNGM <: AbstractModel
    base::NeoclassicalGrowthModel
end

vals_from_policy(_, policy, ::ValueNGM) = (C = policy.C,)
Exogenv(m::ValueNGM) = Exogenv(m.base)
define_endog_grids(m::ValueNGM, stv::AbstractSteadyState) = define_endog_grids(m.base, stv)
bounds(::ValueNGM, p::NamedTuple, t::NamedTuple, v::NamedTuple) =
    ((C = TransitionFunctionIteration.PositiveBounded(),), (K = TransitionFunctionIteration.PositiveBounded(),),
     (V = TransitionFunctionIteration.LowerBounded(-100.0),))

function steadystate(m::ValueNGM, mode::AnalyticalSteadyStateMode)
    s, ok = steadystate(m.base, mode)
    V = TransitionFunctionIteration.crra.u(s.values.C, m.base.γ) / (1 - m.base.β)
    SteadyState(; model = m, exog = s.exog, endog = s.endog, values = (C = s.values.C, V), guess = s.guess), ok
end

function compute_next_period_values(next, m::ValueNGM, mode::EvaluationMode)
    np = compute_next_period_values(next, m.base, mode)
    NextPeriod(; expectations = merge(np.expectations, (EV = next.V,)), continuations = np.continuations)
end

function calc_eqm(state, policy, m::ValueNGM, np::NextPeriod, mode::EvaluationMode = SolveMode())
    out = calc_eqm(state, policy, m.base, np, mode)
    mode isa SolveMode && return out
    V = TransitionFunctionIteration.crra.u(policy.C, m.base.γ) + m.base.β * np.expectations.EV
    mode isa IterateMode && return (; out.fx, out.trans, vals = (C = policy.C, V))
    return out
end

@testset "Fixed-point steady state with equilibrium-value guesses" begin
    vm = ValueNGM(NeoclassicalGrowthModel())
    stv_an, _ = steadystate(vm, AnalyticalSteadyStateMode())

    # V is an extra unknown pinned by V = u(C) + β V
    mode = FixedPointSteadyStateMode((C = 0.9 * stv_an.values.C,), (K = 1.1 * stv_an.endog.K,), (A = 1.0,);
                                     eqmvalguesses = (V = 0.8 * stv_an.values.V,))
    stv_fp, ok, _ = steadystate(vm, mode)
    @test ok
    @test isapprox(stv_fp.endog.K, stv_an.endog.K; rtol = 1e-8)
    @test isapprox(stv_fp.values.V, stv_an.values.V; rtol = 1e-8)

    # main_setup's consistency check finds V among the equilibrium values and passes silently
    @test_logs main_setup(vm)

    # a wrong analytical V is caught by the check
    TransitionFunctionIteration.f_steadystate(stv_an.guess, stv_an.endog, (V = stv_an.values.V,),
        (model = vm, exog = stv_an.exog))[1] |> r -> @test maximum(abs, r) < 1e-10
    r_bad = TransitionFunctionIteration.f_steadystate(stv_an.guess, stv_an.endog, (V = 1.1 * stv_an.values.V,),
        (model = vm, exog = stv_an.exog))[1]
    @test abs(r_bad[end]) > 1e-3

    # guesses that overlap vals_from_policy are rejected
    @test_throws ArgumentError TransitionFunctionIteration.f_steadystate(stv_an.guess, stv_an.endog,
        (C = stv_an.values.C,), (model = vm, exog = stv_an.exog))

    # without guesses, the NGM steady state is unchanged
    stv_ngm, ok_ngm, _ = steadystate(ngm, FixedPointSteadyStateMode((C = 1.0,), (K = 3.0,), (A = 1.0,)))
    @test ok_ngm
    @test isapprox(stv_ngm.endog.K, steadystate(ngm, AnalyticalSteadyStateMode())[1].endog.K; rtol = 1e-8)
end
