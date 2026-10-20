@testset "sample_shocks" begin

    shocks = sample_shocks(10=>5, 42, 0)

    @test size(shocks) == (10, 5)
    @test all(0 .<= shocks .< 1)

    # Same seed → identical draws
    @test sample_shocks(10=>5, 42, 0) == shocks

    # Different seed → different draws (with overwhelming probability)
    @test sample_shocks(10=>5, 99, 0) != shocks

    # Burn adds extra rows (trimmed by simulate, kept here)
    @test size(sample_shocks(10=>5, 42, 3)) == (13, 5)

end

@testset "simulate_markov_chain" begin

    # Identity transition matrix: every state is absorbing
    exog_pts = StructArray([(z=v,) for v in [0.9, 1.0, 1.1]])
    mtrans   = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    exogenv  = Exogenv(pts=exog_pts, mtrans=mtrans, ssidx=2)
    shocks   = fill(0.5, 5, 3)   # 5 periods, 3 runs

    # Output is (states, exog_points); states is (N_periods+1) × N_runs: row 1 is initial state
    states_sc, _ = simulate_markov_chain(exogenv, 1, shocks)
    @test size(states_sc) == (6, 3)

    # Scalar start: all runs start and stay at state 1
    @test all(states_sc .== 1)

    # Vector start: each run stays at its own absorbing state
    states_v, _ = simulate_markov_chain(exogenv, [1, 2, 3], shocks)
    @test states_v[1, :]    == [1, 2, 3]
    @test all(states_v[:,1] .== 1)
    @test all(states_v[:,2] .== 2)
    @test all(states_v[:,3] .== 3)
    @test all(1 .<= states_v .<= 3)

    # sz/seed/burn convenience signature
    states2, _ = simulate_markov_chain(exogenv, 10=>2, 2, 100, 0)
    @test size(states2, 2) == 2     # N_runs
    @test size(states2, 1) == 11    # N_periods + 1

end

@testset "calc_error_stats" begin

    # Synthetic error data: 2 named fields over 100 runs × 50 periods.
    # Using abs(randn()) so all values are positive (as real model errors would be).
    Random.seed!(42)
    errs_data = TFI.SimulationErrorPanel(StructArray(
        [(euler=abs(randn()), K=abs(randn())) for _ in 1:50, _ in 1:100]
    ))

    stats = calc_error_stats(errs_data)

    # Output is a MomentArray with one row per error field, one column per statistic
    @test stats isa TFI.MomentArray
    @test size(stats, 1) == 2                    # :euler and :K
    @test size(stats, 2) == 7                    # mean + 6 default percentiles
    @test row_names(stats) == (:euler, :K)
    @test :mean in col_names(stats)
    @test :p99  in col_names(stats)

    # Basic statistical sanity checks
    @test stats(:euler, :mean) > 0
    @test stats(:K,     :mean) > 0

    # Higher percentiles ≥ lower percentiles
    @test stats(:euler, :p99) >= stats(:euler, :p50)
    @test stats(:K,     :p99) >= stats(:K,     :p50)

    # Custom percentiles list is respected
    stats2 = calc_error_stats(errs_data; percentiles=[50, 95])
    @test size(stats2, 2) == 3    # mean + 2 percentiles

end

# AugModel returns one static (T-length) and one dynamic (T-1-length) extra variable
struct AugModel <: AbstractModel end
function TFI.augment_simulation(run::StructVector, ::AugModel)
    return (; log_x = log.(run.x),   # static: T-length
              dx    = diff(run.x))   # dynamic: (T-1)-length
end

@testset "augment_simulation — default no-op" begin
    run = StructArray([(x = Float64(t), flag = t % 2 == 0) for t in 1:8])
    extra = augment_simulation(run, DummyModel())
    @test isnothing(extra)
end

@testset "augment_simulation — static merge :nan" begin
    N_runs, N_periods = 2, 6
    data = StructArray([(x = Float64(r * t),) for t in 1:N_periods, r in 1:N_runs])

    result = TFI._update_simseries_with_augmented(data, (; model = AugModel()), :nan)

    @test :log_x in propertynames(result)
    @test size(result) == (N_periods, N_runs)
    @test result.log_x ≈ log.(data.x)
end

@testset "augment_simulation — dynamic merge :nan" begin
    N_runs, N_periods = 3, 5
    data = StructArray([(x = Float64(r + t),) for t in 1:N_periods, r in 1:N_runs])

    result = TFI._update_simseries_with_augmented(data, (; model = AugModel()), :nan)

    @test :dx in propertynames(result)
    @test size(result) == (N_periods, N_runs)
    @test all(isnan, result.dx[1, :])                         # first row is NaN padding
    @test result.dx[2, 1] ≈ data.x[2, 1] - data.x[1, 1]
end

@testset "augment_simulation — dynamic merge :drop" begin
    N_runs, N_periods = 3, 5
    data = StructArray([(x = Float64(r + t),) for t in 1:N_periods, r in 1:N_runs])

    result = TFI._update_simseries_with_augmented(data, (; model = AugModel()), :drop)

    @test :dx in propertynames(result)
    @test size(result, 1) == N_periods - 1                    # first period dropped
    @test result.dx[1, 1] ≈ data.x[2, 1] - data.x[1, 1]     # aligns directly, no NaN
end

@testset "describe — filters" begin

    N_runs, N_periods = 4, 20
    Random.seed!(42)
    data = StructArray([(x = randn(), flag = rand() > 0.5)
                        for _ in 1:N_periods, _ in 1:N_runs])
    panel = SimulationPanel(data, (:x,))

    out = describe(panel; filters=(uncond=Returns(true), flagged=row->row.flag))

    # Return is a NamedTuple (; moments, stderrs)
    @test haskey(out, :moments)
    @test haskey(out, :stderrs)

    # Both filters present in moments and stderrs
    @test haskey(out.moments, :uncond)
    @test haskey(out.moments, :flagged)
    @test haskey(out.stderrs, :uncond)
    @test haskey(out.stderrs, :flagged)

    # MomentArrays have correct row names
    @test :x in row_names(out.moments.uncond)
    @test :x in row_names(out.moments.flagged)

    # Unconditional and conditional moments differ
    μ_uncond  = out.moments.uncond(:x, :mean)
    μ_flagged = out.moments.flagged(:x, :mean)
    # Both are finite (not NaN)
    @test isfinite(μ_uncond)
    @test isfinite(μ_flagged)

end

# ── augment_errors ──────────────────────────────────────────────────────────────
# DummyModel (defined in runtests.jl, shared `obj`) has no augment_errors override -- used here
# to confirm the DEFAULT is a true no-op. A SEPARATE model (AugErrModel) below is the only place
# an override is ever defined, precisely so "default" and "with-override" behavior can't be
# confused by method-definition ordering -- once an override exists for a type, every call for
# that type uses it, so keeping them on distinct types sidesteps that entirely.

@testset "augment_errors — default is a true no-op" begin
    state = obj.ssg[1]
    out, trans = TFI.calc_eqm_from_state(state, obj)

    @test augment_errors(obj, state, out, trans) == ((;), (;))

    _, _, fx_nt, transerr = TFI.calc_result(state, obj)
    @test keys(fx_nt) == (:fx1,)              # nothing extra merged in
    @test propertynames(transerr) == (:x,)    # nothing extra merged in
end

# A dedicated model with a state-dependent `v` (unlike DummyModel's constant one) and a real
# augment_errors override -- exercises both the scalar_nt→fx and vector_nt→transerr paths, plus
# the "override replaces an existing fx field" behavior.
struct AugErrModel <: AbstractModel end

function TransitionFunctionIteration.compute_next_period_values(next, ::AugErrModel, ::EvaluationMode)
    expectations = (e = next.v,)
    continuations = (; c = next.x)
    NextPeriod(; expectations, continuations)
end

function TransitionFunctionIteration.calc_eqm(state, policy, ::AugErrModel, np::NextPeriod, mode::EvaluationMode)
    fx = (fx1 = np.expectations.e - policy.p + state.x,)
    mode isa SolveMode && return fx
    trans = (x = 1.0,)
    v = 2.0 * state.x   # state-dependent, unlike the constant Vfct coefficient below
    if mode isa IterateMode
        return (; fx, trans, vals = (; v))
    elseif mode isa SimulateMode
        return (; fx, trans, simvars = (; v))
    end
end

# Vfct's stored coefficient is a CONSTANT 3.0 everywhere, so it disagrees with the
# state-dependent exact v=2*state.x above -- makes the diagnostic's value non-trivial to fake.
function TransitionFunctionIteration.augment_errors(obj::ModelSolution{AugErrModel}, state, out, trans)
    v_interp = obj.Vfct(state).v
    scalar_nt = (; v_check = out.simvars.v - v_interp,
                   fx1 = 1000.0 + out.fx.fx1)   # deliberately collides with the base `fx1` name
    vector_nt = (; x_check = 2.0 .* trans.x)
    return scalar_nt, vector_nt
end

pfct_ec = Interpolant(ssg_ex, (p = 1.0,))
vfct_ec = Interpolant(ssg_ex, (v = 3.0,))
tfct_ec = Interpolant(ssg_ex, (x = SVector(0.9, 1.0, 1.1),))
obj_ec  = ModelSolution(model = AugErrModel(), exogenv = e3, Pfct = pfct_ec, Vfct = vfct_ec, Tfct = tfct_ec, ssg = ssg_ex)

@testset "augment_errors — scalar merges into fx, including a deliberate name collision" begin
    state = obj_ec.ssg[1]
    out, trans = TFI.calc_eqm_from_state(state, obj_ec)
    scalar_nt, vector_nt = augment_errors(obj_ec, state, out, trans)

    @test scalar_nt.v_check ≈ out.simvars.v - obj_ec.Vfct(state).v
    @test scalar_nt.v_check ≈ 2.0 * state.x - 3.0
    @test :x_check in propertynames(vector_nt)

    _, _, fx_nt, transerr = TFI.calc_result(state, obj_ec)
    @test fx_nt.v_check ≈ scalar_nt.v_check
    @test fx_nt.fx1 ≈ 1000.0 + out.fx.fx1        # override WON, not the original fx1
    @test :x_check in propertynames(transerr)
end

@testset "augment_errors — vector merges into transerr and resolves via merge_error" begin
    state = obj_ec.ssg[1]
    _, _, fx_nt, transerr = TFI.calc_result(state, obj_ec)

    @test length(transerr.x_check) == 3
    for next_exst in 1:3
        merged = TFI.merge_error(fx_nt, transerr, next_exst)
        @test merged.x_check ≈ transerr.x_check[next_exst]
    end
end

@testset "augment_errors — flows through simulate() end to end" begin
    simseries_ec, errs_ec = simulate(obj_ec; endog_start = (x = 1.0,), sz = 8 => 3, burn = 0,
                                      report_bounds = false)

    @test :v_check in propertynames(errs_ec)
    @test :x_check in propertynames(errs_ec)
    @test all(isfinite, errs_ec.v_check)
    @test all(isfinite, errs_ec.x_check)
end

# ── describe: dropnan and empty filtered samples ────────────────────────────────

# 3 runs x 6 periods, x = t + 10r, with NaN at (t=2, run 1) and (t=4, run 3); flag = (t > 2).
# (2,1) has flag=false, (4,3) has flag=true, so a flag filter alone removes only the first.
function _nan_panel(; nans=true)
    data = StructArray([(x = Float64(t + 10r), flag = t > 2) for t in 1:6, r in 1:3])
    if nans
        data.x[2, 1] = NaN
        data.x[4, 3] = NaN
    end
    return SimulationPanel(data, (:x,)), data
end

# Expected moment: mean over runs of the per-run mean of x over rows satisfying `keep`
_expected_mean(data, keep) = mean(mean(c.x[keep.(c)]) for c in eachcol(data))

@testset "describe — dropnan excludes rows containing NaN" begin
    panel, data = _nan_panel()
    out = describe(panel; dropnan=true)

    @test out.moments.uncond(:x, :mean) ≈ _expected_mean(data, r -> !isnan(r.x))  rtol=1e-8
    @test out.moments.uncond(:x, :minimum) ≈ 21  rtol=1e-8       # mean of per-run minima 11, 21, 31
    @test isfinite(out.moments.uncond(:x, :std))

    # Without dropnan, the percentile statistics reject NaN input
    @test_throws ArgumentError describe(panel; dropnan=false)
end

@testset "describe — dropnan composes with user filters" begin
    panel, data = _nan_panel()
    filters = (late = row -> row.flag,)

    out = describe(panel; filters, dropnan=true)
    @test out.moments.late(:x, :mean) ≈ _expected_mean(data, r -> r.flag && !isnan(r.x))  rtol=1e-8
    # The unconditional filter is also wrapped
    @test out.moments.uncond(:x, :mean) ≈ _expected_mean(data, r -> !isnan(r.x))  rtol=1e-8

    # Without dropnan, the (unfiltered) uncond sample contains NaN and describe throws
    @test_throws ArgumentError describe(panel; filters, dropnan=false)
end

@testset "describe — empty filtered sample gives NaN moments" begin
    panel, _ = _nan_panel(nans=false)
    out = describe(panel; filters=(never = Returns(false),))   # must not throw

    @test all(isnan, values(out.moments.never(:x, (:mean, :std, :minimum, :maximum, :p50))))
    @test isnan(out.stderrs.never(:x, :mean))
    @test isfinite(out.moments.uncond(:x, :mean))

    # All-NaN panel with dropnan=true leaves every sample empty
    data_nan = StructArray([(x = NaN, flag = true) for _ in 1:5, _ in 1:2])
    out_nan = describe(SimulationPanel(data_nan, (:x,)); dropnan=true)
    @test all(isnan, values(out_nan.moments.uncond(:x, (:mean, :std, :p50))))
end

@testset "simulate — reproducible across calls (threaded over runs)" begin
    kw = (endog_start = (x = 1.0,), sz = 8 => 4, burn = 0, report_bounds = false)
    s1, e1 = simulate(obj_ec; kw..., seed = 7)
    s2, e2 = simulate(obj_ec; kw..., seed = 7)
    s3, _  = simulate(obj_ec; kw..., seed = 8)

    @test all(isequal(getproperty(s1.data, k), getproperty(s2.data, k)) for k in propertynames(s1.data))
    @test all(isequal(getproperty(e1.data, k), getproperty(e2.data, k)) for k in propertynames(e1.data))
    @test s1.data.exst != s3.data.exst          # a different seed changes the draws
end
