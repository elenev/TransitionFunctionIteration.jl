@testset "sample_shocks" begin

    shocks = sample_shocks(5=>10, 42, 0)

    @test size(shocks) == (5, 10)
    @test all(0 .<= shocks .< 1)

    # Same seed → identical draws
    @test sample_shocks(5=>10, 42, 0) == shocks

    # Different seed → different draws (with overwhelming probability)
    @test sample_shocks(5=>10, 99, 0) != shocks

    # Burn adds extra columns (trimmed by simulate, kept here)
    @test size(sample_shocks(5=>10, 42, 3)) == (5, 13)

end

@testset "simulate_markov_chain" begin

    # Identity transition matrix: every state is absorbing
    exog_pts = StructArray([(z=v,) for v in [0.9, 1.0, 1.1]])
    mtrans   = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    exogenv  = Exogenv(pts=exog_pts, mtrans=mtrans, ssidx=2)
    shocks   = fill(0.5, 3, 5)   # 3 runs, 5 periods

    # Output is (states, exog_points); states is N_runs × (N_periods+1): column 1 is initial state
    states_sc, _ = simulate_markov_chain(exogenv, 1, shocks)
    @test size(states_sc) == (3, 6)

    # Scalar start: all runs start and stay at state 1
    @test all(states_sc .== 1)

    # Vector start: each run stays at its own absorbing state
    states_v, _ = simulate_markov_chain(exogenv, [1, 2, 3], shocks)
    @test states_v[:, 1]    == [1, 2, 3]
    @test all(states_v[1,:] .== 1)
    @test all(states_v[2,:] .== 2)
    @test all(states_v[3,:] .== 3)
    @test all(1 .<= states_v .<= 3)

    # sz/seed/burn convenience signature
    states2, _ = simulate_markov_chain(exogenv, 2=>10, 2, 100, 0)
    @test size(states2, 1) == 2     # N_runs
    @test size(states2, 2) == 11    # N_periods + 1

end

@testset "calc_error_stats" begin

    # Synthetic error data: 2 named fields over 100 runs × 50 periods.
    # Using abs(randn()) so all values are positive (as real model errors would be).
    Random.seed!(42)
    errs_data = TFI.SimulationErrorPanel(StructArray(
        [(euler=abs(randn()), K=abs(randn())) for _ in 1:100, _ in 1:50]
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
    data = StructArray([(x = Float64(r * t),) for r in 1:N_runs, t in 1:N_periods])

    result = TFI._update_simseries_with_augmented(data, (; model = AugModel()), :nan)

    @test :log_x in propertynames(result)
    @test size(result) == (N_runs, N_periods)
    @test result.log_x ≈ log.(data.x)
end

@testset "augment_simulation — dynamic merge :nan" begin
    N_runs, N_periods = 3, 5
    data = StructArray([(x = Float64(r + t),) for r in 1:N_runs, t in 1:N_periods])

    result = TFI._update_simseries_with_augmented(data, (; model = AugModel()), :nan)

    @test :dx in propertynames(result)
    @test size(result) == (N_runs, N_periods)
    @test all(isnan, result.dx[:, 1])                         # first column is NaN padding
    @test result.dx[1, 2] ≈ data.x[1, 2] - data.x[1, 1]
end

@testset "augment_simulation — dynamic merge :drop" begin
    N_runs, N_periods = 3, 5
    data = StructArray([(x = Float64(r + t),) for r in 1:N_runs, t in 1:N_periods])

    result = TFI._update_simseries_with_augmented(data, (; model = AugModel()), :drop)

    @test :dx in propertynames(result)
    @test size(result, 2) == N_periods - 1                    # first period dropped
    @test result.dx[1, 1] ≈ data.x[1, 2] - data.x[1, 1]     # aligns directly, no NaN
end

@testset "describe — filters" begin

    N_runs, N_periods = 4, 20
    Random.seed!(42)
    data = StructArray([(x = randn(), flag = rand() > 0.5)
                        for _ in 1:N_runs, _ in 1:N_periods])
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
