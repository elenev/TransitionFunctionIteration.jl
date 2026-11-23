import TransitionFunctionIteration.Transformations: Level, AnnLevel, AnnInvLevel, AnnQuadLevel,
    Pct, AnnPct, AnnInvPct, AnnQuadPct,
    DiffLevel, DiffAnnLevel, DiffAnnInvLevel, DiffAnnQuadLevel,
    DiffPct, DiffAnnPct, DiffAnnInvPct, DiffAnnQuadPct,
    RatioPct, RatioAnnPct,
    DiffInitLevel, DiffInitPct

@testset "Transformations — operators" begin
    x    = [1.0, 2.0, 3.0]
    base = [0.5, 1.0, 1.5]

    @test Level()(x) == x
    @test AnnLevel()(x) ≈ x * 4
    @test AnnInvLevel()(x) ≈ x / 4
    @test AnnQuadLevel()(x) ≈ x * sqrt(4)
    @test Pct()(x) ≈ x * 100
    @test AnnPct()(x) ≈ x * 4 * 100
    @test AnnInvPct()(x) ≈ x / 4
    @test AnnQuadPct()(x) ≈ x * sqrt(4) * 100

    @test DiffLevel()(x; basis=base) ≈ x - base
    @test DiffAnnLevel()(x; basis=base) ≈ (x - base) * 4
    @test DiffAnnInvLevel()(x; basis=base) ≈ (x - base) / 4
    @test DiffAnnQuadLevel()(x; basis=base) ≈ (x - base) * sqrt(4)
    @test DiffPct()(x; basis=base) ≈ (x - base) * 100
    @test DiffAnnPct()(x; basis=base) ≈ (x - base) * 4 * 100
    @test DiffAnnInvPct()(x; basis=base) ≈ (x - base) / 4 * 100
    @test DiffAnnQuadPct()(x; basis=base) ≈ (x - base) * sqrt(4) * 100

    @test RatioPct()(x; basis=ones(3)) ≈ (x - ones(3)) * 100
    @test RatioAnnPct()(x; basis=ones(3)) ≈ (x - ones(3)) * 100

    @test DiffInitLevel()(x) ≈ x .- x[1]
    @test DiffInitPct()(x) ≈ (x .- x[1]) ./ x[1] * 100
end

@testset "Transformations — is_ratio" begin
    @test Transformations.is_ratio(RatioPct)    == true
    @test Transformations.is_ratio(RatioAnnPct) == true
    @test Transformations.is_ratio(Level)       == false
    @test Transformations.is_ratio(DiffPct)     == false
    @test Transformations.is_ratio(AnnPct)      == false
end

@testset "Transformations — labels" begin
    @test Transformations.label(Level)       == ""
end

@testset "Transformations — labels (all remaining types)" begin
    # Every type should return a non-nothing String (even empty is fine for Level)
    all_types = [AnnLevel, AnnInvLevel, AnnQuadLevel,
                 AnnPct, AnnInvPct, AnnQuadPct,
                 DiffAnnLevel, DiffAnnInvLevel, DiffAnnQuadLevel,
                 DiffAnnPct, DiffAnnInvPct, DiffAnnQuadPct,
                 RatioPct, RatioAnnPct]

    for T in all_types
        lbl = Transformations.label(T)
        @test lbl isa String
    end
end

@testset "Transformations — is_ratio (all types)" begin
    ratio_types    = [RatioPct, RatioAnnPct]
    nonratio_types = [Level, AnnLevel, AnnInvLevel, AnnQuadLevel,
                      Pct, AnnPct, AnnInvPct, AnnQuadPct,
                      DiffLevel, DiffAnnLevel, DiffAnnInvLevel, DiffAnnQuadLevel,
                      DiffPct, DiffAnnPct, DiffAnnInvPct, DiffAnnQuadPct,
                      DiffInitLevel, DiffInitPct]

    for T in ratio_types
        @test Transformations.is_ratio(T) == true
    end
    for T in nonratio_types
        @test Transformations.is_ratio(T) == false
    end
end

@testset "Transformations — callable keyword arguments" begin
    x     = [1.0, 2.0, 3.0]
    basis = [0.5, 1.0, 1.5]
    trend = [0.1, 0.2, 0.3]

    # DiffLevel with explicit basis and trend: result = x - basis
    @test DiffLevel()(x; basis=basis) ≈ x .- basis

    # DiffLevel with trend kwarg: trend is accepted but not used by this type
    @test DiffLevel()(x; basis=basis, trend=trend) ≈ x .- basis

    # DiffInitLevel and DiffInitPct with explicit init parameter
    @test DiffInitLevel()(x; init=2.0) ≈ x .- 2.0
    @test DiffInitPct()(x; init=2.0)   ≈ (x .- 2.0) ./ 2.0 .* 100

    # RatioPct with explicit basis (non-unit)
    @test RatioPct()(x; basis=basis) ≈ (x ./ basis .- 1) .* 100
end
