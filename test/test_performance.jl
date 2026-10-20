# ── allocation helper ──────────────────────────────────────────────────────────
# Call f(args...) once to trigger JIT compilation, then measure heap allocations
# on a second call. Uses Base.@allocated — no BenchmarkTools dependency.
# Returns the byte count; zero means fully stack-allocated / register-allocated.
function allocs_after_warmup(f, args...)
    f(args...)                         # compile
    return @allocated f(args...)
end

# ── ForwardDiff.Dual helper ────────────────────────────────────────────────────
# ForwardDiff.Dual{T,V,N}(v, p::Partials{N,V}) requires an explicit Partials object.
# This helper constructs a 1-partial Dual{Nothing,Float64,1} from scalar (value, deriv).
const _D1 = ForwardDiff.Dual{Nothing, Float64, 1}
_d1(v::Float64, dv::Float64=1.0) = _D1(v, ForwardDiff.Partials{1, Float64}((dv,)))

# ── Interpolant evaluation — Float64 ──────────────────────────────────────────

@testset "Performance — Interpolant evaluation (Float64)" begin

    # 1D interpolant with two outputs
    @test allocs_after_warmup(fct1, (; x=3.0))       == 0

    # 2D bilinear interpolant
    @test allocs_after_warmup(fct2, (; x=2.0, y=20.0)) == 0

    # Mixed NoInterp + Linear (exogenous + endogenous state)
    @test allocs_after_warmup(fct_ex, (; exst=1, x=2.0)) == 0

end

# ── Interpolant evaluation — ForwardDiff.Dual ─────────────────────────────────
# The solver differentiates policy functions w.r.t. unknowns using ForwardDiff,
# so Dual-input evaluation must also be allocation-free.

@testset "Performance — Interpolant evaluation (ForwardDiff.Dual)" begin

    x_d = _d1(3.0)
    @test allocs_after_warmup(fct1, (; x=x_d)) == 0

    x2_d = _d1(2.0);  y2_d = _d1(20.0, 0.0)
    @test allocs_after_warmup(fct2, (; x=x2_d, y=y2_d)) == 0

    xe_d = _d1(2.0)    # exst is always an integer, not differentiated
    @test allocs_after_warmup(fct_ex, (; exst=1, x=xe_d)) == 0

end

# ── Interpolant type stability ─────────────────────────────────────────────────
# @inferred checks that Julia's type-inference can determine the return type at
# compile time. A type-inference failure causes dynamic dispatch on every call.

@testset "Performance — Interpolant type stability" begin

    @test @inferred(fct1((; x=3.0)))           isa NamedTuple
    @test @inferred(fct2((; x=2.0, y=20.0)))   isa NamedTuple
    @test @inferred(fct_ex((; exst=1, x=2.0))) isa NamedTuple

end

# ── Economics functions — Float64 ─────────────────────────────────────────────

@testset "Performance — Economics functions (Float64)" begin

    x = SVector(2.0, 3.0);  w = SVector(0.5, 0.5)

    @test allocs_after_warmup(TFI.cobb_douglas,  x, w)               == 0
    @test allocs_after_warmup(TFI.ces_strict,    2.0, 3.0, 0.5, 1.5) == 0
    @test allocs_after_warmup(TFI.crra.u′,       2.0, 2.0)           == 0
    @test allocs_after_warmup(TFI.logutility.u,  2.0)                 == 0

end

# ── Economics functions — ForwardDiff.Dual ────────────────────────────────────
# Inner model functions must be differentiable and allocation-free under Duals
# because the equilibrium solver uses automatic differentiation.

@testset "Performance — Economics functions (ForwardDiff.Dual)" begin

    xd = SVector(_d1(2.0), _d1(3.0, 0.0))
    wd = SVector(0.5, 0.5)                     # weights stay Float64

    @test allocs_after_warmup(TFI.cobb_douglas, xd, wd)                            == 0
    @test allocs_after_warmup(TFI.ces_strict,   _d1(2.0), _d1(3.0,0.0), 0.5, 1.5) == 0
    @test allocs_after_warmup(TFI.crra.u′,      _d1(2.0), 2.0)                     == 0

end

# ── Economics type stability ───────────────────────────────────────────────────

@testset "Performance — Economics type stability" begin

    @test @inferred(TFI.cobb_douglas(SVector(2.0,3.0), SVector(0.5,0.5))) isa Float64
    @test @inferred(TFI.ces_strict(2.0, 3.0, 0.5, 1.5))                  isa Float64
    @test @inferred(TFI.crra.u′(2.0, 2.0))                                isa Float64
    @test @inferred(TFI.logutility.u(2.0))                                 isa Float64

    # Should also be type-stable under Duals (return type is Dual, not Any)
    @test @inferred(TFI.cobb_douglas(SVector(_d1(2.0), _d1(3.0,0.0)), SVector(0.5,0.5))) isa _D1
    @test @inferred(TFI.crra.u′(_d1(2.0), 2.0)) isa _D1

end

# ── Bounded transformations ────────────────────────────────────────────────────
# to_constrained is called in the inner solve loop and must be allocation-free.

@testset "Performance — Bounded transformations (Float64)" begin

    @test allocs_after_warmup(TFI.to_constrained, 1.0, TFI.Unbounded())              == 0
    @test allocs_after_warmup(TFI.to_constrained, 1.0, TFI.LowerBounded(0.0))        == 0
    @test allocs_after_warmup(TFI.to_constrained, 1.0, TFI.UpperBounded(5.0))        == 0
    @test allocs_after_warmup(TFI.to_constrained, 0.5, TFI.IntervalBounded(0.0,1.0)) == 0

    @test @inferred(TFI.to_constrained(1.0, TFI.Unbounded()))              isa Float64
    @test @inferred(TFI.to_constrained(1.0, TFI.LowerBounded(0.0)))        isa Float64
    @test @inferred(TFI.to_constrained(1.0, TFI.UpperBounded(5.0)))        isa Float64
    @test @inferred(TFI.to_constrained(0.5, TFI.IntervalBounded(0.0,1.0))) isa Float64

end

@testset "Performance — Bounded transformations (ForwardDiff.Dual)" begin

    x = _d1(1.0)

    @test allocs_after_warmup(TFI.to_constrained, x, TFI.Unbounded())              == 0
    @test allocs_after_warmup(TFI.to_constrained, x, TFI.LowerBounded(0.0))        == 0
    @test allocs_after_warmup(TFI.to_constrained, x, TFI.UpperBounded(5.0))        == 0
    @test allocs_after_warmup(TFI.to_constrained, _d1(0.0), TFI.IntervalBounded(0.0,1.0)) == 0

end


# ── Hot Loop in main_execute! ─────────────────────────────────

@testset "Performance -- solve_at_point" begin

    old_tfi   = prep_testing_objects(obj, stpfi=false)
    old_stpfi = prep_testing_objects(obj, stpfi=true)

    @test allocs_after_warmup(TFI.solve_at_point, old_tfi,   obj, TFIMethod())   == 0
    @test allocs_after_warmup(TFI.solve_at_point, old_stpfi, obj, STPFIMethod()) == 0

    # @inferred errors if the return type is not fully concrete (e.g. if pfct/vfct/tfct are Any).
    # Destructure to then confirm the field types are statically sized arrays, not abstract.
    res_tfi,   _ = @inferred(TFI.solve_at_point(old_tfi,   obj, TFIMethod()))
    res_stpfi, _ = @inferred(TFI.solve_at_point(old_stpfi, obj, STPFIMethod()))

    @test res_tfi.pfct   isa SVector
    @test res_tfi.vfct   isa SVector
    @test res_tfi.tfct   isa SMatrix

    @test res_stpfi.pfct isa SVector
    @test res_stpfi.vfct isa SVector
    @test res_stpfi.tfct isa SMatrix

end

@testset "Performance — unrolled_any" begin
    wide = NamedTuple{ntuple(i -> Symbol("v", i), 12)}(ntuple(float, 12))
    @test allocs_after_warmup(TFI.unrolled_any, isnan, wide) == 0
end
