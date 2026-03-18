using Test
using TransitionFunctionIteration
using ForwardDiff, StaticArrays, StructArrays, RecipesBase, Statistics, Random
const TFI = TransitionFunctionIteration

# ── shared fixtures ────────────────────────────────────────────────────────────
# All test files are included inside the outer @testset, so these bindings are
# visible to every sub-testset via Julia's scoping rules.

# 1D: y(x) = 2x, z(x) = x²
xgrid  = collect(1.0:5.0)
ssg1   = StateSpaceGrid(; x=xgrid)
fvals1 = StructArray([(y=2.0*xi, z=xi^2) for xi in xgrid])
fct1   = Interpolant(ssg1, fvals1)

# 2D: f(x,y) = x + y  (bilinear → linear interpolation is exact)
xgrid2 = collect(1.0:3.0)
ygrid2 = collect(10.0:10.0:30.0)
ssg2   = StateSpaceGrid(; x=xgrid2, y=ygrid2)
fvals2 = StructArray([(f=xi+yj,) for xi in xgrid2, yj in ygrid2])
fct2   = Interpolant(ssg2, fvals2)

# With exogenous state (NoInterp integer indexing) + endogenous (linear): f(exst,x) = exst*x
ssg_ex   = StateSpaceGrid(; exst=1:3, x=collect(1.0:3.0))
fvals_ex = StructArray([(y=Float64(ei)*xi,) for ei in 1:3, xi in 1.0:3.0])
fct_ex   = Interpolant(ssg_ex, fvals_ex)

# 3-state AR1 exogenous environment (iid, used across exog tests)
e3 = ExogenvAR1(n=3)

# Dummy model 
struct DummyModel <: AbstractModel end

function TransitionFunctionIteration.compute_next_period_values(next, ::DummyModel)
    expectations = (e = next.v,)
    continuations = (; c = next.x)
    NextPeriod(; expectations, continuations)
end

function TransitionFunctionIteration.calc_eqm(state, policy, ::DummyModel, np::NextPeriod, mode::EvaluationMode)
    fx = (fx1 = np.expectations.e - policy.p + state.x,)
    if mode isa SolveMode
        return fx
    elseif mode isa IterateMode
        trans = (x = 1.0,)
        vals = (v = 1.0,)
        return (; fx, trans, vals)
    end
end

pfct = Interpolant(ssg_ex, (p = 1.0,))
vfct = Interpolant(ssg_ex, (v = 2.0,))
tfct = Interpolant(ssg_ex, (x = SVector(0.9, 1.0, 1.1),))
obj = ModelSolution(model = DummyModel(), exogenv = e3, Pfct = pfct, Vfct = vfct, Tfct = tfct, ssg = ssg_ex)

function prep_testing_objects(obj; stpfi=false)
    state = first(points(obj.Tfct.ssg))
    pol = first(coefficients(obj.Pfct))
    val = first(coefficients(obj.Vfct))
    trans = first(coefficients(obj.Tfct))
    guess = deepcopy(pol)
    inittrans = deepcopy(trans)

    fullguess = stpfi ? (guess, inittrans) : guess

    old = TransitionFunctionIteration.FunctionValuesAtPoint(state, fullguess, val, inittrans)
    return old
end


# ── test suites ───────────────────────────────────────────────────────────────

@testset "TransitionFunctionIteration" begin
    include("test_statespace.jl")
    include("test_interpolant.jl")
    include("test_exogenv.jl")
    include("test_bounded.jl")
    include("test_momentarray.jl")
    include("test_economics.jl")
    include("test_simulation.jl")
    include("test_plotting.jl")
    include("test_performance.jl")
end
