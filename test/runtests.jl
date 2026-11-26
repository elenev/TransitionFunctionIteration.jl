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

function TransitionFunctionIteration.compute_next_period_values(next, ::DummyModel, ::EvaluationMode)
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
    elseif mode isa SimulateMode
        # Added for augment_errors/calc_result testing (test_simulation.jl) -- no existing test
        # exercised SimulateMode on DummyModel before this, so this branch is purely additive.
        trans = (x = 1.0,)
        return (; fx, trans, simvars = (v = 1.0,))
    end
end

pfct = Interpolant(ssg_ex, (p = 1.0,))
vfct = Interpolant(ssg_ex, (v = 2.0,))
tfct = Interpolant(ssg_ex, (x = SVector(0.9, 1.0, 1.1),))
obj = ModelSolution(model = DummyModel(), exogenv = e3, Pfct = pfct, Vfct = vfct, Tfct = tfct, ssg = ssg_ex)

function prep_testing_objects(obj; idx=1, stpfi=false)
    state = obj.ssg[idx...]
    pol = coefficients(obj.Pfct)[idx...]
    val = coefficients(obj.Vfct)[idx...]
    trans = coefficients(obj.Tfct)[idx...]
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
    include("test_stubs.jl")
    include("test_template.jl")
    include("test_momentarray.jl")
    include("test_moments.jl")
    include("test_economics.jl")
    include("test_simulation.jl")
    include("test_logging.jl")
    include("test_plotting.jl")
    include("test_performance.jl")
    include("test_irfplotting.jl")
    include("test_execution.jl")
    include("test_utils.jl")
    include("test_transformations.jl")
    include("test_ngm.jl")       # must come first: brings ngm, obj, stv, simseries into global scope
    include("test_irf.jl")       # uses obj, stv, simseries from test_ngm.jl
    include("test_linearize.jl") # uses those fixtures
    include("test_setup.jl")     # uses ngm, obj, stv from test_ngm.jl
    include("test_forecasts.jl") # uses ngm, obj from test_ngm.jl
    include("test_eqmvalguesses.jl") # uses ngm from test_ngm.jl
    include("test_tables.jl")    # uses obj, stv, simseries, errs from test_ngm.jl
    include("test_portable.jl")  # uses obj, stv from test_ngm.jl and ForecastNGM from test_forecasts.jl
    include("test_outputvault.jl") # uses those fixtures and ForecastNGM from test_forecasts.jl
    include("test_Aqua.jl")
end
