using TransitionFunctionIteration
#using StructArrays, StatsBase

function create_test_simvector()
    X = rand(10_000)
    Y = rand(10_000)
    simvector = StructVector((X=X, Y=Y))
    return simvector
end

function create_test_structarray()
    sv1 = create_test_simvector()
    sv2 = create_test_simvector()
    sa = hcat(sv1, sv2)  # This creates a StructArray of StructVectors
    return sa
end

# teststatfcn = (mean = mean, std = std, ac = TransitionFunctionIteration.ac1);
# tss = TransitionFunctionIteration.StatisticsSet(teststatfcn);
# mss = TransitionFunctionIteration.StatisticsSet((:corY,), (TransitionFunctionIteration.MultivariateStatFun{:Y}(cor),));
# sm = create_test_structarray();
# mask = trues(size(sm));

function run_many_times(f, N)
    tmp = ntuple(i -> f(), N)
    return nothing
end

# ds() = TransitionFunctionIteration._describe_sample(runstats, sm, tss, mss, mask);

# runstats = Matrix{Float64}(undef,2,4);
# sv = sm[1,:];
# dr!() = TransitionFunctionIteration._describe_run!(runstats, sv, tss, mss)

# ── STPFIMethod diagnostics ───────────────────────────────────────────────────
# Run after: include("examples/ngm.jl")   (brings ngm, obj into scope)

function test_stpfi_isolated()
    ngm2 = NeoclassicalGrowthModel()
    obj2, _ = main_setup(ngm2)
    main_execute!(obj2; maxiter=200, logging=COMPLETE, solmethod=STPFIMethod())
    result = isapprox(coefficients(obj2.Pfct), coefficients(obj.Pfct), rtol=1e-3)
    println("Isolated: ", result ? "PASS" : "FAIL")
    return result
end

function test_stpfi_blas_threads(n)
    BLAS.set_num_threads(n)
    ngm2 = NeoclassicalGrowthModel()
    obj2, _ = main_setup(ngm2)
    main_execute!(obj2; maxiter=200, logging=COMPLETE, solmethod=STPFIMethod())
    result = isapprox(coefficients(obj2.Pfct), coefficients(obj.Pfct), rtol=1e-3)
    println("BLAS threads=$n: ", result ? "PASS" : "FAIL")
    BLAS.set_num_threads(Sys.CPU_THREADS)
    return result
end

function test_stpfi_threads(use_threads)
    ngm2 = NeoclassicalGrowthModel()
    obj2, _ = main_setup(ngm2)
    main_execute!(obj2; maxiter=200, logging=COMPLETE, solmethod=STPFIMethod(),
                  multithreaded=use_threads)
    result = isapprox(coefficients(obj2.Pfct), coefficients(obj.Pfct), rtol=1e-3)
    println("use_threads=$use_threads: ", result ? "PASS" : "FAIL")
    return result
end

function test_stpfi_after_tfi()
    ngm1 = NeoclassicalGrowthModel()
    obj1, _ = main_setup(ngm1)
    main_execute!(obj1; maxiter=200, logging=COMPLETE)
    ngm2 = NeoclassicalGrowthModel()
    obj2, _ = main_setup(ngm2)
    main_execute!(obj2; maxiter=200, logging=COMPLETE, solmethod=STPFIMethod())
    result = isapprox(coefficients(obj2.Pfct), coefficients(obj1.Pfct), rtol=1e-3)
    println("After TFI solve: ", result ? "PASS" : "FAIL")
    return result
end

function test_stpfi_verbose()
    ngm2 = NeoclassicalGrowthModel()
    obj2, _ = main_setup(ngm2)
    try
        main_execute!(obj2; maxiter=200, logging=ITERATION, solmethod=STPFIMethod())
        result = isapprox(coefficients(obj2.Pfct), coefficients(obj.Pfct), rtol=1e-3)
        println("Converged: ", result ? "match" : "different coefficients")
    catch e
        println("Exception: ", typeof(e))
        showerror(stdout, e)
        println()
    end
end

# ── Dual-number boundary interpolation test ───────────────────────────────────
# Run after: include("examples/ngm.jl")
# Tests whether Vfct (with Flat() extrapolation) handles ForwardDiff Duals
# at the upper K boundary without throwing.
using ForwardDiff
function test_dual_boundary()

    K_max = last(knots(obj.Vfct)[2])   # upper boundary of K grid
    exst  = 1                           # any valid exogenous state index

    println("K_max = ", K_max)

    # 1. Evaluate at exact boundary with a plain Float64 — should work
    try
        v = obj.Vfct((exst=exst, K=K_max))
        println("Float64 at boundary: OK, C = ", v.C)
    catch e
        println("Float64 at boundary: THREW ", typeof(e))
    end

    # 2. Evaluate at exact boundary with a Dual (value=K_max, partial=1) — this is
    #    what ForwardDiff does when differentiating w.r.t. K at the boundary
    K_dual = ForwardDiff.Dual(K_max, 1.0)
    try
        v = obj.Vfct((exst=exst, K=K_dual))
        println("Dual at boundary:    OK, C = ", ForwardDiff.value(v.C))
    catch e
        println("Dual at boundary:    THREW ", typeof(e))
        showerror(stdout, e)
        println()
    end

    # 3. Evaluate slightly above boundary with a Dual — pure extrapolation
    K_above = K_max + 1e-8
    K_dual_above = ForwardDiff.Dual(K_above, 1.0)
    try
        v = obj.Vfct((exst=exst, K=K_dual_above))
        println("Dual above boundary: OK, C = ", ForwardDiff.value(v.C))
    catch e
        println("Dual above boundary: THREW ", typeof(e))
        showerror(stdout, e)
        println()
    end
end