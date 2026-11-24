
"""
    make_point(exst, trans) -> NamedTuple

Merge exogenous state index `exst` with transitions NamedTuple `trans`
into a full state point `(; exst, trans...)`.
"""
@inline make_point(exst, trans) = merge( (; exst,), trans)

"""
    forecast_at_point(i, trans, Vfct, exogenv, forecasts=(;)) -> NamedTuple

Evaluate the value function at exogenous state `i` and transition state `trans`,
returning a merged NamedTuple of state, value outputs and the current-period forecasting
values `forecasts` (the same at every next-period point).
"""
function forecast_at_point(i, trans, Vfct, exogenv, forecasts::NamedTuple=(;))
    point = make_point(i, trans)
    #point = (exst = 2, K = 3.0)
    #point = make_point(2, (K = 3.0,))
    vals = Vfct(point)
    return merge((exst = i,), points(exogenv)[i], trans, vals, forecasts)
end

"""
    forecast_next(trans, obj, forecasts=(;)) -> StructArray

Evaluate `obj.Vfct` at `trans` across all exogenous states in `obj.exogenv`, attaching the
current-period forecasting values `forecasts` (see `forecasts_from_policy`) to every point.
Returns a StructArray used to form expectations over next-period values.
"""
function forecast_next(trans, obj, forecasts::NamedTuple=(;))

    forecast_tuple = ntuple(
        i -> forecast_at_point(i, map( t -> t[i], trans), obj.Vfct, obj.exogenv, forecasts), 
        exnpt(obj.exogenv) )
    return  StaticStructArray( SVector(forecast_tuple) )
end

# Forecasting values as stored in Ffct coefficients (SVector) or `nothing` without an Ffct
@inline _forecast_svec(f::NamedTuple{()}) = nothing
@inline _forecast_svec(f::NamedTuple) = SVector(values(f))
@inline _named_forecasts(::Nothing, obj) = (;)
@inline _named_forecasts(f, obj) = name_values(f, obj.Ffct)



function normalize_next_period_values(np, mdl::AbstractModel, mode::EvaluationMode)
    # Placeholder for any normalization logic needed for expectations and continuations
    # after the initial taking of expecatations. E.g.,
    # the Epstein-Zin SDF normalizes the value function by its certainty equivalent
    # In the TFI method, both are constant over the solver's iterations, so it's efficient
    # to compute them in the compute_expectations call, rather than recompute in calc_eqm
    return np
end

"""
    compute_expectations(next, prnext, model, mode) -> NextPeriod

Aggregate next-period state values into expectations and continuations by calling
`compute_next_period_values` for each exogenous state and probability-weighting the results.
Returns a `NextPeriod` with weighted `expectations` and stacked `continuations`.
"""
function compute_expectations(next::StructArray, prnext, mdl::AbstractModel, mode::EvaluationMode)
    # Ensure probabilities are a statically sized SVector
    N = length(prnext)
    prnext = SVector{N}(prnext)

    # Returns an SVector of Tuples: Tuple(expectations_nt, continuations_nt)
    results = SVector(ntuple(i -> compute_next_period_values(next[i], mdl, mode), N))
    
    # Accumulate Expectations via Functional Fold
    # Initialize with the first element
    # First [1] refers to the first exogenous state, second [1] extracts the expectations NamedTuple
    init_exp = map(x -> x * prnext[1], first(results).expectations)
    
    # foldl forces strict left-to-right accumulation without mutation
    expectations = foldl(2:N; init=init_exp) do acc, i
        map((a, x) -> muladd(prnext[i], x, a), acc, results[i].expectations)
    end

    # Stack Continuations
    # Extract just the continuations into an SVector of NamedTuples and make a StructArray
    continuations = StaticStructArray(map(res -> res.continuations, results))

    np = NextPeriod(
        expectations = expectations,
        continuations = merge_structarrays(next, continuations),
        prnext = prnext,
    )

    return normalize_next_period_values(np, mdl, mode)
end



# function expand_transitions(trans, exogenv::AbstractExogenv)
#     N = Val(exnpt(exogenv))
#     expand_deterministic_transitions(t::Real) = SVector(ntuple(_ -> t, N))
#     expand_deterministic_transitions(t::AbstractVector) = t
#     return map(expand_deterministic_transitions, trans)
# end

function expand_transitions(trans, exogenv::AbstractExogenv)
    T = SVector{exnpt(exogenv)}
    expand_deterministic_transitions(t::Real) = fill(t, T)
    expand_deterministic_transitions(v::AbstractVector) = v
    return map(expand_deterministic_transitions, trans)
end

function _compute_expectations_from_trans(trans::AbstractMatrix, prnext, obj::ModelSolution, forecasts::NamedTuple=(;))
    named_trans = name_values(trans, obj.Tfct)
    return _compute_expectations_from_trans(named_trans, prnext, obj, forecasts)
end

function _compute_expectations_from_trans(named_trans, prnext, obj::ModelSolution, forecasts::NamedTuple=(;))
    next = forecast_next(named_trans, obj, forecasts)
    return compute_expectations(next, prnext, obj.model, SolveMode()), named_trans
end

_fx_to_svec(fx::NamedTuple) = SVector(values(fx))
_fx_to_svec(fx::Tuple) = SVector(fx)
_fx_to_svec(fx::StaticVector) = fx

struct FunctionValuesAtPoint{S, P,V,T,F}
    state::S
    pfct::P
    vfct::V    
    tfct::T
    ffct::F     # forecasting values (SVector), or nothing without an Ffct
end
FunctionValuesAtPoint(state, pfct, vfct, tfct) = FunctionValuesAtPoint(state, pfct, vfct, tfct, nothing)

abstract type EquilibriumSystem end

@inline get_fullstate(sys::EquilibriumSystem) = merge(sys.state, sys.exog)

struct TFISystem{S, E, N, M} <: EquilibriumSystem
    state::S
    exog::E
    np::N
    obj::M
end

function (sys::TFISystem)(u,p::EvaluationMode)
    fullstate = get_fullstate(sys)
    named_policy = name_values(u, sys.obj.Pfct)
    named_policy = map(to_constrained, named_policy, sys.obj.Pfct_bounds)

    if p isa SolveMode
        return _fx_to_svec(calc_eqm(fullstate, named_policy, sys.obj.model, sys.np, SolveMode()))
    else
        out = calc_eqm(fullstate, named_policy, sys.obj.model, sys.np, p)
        #@debug "out.trans: $(out.trans)"
        out_trans = expand_transitions(out.trans, sys.obj.exogenv)
        forecasts = _forecast_svec(forecasts_from_policy(fullstate, named_policy, sys.obj.model))
        return FunctionValuesAtPoint(sys.state, u, SVector(values(out.vals)), hcat(out_trans...)', forecasts)
    end
end

struct STPFISystem{S, E, R, M} <: EquilibriumSystem
    state::S
    exog::E
    prnext::R
    obj::M
end

function (sys::STPFISystem)(u,p::EvaluationMode)
    # Full state
    fullstate = get_fullstate(sys)

    # Construct policy named tuple
    p_nof = nof(sys.obj.Pfct)
    named_policy = name_values(view(u, 1:p_nof), sys.obj.Pfct)
    named_policy = map(to_constrained, named_policy, sys.obj.Pfct_bounds)

    # Construct transition matrix
    u_tfct = SMatrix{ nof(sys.obj.Tfct), exnpt(sys.obj.exogenv)}(view(u,(p_nof+1):length(u)))    
    u_tfct = to_constrained.(u_tfct, values(sys.obj.Tfct_bounds)) # order ensured to be the same by main_setup

    # Name it
    named_trans = name_values(u_tfct, sys.obj.Tfct)
    # forecasting values implied by the policy guess (STPFI: consistent within the solve)
    forecasts = forecasts_from_policy(fullstate, named_policy, sys.obj.model)
    np, _ = _compute_expectations_from_trans(named_trans, sys.prnext, sys.obj, forecasts)
    out = calc_eqm(fullstate, named_policy, sys.obj.model, np, IterateMode())
    out_trans = expand_transitions(out.trans, sys.obj.exogenv)

    transition_residuals = map(-,named_trans, out_trans)

    if p isa SolveMode
        return vcat(_fx_to_svec(out.fx), transition_residuals...)
    else
        return FunctionValuesAtPoint(sys.state, SVector{p_nof}(view(u, 1:p_nof)), SVector(values(out.vals)), hcat(out_trans...)',
                                     _forecast_svec(forecasts))
    end
end

function solve_at_point(old::FunctionValuesAtPoint, obj, solmethod)
    # Unpack initial values
    guess = old.pfct
    init_trans = old.tfct

    # Unpack state
    state = old.state

    exog = points(obj.exogenv)[state.exst]  # exogenous state variables
    prnext = conditional_distro(obj.exogenv, state.exst)

    if solmethod isa TFIMethod
        # Pre-compute expectations using old (read-only) init_trans — avoids WAR aliasing with trans writes below
        # TFI: the forecasting values are the stored ones from the previous iteration
        np, _ = _compute_expectations_from_trans(SMatrix(init_trans), prnext, obj, _named_forecasts(old.ffct, obj))
        sys = TFISystem(state, exog, np, obj)
    end

    if solmethod isa STPFIMethod
        # Add transition function residuals to the system of equations: consistency checking
        sys = STPFISystem(state, exog, prnext, obj) 

        x0p = SVector(guess[1])
        x0t = SArray(guess[2])
        x0t = to_unconstrained.(x0t, values(obj.Tfct_bounds)) # DOES NOT CHECK THAT ORDER IS THE SAME YET
        x0 = vcat(x0p, vec(x0t))
    elseif solmethod isa TFIMethod
        # Use pre-computed expectations, solve only for policy function residuals
        x0 = SVector(guess)
    end

    # Set up the problem
    
    problem = NonlinearProblem{false}(sys, x0, SolveMode())

    # Algorithm sequence
    # alg = NonlinearSolvePolyAlgorithm((
    #     SimpleTrustRegion(;  autodiff = AutoForwardDiff()),   # more robust globalization
    #     TrustRegion(;  autodiff = AutoForwardDiff()),         # robust globalization
    # ))
    # Solver by size of the full system (STPFI adds the transitions to the policies)
    if length(x0) <= 10
        alg = SimpleTrustRegion(;  autodiff = AutoForwardDiff())
    else
        alg = TrustRegion(;  autodiff = AutoForwardDiff())
    end

    # Solve the system of equations
    sol = solve(problem, alg)

    if SciMLBase.successful_retcode(sol)
        solution = sol.u
    else
        @logmsg FAILEDPOINT "" fp_state=string(state) fp_retcode=string(sol.retcode)
        solution = x0 # Revert to original guess if not solved
    end

    new = sys(solution, IterateMode())
    return new, SciMLBase.successful_retcode(sol)
end

function find_neighbors_with_diagonals(i, R, n=1)
    I = R[i]
    Ifirst, Ilast = first(R), last(R)
    I1 = oneunit(Ifirst)*n
    range = max(Ifirst, I-I1):min(Ilast, I+I1)
    return setdiff(range,[I])
end

const THREAD_THRESHOLD = 500  # parallelize when n_points * n_equations exceeds this

@inline _make_guess(guess_idx, _, Pfct_coeffs, Tfct_coeffs, ::TFIMethod)   = Pfct_coeffs[guess_idx]
@inline _make_guess(guess_idx, i, Pfct_coeffs, Tfct_coeffs, ::STPFIMethod) = (Pfct_coeffs[guess_idx], Tfct_coeffs[i])

@inline _guess_indices(i, R, ::Val{1}) = @SVector[R[i]]
@inline _guess_indices(i, R, ::Val)    = find_neighbors_with_diagonals(i, R)

function perform_iteration!(obj, diag, damp, npasses, solmethod, use_threads::Bool, damping::Symbol=:damped_iterate)

    #states = collect(points(obj.ssg))
    Pfct_coeffs = coefficients(obj.Pfct)
    Vfct_coeffs = coefficients(obj.Vfct)
    Tfct_coeffs = coefficients(obj.Tfct)
    Ffct_coeffs = _coefficients_or_nothing(obj.Ffct)

    # Scratch write buffers, pre-allocated once in IterationDiagnostics and
    # reused every iteration (no per-iteration allocation): every point gets
    # written in pass 1 below regardless of success (success is reset to
    # all-false just below), so nothing needs to be copied in first.
    new_Pfct = diag.next_pfct
    new_Tfct = diag.next_tfct
    new_Vfct = diag.next_vfct

    # "Old" reference for the damping recursion, set at the end of the previous iteration
    # (see `damping` in main_execute!):
    # - :damped_iterate: the current damped functions, so new = damp*x_t + (1-damp)*T(x_t)
    #   (standard relaxation; distances are fixed-point residuals T(x_t) - x_t).
    # - :previous_solve: the raw (undamped) solve of the previous iteration, so
    #   new = damp*T(x_{t-1}) + (1-damp)*T(x_t), matching the MATLAB `VF = VFnext`
    #   bookkeeping (main_execute_update.m). A one-step blend of two raw solves does not
    #   damp oscillations of the solution map: a period-2 oscillation in T persists.
    # Vfct_coeffs/Tfct_coeffs (the current, already-damped interpolant) remain the
    # expectations input to solve_at_point below in both cases.
    prev_vfct = diag.prev_vfct
    prev_tfct = diag.prev_tfct
    new_raw_Vfct = diag.next_raw_vfct
    new_raw_Tfct = diag.next_raw_tfct

    distances = diag.value_distances
    distancesT = diag.transition_distances

    success = diag.success
    success .= false  # Reset success flags for this iteration

    R = CartesianIndices(values(dims(obj.ssg)))

    pass = 1
    while true
        unsolved_state_indices = findall(!s for s in success)
        total_unsolved = length(unsolved_state_indices)


        function _solve_loop!(idx, pass)
            i = unsolved_state_indices[idx]
            state = obj.ssg[i]
            guess_indices = _guess_indices(i, R, Val(pass))
            values = Vfct_coeffs[i]
            trans = Tfct_coeffs[i]
            forecasts = _coef_at(Ffct_coeffs, i)

            local new
            for guess_idx in guess_indices
                oldguess = _make_guess(guess_idx, i, Pfct_coeffs, Tfct_coeffs, solmethod)
                old = FunctionValuesAtPoint(state, oldguess, values, trans, forecasts)
                new, success[i] = solve_at_point(old, obj, solmethod)
                if success[i] break end
            end

            # Still unsolved after a later pass, whose guesses come from neighbors: keep the pass-1
            # result, computed from the point's own guess
            if pass > 1 && !success[i]
                return
            end

            new_Pfct[i] = new.pfct
            distances[i] .= new.vfct - prev_vfct[i]
            new_Vfct[i] = @. damp * prev_vfct[i] + (1-damp) * new.vfct
            distancesT[i] .= new.tfct - prev_tfct[i]
            new_Tfct[i] = @. damp * prev_tfct[i] + (1-damp) * new.tfct
            new_raw_Vfct[i] = new.vfct
            new_raw_Tfct[i] = new.tfct
            _update_forecasts_at!(diag, i, new.ffct, damp)
        end

        threaded_foreach(idx -> _solve_loop!(idx, pass), eachindex(unsolved_state_indices); use_threads)

        Pfct_coeffs .= new_Pfct

        total_still_unsolved = sum(!s for s in success)
        solved_this_pass = total_unsolved - total_still_unsolved
        if pass==1 && total_still_unsolved>0
            @logmsg ITERATIONFAILED "" pass=pass n_failed=total_still_unsolved n_total=length(success) revisiting_on=npasses>1
        elseif pass>1
            @logmsg ITERATIONFAILED "" pass=pass n_solved=solved_this_pass n_unsolved=total_still_unsolved
        end

        if all(success)
            # Solved all points, no need for next pass
            break
        elseif pass >= npasses
            # Maximum passes reached, stop
            break
        elseif sum(success) == total_unsolved
            # Pass didn't solve any points, no point for next pass
            pass += 1
        end
        # If we reach here, we need another pass
        pass += 1

    end

    Vfct_coeffs .= new_Vfct
    Tfct_coeffs .= new_Tfct
    if damping === :damped_iterate
        prev_vfct .= new_Vfct
        prev_tfct .= new_Tfct
    else
        prev_vfct .= new_raw_Vfct
        prev_tfct .= new_raw_Tfct
    end
    _finish_forecasts!(Ffct_coeffs, diag, damping)

    # println(coefficients(obj.Vfct)[5,1,1])
    # println(coefficients(obj.Vfct)[5,5,1])
    # println(coefficients(obj.Vfct)[5,10,1])
end

# Forecasting function bookkeeping: no-ops when the model has no Ffct
_coefficients_or_nothing(::Nothing) = nothing
_coefficients_or_nothing(fct) = coefficients(fct)
@inline _coef_at(::Nothing, i) = nothing
@inline _coef_at(coefs, i) = coefs[i]

@inline _update_forecasts_at!(diag, i, ::Nothing, damp) = nothing
@inline function _update_forecasts_at!(diag, i, newf, damp)
    diag.forecast_distances[i] .= newf - diag.prev_ffct[i]
    diag.next_ffct[i] = @. damp * diag.prev_ffct[i] + (1-damp) * newf
    diag.next_raw_ffct[i] = newf
    return nothing
end

_finish_forecasts!(::Nothing, diag, damping) = nothing
function _finish_forecasts!(Ffct_coeffs, diag, damping)
    Ffct_coeffs .= diag.next_ffct
    diag.prev_ffct .= damping === :damped_iterate ? diag.next_ffct : diag.next_raw_ffct
    return nothing
end

function summarize_distances(distances, Fct)
    distances = stack(distances)
    distances .= abs.(distances)

    #@debug "distances size" size(distances)
    idxmax = Tuple(argmax(distances))

    #@debug "idxmax " idxmax
    max_distance = distances[idxmax...]
    dist = log10(max_distance)

    var_idx = idxmax[1]
    state_idx = idxmax[end-ndims(Fct.ssg)+1:end]
    

    worst_var = varnames(Fct)[var_idx]
    worst_point = map((statevar, idx) -> statevar[idx], 
        unigrids(Fct.ssg), 
        NamedTuple{varnames(Fct.ssg)}(state_idx) )

    if length(idxmax) > ndims(Fct.ssg)+1
        exst_idx = idxmax[2]
    else
        exst_idx = nothing
    end

    return dist, worst_var, worst_point, exst_idx
end


struct IterationDiagnostics{B,V,T,PV,PT,NP,NV,NT2,RV,RT,FD,PF,NF,RF}
    success::B
    value_distances::V
    transition_distances::T
    prev_vfct::PV
    prev_tfct::PT
    next_pfct::NP
    next_vfct::NV
    next_tfct::NT2
    next_raw_vfct::RV
    next_raw_tfct::RT
    # forecasting function (all nothing without an Ffct)
    forecast_distances::FD
    prev_ffct::PF
    next_ffct::NF
    next_raw_ffct::RF
end

function IterationDiagnostics(obj)
    similar_to_coef = similar ∘ first ∘ coefficients
    vcoefs = coefficients(obj.Vfct)
    vdistances = [similar_to_coef(obj.Vfct) for _ in vcoefs]

    tcoefs = coefficients(obj.Tfct)
    tdistances = [similar_to_coef(obj.Tfct) for _ in tcoefs]

    pcoefs = coefficients(obj.Pfct)

    success = falses(npt(obj.Vfct))

    # Seeded with the initial guess, matching MATLAB's `VF = obj.Vfct.evaluateAt(gridSt)`
    # before its iteration loop starts (main_execute_update.m).
    prev_vfct = deepcopy(vcoefs)
    prev_tfct = deepcopy(tcoefs)

    # Scratch write buffers for perform_iteration!, pre-allocated once and
    # reused (mutated in place) every iteration for the life of this
    # IterationDiagnostics. Safe uninitialized (`similar`, not `deepcopy`):
    # every point is written to exactly once in pass 1 of every iteration
    # (success resets to all-false each iteration), so no prior contents ever
    # need to be copied in.
    next_pfct = similar(pcoefs)
    next_vfct = similar(vcoefs)
    next_tfct = similar(tcoefs)
    next_raw_vfct = similar(vcoefs)
    next_raw_tfct = similar(tcoefs)

    if isnothing(obj.Ffct)
        fdistances = prev_ffct = next_ffct = next_raw_ffct = nothing
    else
        fcoefs = coefficients(obj.Ffct)
        fdistances = [similar_to_coef(obj.Ffct) for _ in fcoefs]
        prev_ffct = deepcopy(fcoefs)
        next_ffct = similar(fcoefs)
        next_raw_ffct = similar(fcoefs)
    end

    return IterationDiagnostics(success, vdistances, tdistances, prev_vfct, prev_tfct,
        next_pfct, next_vfct, next_tfct, next_raw_vfct, next_raw_tfct,
        fdistances, prev_ffct, next_ffct, next_raw_ffct)
end

# Construct a logger from either a LogLevel (builds TerminalLogger) or a pre-built ExecutionLogger.
_make_logger(level::LogLevel,         obj::AbstractModelSolution) = TerminalLogger(level, obj)
_make_logger(logger::ExecutionLogger, ::AbstractModelSolution)    = logger

# Infer the concrete type of a history snapshot from the interpolant and requested variable names.
function _snap_type(fct::AbstractInterpolant, var_names::Tuple)
    ET = eltype(eltype(coefficients(fct)))
    NamedTuple{var_names, NTuple{length(var_names), Vector{ET}}}
end

"""
    main_execute!(obj; maxiter=100, tol=-5, damp=0.5, npasses=2, solmethod=TFIMethod(),
                  damping=:damped_iterate, logging=ITERATIONFAILED, track_history=nothing,
                  multithreaded=nothing)

Iterate policy and value functions of `obj` to convergence.

Solves equilibrium conditions at every grid point, damps the value function update, and
repeats until one of three exit states, each logged at level `COMPLETE` and returned as `status`:

- `:converged`: the largest change of the value, transition and forecasting functions in the last
  iteration is below `10^tol`;
- `:maxiter`: `maxiter` iterations ran without convergence;
- `:all_failed`: the solver failed at every grid point in the last iteration, so the policy was not
  updated anywhere.

- `tol`: convergence threshold as log₁₀ distance (`-5` means 1e-5)
- `damp`: damping weight in [0,1], applied to value and transition functions
  (0 = full update, 1 = no update)
- `damping`: what `damp` weights against.
  `:damped_iterate` (default): the current damped iterate, `x_{t+1} = damp*x_t + (1-damp)*T(x_t)`
  (standard relaxation; damps oscillations of the solution map, distances are fixed-point
  residuals). `:previous_solve`: the previous iteration's raw solve,
  `x_{t+1} = damp*T(x_{t-1}) + (1-damp)*T(x_t)`, as in the MATLAB implementation
  (main_execute_update.m); a period-2 oscillation of `T` persists under this rule.
- `npasses`: revisit passes for grid points where the solver failed on the first pass
- `solmethod`: `TFIMethod()` (default) or `STPFIMethod()`, controls whether to precompute
  expectations (TFI) or to add transition functions to the system (STPFI, Cao, Luo, and Nie)
- `multithreaded`: `nothing` (default) threads over grid points when more than one thread
  is available and the problem is large enough; `true`/`false` forces the choice
- `logging`: minimum `LogLevel` to display (`COMPLETE` suppresses iteration logs)
- `track_history`: `nothing` (default) or `(fct_sym, var_names)` where `fct_sym` is `:Vfct`
  or `:Tfct` and `var_names` is a tuple of variable name Symbols. When set, returns a
  StructVector of per-iteration snapshots, one entry per iteration, each a NamedTuple of
  full-grid arrays for the requested variables.

Returns a NamedTuple with the settings above, `status`, the number of iterations `iter`, the final
`log₁₀` distances of the value, transition and forecasting functions (`distV`, `distT`, `distF`) with
the variables and grid points where they are largest, and `history` (`nothing` unless
`track_history` is set). Judge the accuracy of a solution by its simulated errors
([`calc_error_stats`](@ref)), not by `status` alone.
"""
function main_execute!(obj; maxiter=100, tol=-5, damp=0.5, npasses=2, solmethod=TFIMethod(),
                        damping::Symbol=:damped_iterate,
                        logging=ITERATIONFAILED, track_history=nothing, multithreaded=nothing)
    damping in (:damped_iterate, :previous_solve) ||
        throw(ArgumentError("damping must be :damped_iterate or :previous_solve, got :$damping"))

    iter = 1

    diag = IterationDiagnostics(obj)
    exec = (; maxiter, tol, damp, npasses, solmethod, damping)
    if isnothing(track_history)
        history_list = nothing
    else
        # History is a NameTuple of tuples indicating what to track 
        # E.g. track_history = (val_fcn_around_ss = (:Vfct, :V, 300), first_Ktrans_fcn = (:Tfct, :K, 1))
        history_list = map(track_history) do track_var
            coefs = coefficients(getproperty(obj, track_var[1]))
            T = eltype(coefs) <: AbstractVector ? eltype(first(coefs)) : typeof(first(eachrow(coefs)))
            return Vector{T}(undef,maxiter)
        end
    end

    #Main.@infiltrate

    # Threading heuristic: parallelize when problem is large enough and threads are available
    n_points = npt(obj.Vfct)
    n_eqs    = nof(obj.Pfct)
    use_threads = multithreaded === nothing ?
        (Threads.nthreads() > 1) && (n_points * n_eqs > THREAD_THRESHOLD) :
        multithreaded::Bool

    with_logger(_make_logger(logging, obj)) do
        while true

            perform_iteration!(obj, diag, damp, npasses, solmethod, use_threads, damping)

            distV, worst_varV, worst_pointV = summarize_distances(diag.value_distances, obj.Vfct)
            distT, worst_varT, worst_pointT, worst_exst = summarize_distances(diag.transition_distances, obj.Tfct)
            # forecasting function: its change counts toward convergence like Vfct's
            distF, worst_varF, worst_pointF = isnothing(obj.Ffct) ? (-Inf, nothing, nothing) :
                summarize_distances(diag.forecast_distances, obj.Ffct)[1:3]

            if isnothing(obj.Ffct)
                @logmsg ITERATION "" iter=iter max_vdist=distV max_tdist=distT worst_vvar=worst_varV worst_tvar=worst_varT worst_vpoint=worst_pointV worst_tpoint=worst_pointT worst_texst=worst_exst
            else
                @logmsg ITERATION "" iter=iter max_vdist=distV max_tdist=distT worst_vvar=worst_varV worst_tvar=worst_varT worst_vpoint=worst_pointV worst_tpoint=worst_pointT worst_texst=worst_exst max_fdist=distF worst_fvar=worst_varF worst_fpoint=worst_pointF
            end

            # Capture per-iteration snapshot if requested
            if !isnothing(track_history)
                foreach(track_history, history_list) do track_var, history_var
                    (fct_sym, var_name, idx) = track_var
                    fct = getproperty(obj, fct_sym)
                    if isa(idx, Integer)
                        coefs = coefficients(fct)
                        snapshot = coefs[idx]  # Extract the requested variable's full grid array at this iteration
                    elseif isa(idx, NamedTuple)
                        snapshot = SVector(fct(idx)...)
                    end

                     # Append to the corresponding history list
                     history_idx = varindex(fct, var_name)
                     history_var[iter] = snapshot[history_idx]
                end
            end

            exec_update = (; iter, distV, worst_varV, worst_pointV, distT, worst_varT, worst_pointT, worst_exst,
                           distF, worst_varF, worst_pointF)

            # Exit states: every point failed (the policy was not updated anywhere, so iterating
            # further cannot help), maximum iterations, or convergence
            if !any(diag.success)
                @logmsg COMPLETE "Solver failed at every grid point."
                exec = merge(exec, exec_update, (; status = :all_failed))
                break
            elseif iter >= maxiter
                @logmsg COMPLETE "Maximum iterations reached without convergence." 
                exec = merge(exec, exec_update, (; status = :maxiter))
                break
            elseif distV < tol && distT < tol && distF < tol
                @logmsg COMPLETE "Converged."
                exec = merge(exec, exec_update, (; status = :converged))
                break
            end
            iter += 1
        end
    end

    if isnothing(history_list) || isempty(history_list)
        history = nothing
    else
        history = StructArray(history_list)
        resize!(history, iter)
    end
    return (; exec..., history)

end