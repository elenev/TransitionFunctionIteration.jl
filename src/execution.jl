
"""
    make_point(exst, trans) -> NamedTuple

Merge exogenous state index `exst` with transitions NamedTuple `trans`
into a full state point `(; exst, trans...)`.
"""
@inline make_point(exst, trans) = merge( (; exst,), trans)

"""
    forecast_at_point(i, trans, Vfct) -> NamedTuple

Evaluate the value function at exogenous state `i` and transition state `trans`,
returning a merged NamedTuple of state and value outputs.
"""
function forecast_at_point(i, trans, Vfct, exogenv)
    point = make_point(i, trans)
    #point = (exst = 2, K = 3.0)
    #point = make_point(2, (K = 3.0,))
    vals = Vfct(point)
    return merge((exst = i,), points(exogenv)[i], trans, vals)
end

"""
    forecast_next(trans, obj) -> StructArray

Evaluate `obj.Vfct` at `trans` across all exogenous states in `obj.exogenv`.
Returns a StructArray used to form expectations over next-period values.
"""
function forecast_next(trans, obj)

    forecast_tuple = ntuple(
        i -> forecast_at_point(i, map( t -> t[i], trans), obj.Vfct, obj.exogenv), 
        exnpt(obj.exogenv) )
    return  StaticStructArray( SVector(forecast_tuple) )
end



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

function _compute_expectations_from_trans(trans::AbstractMatrix, prnext, obj::ModelSolution)
    named_trans = name_values(trans, obj.Tfct)
    return _compute_expectations_from_trans(named_trans, prnext, obj)
end

function _compute_expectations_from_trans(named_trans, prnext, obj::ModelSolution)
    next = forecast_next(named_trans, obj)
    return compute_expectations(next, prnext, obj.model, SolveMode()), named_trans
end

_fx_to_svec(fx::NamedTuple) = SVector(values(fx))
_fx_to_svec(fx::Tuple) = SVector(fx)
_fx_to_svec(fx::StaticVector) = fx

struct FunctionValuesAtPoint{S, P,V,T}
    state::S
    pfct::P
    vfct::V    
    tfct::T
end

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
        return FunctionValuesAtPoint(sys.state, u, SVector(values(out.vals)), hcat(out_trans...)')
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
    np, _ = _compute_expectations_from_trans(named_trans, sys.prnext, sys.obj)
    out = calc_eqm(fullstate, named_policy, sys.obj.model, np, IterateMode())
    out_trans = expand_transitions(out.trans, sys.obj.exogenv)

    transition_residuals = map(-,named_trans, out_trans)

    if p isa SolveMode
        return vcat(_fx_to_svec(out.fx), transition_residuals...)
    else
        return FunctionValuesAtPoint(sys.state, SVector{p_nof}(view(u, 1:p_nof)), SVector(values(out.vals)), hcat(out_trans...)')
    end
end

function solve_at_point(old::FunctionValuesAtPoint, obj, solmethod)
    # Unpack initial values
    guess = old.pfct
    init_trans = old.tfct

    # Unpack state
    state = old.state

    #policy, trans, vals, guess, state,  , init_trans
    p_nof = nof(obj.Pfct)
    exog = points(obj.exogenv)[state.exst]  # exogenous state variables
    prnext = conditional_distro(obj.exogenv, state.exst)

    if solmethod isa TFIMethod
        # Pre-compute expectations using old (read-only) init_trans — avoids WAR aliasing with trans writes below
        np, _ = _compute_expectations_from_trans(SMatrix(init_trans), prnext, obj)
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
    if p_nof <= 10
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

function perform_iteration!(obj, diag, damp, npasses, solmethod, use_threads::Bool)

    #states = collect(points(obj.ssg))
    Pfct_coeffs = coefficients(obj.Pfct)
    Vfct_coeffs = coefficients(obj.Vfct)
    Tfct_coeffs = coefficients(obj.Tfct)

    new_Pfct = deepcopy(Pfct_coeffs)
    new_Tfct = deepcopy(Tfct_coeffs)
    new_Vfct = deepcopy(Vfct_coeffs)    

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

            local new
            for guess_idx in guess_indices
                oldguess = _make_guess(guess_idx, i, Pfct_coeffs, Tfct_coeffs, solmethod)
                old = FunctionValuesAtPoint(state, oldguess, values, trans)
                new, success[i] = solve_at_point(old, obj, solmethod)
                if success[i] break end
            end

            new_Pfct[i] = new.pfct
            distances[i] .= new.vfct - values
            new_Vfct[i] = @. damp * values + (1-damp) * new.vfct
            distancesT[i] .= new.tfct - trans
            new_Tfct[i] = @. damp * trans + (1-damp) * new.tfct
        end

        if use_threads
            @batch per=thread for idx in eachindex(unsolved_state_indices)
                _solve_loop!(idx, pass)
            end
        else
            for idx in eachindex(unsolved_state_indices)
                _solve_loop!(idx, pass)
            end
        end

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

    # println(coefficients(obj.Vfct)[5,1,1])
    # println(coefficients(obj.Vfct)[5,5,1])
    # println(coefficients(obj.Vfct)[5,10,1])
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


struct IterationDiagnostics{B,V,T}
    success::B    
    value_distances::V 
    transition_distances::T
end

function IterationDiagnostics(obj)
    similar_to_coef = similar ∘ first ∘ coefficients
    vcoefs = coefficients(obj.Vfct)
    vdistances = [similar_to_coef(obj.Vfct) for _ in vcoefs]

    tcoefs = coefficients(obj.Tfct)
    tdistances = [similar_to_coef(obj.Tfct) for _ in tcoefs]

    success = falses(npt(obj.Vfct))

    return IterationDiagnostics(success, vdistances, tdistances)
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
    main_execute!(obj; maxiter=100, tol=-5, damp=0.5, npasses=2, logging=ITERATIONFAILED,
                  track_history=nothing)

Iterate policy and value functions of `obj` to convergence.

Solves equilibrium conditions at every grid point, damps the value function update, and
repeats until `log₁₀(max_distance) < tol` or `maxiter` is reached.

- `tol`: convergence threshold as log₁₀ distance (`-5` means 1e-5)
- `damp`: value function damping in [0,1] (1 = full update, 0 = no update)
- `npasses`: revisit passes for grid points where the solver failed on the first pass
- `solmethod`: `TFIMethod()` (default) or `STPFIMethod()`, controls whether to precompute
  expectations (Elenev and Landvoigt) or to add transition functions to the system (Cao et al)
- `logging`: minimum `LogLevel` to display (`COMPLETE` suppresses iteration logs)
- `track_history`: `nothing` (default) or `(fct_sym, var_names)` where `fct_sym` is `:Vfct`
  or `:Tfct` and `var_names` is a tuple of variable name Symbols. When set, returns a
  StructVector of per-iteration snapshots, one entry per iteration, each a NamedTuple of
  full-grid arrays for the requested variables.

Returns a NamedTuple with convergence info and, when `track_history` is set, a `history`
StructVector field.
"""
function main_execute!(obj; maxiter=100, tol=-5, damp=0.5, npasses=2, solmethod=TFIMethod(),
                        logging=ITERATIONFAILED, track_history=nothing, multithreaded=nothing)

    iter = 1
    original_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)  # To avoid thread oversubscription

    diag = IterationDiagnostics(obj)
    exec = (; maxiter, tol, damp, npasses, solmethod)
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

            perform_iteration!(obj, diag, damp, npasses, solmethod, use_threads)

            distV, worst_varV, worst_pointV = summarize_distances(diag.value_distances, obj.Vfct)
            distT, worst_varT, worst_pointT, worst_exst = summarize_distances(diag.transition_distances, obj.Tfct)

            @logmsg ITERATION "" iter=iter max_vdist=distV max_tdist=distT worst_vvar=worst_varV worst_tvar=worst_varT worst_vpoint=worst_pointV worst_tpoint=worst_pointT worst_texst=worst_exst

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

            exec_update = (; iter, distV, worst_varV, worst_pointV, distT, worst_varT, worst_pointT, worst_exst)

            if iter >= maxiter
                @logmsg COMPLETE "Maximum iterations reached without convergence." 
                exec = merge(exec, exec_update)
                break
            elseif distV < tol && distT < tol
                @logmsg COMPLETE "Converged."
                exec = merge(exec, exec_update)
                break
            end
            iter += 1
        end
    end

    BLAS.set_num_threads(original_threads)  # Restore
    if isnothing(history_list) || isempty(history_list)
        history = nothing
    else
        history = StructArray(history_list)
        resize!(history, iter)
    end
    return (; exec..., history)

end