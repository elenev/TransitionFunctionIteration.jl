
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

function fisher_burmeister(a::Real, b::Real)
    return a + b - sqrt(a^2 + b^2)
end


function compute_expectations(next::StructArray, prnext, obj::ModelSolution, mode::EvaluationMode)
    # Ensure probabilities are a statically sized SVector
    N = length(prnext)
    prnext = SVector{N}(prnext)

    # Returns an SVector of Tuples: Tuple(expectations_nt, continuations_nt)
    results = SVector(ntuple(i -> compute_next_period_values(next[i], obj.model), N))
    
    # Accumulate Expectations via Functional Fold
    # Initialize with the first element
    # First [1] refers to the first exogenous state, second [1] extracts the expectations NamedTuple
    init_exp = map(x -> x * prnext[1], results[1][1])
    
    # foldl forces strict left-to-right accumulation without mutation
    expectations = foldl(2:N; init=init_exp) do acc, i
        map((a, x) -> muladd(prnext[i], x, a), acc, results[i][1])
    end

    # Stack Continuations 
    # Extract just the continuations into an SVector of NamedTuples and make a StructArray
    continuations = StaticStructArray(map(res -> res[2], results))

    if mode isa SimulateMode
        addvars = StaticStructArray(map(res -> res[3], results))
    else
        addvars = (;)
    end
    
    return NextPeriod(
        expectations = expectations, 
        continuations = merge_structarrays(next, continuations),
        prnext = prnext,
        addvars = addvars
    )
end


value(x::AbstractFloat) = x 
value(x::ForwardDiff.Dual) = x.value  # Extract value from Dual numbers

function expand_transitions(trans, exogenv::AbstractExogenv)
    N = Val(exnpt(exogenv))
    expand_deterministic_transitions(t::Real) = SVector(ntuple(_ -> t, N))
    expand_deterministic_transitions(t::AbstractVector) = t
    return map(expand_deterministic_transitions, trans)
    #return StructArray(map(expand_deterministic_transitions, trans))
end

function update_transitions!(trans, new)
    new_vals = values(new)
    ntuple(Val(length(new_vals))) do i
        row_vec = new_vals[i]
        for j in axes(trans, 2)
            @inbounds trans[i, j] = row_vec[j]
        end
    end
end

function _compute_expectations_from_trans(trans::AbstractMatrix, prnext, obj::ModelSolution)
    named_trans = name_values(trans, obj.Tfct)
    next = forecast_next(named_trans, obj)
    return compute_expectations(next, prnext, obj, SolveMode()), named_trans
end

function solve_at_point!(policy, trans, vals, guess, state, obj, solmethod, init_trans)
    p_nof = nof(obj.Pfct)
    fullstate = merge(state, points(obj.exogenv)[state.exst])  # Merge in exogenous state variables for calc_eqm
    prnext = conditional_distro(obj.exogenv, state.exst)

    if solmethod isa TFIMethod
        # Pre-compute expectations using old (read-only) init_trans — avoids WAR aliasing with trans writes below
        np, _ = _compute_expectations_from_trans(SMatrix(init_trans), prnext, obj)
        sys = TFISystem(fullstate, np, obj)
    end

    if solmethod isa STPFIMethod
        # Add transition function residuals to the system of equations: consistency checking
        sys = STPFISystem(fullstate, prnext, p_nof, obj) 

        x0p = SVector(guess[1])
        x0t = SArray(guess[2])
        x0 = vcat(x0p, vec(x0t))
    elseif solmethod isa TFIMethod
        # Use pre-computed expectations, solve only for policy function residuals
        x0 = SVector(guess)
    end

    # Set up the problem
    
    problem = NonlinearProblem{false}(sys, x0, nothing)

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
        policy .= view(sol.u, 1:p_nof)
    else
        @failedpoint ("Solver did not converge: $(sol.retcode) at state $state")
    end

    # # Update guess and values with the solution
    out = calc_eqm(fullstate, name_values(policy, obj.Pfct), obj.model, np, IterateMode())


    out_trans = values(expand_transitions(out.trans, obj.exogenv))
    #out_trans = (K = SVector(1.0, 2.0, 3.0, 4.0, 5.0),) # Placeholder for testing, replace
    update_transitions!(trans, out_trans)
    # @debug "trans" typeof(trans)
   #trans .= SMatrix{1,5}(3.1, 3.2, 3.3, 3.4, 3.5) # Placeholder for testing, replace with actual update
    #trans .= values(out.trans)
    vals .= values(out.vals)

    return SciMLBase.successful_retcode(sol)
end

function find_neighbors_with_diagonals(i, R, n=1)
    I = R[i]
    Ifirst, Ilast = first(R), last(R)
    I1 = oneunit(Ifirst)*n
    range = max(Ifirst, I-I1):min(Ilast, I+I1)
    return setdiff(range,[I])
end

function perform_iteration!(obj, damp, npasses, solmethod)

    new_Tfct = deepcopy(coefficients(obj.Tfct))
    new_Vfct = deepcopy(coefficients(obj.Vfct))
    distances = [zeros(nof(obj.Vfct)) for _ in new_Vfct]
    distancesT = [zeros(nof(obj.Tfct), exnpt(obj.exogenv)) for _ in new_Tfct]
    success = falses(npt(obj.Vfct))

    states = collect(points(obj.ssg))
    Pfct_coeffs = coefficients(obj.Pfct)
    Vfct_coeffs = coefficients(obj.Vfct)
    Tfct_coeffs = coefficients(obj.Tfct)

    R = CartesianIndices(values(dims(obj.ssg)))

    guess_selector = i -> @SVector[R[i]]

    pass = 1
    while true
        unsolved_state_indices = findall(!s for s in success)
        total_unsolved = length(unsolved_state_indices)
        new_Pfct = deepcopy(coefficients(obj.Pfct))

        #Threads.@threads :dynamic 
        for idx in eachindex(unsolved_state_indices)
            i = unsolved_state_indices[idx]
            state = states[i]
            guess_indices = guess_selector(i)
            #guess = Pfct_coeffs[i]
            values = Vfct_coeffs[i]
            trans = Tfct_coeffs[i]
            new_trans = new_Tfct[i]
            new_values = new_Vfct[i]
            dist = distances[i]
            distT = distancesT[i]
            
            for guess_idx in guess_indices
                newpolicy = new_Pfct[guess_idx]
                if solmethod isa TFIMethod
                    oldguess = Pfct_coeffs[guess_idx]
                elseif solmethod isa STPFIMethod
                    oldguess = (Pfct_coeffs[guess_idx], Tfct_coeffs[i])
                end
                success[i] = solve_at_point!(newpolicy, new_trans, new_values, oldguess, state, obj, solmethod, trans)
                if success[i] break end
            end

            dist .= new_values - values

            #@debug "trans" trans
            #@debug "new_trans" new_trans
            distT .= new_trans - trans
            new_values .*= damp
            new_values .+= (1 - damp) * values

            new_trans .*= damp
            new_trans .+= (1 - damp) * trans

        end

        Pfct_coeffs .= new_Pfct

        guess_selector = i -> find_neighbors_with_diagonals(i, R)

        total_still_unsolved = sum(!s for s in success)
        solved_this_pass = total_unsolved - total_still_unsolved
        if pass==1 && total_still_unsolved>0
            if npasses > 1
                extramsg = " Revisiting..."
            else
                extramsg = ""
            end
            @iterationfailed "Solver failed at $total_still_unsolved points out of $(length(success)) points. $extramsg"
        elseif pass>1
            @iterationfailed "Pass $pass: Revisiting solved $(solved_this_pass) points. $total_still_unsolved unsolved points remaining."
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

    return success, distances, distancesT
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

"""
    main_execute!(obj; maxiter=100, tol=-5, damp=0.5, npasses=2, logging=ITERATIONFAILED)

Iterate policy and value functions of `obj` to convergence.

Solves equilibrium conditions at every grid point (threaded), damps the value function
update, and repeats until `log₁₀(max_distance) < tol` or `maxiter` is reached.

- `tol`: convergence threshold as log₁₀ distance (`-5` means 1e-5)
- `damp`: value function damping in [0,1] (1 = full update, 0 = no update)
- `npasses`: revisit passes for grid points where the solver failed on the first pass
- `logging`: minimum `LogLevel` to display (`COMPLETE` suppresses iteration logs)
"""
function main_execute!(obj; maxiter=100, tol=-5, damp=0.5, npasses=2, solmethod=TFIMethod(), logging=ITERATIONFAILED)

    iter = 1
    original_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)  # To avoid thread oversubscription

    exec = (; maxiter, tol, damp, npasses, solmethod)

    with_logger(ExecutionLogger(logging)) do

        while true 

            _, distancesV, distancesT = perform_iteration!(obj, damp, npasses, solmethod)

            distV, worst_varV, worst_pointV = summarize_distances(distancesV, obj.Vfct)
            distT, worst_varT, worst_pointT, worst_exst = summarize_distances(distancesT, obj.Tfct)

            @iteration "$iter: max V distance = $distV at point $worst_pointV in variable $worst_varV"
            @iteration "$iter: max T distance = $distT at point $worst_pointT in variable $worst_varT for exog. state $worst_exst"

            if iter >= maxiter
                @complete "Maximum iterations reached without convergence."
                exec = merge(exec, (; iter, distV, worst_varV, worst_pointV, distT, worst_varT, worst_pointT, worst_exst))
                return
            elseif distV < tol && distT < tol
                @complete "Converged."
                exec = merge(exec, (; iter, distV, worst_varV, worst_pointV, distT, worst_varT, worst_pointT, worst_exst))
                return
            end
            iter += 1
        end

    end
    BLAS.set_num_threads(original_threads)  # Restore

    return exec

end