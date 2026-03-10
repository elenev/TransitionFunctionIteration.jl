
function fisher_burmeister(a::Real, b::Real)
    return a + b - sqrt(a^2 + b^2)
end

value(x::AbstractFloat) = x 
value(x::ForwardDiff.Dual) = x.value  # Extract value from Dual numbers

@inline make_point(exst, trans) = merge( (; exst,), trans)

function forecast_at_point(i, trans, Vfct)
    vals = Vfct(make_point(i, trans))
    return merge((exst = i,), trans, vals)
end

function forecast_next(trans, obj)
    exog_iterator = eachindex(obj.exogenv.pts)
    return StructArray( map(i -> forecast_at_point(i, trans, obj.Vfct), exog_iterator))
    #return StructArray(forecast_at_point(i,trans,obj.Vfct) for i in exog_iterator)
end

function solve_at_point!(policy, vals, guess, state, obj)

    # Define the function to solve
    f = (u,p) -> calc_eqm(state, name_values(u, obj.Pfct), obj)

    # Set up the problem
    problem = NonlinearProblem(f, SVector(guess), nothing)

    # Algorithm sequence
    # alg = NonlinearSolvePolyAlgorithm((
    #     SimpleTrustRegion(;  autodiff = AutoForwardDiff()),   # more robust globalization
    #     TrustRegion(;  autodiff = AutoForwardDiff()),         # robust globalization
    # ))
    alg = TrustRegion(;  autodiff = AutoForwardDiff())

    # Solve the system of equations
    sol = solve(problem, alg)

    if SciMLBase.successful_retcode(sol)
        policy .= sol.u
    else
        @failedpoint ("Solver did not converge: $(sol.retcode) at state $state")
    end

    # Update guess and values with the solution
    vals  .= values(calc_eqm(state, name_values(policy, obj.Pfct), obj, mode=IterateMode()).vals)

    return SciMLBase.successful_retcode(sol)
end

function find_neighbors_with_diagonals(i, R, n=1)
    I = R[i]
    Ifirst, Ilast = first(R), last(R)
    I1 = oneunit(Ifirst)*n
    range = max(Ifirst, I-I1):min(Ilast, I+I1)
    return setdiff(range,[I])
end

function perform_iteration!(obj, damp, npasses)

    new_Vfct = deepcopy(coefficients(obj.Vfct))
    distances = [zeros(nof(obj.Vfct)) for _ in new_Vfct]
    success = falses(npt(obj.Vfct))

    states = collect(points(obj.ssg))
    Pfct_coeffs = coefficients(obj.Pfct)
    Vfct_coeffs = coefficients(obj.Vfct)

    R = CartesianIndices(values(dims(obj.ssg)))

    guess_selector = i -> @SVector[R[i]]

    pass = 1
    while true
        unsolved_state_indices = findall(!s for s in success)
        total_unsolved = length(unsolved_state_indices)
        new_Pfct = deepcopy(coefficients(obj.Pfct))

        Threads.@threads :dynamic for idx in eachindex(unsolved_state_indices)
            i = unsolved_state_indices[idx]
            state = states[i]
            guess_indices = guess_selector(i)
            #guess = Pfct_coeffs[i]
            values = Vfct_coeffs[i]
            new_values = new_Vfct[i]
            dist = distances[i]
            
            for guess_idx in guess_indices
                newguess = new_Pfct[guess_idx]
                oldguess = Pfct_coeffs[guess_idx]
                success[i] = solve_at_point!(newguess, new_values, oldguess, state, obj)
                if success[i] break end
            end

            dist .= new_values - values
            new_values .*= damp
            new_values .+= (1 - damp) * values

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

    # println(coefficients(obj.Vfct)[5,1,1])
    # println(coefficients(obj.Vfct)[5,5,1])
    # println(coefficients(obj.Vfct)[5,10,1])

    return success, distances
end

function main_execute!(obj; maxiter=100, tol=-5, damp=0.5, npasses=2, logging=ITERATIONFAILED)

    iter = 1
    original_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)  # To avoid thread oversubscription
    with_logger(ExecutionLogger(logging)) do

        while true 

            _, distances = perform_iteration!(obj, damp, npasses)

            distances = stack(distances)
            distances .= abs.(distances)
            idxmax = Tuple(argmax(distances))
            max_distance = distances[idxmax...]
            worst_var = names(obj.Vfct)[idxmax[1]]
            worst_point = map((statevar, idx) -> statevar[idx], 
                unigrids(obj.ssg), 
                NamedTuple{names(obj.ssg)}(idxmax[2:end]) )

            dist = log10(max_distance)

            @iteration "$iter: max distance = $dist at point $worst_point in variable $worst_var"

            if iter >= maxiter
                @complete "Maximum iterations reached without convergence."
                return
            elseif dist < tol
                @complete "Converged."
                return
            end
            iter += 1
        end

    end
    BLAS.set_num_threads(original_threads)  # Restore

end