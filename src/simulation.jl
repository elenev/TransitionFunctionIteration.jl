"""
    sample_shocks(sz::Pair, seed, burn) -> Matrix{Float64}

Draw a `first(sz) × (last(sz) + burn)` matrix of uniform [0,1) random numbers
using `MersenneTwister(seed)`. `sz` should be `N_runs => N_periods`.
"""
function sample_shocks(sz, seed, burn)
    rng = MersenneTwister(seed)
    return rand(rng, first(sz), last(sz)+burn)
end

"""
    simulate_markov_chain(exogenv, start, shocks) -> Matrix{Integer}

Simulate exogenous state paths from `start` for `size(shocks, 2)` periods.

- `start`: integer (broadcast to all runs) or `Vector{Integer}` (per-run)
- `shocks`: `N_runs × N_periods` matrix of uniform [0,1) draws
- Returns: `N_runs × (N_periods+1)` matrix of 1-indexed state integers
"""
function simulate_markov_chain(exogenv::AbstractExogenv, start::Integer, shocks::AbstractMatrix)
    startvec = fill(start, size(shocks,1))
    return simulate_markov_chain(exogenv, startvec, shocks)
end

function simulate_markov_chain(exogenv::AbstractExogenv, sz, start = ssidx(exogenv), seed=100, burn=0)
    shocks = sample_shocks(sz, seed, burn)
    return simulate_markov_chain(exogenv, start, shocks)
end

function simulate_markov_chain(exogenv::AbstractExogenv, start::AbstractVector{T} where T<:Integer, shocks::AbstractMatrix)
    N_runs, N_periods = size(shocks)
    states = zeros(Int64, N_runs, N_periods+1)
    cumprobs = zeros(N_runs, exnpt(exogenv))
    states[:,1] .= start
    for t in axes(shocks,2)
        cumprobs .= cumsum(exogenv.mtrans[states[:,t],:], dims=2)
        states[:,t+1] .= dropdims(sum(shocks[:,t] .> cumprobs, dims=2), dims=2) .+ 1
    end
    return states, points(exogenv, states)
end

# Returns an MVector type corresponding to the row of an MMatrix
_row_mvec_type(::Type{<:StaticMatrix{R, C, T}}) where {R, C, T} = MVector{C, T}

# Compute simulated results for one period
function calc_result!(fx, transerr, state, obj)
    # Get exogenous state info
    exog = points(obj.exogenv)[state.exst]
    prnext = conditional_distro(obj.exogenv, state.exst)
    fullstate = merge(state, exog)

    # Compute next period's transitions and expectations
    trans = obj.Tfct(state)
    next = forecast_next(trans, obj)
    np = compute_expectations(next, prnext, obj, SimulateMode())
    
    # Compute equilibrium residuals, next period's transitions and simulated variables
    out = calc_eqm(fullstate, obj.Pfct(state), obj.model, np, SimulateMode())

    # Expand transitions
    out_trans = expand_transitions(out.trans, obj.exogenv)
    
    # Fill equilibrium residuals
    fx .= out.fx

    # Compute transition errors
    out_transerr = map(-, components(out_trans), trans)
    fill_transerr!(this_transerr, this_out_transerr) = (this_transerr .= this_out_transerr)
    foreach(fill_transerr!, transerr, out_transerr)

    # Return transitions and simulated variables for this period
    return out_trans, merge(state, exog, out.addvars)
end

# Select transition errors
function select_errors(transerrs, exog_states)
    N, M = size(transerrs)

    # Preallocate
    ScalarTuple = typeof(map(first, first(transerrs)))
    out = StructArray{ScalarTuple}(undef, N, M)

    for CI in CartesianIndices(transerrs)
        n, m = Tuple(CI) # Unpack the spatial coordinates
        
        idx = m < M ? exog_states[n, m+1] : 3
        
        # CI natively indexes into both arrays seamlessly
        out[CI] = map(v -> v[idx], transerrs[CI]) 
    end

    return out
end

"""
    simulate(obj, stv=nothing; sz=12=>100, seed=100, burn=0,
             exog_start=ssidx(obj.exogenv), endog_start=stv.state,
             shocks=sample_shocks(sz, seed, burn))

Simulate the solved model along stochastic paths.

Returns `(simseries, fxerrs)`:
- `simseries`: `N_runs × N_periods` StructArray of state and output variables
- `fxerrs`: matching matrix of equilibrium residual vectors

`sz = N_runs => N_periods`. Burn-in periods are simulated but excluded from output.
"""
function simulate(obj::ModelSolution, stv=nothing;
                                sz::Pair=12=>100,
                                seed=100,
                                burn=0,
                                exog_start=ssidx(obj.exogenv),
                                endog_start=stv.endog,
                                shocks=sample_shocks(sz, seed, burn)
                    )

    exog_states, _ = simulate_markov_chain(obj.exogenv, exog_start, shocks)
    start = make_point(exog_start, endog_start)

    # Simulation size
    N_runs, N_periods = size(exog_states)

    # Preallocate equation errors
    fxerrtype = eltype(coefficients(obj.Pfct))
    fxerrs = [fxerrtype(undef) for _ in 1:N_runs, _ in 1:N_periods]

    # Preallocate transiton errors
    transerrtype = (_row_mvec_type ∘ eltype ∘ coefficients)(obj.Tfct)
    empty_transerr() = NamedTuple{varnames(obj.Tfct)}(ntuple(_ -> transerrtype(undef), nof(obj.Tfct)))
    transerrs = [empty_transerr() for _ in 1:N_runs, _ in 1:N_periods]

    # Identify the size of simulation output to preallocate
    tmp = calc_result!(fxerrs[1,1], transerrs[1,1], start, obj)  # to compile
    simseries = Matrix{typeof(tmp[2])}(undef, N_runs, N_periods)

    # Expand endog_start
    endog_start = expand_transitions(endog_start, obj.exogenv)

    original_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)  # To avoid thread oversubscription
    #Threads.@threads :dynamic 
    for r in axes(exog_states,1)
        trans = endog_start # Start at initial endogenous state
        for t in axes(exog_states,2)
            exog_idx = exog_states[r,t]
            state = make_point(exog_idx, trans[exog_idx])
            trans, simseries[r,t] = calc_result!(fxerrs[r,t],transerrs[r,t], state, obj)
        end
    end
    BLAS.set_num_threads(original_threads)  # Restore

    return StructArray(simseries[:, burn+1:end-1]), fxerrs[:, burn+1:end-1], select_errors(transerrs, exog_states)

end

function abs!(X::AbstractArray)
    X .= abs.(X)
end

"""
    calc_error_stats(fxerrs; percentiles=[50,75,90,95,99,100]) -> StructArray

Summarize absolute equation errors from a simulation.

Returns a StructArray with one entry per equation, each containing `:mean`
and `:pXX` fields for the requested percentiles.
"""
function calc_error_stats(fxerrs; percentiles=[50, 75, 90, 95, 99, 100])

    Nof = length(first(fxerrs))
    allerrs = stack(fxerrs, dims=1)
    abs!(allerrs)

    ntprct = (; (Symbol("p$p") => x -> quantile(x, p/100) for p in percentiles)...)
    statfcns = merge( (mean = mean, ), ntprct)

    function calc_error_stats_for_dim(i)
        this_err = view(allerrs, :, i)
        return map(f -> f(this_err), statfcns)
    end

    return StructArray(calc_error_stats_for_dim(i) for i in 1:Nof)

end

"""
    check_bounds(simseries, obj)

Log whether simulated endogenous variables remain within their grid bounds.
Prints `@info` if within bounds or `@warn` if any variable exceeds the grid.
"""
function check_bounds(simseries, obj)
    endog_vars = varnames(obj.ssg)[2:end]
    for var in endog_vars
        grid = unigrids(obj.ssg)[var]
        mingrid, maxgrid = extrema(grid)
        minsim, maxsim = extrema(getproperty(simseries, var))
        all_good = (minsim >= mingrid) && (maxsim <= maxgrid)
        if all_good
            @info "Variable $var within bounds: [$mingrid, $maxgrid], simulated: [$minsim, $maxsim]"
        else
            @warn "Variable $var out of bounds! Grid: [$mingrid, $maxgrid], simulated: [$minsim, $maxsim]"
        end
    end
end
