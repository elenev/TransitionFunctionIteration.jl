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
    simulate_markov_chain(exogenv, start, shocks) -> Matrix{Int}

Simulate exogenous state paths from `start` for `size(shocks, 2)` periods.

- `start`: integer (broadcast to all runs) or `Vector{Int}` (per-run)
- `shocks`: `N_runs × N_periods` matrix of uniform [0,1) draws
- Returns: `N_runs × (N_periods+1)` matrix of 1-indexed state integers
"""
function simulate_markov_chain(exogenv, start::Int, shocks)
    startvec = fill(start, size(shocks,1))
    return simulate_markov_chain(exogenv, startvec, shocks)
end

function simulate_markov_chain(exogenv::AbstractExogenv, start::AbstractVector{T} where T<:Int, shocks::AbstractMatrix)
    N_runs, N_periods = size(shocks)
    states = zeros(Int64, N_runs, N_periods+1)
    cumprobs = zeros(N_runs, exnpt(exogenv))
    states[:,1] .= start
    for t in axes(shocks,2)
        cumprobs .= cumsum(exogenv.mtrans[states[:,t],:], dims=2)
        states[:,t+1] .= dropdims(sum(shocks[:,t] .> cumprobs, dims=2), dims=2) .+ 1
    end
    return states
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
function simulate(obj::GenericModel, stv=nothing;
                                sz::Pair=12=>100,
                                seed=100,
                                burn=0,
                                exog_start=ssidx(obj.exogenv),
                                endog_start=stv.state,
                                shocks=sample_shocks(sz, seed, burn)
                    )

    exog_states = simulate_markov_chain(obj.exogenv, exog_start, shocks)
    start = make_point(exog_start, endog_start)

    function calc_result!(fx, state)
        out = calc_eqm(state, obj.Pfct(state), obj; mode=SimulateMode())
        exog = points(obj.exogenv)[state.exst]
        fx .= out.fx
        return out.trans, merge(state, exog, out.addvars)
    end

    N_runs, N_periods = size(exog_states)
    errtype = eltype(coefficients(obj.Pfct))
    fxerrs = [errtype(undef) for _ in 1:N_runs, _ in 1:N_periods]
    tmp = calc_result!(fxerrs[1,1], start)  # to compile
    simseries = Matrix{typeof(tmp[2])}(undef, N_runs, N_periods)

    original_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)  # To avoid thread oversubscription
    Threads.@threads :dynamic for r in axes(exog_states,1)
        trans = endog_start
        for t in axes(exog_states,2)
            state = make_point(exog_states[r,t], trans)
            trans, simseries[r,t] = calc_result!(fxerrs[r,t], state)
        end
    end
    BLAS.set_num_threads(original_threads)  # Restore

    return StructArray(simseries[:, burn+1:end]), fxerrs[:, burn+1:end]

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
