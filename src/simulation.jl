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


_fx_to_nt(fx::NamedTuple) = fx
function _fx_to_nt(fx::AbstractVector)
    names = ntuple(i -> Symbol("fx$i"), length(fx))
    return NamedTuple{names}(Tuple(fx))
end

# Next period's exogenous state index; falls back to ssidx at the end of the path.
get_next_exst(t, ssidx, exog_run) = t < length(exog_run) ? exog_run[t+1] : ssidx

# Merge scalar fx errors with transition errors selected for the realized next exog state.
merge_error(fx_nt, transerr, next_exst) = merge(fx_nt, map(v -> v[next_exst], transerr))

# Compute simulated results for one period
function calc_eqm_from_state(state, obj, mode=SimulateMode())
    #@debug "Calculating equilibrium for state: $state"
   # Get exogenous state info
    exog = points(obj.exogenv)[state.exst]
    prnext = conditional_distro(obj.exogenv, state.exst)
    fullstate = merge(state, exog)

    # Compute next period's transitions and expectations
    trans = obj.Tfct(state)
    next = forecast_next(trans, obj)

    np = compute_expectations(next, prnext, obj.model, mode)
    
    # Compute equilibrium residuals, next period's transitions and simulated variables
    policies = map( to_constrained, obj.Pfct(state), obj.Pfct_bounds)
    #@debug "Simulated policys at state $state: " policies
    return calc_eqm(fullstate, policies, obj.model, np, mode), trans
end

function calc_result(state, obj)
    # Get exogenous state info
    exog = points(obj.exogenv)[state.exst]

    # Calculate equilibrium
    out, trans = calc_eqm_from_state(state, obj)

    # Expand transitions and compute errors
    out_trans = expand_transitions(out.trans, obj.exogenv)
    #out_trans = obj.Tfct(state)
    fx_nt = _fx_to_nt(out.fx)
    transerr = map(-, out_trans, trans)   # NamedTuple of SVectors, one per exog state

    return StructArray(trans), merge(state, exog, out.addvars), fx_nt, transerr
end

"""
    simulate(obj, stv=nothing; sz=12=>100, seed=100, burn=0,
             exog_start=ssidx(obj.exogenv), endog_start=stv.endog,
             shocks=sample_shocks(sz, seed, burn))

Simulate the solved model along stochastic paths.

Returns `(simseries, errs)`:
- `simseries`: `N_runs × N_periods` StructArray of state and output variables
- `errs`: matching StructArray of equation and transition errors, named after the
  corresponding equations and endogenous state variables

`sz = N_runs => N_periods`. Burn-in periods are simulated but excluded from output.
"""
function simulate(obj::ModelSolution, stv=nothing;
                                sz::Pair=12=>1000,
                                seed=100,
                                burn=500,
                                exog_start=ssidx(obj.exogenv),
                                endog_start=stv.endog,
                                shocks=sample_shocks(sz, seed, burn),
                                exog_states = simulate_markov_chain(obj.exogenv, exog_start, shocks)[1]
                    )
    start = make_point(exog_start, endog_start)

    # Add another period so that transition errors can be computed
    exog_states = hcat(exog_states, fill(ssidx(obj.exogenv), size(exog_states, 1)))

    # Simulation size
    N_runs, N_periods = size(exog_states)

    # Infer element types from a trial run (also triggers compilation)
    _, _tmp_simvar, _tmp_fxerr, _tmp_transerr = calc_result(start, obj)
    _tmp_errs = merge_error(_tmp_fxerr, _tmp_transerr, exog_start)
    simseries = Matrix{typeof(_tmp_simvar)}(undef, N_runs, N_periods)
    errs      = Matrix{typeof(_tmp_errs)}(undef,   N_runs, N_periods)

    # Expand endog_start
    endog_start = StructArray(expand_transitions(endog_start, obj.exogenv))

    original_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)  # To avoid thread oversubscription
    #Threads.@threads :dynamic
    for r in axes(exog_states, 1)
        exog_run = view(exog_states, r, :)
        trans = endog_start
        for t in axes(exog_states, 2)
            exog_idx = exog_states[r, t]
            state = make_point(exog_idx, trans[exog_idx])
            trans, simseries[r,t], fx_nt, transerr = calc_result(state, obj)
            errs[r,t] = merge_error(fx_nt, transerr, get_next_exst(t, ssidx(obj.exogenv), exog_run))
        end
    end
    BLAS.set_num_threads(original_threads)  # Restore

    t_range = (burn+1):(N_periods-1)
    return SimulationPanel(StructArray(simseries[:, t_range]), propertynames(endog_start)), 
            SimulationErrorPanel(StructArray(errs[:, t_range]))

end

"""
    calc_error_stats(errs; percentiles=[50,75,90,95,99,100]) -> MomentArray

Summarize absolute equation and transition errors from a simulation.

Returns a `MomentArray` with error names as rows and stat names (`:mean`, `:p50`, …)
as columns. Access entries with `getmoment(out, :euler, :mean)` etc.
"""
function calc_error_stats(errs::SimulationErrorPanel; percentiles=[50, 75, 90, 95, 99, 100])

    ntprct = (; (Symbol("p$p") => x -> quantile(x, p/100) for p in percentiles)...)
    statfcns = merge((mean = mean,), ntprct)

    row_names = propertynames(errs)
    col_names = keys(statfcns)

    errvecs = map(c -> vec(abs.(c)), components(errs.data))
    data = [f(v) for v in errvecs, f in values(statfcns)]

    return MomentArray(data, row_names, col_names)

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


const default_percentiles = (5, 25, 50, 75, 95);

const default_stats = merge( (; mean, std, minimum), 
    (; (Symbol("p$k") => x -> quantile(x,k/100) for k in default_percentiles)...),
    (; max = maximum, ac = x -> autocor(x,[1])[1]) )



function _describe(simvector, stat_functions, cor_functions)
    varnames = propertynames(simvector)
    F = promote_type(fieldtypes(eltype(simvector))...)
    stats = [F(f(getproperty(simvector, var))) for var in varnames, f in values(stat_functions)]
    cor_stats = [F(f(getproperty(simvector, var), simvector)) for var in varnames, f in values(cor_functions)]
    #@debug "cor_stats: " cor_stats
    return hcat(stats, cor_stats)
end


function describe(simseries::SimulationPanel, 
    stat_functions = default_stats;
    additional_stats = (;),
    correlation_vars = collect(endog_states(simseries))
    )
    
    if correlation_vars isa Vector{Symbol}
        correlation_vars = [var => 0 for var in correlation_vars]
    end

    stat_functions = merge(stat_functions, additional_stats)
    
    function _cor(v, s, n) 
        if n != 0
            v = view(v, n+1:length(v))
            s = view(s, 1:length(s)-n)
        end
        return cor(v, s)
    end
    
    _corname(s,n) = Symbol("cor_$(s)_lag$(n)")
    cor_functions = (; ( 
        _corname(s,n) => (v,sim) -> _cor(v, getproperty(sim, s), n) 
        for (s,n) in correlation_vars)...)

    first_run, runs_iterator = Iterators.peel(eachrow(simseries))
    run_stats = _describe(first_run, stat_functions, cor_functions)
    
    mean_stats = deepcopy(run_stats)
    new_mean_stats = similar(run_stats)
    stderr_stats = similar(run_stats) 

    fill!(stderr_stats, 0.0)

    for (n, simvector) in enumerate(runs_iterator)
        run_stats .= _describe(simvector, stat_functions, cor_functions)
        new_mean_stats .= mean_stats + (run_stats - mean_stats) / n
        stderr_stats .+= (run_stats - mean_stats) .* (run_stats .- new_mean_stats)
    end

    N = size(simseries, 1)
    stderr_stats .= sqrt.( stderr_stats / (N - 1) ) / sqrt(N)

    stat_names = (keys(stat_functions)..., keys(cor_functions)...)

    return MomentArray(mean_stats,propertynames(simseries), stat_names),
              MomentArray(stderr_stats, propertynames(simseries), stat_names)
end


## Linear approximation from solved model 

@inline _add_primes(names) = map(n -> Symbol("$(n)′"), names)

function _rename_exogenous_points(exogenv, stvexog, demean)
    nextnames = _add_primes(varnames(exogenv))
    expoints = map(demean, components(points(exogenv)), stvexog)
    return StructArray(NamedTuple{nextnames}(values(expoints)))
end

function _regress(states, values, valnames; constant = true)
    states = states isa StructArray ? states : StructArray(states)
    if constant
        states = StructArray(merge((constant = ones(size(states, 1)),), components(states)))
    end
    statesmat = hcat(components(states)...)  # npt x nvar
    #return statesmat
    coefs = statesmat \ values  # nof x nvar
    coefs = SMatrix{size(coefs)...}(coefs)  # Convert to static for easier handling
    return MomentArray(coefs', valnames, keys(components(states)))
end

function linearize(obj::ModelSolution, stv::SteadyState; loglin=true)

    _demean_logs(x,xss) = @. log(x) - log(xss)
    _demean_levels(x,xss) = @. x - xss
    demean = loglin ? _demean_logs : _demean_levels

    function assemble_state(pt)
        exst = pt.exst
        exog = points(obj.exogenv)[exst]
        endog = NamedTuple{filter(n -> n ≠ :exst, keys(pt))}(pt)
        
        exog = map(demean, exog, stv.exog)
        endog = map(demean, endog, stv.endog)
        return merge( exog, endog)
    end

    # Get state variables
    states = vec( StructArray(map(assemble_state, points(obj.ssg))) )

    # Get v functions
    vfct = vec( obj.Vfct(points(obj.ssg)) )  # npt x nof
    vfct = map(demean, components(vfct), stv.values)  # Demean v functions
    vfct = stack(vfct) # npt x nof
    vcoefs = _regress(states, vfct, varnames(obj.Vfct))
    #return vcoefs

    # add future exogenous states
    nextexog = _rename_exogenous_points(obj.exogenv, stv.exog, demean)
    merged_states = StructArray(merge(o,n) for o in states, n in nextexog)
    merged_states = vec(merged_states)

    # Get transition functions
    tfct = vec( obj.Tfct(points(obj.ssg)) )  # npt x nvar
    tfct = map((x,xss) -> demean(stack(x, dims=1), xss), components(tfct), stv.endog)  # Demean v functions
    tfct = stack(map(vec, tfct))
    # tfct = stack(vec(coefficients(obj.Tfct)),dims=1) # npt x nsten x exnpt
    tcoefs = _regress(merged_states, tfct, _add_primes(varnames(obj.Tfct)))


    return vcoefs, tcoefs
end


function linearize(simseries::StructArray{T,2} where T, stv::SteadyState; loglin=true)
    _demean_logs(x,xss) = @. log(x) - log(xss)
    _demean_levels(x,xss) = @. x - xss
    demean = loglin ? _demean_logs : _demean_levels

    statevars = merge(stv.exog, stv.endog)
    stvvars = merge(statevars, stv.values)
    vars = map(demean, NamedTuple{keys(stvvars)}(components(simseries)), stvvars)
    states =  NamedTuple{keys(statevars)}(vars)
    vals = NamedTuple{keys(stv.values)}(vars)
    vals = map(vec, vals)

    vcoefs = _regress(map(vec, states), stack(vals), keys(stv.values))

    next_exog_names = _add_primes(keys(stv.exog))
    next_endog_names = _add_primes(keys(stv.endog))
    nextnames = (next_exog_names..., next_endog_names...)

    nextstates = map(x -> x[:,2:end], states)
    nextstates = NamedTuple{nextnames}(values(nextstates))
    prevstates = map(x -> x[:,1:end-1], states)
    merged_states = merge(prevstates, nextstates)

    x_names = (keys(statevars)..., next_exog_names...)
    X = map(vec, NamedTuple{x_names}(merged_states))
    Y = map(vec, NamedTuple{next_endog_names}(merged_states)) 

    tcoefs = _regress(X, stack(Y), next_endog_names)

    return vcoefs, tcoefs
end