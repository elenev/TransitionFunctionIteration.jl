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

function normalize_fx_errors(::AbstractModel, fx_nt, addvars)
    return fx_nt
end

"""
    augment_simulation(panel, model) -> (static_nt, dynamic_nt)

Compute extra variables to append to the simulation panel after `simulate()` finishes.

Return a 2-tuple `(static_nt, dynamic_nt)` where:
- `static_nt`: `NamedTuple` of `N_runs × N_periods` arrays (same shape as `panel`), or `(;)`
- `dynamic_nt`: `NamedTuple` of `N_runs × (N_periods-1)` arrays (e.g., one-period returns),
  or `(;)` if no dynamic variables are needed

Extend this function for a specific model type to add derived variables. The default returns
empty NamedTuples (no augmentation).
"""
augment_simulation(::SimulationPanel, ::AbstractModel) = ((;), (;))

function _merge_augmentation(simseries::SimulationPanel{T,N,NAMES}, extra_static, extra_dynamic, dynamic_align) where {T,N,NAMES}
    # Merge static extra variables (same length as panel)
    if !isempty(extra_static)
        merged = StructArray(merge(components(simseries.data), extra_static))
        simseries = SimulationPanel(merged, NAMES)
    end

    # Merge dynamic extra variables (one period shorter)
    isempty(extra_dynamic) && return simseries

    if dynamic_align === :nan
        padded = map(arr -> hcat(fill(eltype(arr)(NaN), size(arr, 1), 1), arr), extra_dynamic)
    elseif dynamic_align === :drop
        simseries = simseries[:, 2:end]
        padded = extra_dynamic
    else
        throw(ArgumentError("dynamic_align must be :nan or :drop, got $dynamic_align"))
    end
    merged = StructArray(merge(components(simseries.data), padded))
    return SimulationPanel(merged, NAMES)
end

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

    return StructArray(trans), merge(state, exog, out.vals), fx_nt, transerr
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
                                sz::Pair=12=>10_000,
                                seed=100,
                                burn=500,
                                exog_start=ssidx(obj.exogenv),
                                endog_start=stv.endog,
                                shocks=sample_shocks(sz, seed, burn),
                                exog_states = simulate_markov_chain(obj.exogenv, exog_start, shocks)[1],
                                report_bounds=true,
                                dynamic_align=:nan,
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

    simseries = SimulationPanel(StructArray(@view simseries[:, t_range]), propertynames(endog_start))
    errs = SimulationErrorPanel(StructArray(@view errs[:, t_range]))
    
    if report_bounds
        check_bounds(simseries, obj)
    end

    extra_static, extra_dynamic = augment_simulation(simseries, obj.model)
    simseries = _merge_augmentation(simseries, extra_static, extra_dynamic, dynamic_align)

    return simseries, errs

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
    return hcat(stats, cor_stats)
end


"""
    describe(simseries; filters=(uncond=Returns(true),), ...) -> (; moments, ses)

Compute per-variable summary statistics for the simulation panel.

`filters` is a NamedTuple of callables, each taking one time-period element of the panel
(a NamedTuple of variable values) and returning a `Bool`. Statistics are computed on the
subsample of periods where the filter returns `true`. The default `Returns(true)` passes
every period (unconditional moments).

Returns `(; moments, ses)` where `moments` and `ses` are NamedTuples of `MomentArray`s,
keyed by filter names. For the default `filters=(uncond=Returns(true),)`, access the
unconditional moments via `out.moments.uncond`.
"""
function describe(simseries::SimulationPanel,
    stat_functions = default_stats;
    filters = (uncond = Returns(true),),
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

    var_names  = propertynames(simseries)
    stat_names = (keys(stat_functions)..., keys(cor_functions)...)
    filter_names = keys(filters)
    filter_vals  = Tuple(values(filters))
    n_filters    = length(filter_names)
    n_vars       = length(var_names)
    n_stats      = length(stat_names)

    function _filtered_stats(run, filt)
        mask = map(filt, run)
        sub  = run[mask]
        return isempty(sub) ? fill(NaN, n_vars, n_stats) : _describe(sub, stat_functions, cor_functions)
    end

    N_runs = size(simseries, 1)
    first_run, runs_iterator = Iterators.peel(eachrow(simseries))

    run_stats_f    = [_filtered_stats(first_run, f) for f in filter_vals]
    mean_stats_f   = [copy(r) for r in run_stats_f]
    new_mean_f     = [similar(r) for r in run_stats_f]
    stderr_stats_f = [zeros(eltype(r), size(r)...) for r in run_stats_f]

    for (n, run) in enumerate(runs_iterator)
        for j in 1:n_filters
            run_stats_f[j] .= _filtered_stats(run, filter_vals[j])
            new_mean_f[j]  .= mean_stats_f[j] .+ (run_stats_f[j] .- mean_stats_f[j]) ./ (n + 1)
            stderr_stats_f[j] .+= (run_stats_f[j] .- mean_stats_f[j]) .* (run_stats_f[j] .- new_mean_f[j])
            mean_stats_f[j]   .= new_mean_f[j]
        end
    end

    N_runs > 1 && foreach(se -> se .= sqrt.(se ./ (N_runs - 1)) ./ sqrt(N_runs), stderr_stats_f)

    moments = NamedTuple{filter_names}(Tuple(MomentArray(m, var_names, stat_names) for m in mean_stats_f))
    ses     = NamedTuple{filter_names}(Tuple(MomentArray(m, var_names, stat_names) for m in stderr_stats_f))

    return (; moments, ses)
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


linearize(simseries::SimulationPanel, stv::SteadyState; kwargs...) =
    linearize(simseries.data, stv; kwargs...)

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

## Moment comparison across economies

"""
    ReportedMoment(; var, sample=:uncond, stat=:mean)

Specification for a single entry in a `compare_economies` table.

- `var`: variable name (Symbol) in the `MomentArray`
- `sample`: filter name used in `describe()` (default `:uncond`)
- `stat`: statistic name (default `:mean`)
"""
@kwdef struct ReportedMoment
    var::Symbol
    sample::Symbol = :uncond
    stat::Symbol   = :mean
end

"""
    compare_economies(moment_defs, economies) -> MomentArray

Build a labeled table of moments across multiple economies.

- `moment_defs`: NamedTuple of `ReportedMoment`s (rows of the output table)
- `economies`: NamedTuple of `describe()` outputs (columns of the output table)

Returns a `MomentArray` with moment names as rows and economy names as columns.
Each cell is `economies[econ].moments[moment.sample][moment.var, moment.stat]`.
"""
function compare_economies(moment_defs::NamedTuple, economies::NamedTuple)
    moment_names = keys(moment_defs)
    econ_names   = keys(economies)

    data = [getmoment(economies[econ].moments[moment.sample], moment.var, moment.stat)
            for moment in values(moment_defs), econ in econ_names]

    return MomentArray(data, moment_names, econ_names)
end
