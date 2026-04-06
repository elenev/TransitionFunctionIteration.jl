"""
    sample_shocks(sz::Pair, seed, burn) -> Matrix{Float64}

Draw a `(first(sz) + burn) × last(sz)` matrix of uniform [0,1) random numbers
using `MersenneTwister(seed)`. `sz` should be `N_periods => N_runs`.
"""
function sample_shocks(sz, seed, burn)
    rng = MersenneTwister(seed)
    return rand(rng, first(sz)+burn, last(sz))
end

"""
    simulate_markov_chain(exogenv, start, shocks) -> Matrix{Integer}

Simulate exogenous state paths from `start` for `size(shocks, 1)` periods.

- `start`: integer (broadcast to all runs) or `Vector{Integer}` (per-run)
- `shocks`: `N_periods × N_runs` matrix of uniform [0,1) draws
- Returns: `(N_periods+1) × N_runs` matrix of 1-indexed state integers
"""
function simulate_markov_chain(exogenv::AbstractExogenv, start::Integer, shocks::AbstractMatrix)
    startvec = fill(start, size(shocks,2))
    return simulate_markov_chain(exogenv, startvec, shocks)
end

function simulate_markov_chain(exogenv::AbstractExogenv, sz, start = ssidx(exogenv), seed=100, burn=0)
    shocks = sample_shocks(sz, seed, burn)
    return simulate_markov_chain(exogenv, start, shocks)
end

function simulate_markov_chain(exogenv::AbstractExogenv, start::AbstractVector{T} where T<:Integer, shocks::AbstractMatrix)
    N_periods, N_runs = size(shocks)
    states = zeros(Int64, N_periods+1, N_runs)
    cumprobs = zeros(N_runs, exnpt(exogenv))
    states[1,:] .= start
    for t in axes(shocks,1)
        this_states = view(states,t,:)
        this_shocks = view(shocks,t,:)
        cumprobs .= cumsum(exogenv.mtrans[this_states,:], dims=2)
        states[t+1,:] .= dropdims(sum(this_shocks .> cumprobs, dims=2), dims=2) .+ 1
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

function normalize_fx_errors(::AbstractModel, fx_nt, simvars)
    return fx_nt
end

"""
    augment_simulation(run, model) -> (static_nt, dynamic_nt)

Compute extra variables for one simulation run after `simulate()` finishes.

`run` is a `Vector` of T NamedTuples (one period per element). Return a 2-tuple where:
- `static_nt`: `NamedTuple` of T-element vectors (aligned with `run`), or `(;)`
- `dynamic_nt`: `NamedTuple` of (T-1)-element vectors (e.g., one-period returns computed
  from consecutive pairs), or `(;)` if no dynamic variables are needed

For `dynamic_align=:drop` in `simulate()`, burn is incremented by 1 so that `run` covers one
extra period at the start; static results trim that first element and dynamic results align
directly with the output panel. For `dynamic_align=:nan`, static results are used as-is and
dynamic results receive a leading `NaN`.

The default returns empty NamedTuples (no augmentation).
"""
augment_simulation(::StructVector, ::AbstractModel) = nothing

function _update_simseries_with_augmented(simseries, obj, dynamic_align)
    N_periods = size(simseries, 1)

    # Infer output type. If augment_simulation returns nothing, skip augmentation and return original simseries.
    if augment_simulation(first(eachcol(simseries)), obj.model) === nothing
        return simseries
    end

    # Augment with extra variables and align dynamic variables as needed.
    _pad_aug_var!(v)  = length(v) <  N_periods ? pushfirst!(v, eltype(v)(NaN)) : v
    _trim_aug_var!(v) = length(v) == N_periods ? popfirst!(v)                  : v
    align_aug_var!(v) = dynamic_align === :drop ? _trim_aug_var!(v) : _pad_aug_var!(v)

    # Augment each run, align, and stack
    augseries = let model = obj.model
        stack(
            map(eachcol(simseries)) do runcol
                aug = augment_simulation(runcol, model)
                foreach(align_aug_var!, aug)
                return StructArray(aug)
            end,
            dims=2)
    end

    # Merge
    simseries = dynamic_align === :drop ? view(simseries, 2:N_periods, :) : simseries

    return merge_structarrays(simseries, augseries)
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

    return StructArray(trans), merge(state, exog, out.simvars), fx_nt, transerr
end

"""
    simulate(obj, stv=nothing; sz=10_000=>12, seed=100, burn=500,
             exog_start=ssidx(obj.exogenv), endog_start=stv.endog,
             shocks=sample_shocks(sz, seed, burn))

Simulate the solved model along stochastic paths.

Returns `(simseries, errs)`:
- `simseries`: `N_periods × N_runs` `SimulationPanel` of state and output variables
- `errs`: matching `SimulationErrorPanel` of equation and transition errors

`sz = N_periods => N_runs`. Burn-in periods are simulated but excluded from output.

`dynamic_align` controls how dynamic augmentation variables (length T-1) are aligned with the
output panel: `:drop` drops the first output period (burn incremented by 1 internally);
`:nan` prepends `NaN` to each dynamic series.
"""
function simulate(obj::ModelSolution, stv=nothing;
                                dynamic_align=:drop,
                                sz::Pair=(10_000+(dynamic_align == :drop ? 1 : 0))=>12,
                                seed=100,
                                burn=500,
                                exog_start=ssidx(obj.exogenv),
                                endog_start=stv.endog,
                                shocks=sample_shocks(sz, seed, burn),
                                exog_states = simulate_markov_chain(obj.exogenv, exog_start, shocks)[1],
                                report_bounds=true,
                    )
    dynamic_align ∈ (:drop, :nan) || throw(ArgumentError("dynamic_align must be :drop or :nan, got $dynamic_align"))

    start = make_point(exog_start, endog_start)

    # Add another period so that transition errors can be computed
    exog_states = vcat(exog_states, fill(ssidx(obj.exogenv), 1, size(exog_states, 2)))
    N_periods, N_runs = size(exog_states)

    # Infer element types from a trial run (also triggers compilation)
    _, _tmp_simvar, _tmp_fxerr, _tmp_transerr = calc_result(start, obj)
    _tmp_errs = merge_error(_tmp_fxerr, _tmp_transerr, exog_start)
    sim_raw = Matrix{typeof(_tmp_simvar)}(undef, N_periods, N_runs)
    err_raw = Matrix{typeof(_tmp_errs)}(undef,   N_periods, N_runs)

    # Expand endog_start
    endog_start = StructArray(expand_transitions(endog_start, obj.exogenv))

    # Simulate each run and period
    original_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    for r in axes(exog_states, 2)
        exog_run = view(exog_states, :, r)
        trans = endog_start
        for t in axes(exog_states, 1)
            exog_idx = exog_states[t, r]
            state = make_point(exog_idx, trans[exog_idx])
            trans, sim_raw[t,r], fx_nt, transerr = calc_result(state, obj)
            err_raw[t,r] = merge_error(fx_nt, transerr, get_next_exst(t, ssidx(obj.exogenv), exog_run))
        end
    end
    BLAS.set_num_threads(original_threads)

    # Remove burn-in and package into StructArrays
    t_range = (burn+1):(N_periods-1)
    simseries = StructArray(view(sim_raw, t_range, :))

    # Check bounds
    report_bounds && check_bounds(simseries, obj)

    # Augment with extra variables and align as needed
    simseries = _update_simseries_with_augmented(simseries, obj, dynamic_align)

    # Package into SimulationPanel
    simseries = SimulationPanel(simseries, propertynames(endog_start))

    # Errors
    err_range = (burn+(dynamic_align === :drop ? 1 : 0) + 1):(N_periods-1)
    errs      = SimulationErrorPanel(StructArray(@view err_raw[err_range, :]))

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


const default_percentiles = SVector(5, 25, 50, 75, 95);


ac1(x) = cor(view(x, 2:length(x)), view(x, 1:length(x)-1))


# A multivariate stat function (v, sim) -> scalar that accesses sim field S by name.
# S is a type parameter so getproperty(sim, S) is type-stable; f receives the concrete column.
struct MultivariateStatFun{S, F}
    f::F
end
MultivariateStatFun{S}(f::F) where {S, F} = MultivariateStatFun{S, F}(f)
(msf::MultivariateStatFun{S})(v, sim) where {S} = msf.f(v, getproperty(sim, S))

const default_stats = StatisticsSet(
    (:mean, :std, :minimum, :maximum, :ac, :p5, :p25, :p50, :p75, :p95),
    (mean, std, minimum, maximum, ac1, map(p -> (x -> quantile(x, p/100)), default_percentiles)...)
    )


@inline function _describe_variable!(outrow, v, simvector, unistat_functions, multistat_functions)
    N_uni = length(unistat_functions)
    N_multi = length(multistat_functions)
    N_total = N_uni + N_multi
    outrow[1:N_uni]         .= values(unistat_functions(v))
    outrow[N_uni+1:N_total] .= values(multistat_functions(v, simvector))
end


function _describe_run!(out, simvector, unistat_functions, multistat_functions)
    foreach(eachrow(out), components(simvector)) do outrow, v
       _describe_variable!(outrow, v, simvector, unistat_functions, multistat_functions)
    end
end

# @generated function _describe_run!(out, simvector::StructArray{T,N,C}, unistat_functions, multistat_functions) where {T,N,C}
#     calls = [:(
#         _describe_variable!(view(out, $i, :), getproperty(simvector, $(QuoteNode(name))), simvector, unistat_functions, multistat_functions)
#     ) for (i, name) in enumerate(fieldnames(C))]
#     return quote $(calls...); nothing end
# end

function _preallocate_stats(simvector, stat_functions, cor_functions)
    varnames = propertynames(simvector)
    F = promote_type(fieldtypes(eltype(simvector))...)
    N_stat = length(stat_functions)
    N_cor = length(cor_functions)
    N_total = N_stat + N_cor
    N_var = length(varnames)
    return Matrix{F}(undef, N_var, N_total)
end

function _describe_sample(simseries, stat_functions, cor_functions, mask)
    N_runs = size(simseries, 2)
    (first_run, first_mask), runs_iterator =
        Iterators.peel(zip(
            eachcol(simseries),
            eachcol(mask))
        )

    run_stats    = _preallocate_stats(view(first_run, first_mask), stat_functions, cor_functions)
    _describe_run!(run_stats, view(first_run, first_mask), stat_functions, cor_functions)
    mean_stats   = copy(run_stats)
    delta        = similar(run_stats)
    stderr_stats = zeros(eltype(run_stats), size(run_stats)...)

    for (n, (simvector, maskvector)) in enumerate(runs_iterator)
        filtered_simvector = simvector[maskvector] #view(simvector, maskvector)
        _describe_run!(run_stats, filtered_simvector, stat_functions, cor_functions)
        delta         .= run_stats .- mean_stats
        mean_stats   .+= delta ./ (n + 1)
        stderr_stats .+= delta .* (run_stats .- mean_stats)
    end
    # return run_stats
    N_runs > 1 && (stderr_stats .= sqrt.(stderr_stats ./ (N_runs - 1)) ./ sqrt(N_runs))
    #return mean_stats, stderr_stats 
    all_stats = merge(stat_functions, cor_functions)
    return MomentArray(mean_stats,   simseries, all_stats),
           MomentArray(stderr_stats, simseries, all_stats)

end

function _build_correlation_functions(correlation_vars)
    # If correlation lag isn't specified, default is contemporaneous
    if correlation_vars isa Vector{Symbol}
        correlation_vars = [var => 0 for var in correlation_vars]
    end

    # (Auto)correlation function; col = getproperty(sim, s) resolved by MultivariateStatFun
    function _cor(v, col, n)
        n == 0 && return cor(v, col)
        return cor(view(v, n+1:length(v)), view(col, 1:length(col)-n))
    end

    # Create a map of correlation functions
    _corname(s, n) = Symbol("cor_$(s)_lag$(n)")
    cor_names = Tuple(_corname(s, n) for (s, n) in correlation_vars)
    cor_fcns  = Tuple(MultivariateStatFun{s}((v, col) -> _cor(v, col, n)) for (s, n) in correlation_vars)
    return StatisticsSet(cor_names, cor_fcns)
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
    filters = (;),
    additional_stats = (;),
    correlation_vars = collect(endog_states(simseries))
    )

    simseries = simseries.data

    # Merge user-specified filters  with defaults
    default_filters = (uncond = Returns(true),)
    filters = merge(default_filters, filters)

    # Merge user-specified stats with defaults
    stat_functions = merge(stat_functions, StatisticsSet(additional_stats))

    # Build correlation functions
    cor_functions = _build_correlation_functions(correlation_vars)

    # Iterate over samples
    samples = map(filters) do filter_fcn
        mask = filter_fcn.(simseries)
        #Main.@infiltrate        
        return _describe_sample(simseries, stat_functions, cor_functions, mask)
    end

    #Re-nest 
    moments = map(s -> s[1], samples)
    stderrs = map(s -> s[2], samples)
    
    return (; moments, stderrs)
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
    return MomentArray(coefs', valnames, components(states))
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


"""
    linearize(obj::ModelSolution, stv::SteadyState; loglin=true) -> (MomentArray, MomentArray)
    linearize(simseries, stv::SteadyState; loglin=true)         -> (MomentArray, MomentArray)

Compute linear approximations to the value and transition functions around the steady state `stv`.

The first form linearizes the solution object `obj` at each grid point; the second form runs an
OLS regression over the simulation panel. Returns `(vcoefs, tcoefs)` where each is a `MomentArray`
with rows = variable names and columns = `(constant, state vars...)`.

Set `loglin=false` for level deviations instead of log deviations.
"""
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

    nextstates = map(x -> x[2:end,:], states)
    nextstates = NamedTuple{nextnames}(values(nextstates))
    prevstates = map(x -> x[1:end-1,:], states)
    merged_states = merge(prevstates, nextstates)

    x_names = (keys(statevars)..., next_exog_names...)
    X = map(vec, NamedTuple{x_names}(merged_states))
    Y = map(vec, NamedTuple{next_endog_names}(merged_states)) 

    tcoefs = _regress(X, stack(Y), next_endog_names)

    return vcoefs, tcoefs
end

## Moment comparison across economies


"""
    compare_economies(moment_defs, economies) -> MomentArray

Build a labeled table of moments across multiple economies.

- `moment_defs`: NamedTuple of `ReportedMoment`s (rows of the output table)
- `economies`: NamedTuple of `describe()` outputs (columns of the output table)

Returns a `MomentArray` with moment names as rows and economy names as columns.
Each cell is `economies[econ].moments[moment.sample][moment.var, moment.stat]`.
"""
function compare_economies(moment_defs::NamedTuple, economies::NamedTuple; base::Union{Integer,Symbol}=1)
    if base isa Integer
        base = keys(economies)[base]
    end

    # Promote any non-ReportedMoment entries in moment_defs to ReportedMoment with default specs
    moment_defs  = map(v -> v isa ReportedMoment ? v : ReportedMoment(v), moment_defs)
    moment_names = keys(moment_defs)
    econ_names   = keys(economies)

    data = [getmoment(economies, moment, econ, base)
            for moment in values(moment_defs), econ in econ_names]

    return MomentArray(data, moment_defs, economies)
end
