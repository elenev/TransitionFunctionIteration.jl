# Impulse Responses

"""
    IRFResult

Container for impulse-response function results returned by `simulate_irf`.

- `means`: `NamedTuple` of `StructArray`s (one per shock) of per-period mean paths
- `stdevs`: same structure, per-period standard deviations across simulation draws
- `other_stats`: additional user-requested statistics (see `other_stats` kwarg of `simulate_irf`)
- `paths`: per-draw raw simulation paths (non-empty only when `path_vars` is set)

Access individual shock paths with `irf.means.shock_name.varname`.
"""
struct IRFResult{M,D,S,P}
    means::M
    stdevs::D
    other_stats::S
    paths::P
end

shock_names(irf::IRFResult) = propertynames(irf.means)

function IRFResult(lines, paths)
    outer_keys = propertynames(lines)
    inner_keys = propertynames(first(lines))
    renested_lines = NamedTuple{inner_keys}(
        NamedTuple{outer_keys}(
            lines[outer_key][inner_key] for outer_key in outer_keys
        ) for inner_key in inner_keys
    )
    other_keys = Base.structdiff(renested_lines, (; mean = nothing, std = nothing))
    other_lines = NamedTuple{keys(other_keys)}(renested_lines)
    return IRFResult(renested_lines.mean, renested_lines.std, other_lines, paths)
end

function IRFResult(result)
    lines = map(s -> s[1], result)
    paths = map(s -> s[2], result)
    return IRFResult(lines, paths)
end


"""
    simulate_irf(obj, endogstart; shocks=(;), T=20, N=5_000, seed=100, ...)

Simulate impulse response functions (IRFs) from a solved model `obj`.

- `endogstart`: initial endogenous state, either a `NamedTuple` of values or a `SimulationPanel`
  (in which case the ergodic mean is used)
- `shocks`: `NamedTuple` of named shock specifications `(shock_name = (exst_index, exog_override), ...)`
- `T`: IRF horizon (periods)
- `N`: number of Monte Carlo draws for expectation/fan charts
- `seed`: random seed for reproducibility

Returns an `IRFResult` with mean and standard-deviation paths for each shock and for the
unconditional (`:uncond`) baseline.
"""
function simulate_irf(obj::ModelSolution, endogstart;
    exog_start = obj.exogenv.ssidx,
    shocks = (;), 
    T=20, 
    N=5_000,
    seed = 100,
    path_vars = (),
    other_stats = (;),)

    stats = merge( (; mean, std), other_stats)

    if endogstart isa SimulationPanel
        endogstart = (; (k => mean(getproperty(endogstart, k)) 
            for k in endog_states(endogstart))...)
    end

    if isempty(shocks)
        other_exst = setdiff(1:exnpt(obj.exogenv), exog_start)
        shocks = (; (Symbol("shock_$i") => i for i in other_exst)...)
        shocks = merge( (uncond = Int64[],), shocks)
    end

    # Generate random numbers
    sz = N => T
    rands = sample_shocks(sz, seed, 0)

    # Simulate IRF for each shock
    _simulate_irf_shock(shockvec) = simulate_irf_shock(obj, exog_start, 
        endogstart, shockvec, rands, stats, path_vars)

    result = map(_simulate_irf_shock, shocks)
    return IRFResult(result)
end

function simulate_irf_shock(obj::ModelSolution, 
        start_exog::Integer, 
        start_endog::NamedTuple,
        shock::Union{Integer,AbstractVector{<:Integer}}, 
        rands::AbstractMatrix{<:AbstractFloat}, 
        stat_fcns::NamedTuple, 
        path_vars::NTuple{NP, Symbol} where NP)

    # Simulate exogenus states
    shock = shock isa Integer ? [shock] : shock
    T_fixed = length(shock)
    T = size(rands, 2) + 1
    N = size(rands, 1)
    last_fixed_exog = isempty(shock) ? start_exog : shock[end]
    exog, _ = simulate_markov_chain(obj.exogenv, last_fixed_exog, rands[:, T_fixed+1:end])
    #@debug "exog" exog
    if !isempty(shock)
        exog = hcat(fill(start_exog,N,1), repeat(shock', N, 1), exog[:,2:end])
    else
        exog = hcat(fill(start_exog,N,1), exog[:,2:end])
    end

    # Simulate model
    sim, _ = simulate(obj, burn=0, exog_start=start_exog, 
                endog_start=start_endog, exog_states=exog,report_bounds=false)

    # Summarize
    _f(fcn, var) = map(fcn, eachcol(var))
    _summarize_stat(fcn)= StructArray(map(s -> _f(fcn, s), components(sim.data)))
    lines = map(fcn -> _summarize_stat(fcn), stat_fcns)

    if !isempty(path_vars)
        paths = StructArray(
            NamedTuple{path_vars}(components(sim.data))
        )
    else
        paths = (;)
    end
    
    return lines, paths
end