"""
    main_setup(model; mode=AnalyticalSteadyStateMode(), grids=nothing, check_steadystate=true)
    main_setup(model, initobj; ...)

Assemble a `ModelSolution` for `model`, computing the steady state, building grids, and
constructing initial interpolants. Returns `(obj, stv)` where `obj::ModelSolution` is
passed to `main_execute!` and `stv::SteadyState` holds the steady-state description.

Pass `initobj` to warm-start from a previously solved solution (e.g., when re-solving after
a parameter change).
"""
function main_setup(model::AbstractModel, 
    initobj::Union{Nothing,AbstractModelSolution}=nothing; 
    mode::SteadyStateMode=AnalyticalSteadyStateMode(),
    grids=nothing,
    check_steadystate=true)

    # Create exogenous environment
    exogenv = Exogenv(model)

    # Calculate steady state
    stv, _ = steadystate(model, mode)
    #@debug "Steady state type: " typeof(stv)

    if !(mode isa FixedPointSteadyStateMode) && check_steadystate
        resid, _ = f_steadystate(stv.guess, stv.endog, (model = model, exog = stv.exog))
        if norm(resid) > 1e-8
            @warn "The provided steady state solver is inconsistent with the dynamic steady state."
        end
    end

    # Adjust steady state before using it as initial value 
    stv = initial_guess_from_stv(stv, model)

    # Define endogenous grids
    default_grids = define_endog_grids(model, stv)
    grids = isnothing(grids) ? default_grids : merge(default_grids, grids)

    # Create state space grid
    ssg = StateSpaceGrid(exogenv, grids)

    # Define bounds
    Pfct_bounds, Tfct_bounds = bounds(stv)
    Pfct_bounds = NamedTuple{keys(stv.guess)}(Pfct_bounds)
    Tfct_bounds = NamedTuple{keys(stv.endog)}(Tfct_bounds)

    # Create interpolants
    if isnothing(initobj)
        unconstrained_guess = map(to_unconstrained, stv.guess, Pfct_bounds)
        Pfct = Interpolant(ssg, unconstrained_guess)
        Vfct = Interpolant(ssg, stv.values)

        vectorize_pt(pt) = NamedTuple{keys(stv.endog)}(expand_transitions(pt, exogenv))
        endog = StructArray(vectorize_pt(pt) for pt in points(ssg))
        #@debug "Checking endog" typeof(endog) endog
        Tfct = Interpolant(ssg, endog)
    else
        Pfct = Interpolant(ssg, initobj.Pfct(points(ssg)))
        Vfct = Interpolant(ssg, initobj.Vfct(points(ssg)))
        Tfct = Interpolant(ssg, initobj.Tfct(points(ssg)))
    end

    obj = ModelSolution(; model, exogenv, ssg, Pfct, Vfct, Tfct, Pfct_bounds, Tfct_bounds)
    return obj, stv
end

# Uses user-implemented steadystate_system function to solve for the steady state values given a guess.
# Idea is that the steady state can be represented by smaller, easier-to-solve system than the dynamic system
function steadystate(model::AbstractModel, config::SystemSteadyStateMode)
    x0 = to_unconstrained.(SVector(values(config.ssguess)), values(config.bounds))

    function f_steadystate(x, model, config, mode)
        x = to_constrained.(x, values(config.bounds))
        g = _pack_rows(typeof(config.ssguess), x)
        fx, stv = steadystate_system(g, model, mode)
        return _fx_to_svec(fx), stv
    end

    objfun(x,p) = f_steadystate(x, model, config, SolveMode())[1]
    problem = NonlinearProblem{false}(objfun, x0, nothing)
    alg = SimpleTrustRegion()
    sol = solve(problem, alg)

    _, stv = f_steadystate(sol.u, model, config, IterateMode())

    return stv, SciMLBase.successful_retcode(sol), sol
end

# Default fallback for vals_from_policy, which can be overloaded by users for their specific models
# If Vfct variables include variables that are not in policy, the overload is required to tell 
# the fixed-point steady state solver how to construct the next period values needed.
function vals_from_policy(state, policy, model::AbstractModel, args...)
    return policy
end

# Default initial guess is the steady state 
function initial_guess_from_stv(stv::SteadyState, model::AbstractModel)
    return stv
end

@inline function _guessvec_to_nt(x, nt_bounds)
    xc = to_constrained.(x, values(nt_bounds))
    return _pack_rows(typeof(nt_bounds), xc)
end

@inline function f_steadystate(x::AbstractVector, p)
    policy = _guessvec_to_nt(view(x, 1:p.nof), p.pfct_bounds)
    endog = _guessvec_to_nt(view(x, p.nof+1:length(x)), p.tfct_bounds)
    return f_steadystate(policy, endog, p)
end

function f_steadystate(policy::NamedTuple, endog::NamedTuple, p)
    state = merge(p.exog, endog)
    next_nt = merge(p.exog, endog, vals_from_policy(state, policy, p.model))
    mode = IterateMode()
    np = compute_next_period_values(next_nt, p.model, mode)
    np = normalize_next_period_values(np, p.model, mode)
    out = calc_eqm(state, policy, p.model, np, mode)
    transerr = map(-, out.trans, endog)
    return vcat(_fx_to_svec(out.fx), SVector(values(transerr))), out
end

# Uses calc_eqm with an additional set of equation ensuring the state is a fixed point of the transition function
function steadystate(mdl::AbstractModel, mode::FixedPointSteadyStateMode)

    policy_guess = mode.policies
    transition_guess = mode.transitions
    exogenv = mode.exogenv

    # Get bounds and ensure correct order
    pfct_bounds, tfct_bounds = bounds(mdl, policy_guess, transition_guess)
    pfct_bounds = NamedTuple{keys(policy_guess)}(pfct_bounds)
    tfct_bounds = NamedTuple{keys(transition_guess)}(tfct_bounds)

    @debug "Bounds: Policy function bounds: $pfct_bounds; Transition function bounds: $tfct_bounds"

    payload = (model = mdl, pfct_bounds = pfct_bounds, 
        nof = length(propertynames(policy_guess)),
        tfct_bounds = tfct_bounds,
        exog = merge( (exst = 1,), first(points(exogenv))))

    x0p = to_unconstrained.(SVector(values(policy_guess)), values(pfct_bounds))
    x0t = to_unconstrained.(SVector(values(transition_guess)), values(tfct_bounds))
    x0 = vcat(x0p, x0t)

    objfun(x, p) = f_steadystate(x, p)[1]
    problem = NonlinearProblem{false}(objfun, x0, payload)
    alg = SimpleTrustRegion()

    sol = solve(problem, alg)

    policy_sol = _guessvec_to_nt(view(sol.u, 1:payload.nof), pfct_bounds)
    transition_sol = _guessvec_to_nt(view(sol.u, payload.nof+1:length(sol.u)), tfct_bounds)
    _, out = f_steadystate(sol.u, payload)

    return SteadyState(mdl, policy_sol, transition_sol, first(points(exogenv)), out.vals), SciMLBase.successful_retcode(sol), sol
    
end

# Default is unbounded, but users can define their own bounds by overloading this method for their SteadyState type.
function bounds(stv::SteadyState)
    return map( _ -> Unbounded(), stv.guess), map( _ -> Unbounded(), stv.endog)
end

function bounds(::AbstractModel, policies::NamedTuple, transitions::NamedTuple)
    return map( _ -> Unbounded(), policies), map( _ -> Unbounded(), transitions)
end


# Fallback function for updating the user model with new parameter values
function explicit_parameters(::T) where {T<:AbstractModel}
    return fieldnames(T)
end

function update_model(mdl::T, newparams::NamedTuple) where {T<:AbstractModel}
    oldparams = (; (k => getproperty(mdl, k) for k in explicit_parameters(mdl) )...)
    return T(; merge(oldparams, newparams)...)
end

function update_model(mdl::T; kwargs...) where {T<:AbstractModel}
    return update_model(mdl, (; kwargs...))
end