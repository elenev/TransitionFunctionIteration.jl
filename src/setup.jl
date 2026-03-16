# Methods for configuring the initial guess

function f_steadystate(x, p)
    policy = _pack_rows(typeof(p.policy_guess), view(x,1:p.nof))
    endog = _pack_rows(typeof(p.transition_guess), view(x, (p.nof+1):length(x)))
    state = merge(p.exog, endog)
    next_nt = merge(p.exog, endog, vals_from_policy(state, policy, p.model))
    np = compute_next_period_values(next_nt, p.model)
    out = calc_eqm(state, policy, p.model, np, IterateMode())
    transerr = map(-, out.trans, endog)
    return vcat(out.fx, SVector(values(transerr))), out
end

# Fall back if no steadystate() method is defined
# Uses calc_eqm with an additional set of equation ensuring the state is a fixed point of the transition function
function steadystate_exogenv(exogenv)
    if !(exogenv isa AbstractExogenv) || exnpt(exogenv) > 1
        if exogenv isa AbstractExogenv
            exog_mean = mean(exogenv)
        else
            exog_mean = exogenv
        end

        pts = StaticStructArray(SVector(exog_mean))
        exogenv = Exogenv(pts = pts, mtrans = [1.;;], ssidx=1)
    end
    return exogenv
end

function steadystate(mdl::AbstractModel, policy_guess, transition_guess, exogenv)
    exogenv = steadystate_exogenv(exogenv)

    payload = (model = mdl, policy_guess = policy_guess, 
        nof = length(propertynames(policy_guess)),
        transition_guess = transition_guess, 
        exog = merge( (exst = 1,), first(points(exogenv))), 
        prnext = conditional_distro(exogenv,1))

    x0 = vcat(SVector(values(policy_guess)), SVector(values(transition_guess)))
    objfun(x, p) = f_steadystate(x, p)[1]
    problem = NonlinearProblem{false}(objfun, x0, payload)
    alg = SimpleTrustRegion()

    sol = solve(problem, alg)
    policy_sol = _pack_rows(typeof(policy_guess), sol.u[1:payload.nof])
    transition_sol = _pack_rows(typeof(transition_guess), sol.u[payload.nof+1:end])
    _, out = f_steadystate(sol.u, payload)
    return SteadyState(mdl, policy_sol, transition_sol, first(points(exogenv)), out.vals), SciMLBase.successful_retcode(sol)
    
end


# Main setup 
function main_setup(model::AbstractModel, initobj::Union{Nothing,AbstractModelSolution}=nothing; 
    grids=nothing)

    # Calculate steady state
    stv = steadystate(model)

    # Define endogenous grids
    default_grids = define_endog_grids(model, stv)
    grids = isnothing(grids) ? default_grids : merge(default_grids, grids)

    # Create exogenous environment
    exogenv = Exogenv(model)

    # Create state space grid
    ssg = StateSpaceGrid(exogenv, grids)

    # Create interpolants
    if isnothing(initobj)
        Pfct = Interpolant(ssg, stv.guess)
        Vfct = Interpolant(ssg, stv.values)

        exog_ones = ones(exnpt(exogenv))
        endog = StructArray(
            map(p -> p * exog_ones, Base.structdiff(pt, (exst = nothing,)) ) 
            for pt in points(ssg)
        )
        #@debug "Checking endog" typeof(endog) endog
        Tfct = Interpolant(ssg, endog)

        Pfct = Interpolant(ssg, stv.guess)
        Vfct = Interpolant(ssg, stv.values)
    else
        Pfct = Interpolant(ssg, initobj.Pfct(points(ssg)))
        Vfct = Interpolant(ssg, initobj.Vfct(points(ssg)))
        Tfct = Interpolant(ssg, stv.transitions)
    end

    obj = ModelSolution(; model, exogenv, ssg, Pfct, Vfct, Tfct)
    return obj, stv
end