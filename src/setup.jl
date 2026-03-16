# Methods for configuring the initial guess

# Fall back if no steadystate() method is defined
# Uses calc_eqm with an additional set of equation ensuring the state is a fixed point of the transition function
function steadystate(mdl::AbstractModel)
    
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