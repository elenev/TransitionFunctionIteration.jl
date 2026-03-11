function _make_point(statevars, variablestates, fixedstates)
    return (; (k => (k in keys(fixedstates) ? fixedstates[k] : variablestates[k]) for k in statevars)...)
end

@recipe function f(fct::Interpolant, varname::Symbol, fixedstates::NamedTuple)
    # pop! retrieves transform and removes it before Plots sees an unrecognized attribute
    # (note: `transform` is also a real Plots attribute, so this is mandatory)
    transform = pop!(plotattributes, :transform, identity)

    ssg = statespacegrid(fct)
    statevars = names(ssg)
    grid = (; (k => unigrids(ssg)[k] for k in statevars if k ∉ keys(fixedstates))...)
    Ndims = length(grid)
    fh_make_point = p -> _make_point(statevars, p, fixedstates)

    title  --> string(varname)
    legend --> false

    if Ndims == 1
        xvec = values(grid)[1]
        fvec = transform.([fct(fh_make_point(p))[varname] for p in _points(grid)])
        xlabel --> string(keys(grid)[1])
        marker --> :dot
        xvec, fvec

    elseif Ndims == 2
        x1vec = values(grid)[1]
        x2vec = values(grid)[2]
        fmat  = transform.([fct(fh_make_point(p))[varname] for p in _points(grid)])
        seriestype := :contourf
        levels     --> 20
        colorbar   --> true
        xlabel     --> string(keys(grid)[1])
        ylabel     --> string(keys(grid)[2])
        x1vec, x2vec, fmat

    else
        error("Plotting only implemented in 1 and 2 dimensions.")
    end
end
