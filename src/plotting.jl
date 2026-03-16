function _make_point(statevars, variablestates, fixedstates)
    return (; (k => (k in keys(fixedstates) ? fixedstates[k] : variablestates[k]) for k in statevars)...)
end

"""
    plot(fct::Interpolant, varname::Symbol, fixedstates::NamedTuple; transform=identity)

Plot one output variable of `fct` against the free (non-fixed) state dimensions.

`fixedstates` pins values for all but 1 or 2 free dimensions; the remainder are swept
across the full grid. `transform` is applied element-wise to values before plotting.

Produces a line plot (1 free dimension) or a filled contour plot (2 free dimensions).
Requires a Plots-compatible backend loaded by the caller.
"""
@recipe function f(fct::Interpolant, varname::Union{Symbol, Tuple{Symbol, <:Integer}}, fixedstates::NamedTuple)
    # pop! retrieves transform and removes it before Plots sees an unrecognized attribute
    # (note: `transform` is also a real Plots attribute, so this is mandatory)
    transform = pop!(plotattributes, :transform, identity)

    ssg = statespacegrid(fct)
    statevars = varnames(ssg)
    grid = (; (k => unigrids(ssg)[k] for k in statevars if k ∉ keys(fixedstates))...)
    Ndims = length(grid)
    fh_make_point = p -> _make_point(statevars, p, fixedstates)

    if varname isa Tuple
        # Each output argument is vector valued
        varname, vec_idx = varname
    else
        vec_idx = 1
    end

    title  --> string(varname)
    legend --> false

    if Ndims == 1
        xvec = values(grid)[1]
        fvec = transform.([fct(fh_make_point(p))[varname][vec_idx] for p in _points(grid)])
        xlabel --> string(keys(grid)[1])
        markershape --> :circle
        xvec, fvec

    elseif Ndims == 2
        x1vec = values(grid)[1]
        x2vec = values(grid)[2]
        fmat  = permutedims(transform.([fct(fh_make_point(p))[varname][vec_idx] for p in _points(grid)]))
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
