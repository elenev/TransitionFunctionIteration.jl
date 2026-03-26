"""
    ploterrs(errs::SimulationErrorPanel, sim::SimulationPanel; kwargs...)

Convenience wrapper: equivalent to `plot(errs, sim; kwargs...)`.

Plot simulation errors as a grid: rows = exogenous states, columns = error variables.

# Keyword arguments
- `errnames`:    tuple of error field names to plot (default: all fields of `errs`)
- `exogstates`:  exogenous state values to include as rows (default: all unique values)
- `endogstates`: 1 or 2 endogenous state symbols for the x (and y) axes (default: first two)

# Layout
- 1D (`length(endogstates)==1`): scatter of error vs endog state; shared y-axis limits.
- 2D (`length(endogstates)==2`): scatter colored by absolute error magnitude; shared color limits.
  An extra column is added to hold a single shared colorbar.

Titles appear only on the top row; x-axis labels and tick values only on the bottom row;
y-axis labels and tick values only on the leftmost column.
`legend` is a soft default (`false`) and can be overridden by the caller.
"""
ploterrs(errs::SimulationErrorPanel, sim::SimulationPanel; kw...) = plot(errs, sim; kw...)

@recipe function f(errs::SimulationErrorPanel, sim::SimulationPanel)
    errnames    = pop!(plotattributes, :errnames,    propertynames(errs))
    exogstates  = pop!(plotattributes, :exogstates,  sort(unique(sim.exst)))
    endogstates = pop!(plotattributes, :endogstates, endog_states(sim)[1:2])

    N      = length(exogstates)
    M      = length(errnames)
    n_endog = length(endogstates)
    n_endog in (1, 2) || error("ploterrs only implemented for 1 or 2 endog states")

    _get_exst = exogstates == :all ? (_ -> trues(length(sim.exst))) :
                                     (exst -> sim.exst .== exst)

    # Global error limits; for 2D, color by absolute error magnitude.
    # Use extrema(f, v) to find each panel's range in one pass without allocating.
    _efn = n_endog == 2 ? abs : identity
    err_lims = mapreduce(
        t -> extrema(_efn, getproperty(errs.data[_get_exst(t[1])], t[2])),
        (a, b) -> (min(a[1], b[1]), max(a[2], b[2])),
        Iterators.product(exogstates, errnames)
    )

    # For 2D, add an extra column to hold the shared colorbar
    ncols = n_endog == 2 ? M + 1 : M
    layout := (N, ncols)
    size   --> (ncols * 200, N * 200)

    # Main grid: row r = exogstates[r], col c = errnames[c]
    for r in 1:N
        mask = _get_exst(exogstates[r])
        for c in 1:M
            errname = errnames[c]
            @series begin
                subplot           := (r - 1) * ncols + c
                seriestype        := :scatter
                title             := r == 1 ? string(errname) : ""
                xlabel            := r == N ? string(endogstates[1]) : ""
                xformatter        := r == N ? :auto : _ -> ""
                yformatter        := c == 1 ? :auto : _ -> ""
                markerstrokewidth := 0
                markerstrokecolor := :transparent
                colorbar          := false
                legend            --> false

                if n_endog == 1
                    ylims  := err_lims
                    ylabel := ""
                    getproperty(sim.data[mask], endogstates[1]),
                    getproperty(errs.data[mask], errname)
                else
                    ylabel := c == 1 ? string(endogstates[2]) : ""
                    zcolor := abs.(getproperty(errs.data[mask], errname))
                    clims  := err_lims
                    color  --> :jet
                    getproperty(sim.data[mask], endogstates[1]),
                    getproperty(sim.data[mask], endogstates[2])
                end
            end
        end
    end

    # Shared colorbar column for 2D: one invisible scatter per row, all showing the
    # same colorbar. Together they form a continuous strip spanning the full height.
    # framestyle=:none + explicit xlims/ylims avoids GKS invalid-viewport warnings.
    if n_endog == 2
        for r in 1:N
            @series begin
                subplot           := (r - 1) * ncols + ncols
                seriestype        := :scatter
                ms                := 0
                markerstrokewidth := 0
                colorbar          := true
                clims             := err_lims
                zcolor            := [mean(err_lims)]
                color             --> :jet
                framestyle        := :none
                xlims             := (0, 1)
                ylims             := (0, 1)
                legend            := false
                [0.5], [0.5]
            end
        end
    end
end

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
