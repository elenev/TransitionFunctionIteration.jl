"""
    plot(errs::SimulationErrorPanel, sim:SimulationPanel; kwargs...)

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

    # Shared axis limits across all plotted exog states
    _axis_lims = var -> mapreduce(
        exst -> extrema(getproperty(sim.data[_get_exst(exst)], var)),
        (a, b) -> (min(a[1], b[1]), max(a[2], b[2])),
        exogstates
    )
    x_lims = _axis_lims(endogstates[1])
    y_lims = n_endog == 2 ? _axis_lims(endogstates[2]) : nothing

    # For 2D, add an extra column to hold the shared colorbar
    ncols = n_endog == 2 ? M + 1 : M

    if n_endog == 2
        gridlayout = [(label = Symbol("p$(r*c)"), blank = false) for r=1:N, c=1:M]
        cblayout = (label = :colorbar, width = 0.12, height = :auto)
        combolayout = Matrix{Any}(undef,1,2)
        combolayout[1,1] = gridlayout
        combolayout[1,2] = cblayout

        layout  := combolayout
    else
        layout := (N, ncols)
    end
    #layout := @layout [grid(N, ncols) tall_col{0.25w}]
    cbwidth = n_endog == 2 ? 80 : 0
    size   --> (M * 200 + cbwidth, N * 200)

    # Main grid: row r = exogstates[r], col c = errnames[c]
    for r in 1:N
        mask = _get_exst(exogstates[r])
        for c in 1:M
            errname = errnames[c]
            @series begin
                subplot           := (r - 1) * M + c
                seriestype        := :scatter
                title             := r == 1 ? string(errname) : ""
                xlabel            := r == N ? string(endogstates[1]) : ""
                xformatter        := r == N ? :auto : _ -> ""
                xrotation         := 45
                yformatter        := c == 1 ? :auto : _ -> ""
                markerstrokewidth := 0
                markerstrokecolor := :transparent
                colorbar          := false
                legend            --> false

                xlims := x_lims

                if n_endog == 1
                    ylims  := err_lims
                    ylabel := ""
                    getproperty(sim.data[mask], endogstates[1]),
                    getproperty(errs.data[mask], errname)
                else
                    ylims  := y_lims
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
        #for r in 1:N
            @series begin
                subplot           := N * M + 1
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
        #end
    end
end

# The recipe above implements something like this function. We can use it for prototyping
# since recipes are kinda hard to debug

# function plot_errs(errs::SimulationErrorPanel, sim::SimulationPanel;
#                     errnames=propertynames(errs),
#                     exogstates=sort(unique(sim.exst)),
#                     endogstates=endog_states(sim)[1:2])

#     N = length(exogstates)
#     M = length(errnames)
#     n_endog = length(endogstates)

#     if exogstates == :all
#         _get_exst = _ -> trues(length(sim.exst))
#     else
#         _get_exst = exst -> sim.exst .== exst
#     end

#     # Global error limits; for 2D, color by absolute error magnitude.
#     # Use extrema(f, v) to find each panel's range in one pass without allocating.
#     _efn = n_endog == 2 ? abs : identity
#     err_lims = mapreduce(
#         t -> extrema(_efn, getproperty(errs.data[_get_exst(t[1])], t[2])),
#         (a, b) -> (min(a[1], b[1]), max(a[2], b[2])),
#         Iterators.product(exogstates, errnames)
#     )

#     # Shared axis limits across all plotted exog states
#     _axis_lims = var -> mapreduce(
#         exst -> extrema(getproperty(sim.data[_get_exst(exst)], var)),
#         (a, b) -> (min(a[1], b[1]), max(a[2], b[2])),
#         exogstates
#     )
#     x_lims = _axis_lims(endogstates[1])
#     y_lims = n_endog == 2 ? _axis_lims(endogstates[2]) : nothing

#     function _plot_err_1D(mask, errname, statename; kw...)
#         scatter(getproperty(sim.data[mask], statename), getproperty(errs.data[mask], errname);
#             xlims=x_lims, ylims=err_lims, markerstrokewidth=0, markerstrokecolor=:transparent, kw...)
#     end

#     function _plot_err_2D(mask, errname, statename1, statename2; kw...)
#         abs_errs = abs.(getproperty(errs.data[mask], errname))
#         scatter(getproperty(sim.data[mask], statename1), getproperty(sim.data[mask], statename2);
#             xlims=x_lims, ylims=y_lims, zcolor=abs_errs, colormap=:jet, clims=err_lims, colorbar=false,
#             markerstrokewidth=0, markerstrokecolor=:transparent, kw...)
#     end

#     if n_endog == 1
#         _plot = (mask, errname; kw...) -> _plot_err_1D(mask, errname, endogstates[1]; kw...)
#     elseif n_endog == 2
#         _plot = (mask, errname; kw...) -> _plot_err_2D(mask, errname, endogstates[1], endogstates[2]; kw...)
#     else
#         error("Plotting only implemented for 1 or 2 endog states")
#     end

#     # Build as a 1D vector in row-major order so Plots.jl places them correctly:
#     # row r = exogstates[r], col c = errnames[c]; title on top row, xlabel on bottom, ylabel on left
#     plots = [_plot(_get_exst(exogstates[r]), errnames[c];
#                    title      = r == 1 ? string(errnames[c]) : "",
#                    xlabel     = r == N ? string(endogstates[1]) : "",
#                    ylabel     = c == 1 && n_endog == 2 ? string(endogstates[2]) : "",
#                    xformatter = r == N ? :auto : _ -> "",
#                    yformatter = c == 1 ? :auto : _ -> "",
#                    colorbar   = false,
#                    legend     = false)
#              for r in 1:N for c in 1:M]

#     min_width = 200
#     min_height = 200
#     total_size = (M * min_width, N * min_height)

#     if n_endog == 2
#         # Invisible scatter panel: no visible body, colorbar strip only
#         cbar_plot = scatter([0.5], [0.5]; zcolor=[mean(err_lims)], clims=err_lims,
#                             color=:jet, ms=0, markerstrokewidth=0, colorbar=true,
#                             framestyle=:none, xlims=(0,1), ylims=(0,1), label="")
#         l = @layout [grid(N, M) a{0.12w}]
#         plot(plots..., cbar_plot, layout=l, size=(total_size[1] + 80, total_size[2]))
#     else
#         plot(plots..., layout=(N, M), size=total_size)
#     end
# end

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
