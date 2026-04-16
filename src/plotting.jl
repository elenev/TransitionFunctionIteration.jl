"""
    plot(errs::SimulationErrorPanel, sim:SimulationPanel; kwargs...)

Plot simulation errors as a grid: rows = exogenous states, columns = error variables.

# Keyword arguments
- `errnames`:    tuple of error field names to plot (default: all fields of `errs`)
- `exogstates`:  exogenous state values to include as rows (default: all unique values)
- `endogstates`: 1 or 2 endogenous state symbols for the x (and y) axes (default: first two)
- `stategrid`:   an `AbstractStateSpaceGrid`; if supplied, overlays thin grey lines at each
                 grid node for the selected endogenous states (within the data axis bounds)
- `maxpoints`:   if set, randomly subsample each panel to at most this many points before
                 plotting (default: `nothing` — plot all). Axis limits are always computed
                 from the full data. The same sample is used across all error columns for a
                 given exogenous state row.

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
    endogstates = pop!(plotattributes, :endogstates, endog_states(sim)[1:min(2, length(endog_states(sim)))])
    stategrid   = pop!(plotattributes, :stategrid,   nothing)
    maxpoints   = pop!(plotattributes, :maxpoints,   nothing)

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

    # Grid points filtered to data bounds (so axis limits are unaffected)
    if !isnothing(stategrid)
        ug = unigrids(stategrid)
        xgrid_pts = collect(ug[endogstates[1]])
        # Combined axis limits: union of simulation data and full grid range
        x_lims = (min(x_lims[1], minimum(xgrid_pts)), max(x_lims[2], maximum(xgrid_pts)))
        if n_endog == 2
            ygrid_pts = collect(ug[endogstates[2]])
            y_lims    = (min(y_lims[1], minimum(ygrid_pts)), max(y_lims[2], maximum(ygrid_pts)))
        else
            ygrid_pts = nothing
        end
    end

    # Main grid: row r = exogstates[r], col c = errnames[c]
    for r in 1:N
        mask   = _get_exst(exogstates[r])
        n_mask = sum(mask)
        row_idx = (!isnothing(maxpoints) && n_mask > maxpoints) ?
                      sort(randperm(n_mask)[1:maxpoints]) : Colon()
        for c in 1:M
            errname = errnames[c]
            sp = (r - 1) * M + c

            # Overlay state space grid lines
            if !isnothing(stategrid)
                @series begin
                    subplot    := sp
                    seriestype := :vline
                    linecolor  := :black
                    linewidth  := 0.5
                    label      := ""
                    colorbar   := false
                    legend     := false
                    xlims      := x_lims
                    ylims      := n_endog == 1 ? err_lims : y_lims
                    grid       --> false                    
                    xgrid_pts
                end
                if n_endog == 2
                    @series begin
                        subplot    := sp
                        seriestype := :hline
                        linecolor  := :black
                        linewidth  := 0.5
                        label      := ""
                        colorbar   := false
                        legend     := false
                        xlims      := x_lims
                        ylims      := y_lims
                        grid       --> false
                        ygrid_pts
                    end
                end
            end

            @series begin
                subplot           := sp
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
                    getproperty(sim.data[mask][row_idx],  endogstates[1]),
                    getproperty(errs.data[mask][row_idx], errname)
                else
                    ylims  := y_lims
                    ylabel := c == 1 ? string(endogstates[2]) : ""
                    zcolor := abs.(getproperty(errs.data[mask][row_idx], errname))
                    clims  := err_lims
                    color  --> :jet
                    getproperty(sim.data[mask][row_idx], endogstates[1]),
                    getproperty(sim.data[mask][row_idx], endogstates[2])
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

"""
    plot(sim::SimulationPanel; kwargs...)

Plot a P×P scatterplot matrix of endogenous state (or other) variables.

# Keyword arguments
- `varnames`:  variables to plot (default: all endogenous states)
- `exst`:      if set, filter to rows where `sim.exst == exst` before plotting (default: all)
- `maxpoints`:  if set, randomly subsample to at most this many points for scatter panels
               (default: `nothing` — plot all). Histogram panels always use the full data.
- `stategrid`:  an `AbstractStateSpaceGrid`; if supplied, overlays thin grey lines at each
               grid node for any plotted variable that is an endogenous state (within data
               axis bounds). Scatter panels get vlines and hlines; histogram panels get vlines.

# Layout
Diagonal panels: histogram of that variable (no bar outlines).
Off-diagonal panels: scatter of the column variable (x) against the row variable (y).
All panels show axis labels and tick numerics. Histogram y-axis is labeled "percent".
"""
@recipe function f(sim::SimulationPanel)
    varnames_kw = pop!(plotattributes, :varnames,   endog_states(sim))
    exst_filt   = pop!(plotattributes, :exst,       nothing)
    maxpoints   = pop!(plotattributes, :maxpoints,  nothing)
    stategrid   = pop!(plotattributes, :stategrid,  nothing)

    varnames_kw = collect(varnames_kw)

    mask = isnothing(exst_filt) ? trues(length(sim.exst)) : (sim.exst .== exst_filt)
    data = sim.data[mask]
    P    = length(varnames_kw)

    layout := (P, P)
    size   --> (P * 300, P * 300)
    legend --> false

    # Single subsample index shared across all scatter panels
    n_data   = length(data)
    scat_idx = (!isnothing(maxpoints) && n_data > maxpoints) ?
                   sort(randperm(n_data)[1:maxpoints]) : Colon()

    # Data bounds per variable — used for axis limits and grid-point filtering
    _bounds = Dict(v => extrema(vec(getproperty(data, v))) for v in varnames_kw)

    # Grid points: interior (within data bounds) + one immediate neighbor on each side
    _gridpts = if !isnothing(stategrid)
        ug       = unigrids(stategrid)
        gridkeys = keys(ug)
        Dict(v => begin
            lo, hi = _bounds[v]
            pts    = ug[v]
            below  = filter(x -> x < lo, pts)
            above  = filter(x -> x > hi, pts)
            nb     = isempty(below) ? eltype(pts)[] : [maximum(below)]
            na     = isempty(above) ? eltype(pts)[] : [minimum(above)]
            vcat(nb, filter(x -> lo <= x <= hi, pts), na)
        end for v in varnames_kw if v in gridkeys)
    else
        Dict{Symbol, Vector}()
    end

    for r in 1:P, c in 1:P
        sp   = (r - 1) * P + c
        varc = varnames_kw[c]   # x-axis / histogram variable
        varr = varnames_kw[r]   # y-axis variable (off-diagonal only)

        if r == c
            # 1. data
            @series begin
                subplot    := sp
                seriestype := :histogram
                normalize  := :probability
                linewidth  := 0
                linecolor  := :match
                legend     := false
                vec(getproperty(data, varc))
            end
            # 2. grid lines
            if haskey(_gridpts, varc)
                @series begin
                    subplot    := sp
                    seriestype := :vline
                    linecolor  := :red
                    linewidth  := 0.5
                    primary    := false
                    colorbar   := false
                    legend     := false
                    _gridpts[varc]
                end
            end
            # 3. axis config (last word on subplot attributes)
            @series begin
                subplot    := sp
                seriestype := :scatter
                primary    := false
                legend     := false
                title      := string(varc)
                xlabel     := string(varc)
                ylabel     := ""
                xrotation  := 45
                xgrid      --> !haskey(_gridpts, varc)
                Float64[], Float64[]
            end
        else
            # 1. data
            @series begin
                subplot           := sp
                seriestype        := :scatter
                markerstrokewidth := 0
                markerstrokecolor := :transparent
                legend            := false
                vec(getproperty(data[scat_idx], varc)),
                vec(getproperty(data[scat_idx], varr))
            end
            # 2. grid lines
            if haskey(_gridpts, varc)
                @series begin
                    subplot    := sp
                    seriestype := :vline
                    linecolor  := :red
                    linewidth  := 0.5
                    primary    := false
                    colorbar   := false
                    legend     := false
                    _gridpts[varc]
                end
            end
            if haskey(_gridpts, varr)
                @series begin
                    subplot    := sp
                    seriestype := :hline
                    linecolor  := :red
                    linewidth  := 0.5
                    primary    := false
                    colorbar   := false
                    legend     := false
                    _gridpts[varr]
                end
            end
            # 3. axis config (last word on subplot attributes)
            @series begin
                subplot    := sp
                seriestype := :scatter
                primary    := false
                legend     := false
                title      := string(varc)
                xlabel     := string(varc)
                ylabel     := string(varr)
                xrotation  := 45
                xgrid      --> !haskey(_gridpts, varc)
                ygrid      --> !haskey(_gridpts, varr)
                Float64[], Float64[]
            end
        end
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
