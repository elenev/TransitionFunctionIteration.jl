function _make_point(statevars, variablestates, fixedstates)
    return (; (k => (k in keys(fixedstates) ? fixedstates[k] : variablestates[k]) for k in statevars)... )
end

function _plotfct2D(fct, varname::Symbol, grid, fh_make_point, transform)
    pts = _points(grid)
    xvec = values(grid)[1]
    fvec = [fct(fh_make_point(p))[varname] for p in pts]
    fvec .= transform.(fvec)
    plot(xvec, fvec, xlabel=string(keys(grid)[1]), title=string(varname), marker=:dot, legend=false)
end

function _plotfct3D(fct, varname::Symbol, grid, fh_make_point, transform)
    x1vec = grid[1]
    x2vec = grid[2]
    pts =  _points(grid)
    fvec = [fct(fh_make_point(p)...)[varname] for p in pts]
    fvec .= transform.(fvec)
    contourf(x1vec, x2vec, fvec; levels=20,  colorbar=true, 
        xlabel=string(keys(grid)[1]), 
        ylabel=string(keys(grid)[2]), 
        title=string(varname))
    end

function Plots.plot(fct::Interpolant, varname::Symbol, fixedstates::NamedTuple; transform = x -> x)
    ssg = statespacegrid(fct)
    statevars = names(ssg)
    grid = (; (k => unigrids(ssg)[k] for k in statevars if k ∉ keys(fixedstates))... )
    Ndims = length(grid)

    fh_make_point = p -> _make_point(statevars, p, fixedstates)

    if Ndims == 1
        _plotfct2D(fct, varname, grid, fh_make_point, transform)
    elseif Ndims == 2
        _plotfct3D(fct, varname, grid, fh_make_point, transform)
    else
        error("Plotting only implemented in 1 and 2 dimensions.")
    end
end