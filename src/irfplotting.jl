using .Transformations

"""
    ReportedVariable(varname; id=varname, label=string(varname), transformation=Level())

Specification for one panel in an IRF plot. Wraps `varname` with an optional display `label`
and a `Transformations` transformation (e.g., `DiffLevel()`, `RatioPct()`).

The `id` type parameter is used as the dispatch key when looking up the variable in recipes.
"""
struct ReportedVariable{ID, T<:Transformations.AbstractTransformation}
    varname::Symbol
    label::String
    transformation::T
end

function ReportedVariable(varname;
    id = varname,
    label = string(varname), 
    transformation::T=Transformations.Level()) where T<:Transformations.AbstractTransformation
    return ReportedVariable{id, T}(varname, label, transformation)
end
ReportedVariable(rv::ReportedVariable) = rv
id(rv::ReportedVariable{ID}) where ID = ID

"""
    IRFPage(vars)

An `AbstractArray` of `ReportedVariable`s that defines one page (subplot group) in an IRF figure.
Construct from a vector or matrix of `ReportedVariable` specs; each element becomes one subplot.

    page = IRFPage([ReportedVariable(:C), ReportedVariable(:K; transformation=DiffLevel())])
"""
struct IRFPage{T<:ReportedVariable,N,A<:AbstractVecOrMat{T}} <: AbstractArray{T,N}
    data::A
end
function IRFPage(data)
    data = ReportedVariable.(data)
    IRFPage{eltype(data), ndims(data), typeof(data)}(data)
end

# --- AbstractMatrix Interface ---
Base.size(sp::IRFPage) = size(sp.data)
Base.getindex(sp::IRFPage, I::Int...) = getindex(sp.data, I...)
Base.IndexStyle(::Type{<:IRFPage}) = IndexLinear()

# --- REPL Visualization ---
function Base.summary(io::IO, sp::IRFPage{T}) where T
    print(io, join(size(sp), "x"), " IRFPage{$T}")
end

# --- Recipes ---

@recipe function f(irf::IRFResult, rv::Union{Symbol,ReportedVariable}, shock::Symbol;
    pers      = nothing,
    basis     = nothing,
    trend     = nothing,
    init      = nothing,
    sd        = false,
    fan       = false,   # true or integer n_bands to plot path fan chart
    n_bands   = 10,      # default number of fan bands
    tail_pct  = 1.0,     # lower tail percentile cutoff (1st–99th by default)
    fan_alpha = 0.8)     # target max cumulative opacity at the center band

    rv = rv isa Symbol ? ReportedVariable(rv) : rv
    basis = isnothing(basis) ? first(keys(irf.means)) : basis

    data  = getproperty(irf.means, shock)
    x     = getproperty(data, rv.varname)
    T     = length(x)
    _pers = isnothing(pers) ? (0:T-1) : pers

    # Build kwargs for transformation -- only pass non-nothing, using correct sources
    kw_basis = isnothing(basis) ? (;) : (; basis = getproperty(getproperty(irf.means, basis), rv.varname))
    kw_trend = isnothing(trend) ? (;) : (; trend = getproperty(data, trend))
    kw_init  = isnothing(init)  ? (;) : (; init  = getproperty(data, init)[1])
    kw = merge(kw_basis, kw_trend, kw_init)

    y = rv.transformation(x; kw...)[_pers .+ 1]

    title     --> rv.label
    yguide    --> Transformations.label(typeof(rv.transformation))
    legend    --> false
    linewidth --> 2

    if fan !== false
        n   = fan === true ? n_bands : Int(fan)
        n_t = length(_pers)

        paths_data = getproperty(irf.paths, shock)
        paths_mat  = getproperty(paths_data, rv.varname)   # N × T_full
        N          = size(paths_mat, 1)

        # Apply transformation to each simulation path
        trans_paths = Matrix{Float64}(undef, N, n_t)
        for i in 1:N
            trans_paths[i, :] = rv.transformation(paths_mat[i, 1:T]; kw...)[_pers .+ 1]
        end

        pct_step = (50.0 - tail_pct) / n

        for k in 1:n
            p_lo = (tail_pct + (k - 1) * pct_step) / 100
            p_hi = 1.0 - p_lo
            y_lo = [quantile(view(trans_paths, :, t), p_lo) for t in 1:n_t]
            y_hi = [quantile(view(trans_paths, :, t), p_hi) for t in 1:n_t]

            @series begin
                primary    := (k == 1)  # band 1 claims the color cycle slot
                seriestype := :path
                fillrange  := y_hi
                # Per-band alpha chosen so cumulative effective alpha grows linearly
                # from 0 to fan_alpha as k goes from 0 to n.
                fillalpha  := 1 - (1 - fan_alpha * k / n) / (1 - fan_alpha * (k - 1) / n)
                linealpha  := 0        # hide border line without affecting fillcolor
                label      := ""
                _pers, y_lo
            end
        end

        # Mean line on top; primary=false so it shares band 1's color slot
        @series begin
            primary   := false
            linewidth := 2
            label     := ""
            _pers, y
        end

    else
        sd_mult = sd === true ? 1 : sd
        if sd_mult != 0
            sd_data   = getproperty(irf.stdevs, shock)
            sd_x      = getproperty(sd_data, rv.varname) .* sd_mult
            y_upper   = rv.transformation(x .+ sd_x; kw...)[_pers .+ 1]
            y_lower   = rv.transformation(x .- sd_x; kw...)[_pers .+ 1]
            ribbon    := (y .- y_lower, y_upper .- y)
            fillalpha := 0.3
        end

        _pers, y
    end
end

@recipe function f(irf::IRFResult, page::IRFPage, shock::Symbol;
    pers      = nothing,
    basis     = nothing,
    trend     = nothing,
    init      = nothing,
    sd        = false,
    fan       = false,
    n_bands   = 10,
    tail_pct  = 1.0,
    fan_alpha = 0.8,
    legendon  = false)   # which subplot gets the legend (false=none, true=1, integer=that subplot)

    legend_subplot = legendon === true ? 1 : legendon

    layout --> size(page)
    for (i, rv) in enumerate(page)
        @series begin
            subplot   := i
            legend    := (i == legend_subplot)
            pers      := pers
            basis     := basis
            trend     := trend
            init      := init
            sd        := sd
            fan       := fan
            n_bands   := n_bands
            tail_pct  := tail_pct
            fan_alpha := fan_alpha
            irf, rv, shock
        end
    end
end

