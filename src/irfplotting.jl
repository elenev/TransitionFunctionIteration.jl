module Transformations
## Structs

abstract type AbstractTransformation end

Base.broadcastable(t::Type{<:AbstractTransformation}) = Ref(t)

# Levels
struct Level        <: AbstractTransformation end
struct AnnLevel     <: AbstractTransformation end
struct AnnInvLevel  <: AbstractTransformation end
struct AnnQuadLevel <: AbstractTransformation end

# Percentages
struct Pct          <: AbstractTransformation end
struct AnnPct       <: AbstractTransformation end
struct AnnInvPct    <: AbstractTransformation end
struct AnnQuadPct   <: AbstractTransformation end

# Differences from (e.g. unconditional) baseline
struct DiffLevel        <: AbstractTransformation end
struct DiffAnnLevel     <: AbstractTransformation end
struct DiffAnnInvLevel  <: AbstractTransformation end
struct DiffAnnQuadLevel <: AbstractTransformation end

struct DiffPct          <: AbstractTransformation end
struct DiffAnnPct       <: AbstractTransformation end
struct DiffAnnInvPct    <: AbstractTransformation end
struct DiffAnnQuadPct   <: AbstractTransformation end

# Ratio differences from (e.g. unconditional) baseline
struct RatioPct       <: AbstractTransformation end
struct RatioAnnPct    <: AbstractTransformation end

# Differences from initial level
struct DiffInitLevel        <: AbstractTransformation end
struct DiffInitPct          <: AbstractTransformation end

## Functions

const ANN = 4
const PCT = 100

is_ratio(::Type{<:Union{RatioPct, RatioAnnPct}}) = true
is_ratio(::Type{<:AbstractTransformation}) = false

function (t::AbstractTransformation)(x;
    basis      = is_ratio(typeof(t)) ? 1.0 : 0.0,
    trend      = zeros(size(x)),
    basistrend = zeros(size(x)),
    init       = x[1])
    return t(x, basis, trend, basistrend, init)
end

(t::Level)(x, basis, trend, basistrend, init) = x
(t::AnnLevel)(x, basis, trend, basistrend, init) = x * ANN
(t::AnnInvLevel)(x, basis, trend, basistrend, init) = x / ANN
(t::AnnQuadLevel)(x, basis, trend, basistrend, init) = x * sqrt(ANN)
(t::Pct)(x, basis, trend, basistrend, init) = x * PCT
(t::AnnPct)(x, basis, trend, basistrend, init) = x * ANN * PCT
(t::AnnInvPct)(x, basis, trend, basistrend, init) = x / ANN
(t::AnnQuadPct)(x, basis, trend, basistrend, init) = x * sqrt(ANN) * PCT
(t::DiffLevel)(x, basis, trend, basistrend, init) = x - basis
(t::DiffAnnLevel)(x, basis, trend, basistrend, init) = (x - basis) * ANN
(t::DiffAnnInvLevel)(x, basis, trend, basistrend, init) = (x - basis) / ANN
(t::DiffAnnQuadLevel)(x, basis, trend, basistrend, init) = (x - basis) * sqrt(ANN)
(t::DiffPct)(x, basis, trend, basistrend, init) = (x - basis) * PCT
(t::DiffAnnPct)(x, basis, trend, basistrend, init) = (x - basis) * ANN * PCT
(t::DiffAnnInvPct)(x, basis, trend, basistrend, init) = (x - basis) / ANN * PCT
(t::DiffAnnQuadPct)(x, basis, trend, basistrend, init) = (x - basis)  * sqrt(ANN) * PCT
(t::RatioPct)(x, basis, trend, basistrend, init) = (x / basis - 1) * PCT
(t::RatioAnnPct)(x, basis, trend, basistrend, init) = (x / basis - 1) * PCT
(t::DiffInitLevel)(x, basis, trend, basistrend, init) = x - init
(t::DiffInitPct)(x, basis, trend, basistrend, init) = (x - init) / init * PCT

const ANNLABEL = "Annual"
const PCTLABEL = "%"
const DIFFLABEL = " Diff. from BGP"
const INITLABEL = "Diff. from Initial"

label(::Type{<:AbstractTransformation}) = ""
label(::Type{<:Level}) = ""
label(::Type{<:AnnLevel}) = ANNLABEL
label(::Type{<:AnnInvLevel}) = ANNLABEL
label(::Type{<:AnnQuadLevel}) = ANNLABEL
label(::Type{<:Pct}) = PCTLABEL
label(::Type{<:AnnPct}) = "$ANNLABEL $PCTLABEL"
label(::Type{<:AnnInvPct}) = "$ANNLABEL $PCTLABEL"
label(::Type{<:AnnQuadPct}) = "$ANNLABEL $PCTLABEL"
label(::Type{<:DiffLevel}) = DIFFLABEL
label(::Type{<:DiffAnnLevel}) = "$DIFFLABEL, $ANNLABEL"
label(::Type{<:DiffAnnInvLevel}) = "$DIFFLABEL, $ANNLABEL"
label(::Type{<:DiffAnnQuadLevel}) = "$DIFFLABEL, $ANNLABEL"
label(::Type{<:DiffPct}) = "$DIFFLABEL ($PCTLABEL)"
label(::Type{<:DiffAnnPct}) = "$DIFFLABEL, $ANNLABEL ($PCTLABEL)"
label(::Type{<:DiffAnnInvPct}) = "$DIFFLABEL, $ANNLABEL ($PCTLABEL)"
label(::Type{<:DiffAnnQuadPct}) = "$DIFFLABEL, $ANNLABEL ($PCTLABEL)"
label(::Type{<:RatioPct}) = "$DIFFLABEL ($PCTLABEL)"
label(::Type{<:RatioAnnPct}) = "$DIFFLABEL, $ANNLABEL ($PCTLABEL)"
label(::Type{<:DiffInitLevel}) = "$INITLABEL"
label(::Type{<:DiffInitPct}) = "$INITLABEL ($PCTLABEL)"


end

using .Transformations

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

    y = rv.transformation.(x; kw...)[_pers .+ 1]

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

