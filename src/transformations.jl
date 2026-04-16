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
    basis      = is_ratio(typeof(t)) ? one.(x) : zero(x),
    trend      = zero(x),
    basistrend = zero(x),
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
(t::RatioPct)(x, basis, trend, basistrend, init) = (x ./ basis .- 1) .* PCT
(t::RatioAnnPct)(x, basis, trend, basistrend, init) = (x ./ basis .- 1) .* PCT
(t::DiffInitLevel)(x, basis, trend, basistrend, init) = x .- init
(t::DiffInitPct)(x, basis, trend, basistrend, init) = (x .- init) ./ init .* PCT

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