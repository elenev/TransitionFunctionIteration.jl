abstract type AbstractStateSpaceGrid end

@generated singleton_keys(::Type{NamedTuple{N,T}}) where {N,T} =
    tuple((N[i] for i in 1:length(N) if fieldtype(T, i) <: Number)...)

struct StateSpaceGrid{N,T,S} <: AbstractStateSpaceGrid
    data::NamedTuple{N,T}

    function StateSpaceGrid(data::NamedTuple{N,T}) where {N,T}
        S = singleton_keys(typeof(data))
        new{N,T,S}(data)
    end
end

singletons(::StateSpaceGrid{N,T,S}) where {N,T,S} = S

StateSpaceGrid(; kwargs...) = StateSpaceGrid(NamedTuple{keys(kwargs)}(values(kwargs)))
@inline unigrids(ssg::AbstractStateSpaceGrid) = ssg.data
@inline npt(ssg::AbstractStateSpaceGrid) = prod(dims(ssg))
@inline ndims(ssg::AbstractStateSpaceGrid) = length(keys(unigrids(ssg)))
@inline dims(ssg::AbstractStateSpaceGrid) = map(length, unigrids(ssg))
@inline varnames(ssg::AbstractStateSpaceGrid) = keys(unigrids(ssg))

_points(grids::NamedTuple{NAMES,T}) where {NAMES,T} = (NamedTuple{NAMES}(x) for x in Iterators.product(grids...))
points(ssg::AbstractStateSpaceGrid) = _points(unigrids(ssg))

abstract type AbstractInterpolant end
struct Interpolant{SSG, ITP<:Interpolations.AbstractInterpolation, NAMES} <: AbstractInterpolant
    ssg::SSG
    itp::ITP
    function Interpolant(ssg::SSG, itp::ITP, ::NamedTuple{NAMES,T}) where {SSG, ITP, NAMES,T}
        return new{SSG,ITP,NAMES}(ssg, itp)
    end
end
@inline varnames(::Interpolant{SSG,ITP,NAMES}) where {SSG, ITP, NAMES} = NAMES

function Interpolant(ssg::AbstractStateSpaceGrid, returnvals::NamedTuple)
    Y = StructArray(map(x -> deepcopy(returnvals), points(ssg)))
    return Interpolant(ssg, Y)
end

function Interpolant(ssg::AbstractStateSpaceGrid, returnvals::StructArray)
    Y = [MVector(values(x)) for x in returnvals]

    function interpolation_rule(statevar)
        if statevar==:exst
            return NoInterp()           # No interpolation in exog state
        else
            return Gridded(Linear())    # Linear interpolation for other dimensions
        end
    end

    function extrapolation_rule(statevar)
        if statevar==:exst
            return Throw()           # No extrapolation in exog state
        else
            return Flat()            # Constant extrapolation for other dimensions
        end
    end

    ugrids = Base.structdiff(unigrids(ssg), NamedTuple{singletons(ssg)})
    interp_rules = map(interpolation_rule, keys(ugrids))
    extrap_rules = map(extrapolation_rule, keys(ugrids))

    itp = extrapolate(interpolate(values(ugrids), Y, values(interp_rules)), values(extrap_rules) )
    example_returnvals = first(returnvals)
    return Interpolant(ssg, itp, example_returnvals)
end

@inline function (fct::Interpolant)(; x...)
    return fct((;x...))
end

function (fct::Interpolant)(statevars::NamedTuple)
    statevars = Base.structdiff(statevars, NamedTuple{singletons(statespacegrid(fct))})
    return _evaluate(fct, statevars)
end

function (fct::Interpolant)(pts::Union{AbstractArray,Base.Generator})
    return StructArray(_evaluate(fct,x) for x in pts)
end

function _evaluate(fct, x)
    out = fct.itp(x...)
    return name_values(out, fct)
end

@inline statespacegrid(fct::AbstractInterpolant) = fct.ssg

@inline coefficients(itp::Interpolations.GriddedInterpolation) = itp.coefs
@inline coefficients(itp::Interpolations.Extrapolation) = coefficients(itp.itp)
@inline coefficients(fct::AbstractInterpolant) = reshape(coefficients(fct.itp), values(dims(fct.ssg)))

@inline knots(itp::Interpolations.GriddedInterpolation) = itp.knots
@inline knots(itp::Interpolations.Extrapolation) = knots(itp.itp)
@inline knots(fct::AbstractInterpolant) = knots(fct.itp)

@inline npt(fct::AbstractInterpolant) = npt(statespacegrid(fct))
@inline nof(fct::AbstractInterpolant) = length(varnames(fct))
@inline name_values(x, fct::AbstractInterpolant) = NamedTuple{varnames(fct)}(x)

abstract type AbstractExogenv end

@kwdef struct Exogenv{V<:StructArray,M} <: AbstractExogenv
    pts::V
    mtrans::M
    ssidx::Int
end

@inline nex(exogenv::AbstractExogenv) = length(exogenv.pts[1])
@inline exnpt(exogenv::AbstractExogenv) = size(exogenv.pts,1)
@inline points(exogenv::AbstractExogenv) = exogenv.pts
@inline transition_matrix(exogenv::AbstractExogenv) = exogenv.mtrans
@inline ssidx(exogenv::AbstractExogenv) = exogenv.ssidx

# Interface: model-specific economics type; must extend calc_eqm
abstract type AbstractEconomics end
function calc_eqm end  # interface stub; model packages add methods dispatching on GenericModel{M}

abstract type AbstractModel end
@kwdef struct GenericModel{M<:AbstractEconomics,
            E<:AbstractExogenv,
            S<:AbstractStateSpaceGrid,
            PF<:AbstractInterpolant,
            VF<:AbstractInterpolant} <: AbstractModel
    econ::M
    exogenv::E
    ssg::S
    Pfct::PF
    Vfct::VF
end

abstract type EvaluationMode end
struct SolveMode <: EvaluationMode end
struct IterateMode <: EvaluationMode end
struct SimulationMode <: EvaluationMode end
struct SimulateMode <: EvaluationMode end

function StateSpaceGrid(exogenv::AbstractExogenv, endog_grids::NamedTuple)
    return StateSpaceGrid(; exst = eachindex(exogenv.pts), endog_grids...)
end
