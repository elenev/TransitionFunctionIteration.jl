"Abstract base for multidimensional state space grids."
abstract type AbstractStateSpaceGrid end

@generated singleton_keys(::Type{NamedTuple{N,T}}) where {N,T} =
    tuple((N[i] for i in 1:length(N) if fieldtype(T, i) <: Number)...)

"""
    StateSpaceGrid(; name=grid, ...)
    StateSpaceGrid(exogenv, endog_grids::NamedTuple)

Multidimensional state space grid with named axes.

Each keyword argument defines one axis. Scalar-valued fields are marked as
*singletons* — they are stored in the grid but excluded from interpolation axes.
The second constructor prepends the exogenous state indices from `exogenv` as the `:exst` axis.
"""
struct StateSpaceGrid{N,T,S} <: AbstractStateSpaceGrid
    data::NamedTuple{N,T}

    function StateSpaceGrid(data::NamedTuple{N,T}) where {N,T}
        S = singleton_keys(typeof(data))
        new{N,T,S}(data)
    end
end

"Return the names of scalar-valued (singleton) fields in `ssg`. Singletons are excluded from interpolation but appear in state NamedTuples."
singletons(::StateSpaceGrid{N,T,S}) where {N,T,S} = S

StateSpaceGrid(; kwargs...) = StateSpaceGrid(NamedTuple{keys(kwargs)}(values(kwargs)))

"Return the NamedTuple of univariate grid vectors keyed by state variable name."
@inline unigrids(ssg::AbstractStateSpaceGrid) = ssg.data

"Return the total number of grid points (product of all dimension lengths)."
@inline npt(ssg::AbstractStateSpaceGrid) = prod(dims(ssg))

@inline ndims(ssg::AbstractStateSpaceGrid) = length(keys(unigrids(ssg)))

"Return the number of grid points per dimension as a NamedTuple."
@inline dims(ssg::AbstractStateSpaceGrid) = map(length, unigrids(ssg))

"Return the named state dimensions of a grid, or the named output variables of an interpolant."
@inline varnames(ssg::AbstractStateSpaceGrid) = keys(unigrids(ssg))

_points(grids::NamedTuple{NAMES,T}) where {NAMES,T} = (NamedTuple{NAMES}(x) for x in Iterators.product(grids...))

"""
    points(ssg) -> Generator
    points(exogenv) -> StructArray

For a state space grid: lazy generator of NamedTuples over the full Cartesian product
(first dimension varies fastest). For an exogenous environment: the array of exogenous state vectors.
"""
points(ssg::AbstractStateSpaceGrid) = _points(unigrids(ssg))

"Abstract base for interpolated functions over a state space."
abstract type AbstractInterpolant end

"""
    Interpolant(ssg, returnvals::NamedTuple)
    Interpolant(ssg, returnvals::StructArray)

Gridded interpolation of multiple functions over a `StateSpaceGrid`.

The NamedTuple constructor fills every grid point with the same constant values.
The StructArray constructor uses per-point values; its shape must match `dims(ssg)`.

Uses `Gridded(Linear())` for endogenous dimensions and `NoInterp()` for `:exst`.
Extrapolation is `Flat()` for endogenous axes and `Throw()` for `:exst`.
"""
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
    first_values = values(first(returnvals))
    Nof = length(first_values)

    if eltype(first_values) <: AbstractArray 
        # @debug eltype(first_values)
        # Output is vector-valued.
        # Must be the same length for each variable.
        all_lengths = map(length, first_values)
        if length(unique(all_lengths)) != 1
            throw(ArgumentError("All output vectors must have the same length; got lengths: $all_lengths"))
        end
        Nvals = length(first(first_values))
        @debug "Output is vector-valued with length $Nvals; constructing MMatrix of size ($Nof, $Nvals) for each grid point."
        Y = [MMatrix{Nof,Nvals}(reduce(vcat, x)) for x in returnvals]
    else
        Y = [MVector(values(x)) for x in returnvals]
    end

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

    #@debug "Knots: " typeof(values(ugrids)) map(size, values(ugrids))
    #@debug "Values: " typeof(Y) size(Y) Y[1]
    #@debug "Interpolation rules: " interp_rules
    itp = extrapolate(interpolate(values(ugrids), Y, values(interp_rules)), values(extrap_rules) )
    example_returnvals = first(returnvals)
    return Interpolant(ssg, itp, example_returnvals)
end

"""
    (fct::Interpolant)(state::NamedTuple) -> NamedTuple
    (fct::Interpolant)(; kwargs...) -> NamedTuple
    (fct::Interpolant)(pts) -> StructArray

Evaluate the interpolant at one state point or a collection of points.
Singleton dimensions in the state are silently dropped before evaluation.
Calling on a collection returns a StructArray.
"""
@inline function (fct::Interpolant)(; x...)
    return fct((;x...))
end

@inline function (fct::Interpolant)(statevars::NamedTuple)
    statevars = Base.structdiff(statevars, NamedTuple{singletons(statespacegrid(fct))})
    return _evaluate(fct, statevars)
end

@inline function (fct::Interpolant)(pts::Union{AbstractArray,Base.Generator})
    return StructArray(_evaluate(fct,x) for x in pts)
end

@inline function _evaluate(fct, x)
    out = fct.itp(x...)
    return name_values(out, fct)
end

"Return the `StateSpaceGrid` underlying the interpolant."
@inline statespacegrid(fct::AbstractInterpolant) = fct.ssg

@inline coefficients(itp::Interpolations.GriddedInterpolation) = itp.coefs
@inline coefficients(itp::Interpolations.Extrapolation) = coefficients(itp.itp)

"Return interpolation coefficients shaped as an array matching `dims(statespacegrid(fct))`. Each element is a vector of length `nof(fct)`."
@inline coefficients(fct::AbstractInterpolant) = reshape(coefficients(fct.itp), values(dims(fct.ssg)))

@inline knots(itp::Interpolations.GriddedInterpolation) = itp.knots
@inline knots(itp::Interpolations.Extrapolation) = knots(itp.itp)

"Return the grid knots of the interpolant as a tuple of univariate vectors."
@inline knots(fct::AbstractInterpolant) = knots(fct.itp)

@inline npt(fct::AbstractInterpolant) = npt(statespacegrid(fct))

"Return the number of output functions (components) of the interpolant."
@inline nof(fct::AbstractInterpolant) = length(varnames(fct))

@inline elsize(fct::AbstractInterpolant) = (size ∘ first ∘ coefficients)(fct)

# 1. Helper to uniformly extract scalars from vectors and SVectors from matrices
@inline _extract_row(y::AbstractVector, i) = @inbounds y[i]
@inline _extract_row(y::AbstractMatrix, i) = @inbounds SVector{size(y, 2)}(y[i, :])

# 2. The zero-allocation constructor
@inline function _pack_rows(::Type{NamedTuple{K}}, y) where {K}
    # Val(length(K)) forces the compiler to completely unroll the tuple construction
    vals = ntuple(i -> _extract_row(y, i), Val(length(K)))
    return NamedTuple{K}(vals)
end

"Wrap values `x` into a NamedTuple keyed by the output names of `fct`."
@inline name_values(y, fct::AbstractInterpolant) = _pack_rows(NamedTuple{varnames(fct)}, y)

"Abstract base for exogenous environments (Markov chains)."
abstract type AbstractExogenv end

"""
    Exogenv(; pts, mtrans, ssidx)

Exogenous environment for a discrete Markov chain.

- `pts`: StructArray of exogenous state vectors (one row per state)
- `mtrans`: transition probability matrix (rows sum to 1)
- `ssidx`: index of the steady-state (default starting) state
"""
@kwdef struct Exogenv{V<:StructArray,M} <: AbstractExogenv
    pts::V
    mtrans::M
    ssidx::Int
    function Exogenv(pts::V, mtrans::M, ssidx::Int) where {V<:StructArray,M} 
        pts = StaticStructArray(pts) # Convert to SVector for performance
        if length(pts) <= 10
            mtrans = SMatrix{size(mtrans)...}(mtrans) # Convert to SMatrix for performance
            # not worth it for large exnpt
        end
        return new{typeof(pts),typeof(mtrans)}(pts, mtrans, ssidx)
    end
end

"Return the number of variables in each exogenous state vector."
@inline nex(exogenv::Exogenv) = length(exogenv.pts[1])

"Return the number of discrete exogenous states."
@inline exnpt(exogenv::Exogenv) = length(exogenv.pts)

@inline points(exogenv::Exogenv) = exogenv.pts
@inline transition_matrix(exogenv::Exogenv) = exogenv.mtrans
@inline conditional_distro(exogenv::AbstractExogenv, state_idx) = view(transition_matrix(exogenv), state_idx, :)
#@inline conditional_distro(exogenv::Exogenv, state_idx) = transition_matrix(exogenv)[state_idx, :]


"Return the steady-state index of the exogenous environment."
@inline ssidx(exogenv::Exogenv) = exogenv.ssidx

"Return a StructArray of exogenous points given an AbstractArray of exogenous state integers."
function points(exogenv::AbstractExogenv, states::AbstractArray{<:Integer})
    return StructArray(points(exogenv)[states])
end

@kwdef struct NextPeriod{E,C,A,V<:AbstractVector}
    expectations::E
    continuations::C
    prnext::V
    addvars::A
end

E(x, np::NextPeriod) = dot(np.prnext, x)

"Abstract base for user-supplied economic model logic. Implement `calc_eqm` for each subtype."
abstract type AbstractModel end

"Abstract base for assembled solution objects."
abstract type AbstractModelSolution end

"""
    ModelSolution(; model, exogenv, ssg, Pfct, Vfct)

Assembled solution container for an economic model.

- `model`: user-supplied `AbstractModel` instance (economics logic)
- `exogenv`: `Exogenv` for the Markov exogenous process
- `ssg`: combined `StateSpaceGrid` (exogenous + endogenous axes)
- `Pfct`: policy function `Interpolant`
- `Vfct`: value function `Interpolant`
"""
@kwdef struct ModelSolution{M<:AbstractModel,
            E<:AbstractExogenv,
            S<:AbstractStateSpaceGrid,
            PF<:AbstractInterpolant,
            VF<:AbstractInterpolant,
            TF<:AbstractInterpolant} <: AbstractModelSolution
    model::M
    exogenv::E
    ssg::S
    Pfct::PF
    Vfct::VF
    Tfct::TF
end


abstract type AbstractSteadyState end

@kwdef struct SteadyState{M<:AbstractModel, GNAMES, GTYPES, XNAMES, XTYPES, VNAMES, VTYPES, NNAMES, NTYPES} <: AbstractSteadyState
    model::M
    guess::NamedTuple{GNAMES,GTYPES}
    endog::NamedTuple{NNAMES,NTYPES}
    exog::NamedTuple{XNAMES,XTYPES}
    values::NamedTuple{VNAMES,VTYPES}
end

state(ss::SteadyState) = merge(ss.exog, ss.endog)


"Abstract dispatch tag for evaluation context."
abstract type EvaluationMode end

"Evaluation context: Solving for equilibrium policy functions."
struct SolveMode <: EvaluationMode end

"Evaluation context: Iterating (updating) the value function."
struct IterateMode <: EvaluationMode end

"Evaluation context: Simulating the model."
struct SimulateMode <: EvaluationMode end

"Abstract dispatch tag for solution method."
abstract type SolutionMethod end

"Solution method: Transition Function Iteration (TFI) of Elenev and Landvoigt."
struct TFIMethod <: SolutionMethod end

"Solution method: Simultaneous Transition and Policy Function Iteration (STPFI) of Cao, Luo, and Nie."
struct STPFIMethod <: SolutionMethod end

function StateSpaceGrid(exogenv::AbstractExogenv, endog_grids::NamedTuple)
    return StateSpaceGrid(; exst = eachindex(exogenv.pts), endog_grids...)
end

# NonlinearProblem functors
abstract type EquilibriumSystem end

struct TFISystem{S, N, P, M} <: EquilibriumSystem
    fullstate::S
    np::N
    pfct::P
    model::M
end
TFISystem(fullstate, np, obj) = TFISystem(fullstate, np, obj.Pfct, obj.model)

function (sys::TFISystem)(u,p)
    named_policy = name_values(u, sys.pfct)
    return calc_eqm(sys.fullstate, named_policy, sys.model, sys.np, SolveMode())
end

struct STPFISystem{S, R, P, M} <: EquilibriumSystem
    fullstate::S
    prnext::R
    p_nof::Int
    pfct::P
    model::M
end
STPFISystem(fullstate, prnext, p_nof, obj) = STPFISystem(fullstate, prnext, p_nof, obj.Pfct, obj.model)

function (sys::STPFISystem)(u,p)
    named_policy = name_values(u[1:sys.p_nof], sys.pfct)
    np, named_transitions = _compute_expectations_from_trans(u[sys.p_nof+1:end], sys.prnext, sys.obj)
    out = calc_eqm(sys.fullstate, named_policy, sys.model, np, IterateMode())
    transition_residuals = map(-,named_transitions, out.transitions)
    return vcat(out.fx, transition_residuals)
end

## Stubs

# 1. Define the custom exception type
struct MethodNotImplementedError <: Exception
    func_name::Symbol
    type_missing::Type
end

# 2. Overload showerror for clean REPL output
function Base.showerror(io::IO, e::MethodNotImplementedError)
    print(io, "MethodNotImplementedError: Solver requires that `", e.func_name, 
              "` be implemented for `", e.type_missing, "`.")
end

"""
    calc_eqm(state, policy, obj; mode=SolveMode())

Compute equilibrium conditions at `state` given `policy`. Interface stub —
model packages must add methods dispatching on `ModelSolution{<:MyModel}`.

Expected return: a NamedTuple with fields `vals` (value outputs), `fx` (residuals),
`trans` (next-period transient state), and `addvars` (extra simulation variables).
"""
function calc_eqm(state, policy, mdl::AbstractModel, np::NextPeriod, mode::EvaluationMode) 
    # Not Implemented Error 
    throw(MethodNotImplementedError(:calc_eqm, typeof(mdl)))
end

function compute_next_period_values(next, mdl::AbstractModel)
    throw(MethodNotImplementedError(:compute_next_period_values, typeof(mdl)))
end

function Exogenv(mdl::AbstractModel)
    throw(MethodNotImplementedError(:Exogenv, typeof(mdl)))
end

function define_endog_grids(mdl::AbstractModel, ::AbstractSteadyState)
    throw(MethodNotImplementedError(:define_endog_grids, typeof(mdl)))
end