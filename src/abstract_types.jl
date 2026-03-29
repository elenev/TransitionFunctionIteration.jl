"Abstract base for multidimensional state space grids."
abstract type AbstractStateSpaceGrid{T, N} <: AbstractArray{T, N} end

@generated singleton_keys(::Type{NamedTuple{N,T}}) where {N,T} =
    tuple((N[i] for i in 1:length(N) if fieldtype(T, i) <: Number)...)

"Abstract base for interpolated functions over a state space."
abstract type AbstractInterpolant end

"Abstract base for exogenous environments (Markov chains)."
abstract type AbstractExogenv end

"Abstract base for user-supplied economic model logic. Implement `calc_eqm` for each subtype."
abstract type AbstractModel end

"Abstract base for assembled solution objects."
abstract type AbstractModelSolution end

"Abstract base for steady-state descriptions returned by `steadystate`."
abstract type AbstractSteadyState end

"Abstract dispatch tag for evaluation context."
abstract type EvaluationMode end

"Abstract dispatch tag for solution method."
abstract type SolutionMethod end

# Steady state solver modes
abstract type SteadyStateMode end

# 1. Define the custom exception type
struct MethodNotImplementedError <: Exception
    func_name::Symbol
    type_missing::Type
    msg::String
end
function MethodNotImplementedError(func_name::Symbol, type_missing::Type)
    MethodNotImplementedError(func_name, type_missing, ".")
end


# 2. Overload showerror for clean REPL output
function Base.showerror(io::IO, e::MethodNotImplementedError)
    print(io, "MethodNotImplementedError: Solver requires that `", e.func_name,
              "` be implemented for `", e.type_missing, "`", e.msg)
end

"""
    define_endog_grids(model, stv) -> NamedTuple

Return a `NamedTuple` of endogenous grid vectors (e.g., `(K = range(0.5, 10.0, 50),)`).
Called by `main_setup` to build the `StateSpaceGrid`. Interface stub -- implement for your
`AbstractModel` subtype when the default steady-state-proportional grid is not appropriate.
"""
function define_endog_grids(mdl::AbstractModel, ::AbstractSteadyState)
    throw(MethodNotImplementedError(:define_endog_grids, typeof(mdl)))
end

# This is only necessary if steady state is not implemented. Must return a tuple (policy_guess, endog_guess, exogenv)
# of NamedTuples keyed by the policy variable names and endogenous state variable names, respectively.
# Exogenv should contain the 1-state MC in deterministic steady state.
# It can also be nothing, in which case the default Exogenv(mdl) will be used to construct a constant MC at the ergodic mean.
function steadystate_guess(mdl::AbstractModel)
    throw(MethodNotImplementedError(:steadystate_guess, typeof(mdl)))
end

function steadystate_system(::NamedTuple, mdl::AbstractModel, mode::EvaluationMode)
   throw(MethodNotImplementedError(:steadystate_system, typeof(mdl)))
end
