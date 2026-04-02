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

"Abstract dispatch tag for steady-state solution method."
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



