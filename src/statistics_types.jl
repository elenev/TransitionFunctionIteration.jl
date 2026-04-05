using .Transformations

## Statistics Function Sets 
struct StatisticsSet{NAMES, FUNCS<:Tuple}
    funcs::FUNCS
end

# Outer constructor
function StatisticsSet(names::NTuple{N, Symbol}, funcs::Tuple) where N
    return StatisticsSet{names, typeof(funcs)}(funcs)
end
function StatisticsSet(nt::NamedTuple{NAMES, FUNCS}) where {NAMES, FUNCS<:Tuple}
    return StatisticsSet{NAMES, FUNCS}(values(nt))
end

StatisticsSet(ss::StatisticsSet) = ss  # Identity for already-constructed sets

# The "Call" method
function (statset::StatisticsSet{NAMES})(args...) where NAMES
    nested = fmap(statset.funcs, args...)
    flattened = (Iterators.flatten(nested)...,)
    return NamedTuple{NAMES}(flattened)
end

# Length methods
Base.length(::Type{<:StatisticsSet{NAMES}}) where NAMES = length(NAMES)
Base.length(statset::StatisticsSet) = length(typeof(statset))
Base.merge(s1::StatisticsSet{NAMES1}, s2::StatisticsSet{NAMES2}) where {NAMES1, NAMES2} = 
    StatisticsSet( (NAMES1..., NAMES2...), (s1.funcs..., s2.funcs...) )
varnames(::StatisticsSet{NAMES}) where NAMES = NAMES

## Simulation Stats

"""
    MomentArray(data, rows, cols)

A labeled matrix of simulation moments. `rows` and `cols` are tuples of unique
symbols baked into the type, enabling compile-time lookup via `getmoment`.
"""
struct MomentArray{ROWS, COLS, T, M <: AbstractMatrix{T}} <: AbstractMatrix{T}
    data::M

    function MomentArray{ROWS, COLS}(data::M) where {ROWS, COLS, M<:AbstractMatrix}
        allunique(ROWS) || throw(ArgumentError("Row names must be unique: $ROWS"))
        allunique(COLS) || throw(ArgumentError("Column names must be unique: $COLS"))
        size(data, 1) == length(ROWS) || throw(DimensionMismatch(
            "data has $(size(data,1)) rows but $(length(ROWS)) row names"))
        size(data, 2) == length(COLS) || throw(DimensionMismatch(
            "data has $(size(data,2)) columns but $(length(COLS)) column names"))
        new{ROWS, COLS, eltype(M), M}(data)
    end
end

# MomentArray(data, row_names::Tuple{Vararg{Symbol}}, col_names::Tuple{Vararg{Symbol}}) =
#     MomentArray(data, Val(row_names), Val(col_names))

# Convenience constructors for types whose names are statically known
const _HasVarnames = Union{StructArray, NamedTuple, StatisticsSet, Tuple{Vararg{Symbol}}, Val{N} where N}
_varnames(x::StructArray)               = propertynames(x)
_varnames(::NamedTuple{N})    where N   = N
_varnames(::StatisticsSet{N}) where N   = N
_varnames(::Val{N})           where N   = N
_varnames(names::Tuple{Vararg{Symbol}}) = names
MomentArray(data, rows::_HasVarnames, cols::_HasVarnames) =
    MomentArray{_varnames(rows), _varnames(cols)}(data)

macro MomentArray(data, rows, cols)
    return :( MomentArray{$(esc(rows)), $(esc(cols))}($(esc(data))) )
end

# Names are free: extracted directly from type parameters
"Return the row names (variable names) of a `MomentArray` as a tuple of Symbols."
@inline row_names(::MomentArray{ROWS}) where {ROWS} = ROWS
"Return the column names (statistic names) of a `MomentArray` as a tuple of Symbols."
@inline col_names(::MomentArray{ROWS, COLS}) where {ROWS, COLS} = COLS

# --- AbstractArray Interface ---
Base.size(m::MomentArray) = size(m.data)
Base.IndexStyle(::Type{<:MomentArray}) = IndexCartesian()

# Primary scalar access (required by IndexCartesian)
@inline Base.getindex(A::MomentArray, i::Int, j::Int) = A.data[i, j]
@inline Base.setindex!(A::MomentArray, v, i::Int, j::Int) = setindex!(A.data, v, i, j)


# function Base.similar(A::MomentArray{ROWS, COLS}, ::Type{T}, dims::Dims) where {ROWS, COLS, T}
#     return MomentArray{ROWS, COLS}(similar(A.data, T, dims))
# end

Base.BroadcastStyle(::Type{<:MomentArray}) = Broadcast.ArrayStyle{MomentArray}()

find_ma(bc::Broadcast.Broadcasted) = find_ma(bc.args)
find_ma(args::Tuple) = find_ma(args[1], Base.tail(args))
find_ma(x) = x
find_ma(::Tuple{}) = nothing
find_ma(a::MomentArray, rest) = a
find_ma(::Any, rest) = find_ma(rest)

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{MomentArray}}, ::Type{T}) where {T}
    ma = find_ma(bc)
    return MomentArray(similar(Array{T}, axes(bc)), row_names(ma), col_names(ma))
end

# --- REPL Printing ---
function Base.show(io::IO, ::MIME"text/plain", m::MomentArray)
    println(io, "MomentArray ($(size(m, 1)) × $(size(m, 2))):")
    pretty_table(io, m.data; column_labels = collect(col_names(m)), row_labels = collect(row_names(m)))
end

function Base.show(io::IO, m::MomentArray)
    print(io, "MomentArray{$(eltype(m))}($(size(m,1)) × $(size(m,2)))")
end

# --- Named Access ---
# findfirst over a type-parameter tuple is constant-folded when the symbol is known at compile time
@inline function row_index(::MomentArray{ROWS}, var::Symbol) where {ROWS}
    idx = findfirst(==(var), ROWS)
    idx === nothing && throw(KeyError(var))
    return idx
end

@inline function col_index(::MomentArray{ROWS, COLS}, stat::Symbol) where {ROWS, COLS}
    idx = findfirst(==(stat), COLS)
    idx === nothing && throw(KeyError(stat))
    return idx
end

_expand_colon(needle, haystack) = needle isa Colon ? haystack : needle

# Scalar / Scalar -> Number
@inline getmoment(m::MomentArray, var::Symbol, stat::Symbol) =
    m.data[row_index(m, var), col_index(m, stat)]

# Collection / Scalar -> mapped collection
function getmoment(m::MomentArray, vars::Union{Tuple, AbstractVector, Colon}, stat::Symbol)
    vars = _expand_colon(vars, row_names(m))
    return NamedTuple{vars}(map(v -> getmoment(m, v, stat), vars))
end

# Scalar / Collection -> NamedTuple
function getmoment(m::MomentArray, var::Symbol, stats::Union{Tuple, AbstractVector, Colon})
    stats = _expand_colon(stats, col_names(m))
    #keys_tup = Tuple(stats)
    return NamedTuple{stats}(map(s -> getmoment(m, var, s), stats))
end

# Collection / Collection -> MomentArray (Val-wrapped so requested names are compile-time)
@inline function getmoment(m::MomentArray, ::Val{VARS}, ::Val{STATS}) where {VARS, STATS}
    row_inds = SVector(map(v -> row_index(m, v), VARS))
    col_inds = SVector(map(s -> col_index(m, s), STATS))
    return MomentArray(view(m.data, row_inds, col_inds), VARS, STATS)
end

@inline function getmoment(m::MomentArray, vars::Union{Tuple, AbstractVector, Colon}, stats::Union{Tuple, AbstractVector, Colon})
    vars = _expand_colon(vars, row_names(m))
    stats = _expand_colon(stats, col_names(m))
    return getmoment(m, Val(vars), Val(stats))
end

(m::MomentArray)(vars,stats) = getmoment(m, vars, stats)


## Equation Errors
"""
    SimulationErrorPanel(data::StructArray)

Thin wrapper around a 2-D `StructArray` (periods × runs) of equation residuals produced by
`simulate`. Each named field is one residual series (e.g., `:euler`). Passed to
`calc_error_stats` for summary statistics.
"""
struct SimulationErrorPanel{T, N, A<:StructArray{T, N}} <: AbstractArray{T,N}
    data::A
end

# --- AbstractMatrix Interface ---
Base.size(sp::SimulationErrorPanel) = size(sp.data)
Base.getindex(sp::SimulationErrorPanel, I::Int...) = getindex(sp.data, I...)
Base.getindex(sp::SimulationErrorPanel, I...) = SimulationErrorPanel(sp.data[I...])
Base.view(sp::SimulationErrorPanel, I...) = SimulationErrorPanel(view(sp.data, I...))
Base.IndexStyle(::Type{<:SimulationErrorPanel}) = IndexLinear()

# --- Property Access Delegation ---
Base.propertynames(sp::SimulationErrorPanel) = propertynames(getfield(sp, :data))
Base.getproperty(sp::SimulationErrorPanel, s::Symbol) = (s === :data) ? getfield(sp, :data) : getproperty(getfield(sp, :data), s)

# --- REPL Visualization ---
function Base.summary(io::IO, sp::SimulationErrorPanel{T,N}) where {T,N}
    print(io, join(size(sp), "x"), " SimulationErrorPanel{$T,$N} with fields: (", join(propertynames(sp), ", "), ")")
end


### Simulated output
"""
    SimulationPanel(data::StructArray, endog_names)

Thin wrapper around a 2-D `StructArray` (periods × runs) of simulated variables. The
`NAMES` type parameter records which fields are endogenous states (used by `linearize` and
`endog_states`). All fields are accessible via `panel.varname`.
"""
struct SimulationPanel{T, N, NAMES, A<:StructArray{T, N}} <: AbstractArray{T,N}
    data::A
end
SimulationPanel(data::A, endog_names) where {T, N, A<:StructArray{T, N}} = SimulationPanel{T, N, Tuple(endog_names), A}(data)

# --- AbstractMatrix Interface ---
Base.size(sp::SimulationPanel) = size(sp.data)
Base.getindex(sp::SimulationPanel, I::Int...) = getindex(sp.data, I...)
# Slice indexing: delegate to StructArray (returns StructArray) then rewrap
Base.getindex(sp::SimulationPanel{T,N,NAMES}, I...) where {T,N,NAMES} =
    SimulationPanel(sp.data[I...], NAMES)
Base.view(sp::SimulationPanel{T,N,NAMES}, I...) where {T,N,NAMES} =
    SimulationPanel(view(sp.data, I...), NAMES)
Base.IndexStyle(::Type{<:SimulationPanel}) = IndexLinear()

# --- Property Access Delegation ---
Base.propertynames(sp::SimulationPanel) = propertynames(getfield(sp, :data))
Base.getproperty(sp::SimulationPanel, s::Symbol) = (s === :data) ? getfield(sp, :data) : getproperty(getfield(sp, :data), s)

# --- REPL Visualization ---
function Base.summary(io::IO, sp::SimulationPanel{T,N}) where {T,N}
    print(io, join(size(sp), "x"), " SimulationPanel{$T,$N} with fields: (", join(propertynames(sp), ", "), ")")
end

"Return the tuple of endogenous state variable names encoded in the `SimulationPanel` type."
endog_states(::SimulationPanel{T, N, NAMES, A}) where {T, N, NAMES, A} = NAMES
"Return `(:exst, endog_states(sp)...)` -- all state variable names including the exogenous index."
states(sp::SimulationPanel) = (:exst, endog_states(sp)...)


"""
    ReportedMoment(; var, sample=:uncond, stat=:mean, transformation=Transformations.Level)

Specification for a single entry in a `compare_economies` table.

- `var`: variable name (Symbol) in the `MomentArray`
- `sample`: filter name used in `describe()` (default `:uncond`)
- `stat`: statistic name (default `:mean`)
- `transformation`: optional T<:AbstractTransformation to apply to the raw moment (e.g., for reporting in levels or percentages)
"""
@kwdef struct ReportedMoment{T<:Transformations.AbstractTransformation, CT<:Transformations.AbstractTransformation}
    var::Symbol
    sample::Symbol = :uncond
    stat::Symbol   = :mean
    transformation::T = Transformations.Level()
    compare_transformation::CT = Transformations.Level()
end
ReportedMoment(var::Symbol; kwargs...) = ReportedMoment(; var, kwargs...)

const SamplesTuple{SAMPLES, N, VARIABLES, STATS} = NamedTuple{SAMPLES, T} where {T <: NTuple{N, MomentArray{VARIABLES, STATS}}} 
@inline function getmoment(s::SamplesTuple{SAMPLES, N, VARIABLES, STATS}, rm::ReportedMoment) where {SAMPLES, N, VARIABLES, STATS}
    rm.sample ∈ SAMPLES   || throw(KeyError("Sample name $(rm.sample) not found in describe() output"))
    rm.var    ∈ VARIABLES || throw(KeyError("Variable name $(rm.var) not found in moments for sample $(rm.sample)"))
    rm.stat   ∈ STATS     || throw(KeyError("Statistic name $(rm.stat) not found in moments for sample $(rm.sample)"))
    return rm.transformation(getmoment(s[rm.sample], rm.var, rm.stat))
end

const EconomiesTuple{ECONS, N, SAMPLES} = NamedTuple{ECONS, NTuple{N, SAMPLES}}
@inline function getmoment(nt::EconomiesTuple, rm::ReportedMoment, economy::Symbol, base::Union{Symbol, Integer}=1)
    base = base isa Integer ? keys(nt)[base] : base
    economy ∈ keys(nt) || throw(KeyError("Economy name $economy not found in economies tuple"))
    val     = getmoment(nt[economy], rm)
    valbase = getmoment(nt[base], rm)
    if economy == base 
        return rm.transformation(val, basis=valbase)
    else
        return rm.compare_transformation(val, basis=valbase)
    end
end