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

MomentArray(data, rownames, colnames) = MomentArray{tuple(rownames...), tuple(colnames...)}(data)

# Names are free: extracted directly from type parameters
@inline row_names(::MomentArray{ROWS}) where {ROWS} = ROWS
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

endog_states(::SimulationPanel{T, N, NAMES, A}) where {T, N, NAMES, A} = NAMES
states(sp::SimulationPanel) = (:exst, endog_states(sp)...)
