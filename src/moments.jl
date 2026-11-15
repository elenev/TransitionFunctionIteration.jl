# Moments: specifications of one number per economy and how to read them from results.
#
# Maintainer notes
# - A moment is a lazy specification (`AbstractMoment`); it holds no data. `getmoment(source, moment)`
#   reads it from a source of results, and dispatch on the pair (source, moment) picks the right
#   property: `SimulationMoment`s read the `describe` statistics, `ImpulseMoment`s the IRFs.
# - `momentvalues(source, moment)` is the single reading function: the value before the moment's
#   transformation and the basis the transformation uses (`nothing`: its default). Every `getmoment`
#   method is built on it, so a new source needs only a `momentvalues` method.
# - Sources: `SamplesTuple` (describe output) for simulation moments, `IRFResult` for impulse moments,
#   `EconomyResults` (both), a NamedTuple of economies, and a `MomentTable` (long table, below).
# - `moment_keys(moment)` lists the table entries a moment needs; the DuckDB query in OutputVault and
#   `MomentTable` both use it, so the two cannot drift apart.

"""
    AbstractMoment

A specification of one number per economy: a statistic of a simulated variable
([`SimulationMoment`](@ref)) or a point on an impulse response ([`ImpulseMoment`](@ref)).
Read it with [`getmoment`](@ref), tabulate it with [`compare_economies`](@ref).
"""
abstract type AbstractMoment end

"""
    SimulationMoment(; var, sample=:uncond, stat=:mean, transformation=Transformations.Level,
                       compare_transformation=Transformations.Level)

A statistic of a simulated variable, as in a `compare_economies` table.

- `var`: variable name (Symbol) in the `MomentArray`
- `sample`: filter name used in `describe()` (default `:uncond`)
- `stat`: statistic name (default `:mean`)
- `transformation`: applied to the statistic when it is reported (e.g. levels or percentages)
- `compare_transformation`: applied to a non-base economy's statistic with the base economy's as
  the basis (e.g. `DiffPct`)
"""
@kwdef struct SimulationMoment{T<:Transformations.AbstractTransformation, CT<:Transformations.AbstractTransformation} <: AbstractMoment
    var::Symbol
    sample::Symbol = :uncond
    stat::Symbol   = :mean
    transformation::T = Transformations.Level()
    compare_transformation::CT = Transformations.Level()
end
SimulationMoment(var::Symbol; kwargs...) = SimulationMoment(; var, kwargs...)

"""
    ImpulseMoment(; var, shock, period, stat=:mean, transformation=Transformations.Level,
                    compare_transformation=Transformations.Level)

A point on an impulse response: `stat` of `var` under `shock`, `period` periods into the IRF
(`period` indexes the IRF lines, so it is the `__period` of the IRF table; 1 is the starting period).

The `transformation` is applied with the same statistic of the unconditional path, `uncond`, at the
same period as its basis, so `DiffLevel` is the response and `RatioPct` the percentage response. The
IRF is assumed to start at the stochastic steady state, where `uncond` is the no-shock baseline.
`compare_transformation` is applied to the shock path's values in two economies, as for a
[`SimulationMoment`](@ref).
"""
@kwdef struct ImpulseMoment{T<:Transformations.AbstractTransformation, CT<:Transformations.AbstractTransformation} <: AbstractMoment
    var::Symbol
    shock::Symbol
    period::Int
    stat::Symbol = :mean
    transformation::T = Transformations.Level()
    compare_transformation::CT = Transformations.Level()
end
ImpulseMoment(var::Symbol, shock::Symbol, period::Integer; kwargs...) = ImpulseMoment(; var, shock, period, kwargs...)

"""
    EconomyResults(stats, irf=nothing)

The results of one economy that moments are read from: `stats` is `describe(...).moments` (a
`SamplesTuple`) and `irf` an [`IRFResult`](@ref) (or `nothing`, if no `ImpulseMoment` is read).
"""
struct EconomyResults{S,I}
    stats::S
    irf::I
end
EconomyResults(stats) = EconomyResults(stats, nothing)
EconomyResults(; stats, irf = nothing) = EconomyResults(stats, irf)

## Table entries of a moment

"""
    moment_keys(moment) -> Tuple{Vararg{NamedTuple}}

The entries of a moment table needed to evaluate `moment`, each with fields `kind`
(`:simulation` or `:impulse`), `group` (the sample or the shock), `variable`, `stat` and `period`
(0 for simulation moments). The first entry is the moment's value; an `ImpulseMoment` has a second,
the `uncond` baseline that its transformation uses as the basis.
"""
moment_keys(m::SimulationMoment) = ((; kind = :simulation, group = m.sample, variable = m.var, stat = m.stat, period = 0),)
moment_keys(m::ImpulseMoment) = (
    (; kind = :impulse, group = m.shock, variable = m.var, stat = m.stat, period = m.period),
    (; kind = :impulse, group = :uncond, variable = m.var, stat = m.stat, period = m.period),
)

## Reading a moment: value and basis

"""
    momentvalues(source, moment) -> (value, basis)

The moment's value in `source` before its transformation, and the basis the transformation uses
(`nothing` for the transformation's default).
"""
function momentvalues end

function momentvalues(s::SamplesTuple{SAMPLES}, m::SimulationMoment) where {SAMPLES}
    m.sample ∈ SAMPLES || throw(KeyError("Sample name $(m.sample) not found in describe() output"))
    m_array = s[m.sample]
    m.var  ∈ _get_vars(m_array)  || throw(KeyError("Variable name $(m.var) not found in moments for sample $(m.sample)"))
    m.stat ∈ _get_stats(m_array) || throw(KeyError("Statistic name $(m.stat) not found in moments for sample $(m.sample)"))
    return getmoment(m_array, m.var, m.stat), nothing
end

# the line of a statistic ("mean", "std" or an other_stats name) of one shock, as a StructArray over periods
function _irf_line(irf::IRFResult, stat::Symbol, shock::Symbol)
    lines = merge((mean = irf.means, std = irf.stdevs), irf.other_stats)
    stat ∈ keys(lines) || throw(KeyError("IRF statistic $stat not found; available: $(keys(lines))"))
    shocks = getproperty(lines, stat)
    shock ∈ keys(shocks) || throw(KeyError("IRF shock $shock not found; available: $(keys(shocks))"))
    return getproperty(shocks, shock)
end

function momentvalues(irf::IRFResult, m::ImpulseMoment)
    path(shock) = begin
        line = _irf_line(irf, m.stat, shock)
        m.var ∈ propertynames(line) || throw(KeyError("Variable name $(m.var) not found in the IRF"))
        1 <= m.period <= length(line) || throw(BoundsError(line, m.period))
        getproperty(line, m.var)[m.period]
    end
    return path(m.shock), path(:uncond)
end

momentvalues(r::EconomyResults, m::SimulationMoment) = momentvalues(r.stats, m)
function momentvalues(r::EconomyResults, m::ImpulseMoment)
    isnothing(r.irf) && throw(ArgumentError("moment on shock $(m.shock) needs an IRF, but these EconomyResults have none"))
    return momentvalues(r.irf, m)
end

_transform(t, value, ::Nothing) = t(value)
_transform(t, value, basis) = t(value; basis)

"""
    getmoment(source, moment::AbstractMoment)
    getmoment(economies, moment, economy, base=1)

The moment's value in `source` (a `describe` output for a [`SimulationMoment`](@ref), an
[`IRFResult`](@ref) for an [`ImpulseMoment`](@ref), or [`EconomyResults`](@ref) for either), with
its `transformation` applied.

With `economies` (a NamedTuple of sources, or a [`MomentTable`](@ref)) the moment is read for
`economy`: its `transformation` if `economy` is the `base` (a name or position), otherwise its
`compare_transformation` of the economy's value with the base's as the basis.
"""
getmoment(source, m::AbstractMoment) = _transform(m.transformation, momentvalues(source, m)...)

## Moment tables: the long form of saved results

"""
    MomentTable(table)

A long table of statistics, one row per (economy, statistic), as read from saved results (see
OutputVault's `assemble_moments`). `table` is any Tables.jl source with the columns `economy`,
`kind`, `group`, `variable`, `stat`, `period` and `value`, whose rows are the entries listed by
[`moment_keys`](@ref): a simulation statistic is (`:simulation`, sample, variable, statistic, 0),
an impulse response point (`:impulse`, shock, variable, statistic, period). Names may be Strings or
Symbols. `getmoment(table, moment, economy, base)` reads a moment from it; `table[economy]` is a
source for one economy.
"""
struct MomentTable
    index::Dict{NTuple{6,Any}, Float64}
    economies::Vector{Symbol}
end

const MOMENT_TABLE_COLUMNS = (:economy, :kind, :group, :variable, :stat, :period, :value)

function MomentTable(table)
    Tables.istable(table) || throw(ArgumentError("a MomentTable needs a Tables.jl table"))
    cols = Tables.columns(table)
    missing_cols = setdiff(MOMENT_TABLE_COLUMNS, Tables.columnnames(cols))
    isempty(missing_cols) || throw(ArgumentError("moment table lacks columns $(Tuple(missing_cols))"))
    index = Dict{NTuple{6,Any}, Float64}()
    economies = Symbol[]
    for (e, k, g, v, s, p, x) in zip(map(c -> Tables.getcolumn(cols, c), MOMENT_TABLE_COLUMNS)...)
        key = (String(e), String(k), String(g), String(v), String(s), Int(p))
        haskey(index, key) && throw(ArgumentError("moment table has duplicate rows for $key"))
        index[key] = x
        Symbol(e) ∈ economies || push!(economies, Symbol(e))
    end
    return MomentTable(index, economies)
end
MomentTable(t::MomentTable) = t

Base.keys(t::MomentTable) = Tuple(t.economies)

struct TableEconomy
    table::MomentTable
    economy::Symbol
end
Base.getindex(t::MomentTable, economy::Symbol) = TableEconomy(t, economy)

function momentvalues(e::TableEconomy, m::AbstractMoment)
    values = map(moment_keys(m)) do k
        key = (String(e.economy), String(k.kind), String(k.group), String(k.variable), String(k.stat), k.period)
        haskey(e.table.index, key) || throw(KeyError("moment table has no entry $key"))
        e.table.index[key]
    end
    return first(values), length(values) == 2 ? values[2] : nothing
end

## Moments across economies

function getmoment(economies::Union{NamedTuple, MomentTable}, m::AbstractMoment, economy::Symbol, base::Union{Symbol, Integer} = 1)
    base = base isa Integer ? keys(economies)[base] : base
    economy ∈ keys(economies) || throw(KeyError("Economy name $economy not found in economies"))
    value, basis = momentvalues(economies[economy], m)
    economy == base && return _transform(m.transformation, value, basis)
    basevalue, _ = momentvalues(economies[base], m)
    return m.compare_transformation(value; basis = basevalue)
end
