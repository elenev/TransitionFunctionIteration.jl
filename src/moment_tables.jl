# Moments as a long table: the numbers that moments need, one row per economy and number, the same
# whether the results are in memory or saved with OutputVault. `MomentRows` reads moments from such
# a table, so `compare_economies` works on it unchanged.
#
# Maintainer notes
# - A number is identified by a MomentRow (shape, group, variable, stat, period). `shape` names the
#   saved table the number comes from, as in `result_shape`: :stats (group: the sample, period 0),
#   :irf (group: the shock) or :params (group "", stat "value", period 0). A leaf moment needs one
#   number, or two for an impulse moment whose transformation needs the basis path.
# - The OutputVault method of assemble_moments (in the extension) queries the same numbers and must
#   return the same table; a test checks it.

const MomentRow = NamedTuple{(:shape, :group, :variable, :stat, :period), Tuple{Symbol, String, String, String, Int}}

# a row of the moment table: the run, then the number's MomentRow fields and its value
const MomentRecord = NamedTuple{(:branch, :experdef, :variant, :shape, :group, :variable, :stat, :period, :value),
                                Tuple{String, String, String, Symbol, String, String, String, Int, Float64}}

const MOMENT_TABLE_COLUMNS = fieldnames(MomentRecord)

_record(branch, experdef, variant, r::MomentRow, value) = MomentRecord((branch, experdef, variant, r..., value))

# the numbers a moment needs
_moment_rows(m::SimulationMoment) = [MomentRow((:stats, string(m.sample), string(m.var), string(m.stat), 0))]
function _moment_rows(m::ImpulseMoment)
    rows = [MomentRow((:irf, string(m.shock), string(m.var), string(m.stat), m.period))]
    if Transformations.requires_basis(transformation(m))
        push!(rows, MomentRow((:irf, string(m.basis), string(m.var), string(m.stat), m.period)))
    end
    return rows
end
_moment_rows(m::ParameterMoment) = [MomentRow((:params, "", string(m.name), "value", 0))]
_moment_rows(m::CalculatedMoment) = collect(MomentRow, Iterators.flatmap(_moment_rows, values(m.inputs)))

# a NamedTuple of moments, bare names read as SimulationMoments
_as_moments(moments::NamedTuple) = map(v -> v isa AbstractMoment ? v : SimulationMoment(v), moments)

# the numbers all moments need, each once
_requested_rows(moments) = collect(Set(r for m in values(_as_moments(moments)) for r in _moment_rows(m)))

## Results in memory

"""
    assemble_moments(moments, economies::NamedTuple; branch = "memory", experdef = "default") -> NamedTuple
    assemble_moments(moments, runs)

The numbers that `moments` (a NamedTuple of [`AbstractMoment`](@ref)s or variable names) need, as
a long column table: `branch`, `experdef`, `variant`, `shape` (the saved table the number comes
from: `:stats`, `:irf` or `:params`, as in [`result_shape`](@ref)), `group` (the sample, the shock,
or `""` for a parameter), `variable`, `stat`, `period` (of an IRF; 0 otherwise) and `value`.

With `economies` (the results of each economy in memory, as [`compare_economies`](@ref) takes
them), every row gets `branch` and `experdef`, and the variant is the economy's name. With `runs`
(OutputVault `RunHandle`s or a `RunPattern`; needs `using OutputVault`) the numbers are read from
the saved runs, IRFs from the starting point `start` (keyword, default `:ss`). The same results
give the same rows, though not necessarily in the same order. A number that a result lacks has no
row.
"""
function assemble_moments(moments::NamedTuple, economies::NamedTuple; branch::AbstractString = "memory",
                          experdef::AbstractString = "default")
    rows = _requested_rows(moments)
    found = ((name, r, _row_value(source, r)) for (name, source) in pairs(economies) for r in rows)
    return Tables.columntable(MomentRecord[_record(branch, experdef, string(name), r, v) for (name, r, v) in found if !isnothing(v)])
end

# the value of a number in a source of results, or `nothing` if the source lacks it
function _row_value(r::EconomyResult, row)
    part = row.shape === :stats ? r.stats : row.shape === :irf ? r.irf : r.params
    return isnothing(part) ? nothing : _row_value(part, row)
end
function _row_value(s::SamplesTuple, row)
    if row.shape !== :stats || !haskey(s, Symbol(row.group))
        return nothing
    end
    return _row_value(s[Symbol(row.group)], row)
end
function _row_value(m::MomentArray, row)
    var, stat = Symbol(row.variable), Symbol(row.stat)
    if !(var in _get_vars(m) && stat in _get_stats(m))
        return nothing
    end
    return Float64(getmoment(m, var, stat))
end
function _row_value(irf::IRFResult, row)
    if row.shape !== :irf
        return nothing
    end
    lines = merge((mean = irf.means, std = irf.stdevs), irf.other_stats)
    stat, shock, var = Symbol(row.stat), Symbol(row.group), Symbol(row.variable)
    if !(stat in keys(lines) && shock in keys(lines[stat]))
        return nothing
    end
    line = lines[stat][shock]
    if !(var in propertynames(line) && 1 <= row.period <= length(line))
        return nothing
    end
    return Float64(getproperty(line, var)[row.period])
end
_row_value(params::NamedTuple, row) = row.shape === :params && haskey(params, Symbol(row.variable)) ? Float64(params[Symbol(row.variable)]) : nothing
_row_value(source, row) = nothing

## Reading moments from a long table

"""
    MomentRows

The numbers of one economy from a long table of [`assemble_moments`](@ref), a source that
[`getmoment`](@ref) and [`compare_economies`](@ref) read every moment from.
"""
struct MomentRows
    values::Dict{MomentRow, Float64}
end

function momentvalues(d::MomentRows, m::Union{SimulationMoment, ImpulseMoment, ParameterMoment})
    values = map(_moment_rows(m)) do r
        if !haskey(d.values, r)
            throw(KeyError("the results have no $(r.shape) number $(r.group), $(r.variable), $(r.stat), period $(r.period)"))
        end
        d.values[r]
    end
    return first(values), length(values) == 2 ? values[2] : nothing
end

"""
    split_by_economy(table; order = nothing) -> NamedTuple

The long table of [`assemble_moments`](@ref) as one [`MomentRows`](@ref) per economy, keyed by
variant, the `economies` that [`compare_economies`](@ref) takes. Variants that repeat across
branches or experdefs are an error: select fewer runs. `order` lists the variants to keep, in that
order; without it, variants come in the order of the table.
"""
function split_by_economy(table; order = nothing)
    runs = Dict{String, Set{Tuple{String, String}}}()           # variant => its (branch, experdef)s
    for (b, e, v) in zip(table.branch, table.experdef, table.variant)
        push!(get!(Set{Tuple{String, String}}, runs, v), (b, e))
    end
    clashes = sort!([v for (v, rs) in runs if length(rs) > 1])
    if !isempty(clashes)
        throw(ArgumentError("variants $(join(clashes, ", ")) appear in several runs; select runs with a narrower pattern"))
    end
    dicts = Dict{String, Dict{MomentRow, Float64}}()
    for j in eachindex(table.value)
        row = MomentRow((table.shape[j], table.group[j], table.variable[j], table.stat[j], table.period[j]))
        d = get!(Dict{MomentRow, Float64}, dicts, table.variant[j])
        if haskey(d, row)
            throw(ArgumentError("the table has two values of $row for variant $(table.variant[j])"))
        end
        d[row] = table.value[j]
    end
    names = isnothing(order) ? unique(table.variant) : [string(n) for n in order]
    return NamedTuple{Tuple(Symbol(n) for n in names)}(Tuple(MomentRows(get(dicts, n, Dict{MomentRow, Float64}())) for n in names))
end
