# Moments: specifications of one number per economy and how to read them from results.
#
# Maintainer notes
# - A moment is a lazy specification (`AbstractMoment`); it holds no data. `getmoment(source, moment)`
#   reads it from a source of results, and dispatch on the pair (source, moment) picks the right
#   property: `SimulationMoment`s read the `describe` statistics, `ImpulseMoment`s the IRFs and
#   `ParameterMoment`s the model parameters.
# - `momentvalues(source, moment)` is the single reading function: the value before the moment's
#   transformation and the basis the transformation uses (`nothing`: its default). Every `getmoment`
#   method is built on it, so a new source needs only a `momentvalues` method for each leaf moment
#   (`SimulationMoment`, `ImpulseMoment`, `ParameterMoment`); a `CalculatedMoment` reads its inputs
#   from the same source and needs nothing more.
# - Sources in this package: `SamplesTuple` (describe output) for simulation moments, `IRFResult` for
#   impulse moments, `EconomyResult` (all three kinds). A NamedTuple of sources, keyed by economy, is
#   what the comparisons across economies read. Sources defined outside this package (over tables
#   read from files, say) need only the `momentvalues` methods.
# - `@CalculatedMoment` writes a `CalculatedMoment` as an expression: every function call is
#   replaced by `lifted(f, args...)`, which is the plain call unless an argument is a moment, so any
#   function of numbers works on moments with no list of supported functions.

"""
    AbstractMoment

A specification of one number per economy: a statistic of a simulated variable
([`SimulationMoment`](@ref)), a point on an impulse response ([`ImpulseMoment`](@ref)), a model
parameter ([`ParameterMoment`](@ref)) or a function of other moments ([`CalculatedMoment`](@ref)).
Read it with [`getmoment`](@ref), tabulate it with [`compare_economies`](@ref).

The interface of a moment type: a `momentvalues` method for each kind of source it can be read from,
and the accessors `transformation` and `compare_transformation` (by default the fields of those
names).
"""
abstract type AbstractMoment end

"""
    transformation(moment)

The transformation applied when `moment` is reported. Defaults to the field of the same name.
"""
transformation(m::AbstractMoment) = m.transformation

"""
    compare_transformation(moment)

The transformation that compares a non-base economy with the base: applied to the economy's
reported value (its [`transformation`](@ref) applied) with the base economy's reported value as the
basis. Defaults to the field of the same name.
"""
compare_transformation(m::AbstractMoment) = m.compare_transformation

"""
    SimulationMoment(; var, sample=:uncond, stat=:mean, transformation=Transformations.Level,
                       compare_transformation=Transformations.Level)

A statistic of a simulated variable, as in a `compare_economies` table.

- `var`: variable name (Symbol) in the `MomentArray`
- `sample`: filter name used in `describe()` (default `:uncond`)
- `stat`: statistic name (default `:mean`)
- `transformation`: applied to the statistic when it is reported (e.g. levels or percentages)
- `compare_transformation`: applied to a non-base economy's reported statistic with the base
  economy's as the basis (e.g. `DiffPct`)
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
    ImpulseMoment(; var, shock, period, stat=:mean, basis=:uncond, transformation=Transformations.Level,
                    compare_transformation=Transformations.Level)

A point on an impulse response: `stat` of `var` under `shock`, `period` periods into the IRF
(`period` indexes the IRF lines, so it is the `__period` of the IRF table; 1 is the starting period).

For a transformation that needs a basis, `basis` names the reference shock path (default `:uncond`),
read at the same variable, statistic and period. Thus `DiffLevel` is the response and `RatioPct`
the percentage response relative to that path. A transformation that ignores the basis, such as
`Level`, reads only `shock`; the reference path need not exist. See
[`Transformations.requires_basis`](@ref).
`compare_transformation` compares the reported values (the responses, for `DiffLevel`) of two
economies, as for a [`SimulationMoment`](@ref).
"""
@kwdef struct ImpulseMoment{T<:Transformations.AbstractTransformation, CT<:Transformations.AbstractTransformation} <: AbstractMoment
    var::Symbol
    shock::Symbol
    period::Int
    stat::Symbol = :mean
    basis::Symbol = :uncond
    transformation::T = Transformations.Level()
    compare_transformation::CT = Transformations.Level()
end
ImpulseMoment(var::Symbol, shock::Symbol, period::Integer; kwargs...) = ImpulseMoment(; var, shock, period, kwargs...)

"""
    ParameterMoment(; name, transformation=Transformations.Level, compare_transformation=Transformations.Level)

A model parameter, as saved by [`collect_parameters`](@ref) (the `params` of an [`EconomyResult`](@ref)).
Its value is read with the `transformation` and no basis.
"""
@kwdef struct ParameterMoment{T<:Transformations.AbstractTransformation, CT<:Transformations.AbstractTransformation} <: AbstractMoment
    name::Symbol
    transformation::T = Transformations.Level()
    compare_transformation::CT = Transformations.Level()
end
ParameterMoment(name::Symbol; kwargs...) = ParameterMoment(; name, kwargs...)

"""
    CalculatedMoment(f, inputs; transformation=Transformations.Level(),
                     compare_transformation=Transformations.Level())

A function of other moments: `f(nt)` where `nt` is a `NamedTuple` with the keys of `inputs` and, as
values, the inputs read from the same economy with their own `transformation`s (their
`compare_transformation`s are ignored: a comparison across economies is made once, on the result).
`inputs` is a `NamedTuple` of moments, which may themselves be `CalculatedMoment`s.

```julia
CalculatedMoment(nt -> nt.irf_C / nt.irf_Y, (irf_C = ImpulseMoment(:C, :hike, 3), irf_Y = ImpulseMoment(:Y, :hike, 3)))
```

[`@CalculatedMoment`](@ref) writes the same thing as an expression. The `transformation` is applied
to the result without a basis. The type is inferable when `f` is, and a lambda is a type of its own,
so two specs with equal lambdas are not `==` and a spec cannot be stored by content; use a named
function or a callable struct for that.
"""
struct CalculatedMoment{F, M<:NamedTuple, T<:Transformations.AbstractTransformation, CT<:Transformations.AbstractTransformation} <: AbstractMoment
    f::F
    inputs::M
    transformation::T
    compare_transformation::CT
end
function CalculatedMoment(f, inputs::NamedTuple; transformation = Transformations.Level(), compare_transformation = Transformations.Level())
    all(m -> m isa AbstractMoment, values(inputs)) || throw(ArgumentError("the inputs of a CalculatedMoment must all be moments"))
    return CalculatedMoment(f, inputs, transformation, compare_transformation)
end

## Calculated moments as expressions

# a call `f(args...; kwargs...)` with moments among the arguments, evaluated when the moment is read:
# `template` is the argument list with a `Slot{k}` where the kth moment's value goes
struct Slot{K} end
struct Applied{F, A<:Tuple, K<:NamedTuple}
    f::F
    template::A
    kwargs::K
end
_fill(::Slot{K}, inputs) where {K} = inputs[K]
_fill(x, inputs) = x
(a::Applied)(inputs::NamedTuple) = a.f(map(x -> _fill(x, inputs), a.template)...; a.kwargs...)

"""
    lifted(f, args...; kwargs...)

`f(args...; kwargs...)` if no argument is a moment, otherwise the [`CalculatedMoment`](@ref) of the
moment arguments (`kwargs` are carried along as constants). What [`@CalculatedMoment`](@ref) turns
every call into.
"""
function lifted(f, args...; kwargs...)
    is_moment = map(a -> a isa AbstractMoment, args)
    any(is_moment) || return f(args...; kwargs...)
    positions = [i for i in eachindex(args) if is_moment[i]]
    template = ntuple(i -> is_moment[i] ? Slot{count(j -> is_moment[j], 1:i)}() : args[i], length(args))
    inputs = NamedTuple{ntuple(k -> Symbol(:m, k), length(positions))}(ntuple(k -> args[positions[k]], length(positions)))
    return CalculatedMoment(Applied(f, template, NamedTuple(kwargs)), inputs)
end

# a bare name in an expression with named moments: a variable in scope, else a named moment
function named_moment(names, key::Symbol)
    haskey(names, key) || throw(ArgumentError("$key is neither a variable in scope nor one of the named moments"))
    return names[key]
end

# the expression with every call f(args...) turned into lifted(f, args...) and, if `names` is a
# variable name, every bare name x into `@isdefined(x) ? x : named_moment(names, :x)`
function lift_calls(ex, names)
    lifted_ref, named_ref = GlobalRef(@__MODULE__, :lifted), GlobalRef(@__MODULE__, :named_moment)
    lift(x) = lift_calls(x, names)
    if ex isa Symbol
        isnothing(names) && return ex
        return :(@isdefined($ex) ? $ex : $named_ref($names, $(QuoteNode(ex))))
    elseif !(ex isa Expr)
        return ex
    elseif ex.head === :call
        f, args = ex.args[1], ex.args[2:end]
        is_parameters(a) = a isa Expr && a.head === :parameters      # keywords after `;` come first in a call
        return Expr(:call, lifted_ref, map(lift, filter(is_parameters, args))..., f, map(lift, filter(!is_parameters, args))...)
    elseif ex.head === :parameters
        return Expr(:parameters, map(lift, ex.args)...)
    elseif ex.head === :kw
        return Expr(:kw, ex.args[1], lift(ex.args[2]))
    elseif ex.head === :.                   # a field access: only the object is an expression
        return Expr(:., lift(ex.args[1]), ex.args[2])
    elseif ex.head in (:->, :quote, :macrocall, :function, :(=))
        return ex
    else
        return Expr(ex.head, map(lift, ex.args)...)
    end
end

"""
    @CalculatedMoment expression
    @CalculatedMoment names expression

A [`CalculatedMoment`](@ref) written as an expression of moments. Every function call in the
expression is evaluated as `lifted(f, args...)`: the plain call if no argument is a moment (so
`SimulationMoment(:C)` is still a moment), otherwise a `CalculatedMoment` of the moment arguments, with
the other arguments as constants. Any function of numbers works, operators included, and the result
may be used in further expressions.

```julia
ols = @CalculatedMoment cor_C_Y * sd_C / sd_Y    # three moments in scope
@CalculatedMoment sqrt(SimulationMoment(:K; stat = :std)^2 + 2.0)
```

With `names` (a NamedTuple or Dict of moments), a bare name that is not a variable in scope is read
as `names[:name]`, so `@CalculatedMoment moments mean_C / mean_Y` refers to named moments directly; a
name that is neither is an error. Only calls are lifted: `if`, `&&`, broadcasting and anonymous
functions are left as they are, and a moment passed as a keyword argument is not lifted.
"""
macro CalculatedMoment(ex)
    return esc(lift_calls(ex, nothing))
end
macro CalculatedMoment(names, ex)
    key = gensym(:names)
    return esc(Expr(:let, Expr(:block, Expr(:(=), key, names)), lift_calls(ex, key)))
end

## Model parameters

"""
    collect_parameters(obj_or_model) -> MomentArray

The model's parameters as a table, one row per parameter and one column `value`, as the `params` of
an [`EconomyResult`](@ref) expects them. The default takes
every scalar real or Boolean property of the model (a Boolean is 0 or 1), so a model that forwards the
properties of its parts is covered; a project can add a method for its model type to choose differently.
"""
function collect_parameters(mdl::AbstractModel)
    names = Tuple(p for p in propertynames(mdl) if getproperty(mdl, p) isa Real)
    values = Float64[getproperty(mdl, p) for p in names]
    return MomentArray(reshape(values, :, 1), names, (:value,); rowkey = :parameter)
end
collect_parameters(obj::ModelSolution) = collect_parameters(obj.model)

## Results of one economy

"""
    EconomyResult(; stats=nothing, irf=nothing, params=nothing)

The results of one economy that moments are read from, each optional: `stats` is
`describe(...).moments` (a `SamplesTuple`), `irf` an [`IRFResult`](@ref) and `params` the model's
parameters, as [`collect_parameters`](@ref) returns them. A moment that needs a part the result lacks
is an `ArgumentError` that names it.
"""
@kwdef struct EconomyResult{S, I, P}
    stats::S  = nothing
    irf::I    = nothing
    params::P = nothing
end

## Reading a moment: value and basis

"""
    momentvalues(source, moment) -> (value, basis)

The moment's value in `source` before its transformation, and the basis the transformation uses
(`nothing` for the transformation's default).
"""
function momentvalues end

momentvalues(source, m::AbstractMoment) =
    throw(ArgumentError("a $(nameof(typeof(m))) cannot be read from a $(nameof(typeof(source)))"))

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
    return path(m.shock), Transformations.requires_basis(transformation(m)) ? path(m.basis) : nothing
end

function _parameter(params::MomentArray, name::Symbol)
    name ∈ _get_vars(params) || throw(KeyError("Parameter name $name not found in the parameters"))
    return getmoment(params, name, :value)
end
function _parameter(params::NamedTuple, name::Symbol)
    haskey(params, name) || throw(KeyError("Parameter name $name not found in the parameters"))
    return params[name]
end

function momentvalues(r::EconomyResult, m::SimulationMoment)
    isnothing(r.stats) && throw(ArgumentError("moment on variable $(m.var) needs simulation statistics, but this EconomyResult has none"))
    return momentvalues(r.stats, m)
end
function momentvalues(r::EconomyResult, m::ImpulseMoment)
    isnothing(r.irf) && throw(ArgumentError("moment on shock $(m.shock) needs an IRF, but this EconomyResult has none"))
    return momentvalues(r.irf, m)
end
function momentvalues(r::EconomyResult, m::ParameterMoment)
    isnothing(r.params) && throw(ArgumentError("moment on parameter $(m.name) needs the parameters, but this EconomyResult has none"))
    return _parameter(r.params, m.name), nothing
end

# a calculated moment reads its inputs from the same source; this also nests
momentvalues(source, m::CalculatedMoment) = (m.f(map(c -> getmoment(source, c), m.inputs)), nothing)

_transform(t, value, ::Nothing) = t(value)
_transform(t, value, basis) = t(value; basis)

"""
    getmoment(source, moment::AbstractMoment)
    getmoment(economies, moment, economy, base=1)

The moment's value in `source` (a `describe` output for a [`SimulationMoment`](@ref), an
[`IRFResult`](@ref) for an [`ImpulseMoment`](@ref), or an [`EconomyResult`](@ref) for any), with its
`transformation` applied.

With `economies` (a NamedTuple of sources, keyed by economy name) the moment is read for
`economy`: its reported value (`transformation` applied) if `economy` is the `base` (a name or
position), otherwise `compare_transformation` of the economy's reported value with the base's
reported value as the basis. With the default `Level` comparison, every economy shows its own
reported value.
"""
getmoment(source, m::AbstractMoment) = _transform(transformation(m), momentvalues(source, m)...)

## Moments across economies

function getmoment(economies::NamedTuple, m::AbstractMoment, economy::Symbol, base::Union{Symbol, Integer} = 1)
    base = base isa Integer ? keys(economies)[base] : base
    economy ∈ keys(economies) || throw(KeyError("Economy name $economy not found in economies"))
    value = getmoment(economies[economy], m)
    if economy == base
        return value
    end
    return compare_transformation(m)(value; basis = getmoment(economies[base], m))
end

## Moments as a long table
#
# The numbers that moments need, one row per economy and number, the same whether the results are
# in memory or saved with OutputVault. `MomentRows` reads moments from such a table, so
# `compare_economies` works on it unchanged.
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
