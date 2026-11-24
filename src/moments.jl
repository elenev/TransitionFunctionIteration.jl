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

The transformation applied to a non-base economy's value with the base economy's as the basis.
Defaults to the field of the same name.
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
    return path(m.shock), path(:uncond)
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
`economy`: its `transformation` if `economy` is the `base` (a name or position), otherwise its
`compare_transformation` of the economy's value with the base's as the basis.
"""
getmoment(source, m::AbstractMoment) = _transform(transformation(m), momentvalues(source, m)...)

## Moments across economies

function getmoment(economies::NamedTuple, m::AbstractMoment, economy::Symbol, base::Union{Symbol, Integer} = 1)
    base = base isa Integer ? keys(economies)[base] : base
    economy ∈ keys(economies) || throw(KeyError("Economy name $economy not found in economies"))
    value, basis = momentvalues(economies[economy], m)
    economy == base && return _transform(transformation(m), value, basis)
    basevalue, _ = momentvalues(economies[base], m)
    return compare_transformation(m)(value; basis = basevalue)
end
