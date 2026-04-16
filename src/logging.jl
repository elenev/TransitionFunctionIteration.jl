using Logging
using Printf

"""Log level: solver converged or iteration limit reached."""
const COMPLETE        = LogLevel(200)
"""Log level: per-iteration convergence summary."""
const ITERATION       = LogLevel(100)
"""Log level: some grid points did not converge this iteration."""
const ITERATIONFAILED = LogLevel(50)
"""Log level: a single grid point failed to converge."""
const FAILEDPOINT     = LogLevel(-100)

const CSV_COLUMNS = (:level, :iter, :max_vdist, :max_tdist,
                     :worst_vvar, :worst_tvar, :worst_vpoint, :worst_tpoint, :worst_texst,
                     :n_failed, :n_total, :pass, :n_solved, :n_unsolved,
                     :fp_state, :fp_retcode, :msg)
const CSV_HEADER  = join(CSV_COLUMNS, ",")


# -- Appearance -------
_level_color(level) =
    level == ITERATION       ? :cyan   :
    level == COMPLETE        ? :green  :
    level == ITERATIONFAILED ? :yellow :
    level == FAILEDPOINT     ? :red    : :normal

_level_prefix(level) =
    level == COMPLETE        ? "✅ COMPLETE"      :
    level == ITERATION       ? "🔄 ITERATION"     :
    level == ITERATIONFAILED ? "⚠️ FAILED POINTS" :
    level == FAILEDPOINT     ? "❌ FAILED POINT"  : "" #"[$(Int(level))]"

_level_name(level) =
    level == COMPLETE        ? "COMPLETE"        :
    level == ITERATION       ? "ITERATION"       :
    level == ITERATIONFAILED ? "ITERATIONFAILED" :
    level == FAILEDPOINT     ? "FAILEDPOINT"     : string(Int(level))


# ── Types ─────────────────────────────────────────────────────────────────────

"""
    ExecutionLogger

Abstract supertype for TFI solver loggers. Concrete subtypes:
[`TerminalLogger`](@ref) (colored stderr) and [`CSVLogger`](@ref) (CSV file).
"""
abstract type ExecutionLogger <: AbstractLogger end

Logging.min_enabled_level(logger::ExecutionLogger) = logger.min_level
Logging.shouldlog(logger::ExecutionLogger, level, ::Any, ::Any, ::Any) = level >= logger.min_level

"""
    TerminalLogger(min_level; [io=stderr, varname_width=0, gridpoint_widths=(;)])
    TerminalLogger(min_level, obj::AbstractModelSolution; [io=stderr])

Colored iteration-progress logger. The second form infers all formatting from the model solution.
"""
struct TerminalLogger{S <: IO, GW <: NamedTuple} <: ExecutionLogger
    min_level        :: LogLevel
    io               :: S
    varname_width    :: Int
    gridpoint_widths :: GW
end

TerminalLogger(min_level::LogLevel = ITERATION;
        io = stderr, varname_width::Int = 0, gridpoint_widths::NamedTuple = (;)) =
    TerminalLogger(min_level, io, varname_width, gridpoint_widths)

function TerminalLogger(min_level::LogLevel, obj::AbstractModelSolution; io = stderr)
    w_v = maximum(v -> length(":$v"),       varnames(obj.Vfct); init=1)
    w_t = maximum(v -> length("(:$v, 99)"), varnames(obj.Tfct); init=1)
    TerminalLogger(min_level, io, max(w_v, w_t), _gridpoint_widths(obj.ssg))
end

Base.close(::TerminalLogger) = nothing

"""
    CSVLogger(filename; [min_level=ITERATION])

Logger that writes structured CSV rows to `filename`. Header written on construction.
"""
struct CSVLogger{S <: IO} <: ExecutionLogger
    min_level :: LogLevel
    io        :: S
end

function CSVLogger(filename::AbstractString; min_level::LogLevel = ITERATION)
    io = open(filename, "w")
    println(io, CSV_HEADER)
    CSVLogger(min_level, io)
end

Base.close(logger::CSVLogger) = close(logger.io)

# ── Grid-width helpers ────────────────────────────────────────────────────────

_fmt_width(grid::AbstractArray{<:Integer}) = isempty(grid) ? 1 : maximum(ndigits, grid)
function _fmt_width(grid::AbstractArray{<:Real})
    isempty(grid) && return 5
    int_digits = ndigits(floor(Int, max(1.0, maximum(abs, grid))))
    int_digits + (any(<(0), grid) ? 1 : 0) + 4
end

_gridpoint_widths(ssg) = map(maximum ∘ _fmt_width, unigrids(ssg))

# ── Formatters ────────────────────────────────────────────────────────────────

# Distances
_fmt_val(v::AbstractFloat) = @sprintf("%.4g", v)
_fmt_val(v)                = string(v)

# Variables
_fmt_invar(v::Symbol, ::Nothing)     = ":$v"
_fmt_invar(v::Symbol, exst::Integer) = "(:$v, $exst)"

# Grid point elements
_fmt_gridval(v::Integer,       w::Int) = w > 0 ? lpad(string(v), w) : string(v)
_fmt_gridval(v::AbstractFloat, w::Int) = lpad(@sprintf("%.3f", v), max(w, 5))

# Grid point
_fmt_gridpoint(point::NamedTuple, widths::NamedTuple) = 
    replace(string(map(_fmt_gridval, point, widths)), '"' => "")

# ── Renderers ─────────────────────────────────────────────────────────────────

function _render_csv!(logger::CSVLogger, level, message, kw)
    d = Dict{Symbol,String}(k => _fmt_val(v) for (k,v) in kw)
    d[:level] = _level_name(level)
    haskey(d, :msg) || isempty(message) || (d[:msg] = message)
    println(logger.io, join((get(d, col, "") for col in CSV_COLUMNS), ","))
    flush(logger.io)
end

function _prepare_message(logger::TerminalLogger, kw, ::Val{ITERATION})
    haskey(kw, :max_vdist) && haskey(kw, :max_tdist) || return
    iter         = get(kw, :iter, 0)
    worst_vvar   = get(kw, :worst_vvar,   :?)
    worst_tvar   = get(kw, :worst_tvar,   :?)
    worst_vpoint = get(kw, :worst_vpoint, (;))
    worst_tpoint = get(kw, :worst_tpoint, (;))
    worst_texst  = get(kw, :worst_texst,  nothing)
    W  = logger.varname_width
    gw = logger.gridpoint_widths
    vname = rpad(_fmt_invar(worst_vvar, nothing),     W)
    tname = rpad(_fmt_invar(worst_tvar, worst_texst), W)
    
    vmsg = @sprintf("%4d -- Max V Dist = %7.3f in %s at %s",
            iter, kw[:max_vdist], vname, _fmt_gridpoint(worst_vpoint, gw))
    tmsg = @sprintf("%4d -- Max T Dist = %7.3f in %s at %s",
            iter, kw[:max_tdist], tname, _fmt_gridpoint(worst_tpoint, gw))
    return (vmsg, tmsg)
end

function _prepare_message(::TerminalLogger, kw, ::Val{COMPLETE})
    msg = get(kw, :msg, "")
    return (msg,)
end

function _prepare_message(::TerminalLogger, kw, ::Val{ITERATIONFAILED})
    pass = get(kw, :pass, 0)
    n_failed = get(kw, :n_failed, 0)
    n_total = get(kw, :n_total, 0)
    n_solved = get(kw, :n_solved, 0)
    n_unsolved = get(kw, :n_unsolved, 0)
    revisiting_on = get(kw, :revisiting_on, false)

    if pass == 1
        extramsg = revisiting_on ? " Revisiting..." : ""
        msg = @sprintf("Solver failed at %d points out of %d.%s",n_failed, n_total, extramsg)
    else
        msg = @sprintf("Pass %2d: Revisiting solved %d points. %d remain unsolved.",
            pass, n_solved, n_unsolved)
    end
    return (msg,)
end

function _prepare_message(logger::TerminalLogger, kw, ::Val{FAILEDPOINT})
    fp_state = get(kw, :fp_state, "")
    fp_retcode = get(kw, :fp_retcode, "")

    fp_state = _fmt_gridpoint(fp_state, logger.gridpoint_widths)
    msg = "Solver did not converge: $(fp_retcode) at state $fp_state"
    return (msg,)
end

function _prepare_message(::TerminalLogger, kw, ::Any)
    text = isempty(kw) ? message : join(("$k=$(_fmt_val(v))" for (k,v) in kw), "  ")
    return (text,)
end

# ── handle_message ────────────────────────────────────────────────────────────

function Logging.handle_message(logger::TerminalLogger, level, message,
        _module, group, id, file, line; 
        msg=message, kwargs...)
    level < logger.min_level && return nothing
    kw = (; msg, kwargs...)
    color = _level_color(level)

    msg = _prepare_message(logger, kw, Val(level))
    prefix = _level_prefix(level)
    foreach(t -> printstyled(logger.io, "$prefix: $t\n"; color), msg)

    return nothing
end

function Logging.handle_message(logger::CSVLogger, level, message,
        _module, group, id, file, line; kwargs...)
    level < logger.min_level && return nothing
    _render_csv!(logger, level, message, pairs(kwargs))
    return nothing
end