using Logging
using Printf

# Custom log levels for solver progress
"""Log level: solver converged or iteration limit reached."""
const COMPLETE = LogLevel(200)

"""Log level: per-iteration convergence summary."""
const ITERATION = LogLevel(100)

"""Log level: some grid points did not converge this iteration."""
const ITERATIONFAILED = LogLevel(50)

"""Log level: a single grid point failed to converge."""
const FAILEDPOINT = LogLevel(-100)

# Mappings
const DEFAULT_PREFIX_MAP = Dict(
    COMPLETE         => "✅ COMPLETE",
    ITERATION        => "▶️  ITERATION",
    ITERATIONFAILED  => "⚠️  FAILED POINTS",
    FAILEDPOINT      => "❌ FAILED POINT"
)

const DEFAULT_COLOR_MAP = Dict(
    ITERATION        => :cyan,
    FAILEDPOINT      => :red,
    ITERATIONFAILED  => :yellow,
    COMPLETE         => :green,
)

# --- CSV column schema (fixed-width rows; missing fields left blank) ---
const CSV_COLUMNS = (:level, :iter, :max_vdist, :max_tdist,
                     :worst_vvar, :worst_tvar, :worst_vpoint, :worst_tpoint, :worst_texst,
                     :n_failed, :n_total, :pass, :n_solved, :n_unsolved,
                     :fp_state, :fp_retcode, :msg)
const CSV_HEADER = join(CSV_COLUMNS, ",")

"""
    ExecutionLogger(min_level; [io=stderr])
    ExecutionLogger(filename::AbstractString; [min_level=ITERATION])

Colored, prefixed logger for solver progress. In the first form, logs to `io`
(default `stderr`) with ANSI color. In the second form, opens `filename` and
writes comma-separated rows for machine parsing (no color).

Use inside `Logging.with_logger(ExecutionLogger(...)) do ... end` or pass via
the `logging` keyword of `main_execute!`.

Note: the current implementation logs to `stderr` (not `stdout`) for the colored
mode, matching Julia's Logging convention.
"""
struct ExecutionLogger <: AbstractLogger
    min_level::LogLevel
    io::IO
    csv_mode::Bool
    csv_header_written::Ref{Bool}
    prefix_map::Dict{LogLevel, String}
    color_map::Dict{LogLevel, Symbol}
end

function ExecutionLogger(min_level::LogLevel = ITERATION;
        io::IO = stderr,
        prefix_map = DEFAULT_PREFIX_MAP,
        color_map  = DEFAULT_COLOR_MAP)
    ExecutionLogger(min_level, io, false, Ref(false), prefix_map, color_map)
end

function ExecutionLogger(filename::AbstractString;
        min_level::LogLevel = ITERATION,
        prefix_map = DEFAULT_PREFIX_MAP,
        color_map  = DEFAULT_COLOR_MAP)
    io = open(filename, "w")
    ExecutionLogger(min_level, io, true, Ref(false), prefix_map, color_map)
end

Base.close(logger::ExecutionLogger) = logger.csv_mode ? close(logger.io) : nothing

Logging.min_enabled_level(logger::ExecutionLogger) = logger.min_level

function Logging.shouldlog(logger::ExecutionLogger, level, _module, group, id)
    return level >= logger.min_level
end

# Format a single value: floats get 4 significant figures, others are shown as-is
_fmt_val(v::AbstractFloat) = @sprintf("%.4g", v)
_fmt_val(v) = string(v)

function Logging.handle_message(logger::ExecutionLogger, level, message,
        _module, group, id, file, line; kwargs...)

    level < logger.min_level && return nothing
    kw = pairs(kwargs)

    if logger.csv_mode
        # --- CSV mode ---
        if !logger.csv_header_written[]
            println(logger.io, CSV_HEADER)
            logger.csv_header_written[] = true
        end
        # Build a dict for fast lookup of provided kwargs
        kw_dict = Dict(string(k) => _fmt_val(v) for (k, v) in kw)
        level_str = get(Dict(
            ITERATION       => "ITERATION",
            ITERATIONFAILED => "ITERATIONFAILED",
            FAILEDPOINT     => "FAILEDPOINT",
            COMPLETE        => "COMPLETE",
        ), level, string(level))
        kw_dict["level"] = level_str
        if !isempty(message)
            kw_dict["msg"] = message
        end
        row = join((get(kw_dict, string(col), "") for col in CSV_COLUMNS), ",")
        println(logger.io, row)
        flush(logger.io)
    else
        # --- Colored terminal mode (stderr) ---
        prefix = get(logger.prefix_map, level, "[$(string(level))]")
        color  = get(logger.color_map,  level, :normal)

        if isempty(kw)
            # Fall back to raw message string
            printstyled(logger.io, "$prefix: $message\n"; color)
        else
            # Format kwargs as "key=value" pairs, floats rounded to 4 sig figs
            parts = ["$k=$(_fmt_val(v))" for (k, v) in kw]
            printstyled(logger.io, "$prefix: $(join(parts, "  "))\n"; color)
        end
    end
    return nothing
end

# Custom macros for each log level.
# Accept Julia @logmsg-style keyword arguments; message is optional.
"Log structured data at the `COMPLETE` level."
macro complete(exs...)
    :($Logging.@logmsg($COMPLETE, "", $(exs...)))
end

"Log structured data at the `ITERATION` level."
macro iteration(exs...)
    :($Logging.@logmsg($ITERATION, "", $(exs...)))
end

"Log structured data at the `ITERATIONFAILED` level."
macro iterationfailed(exs...)
    :($Logging.@logmsg($ITERATIONFAILED, "", $(exs...)))
end

"Log structured data at the `FAILEDPOINT` level."
macro failedpoint(exs...)
    :($Logging.@logmsg($FAILEDPOINT, "", $(exs...)))
end
