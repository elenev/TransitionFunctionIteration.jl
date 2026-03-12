using Logging

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
    COMPLETE => "✅ COMPLETE",
    ITERATION => "▶️ ITERATION",
    ITERATIONFAILED => "▶️ FAILED POINTS",
    FAILEDPOINT => "❌ FAILED POINT"
)

const DEFAULT_COLOR_MAP = Dict(
    ITERATION => :cyan,
    FAILEDPOINT => :red,
    ITERATIONFAILED => :yellow,
    COMPLETE => :green,
)

"""
    ExecutionLogger(min_level::LogLevel)

Colored, prefixed logger for solver progress. Messages below `min_level` are suppressed.
Use with `Logging.with_logger(ExecutionLogger(ITERATION)) do ... end`.
"""
struct ExecutionLogger <: AbstractLogger
    min_level::LogLevel
    prefix_map::Dict{LogLevel, String}
    color_map::Dict{LogLevel, Symbol}
end
ExecutionLogger(min_level::LogLevel) = ExecutionLogger(min_level, DEFAULT_PREFIX_MAP, DEFAULT_COLOR_MAP)
Logging.min_enabled_level(logger::ExecutionLogger) = logger.min_level

# It checks if a message should be logged based on its level.
function Logging.shouldlog(logger::ExecutionLogger, level, _module, group, id)
    return level >= logger.min_level
end

# --- 3. Implement the message handling function ---
# This function defines what happens when a log message is sent.
# We ignore all the other arguments (file, line, module) to keep the output clean.
function Logging.handle_message(
    logger::ExecutionLogger,
    level,
    message,
    _module,
    group,
    id,
    file,
    line;
    kwargs...
)
    # Check if the message's level is at or above the logger's minimum level
    if level >= logger.min_level
        # Get the prefix and color from the maps, or use defaults
        prefix = get(logger.prefix_map, level, "[$(string(level))]")
        color = get(logger.color_map, level, :normal)
        
        # Use printstyled to apply the color
        # This will print the prefix and message in a single line
        printstyled(stderr, "$prefix: $message\n"; color=color)
    end
    nothing
end

# Custom macros for each log level
"Log a message at the `COMPLETE` level."
macro complete(ex)
    @nospecialize
    level = COMPLETE
    _module = __module__
    :($Logging.with_logger($Logging.current_logger()) do
        $Logging.@logmsg($level, $(esc(ex)))
    end)
end

"Log a message at the `ITERATION` level."
macro iteration(ex)
    @nospecialize
    level = ITERATION
    _module = __module__
    :($Logging.with_logger($Logging.current_logger()) do
        $Logging.@logmsg($level, $(esc(ex)))
    end)
end

"Log a message at the `ITERATIONFAILED` level."
macro iterationfailed(ex)
    @nospecialize
    level = ITERATIONFAILED
    _module = __module__
    :($Logging.with_logger($Logging.current_logger()) do
        $Logging.@logmsg($level, $(esc(ex)))
    end)
end

"Log a message at the `FAILEDPOINT` level."
macro failedpoint(ex)
    @nospecialize
    level = FAILEDPOINT
    _module = __module__
    :($Logging.with_logger($Logging.current_logger()) do
        $Logging.@logmsg($level, $(esc(ex)))
    end)
end
