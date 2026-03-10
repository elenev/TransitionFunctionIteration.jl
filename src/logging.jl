using Logging

# --- 1. Define the custom log levels ---
# We'll use LogLevel constants to make them easy to use.
# LogLevel(-100) is a value between Debug and Info.
const COMPLETE = LogLevel(200)
const ITERATION = LogLevel(100)
const ITERATIONFAILED = LogLevel(50)
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

# --- 2. Define the custom logger struct ---
# This logger will hold the minimum level it should print.
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

# --- 4. Define the custom macros ---
# These macros provide a clean way to use our new log levels.
macro complete(ex)
    @nospecialize
    level = COMPLETE
    _module = __module__
    :($Logging.with_logger($Logging.current_logger()) do
        $Logging.@logmsg($level, $(esc(ex)))
    end)
end

macro iteration(ex)
    @nospecialize
    level = ITERATION
    _module = __module__
    :($Logging.with_logger($Logging.current_logger()) do
        $Logging.@logmsg($level, $(esc(ex)))
    end)
end

macro iterationfailed(ex)
    @nospecialize
    level = ITERATIONFAILED
    _module = __module__
    :($Logging.with_logger($Logging.current_logger()) do
        $Logging.@logmsg($level, $(esc(ex)))
    end)
end

macro failedpoint(ex)
    @nospecialize
    level = FAILEDPOINT
    _module = __module__
    :($Logging.with_logger($Logging.current_logger()) do
        $Logging.@logmsg($level, $(esc(ex)))
    end)
end
