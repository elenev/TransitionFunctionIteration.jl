using Documenter
using TransitionFunctionIteration

ENV["GKSwstype"] = "100"   # headless GR for the plots of the @example blocks

makedocs(
    sitename = "TransitionFunctionIteration.jl",
    modules = [TransitionFunctionIteration],
    pages = [
        "Home" => "index.md",
        "Guide" => [
            "guide/quickstart.md",
            "guide/model.md",
            "guide/simulation.md",
            "guide/moments.md",
            "guide/tables.md",
            "guide/outputvault.md",
        ],
        "API Reference" => "api.md",
    ]
)

# Host docs on GitHub Pages
deploydocs(
    repo = "github.com/elenev/TransitionFunctionIteration.jl.git",
)