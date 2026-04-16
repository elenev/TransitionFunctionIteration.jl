using Documenter
using TransitionFunctionIteration

makedocs(
    sitename = "TransitionFunctionIteration.jl",
    modules = [TransitionFunctionIteration],
    pages = [
        "Home" => "index.md",
        "Guide" => [
            "guide/quickstart.md",
            "guide/model.md",
            "guide/simulation.md",
        ],
        "API Reference" => "api.md",
    ]
)

# Host docs on GitHub Pages
deploydocs(
    repo = "github.com/elenev/TransitionFunctionIteration.jl.git",
)