using Documenter, NumericalMethods

makedocs(
    sitename="NumericalMethods.jl",
    format = Documenter.HTML(; prettyurls = true),
    pages=[
        "index.md",
        "deriv.md",
        "interp.md",
        "min.md",
        "roots.md"
    ]
)

deploydocs(
    repo = "github.com/michaelbennett99/NumericalMethods.jl",
    devurl = "current"
)
