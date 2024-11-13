using Documenter: Documenter, makedocs, deploydocs
using TIFFDatasets: TIFFDatasets
import Literate

Literate.markdown(
    joinpath(@__DIR__, "..", "examples", "flood_example.jl"),
    joinpath(@__DIR__, "src", "examples"),
    execute = true,
    documenter = true,
    # We add the credit to Literate.jl the footer
    credit = false,
)

makedocs(;
    modules=[TIFFDatasets],
    repo="https://github.com/Alexander-Barth/TIFFDatasets.jl/blob/{commit}{path}#{line}",
    sitename="TIFFDatasets.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://alexander-barth.github.io/TIFFDatasets.jl",
        assets=String[],
        footer = "Powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl), [Literate.jl](https://github.com/fredrikekre/Literate.jl) and the [Julia Programming Language](https://julialang.org/)"
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => [
            "MODIS/Aqua+Terra Global Flood Product" => "examples/flood_example.md",
            ]
    ],
)

deploydocs(;
    repo="github.com/Alexander-Barth/TIFFDatasets.jl",
)
