using Documenter: Documenter, makedocs, deploydocs
using TIFFDatasets: TIFFDatasets

makedocs(;
    modules=[TIFFDatasets],
    repo="https://github.com/Alexander-Barth/TIFFDatasets.jl/blob/{commit}{path}#{line}",
    sitename="TIFFDatasets.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://alexander-barth.github.io/TIFFDatasets.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Alexander-Barth/TIFFDatasets.jl",
)
