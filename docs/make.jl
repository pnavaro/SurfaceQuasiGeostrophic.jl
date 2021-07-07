using SurfaceQuasiGeostrophic
using Documenter
using Plots

ENV["GKSwstype"] = "100"

DocMeta.setdocmeta!(SurfaceQuasiGeostrophic, :DocTestSetup, :(using SurfaceQuasiGeostrophic); recursive=true)

makedocs(;
    modules=[SurfaceQuasiGeostrophic],
    authors="Mingus Team",
    repo="https://github.com/pnavaro/SurfaceQuasiGeostrophic.jl/blob/{commit}{path}#{line}",
    sitename="SurfaceQuasiGeostrophic.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pnavaro.github.io/SurfaceQuasiGeostrophic.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Buoyancy" => "buoyancy.md"
    ],
)

deploydocs(;
    branch = "gh-pages",
    repo="github.com/pnavaro/SurfaceQuasiGeostrophic.jl",
    devbranch = "main"
)
