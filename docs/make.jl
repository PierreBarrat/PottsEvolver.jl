using PottsEvolver
using Documenter

DocMeta.setdocmeta!(PottsEvolver, :DocTestSetup, :(using PottsEvolver); recursive=true)

makedocs(;
    modules=[PottsEvolver],
    authors="Pierre Barrat-Charlaix",
    sitename="PottsEvolver.jl",
    format=Documenter.HTML(;
        canonical="https://pierrebarrat.github.io/PottsEvolver.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=["Quickstart" => "index.md", "Reference" => "reference.md"],
    checkdocs=:exports,
)

deploydocs(; repo="github.com/PierreBarrat/PottsEvolver.jl.git", devbranch="master")
