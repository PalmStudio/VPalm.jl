using VPalm
using Documenter

DocMeta.setdocmeta!(VPalm, :DocTestSetup, :(using VPalm); recursive=true)

makedocs(;
    modules=[VPalm],
    authors="Rémi Vezy <VEZY@users.noreply.github.com> and contributors",
    sitename="VPalm.jl",
    format=Documenter.HTML(;
        canonical="https://PalmStudio.github.io/VPalm.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/PalmStudio/VPalm.jl",
    devbranch="main",
)
