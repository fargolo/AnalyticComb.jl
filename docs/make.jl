using AnalyticComb
using Documenter

DocMeta.setdocmeta!(AnalyticComb, :DocTestSetup, :(using AnalyticComb); recursive=true)

makedocs(;
    modules=[AnalyticComb],
    authors="fargolo <felipe.c.argolo@protonmail.com> and contributors",
    repo="https://github.com/fargolo/AnalyticComb.jl/blob/{commit}{path}#{line}",
    sitename="AnalyticComb.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fargolo.github.io/AnalyticComb.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Introduction" => "intro.md" ,
        "Function definitions" => "index.md",
        
    ],
)

deploydocs(;
    repo="github.com/fargolo/AnalyticComb.jl",
    devbranch="main",
)
