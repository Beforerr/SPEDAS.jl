using SPEDAS
using SPEDAS.SpaceDataModel
using SPEDAS.GeoCotrans
using Documenter
using DocumenterCitations
# using DemoCards

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

# demopage, postprocess_cb, gallery_assets = makedemos("gallery")

md_filter(x) = filter(endswith(".md"), x)
list_pages(dir) = ["$dir/$f" for f in readdir(joinpath(@__DIR__, "src", dir))]

makedocs(
    sitename = "SPEDAS.jl",
    pages = [
        "Home" => "index.md",
        "Tutorials" => list_pages("tutorials"),
        "Examples" => [
            "examples/index.md",
            "examples/speasy.md",
            "examples/tplot.md",
            "examples/interactive.md",
            "examples/interactive_speasy.md"
        ],
        "Explanation" => [
            "explanations/data.md",
            "explanations/data_model.md",
            "explanations/tplot.md",
            "explanations/coords.md",
            "explanations/multispacecraft.md",
            "explanations/resampling.md",
            "explanations/waves.md",
            "explanations/analysis.md"
        ],
        "Observatories" => list_pages("observatory"),
        "Validation" => list_pages("validation"),
        "API" => "api.md"
    ],
    format = Documenter.HTML(size_threshold = nothing),
    modules = [SPEDAS, SPEDAS.SpaceDataModel, GeoCotrans],
    warnonly = Documenter.except(:doctest),
    plugins = [bib],
    doctest = true
)

# postprocess_cb() # redirect url for DemoCards generated files

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(repo = "github.com/Beforerr/SPEDAS.jl", push_preview = true)
