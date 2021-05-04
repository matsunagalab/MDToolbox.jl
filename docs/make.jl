using Documenter, MDToolbox

makedocs(sitename="MDToolbox.jl",
        pages = [
        "MDToolbox.jl" => "index.md", 
        "Getting started" => [
            "getting_started01.md", 
            "getting_started02.md",
            "getting_started03.md"
        ], 
        "Examples" => [
        ], 
        "References" => [
            "fileio.md",
            "structure.md",
            "reduction.md",
            "wham.md",
            "mbar.md",
            "msm.md"],
        "License" => "license.md"
            ]
)

deploydocs(
    repo = "github.com/matsunagalab/MDToolbox.jl.git",
)
