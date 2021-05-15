using Documenter, MDToolbox

makedocs(sitename="MDToolbox.jl",
        pages = [
        "MDToolbox.jl" => "index.md", 
        "Installation" => "installation.md", 
        "Getting started" => [
            "getting_started01.md", 
            "getting_started02.md",
            "getting_started03.md"
        ], 
        "Examples" => [
            "superimpose_rmsd.md",
            "free_energy_surface.md",
            "wham_example.md"
        ], 
        "References" => [
            "fileio.md",
            "structure.md",
            "reduction.md",
            "wham.md",
            "mbar.md",
            "msm.md"],
        "Workflow for developers" => "workflow_for_developers.md"
            ]
)

deploydocs(
    repo = "github.com/matsunagalab/MDToolbox.jl.git",
)
