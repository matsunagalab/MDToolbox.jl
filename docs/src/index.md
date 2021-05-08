# MDToolbox.jl

## What is MDToolbox.jl?

MDToolbox.jl is a Julia package for molecular dynamics (MD) trajectories analysis and modeling of biomolecules. The package contains functions covering the following types of scientific computations:

- I/O for trajectory, coordinate, and topology files used for MD simulations
- Unique type (TrjArray) for storing and processing MD trajectory data
- Flexible atom selections
- Various structural calculations
- Least-squares fitting of structures
- Potential of mean force (PMF) or free energy profile from unbiased MD trajectories
- Statistical estimates (WHAM and MBAR methods) from biased MD trajectories
- Dimensional reductions (Principal Component Analysis, and others)
- Markov state model, etc. 

Some functions have been ported into this package from our old MATLAB toolbox https://github.com/ymatsunaga/mdtoolbox

## Citations

In preparation

## Contents

```@contents
Pages = ["installation.md",
        "getting_started01.md", 
        "getting_started02.md", 
        "getting_started03.md", 
        "superimpose_rmsd.md",
        "free_energy_surface.md",
        "wham_example.md",
        "fileio.md", 
        "structure.md", 
        "reduction.md", 
        "wham.md", 
        "mbar.md", 
        "msm.md", 
        "workflow_for_developers.md", 
        "license.md"]
Depth = 2
```


