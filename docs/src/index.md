# MDToolbox.jl

MDToolbox.jl is a Julia package for molecular dynamics (MD) trajectories analysis and modeling of biomolecules. The package contains functions covering the following types of scientific computations:

- I/O for trajectory, coordinate, and topology files used for MD simulations
- Unique type (TrjArray) for MD trajecory data
- Flexible atom selections
- Least-squares fitting of structures
- Potential of mean force (PMF) or free energy profile from unbiased MD trajectories
- Statistical estimates (WHAM and MBAR methods) from biased MD trajectories
- Dimensional reductions (Principal Component Analysis, and others)
- Markov state model, etc. 

Some functions are ported into Julia from our old MATLAB toolbox https://github.com/ymatsunaga/mdtoolbox

## Installation

Julia version 1.0 or later is required.

```julia
]add https://github.com/matsunagalab/MDToolbox.jl.git
julia> using MDToolbox
```

## Citations

In preparation

## License



## Contents

```@contents
Pages = ["fileio.md", "structure.md", "reduction.md", "wham.md", "mbar.md", "msm.md", "license.md"]
Depth = 2
```


