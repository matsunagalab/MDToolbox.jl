# MDToolbox.jl

MDToolbox.jl is a Julia package for molecular dynamics (MD) trajectories analysis and modeling of biomolecules. It consists of a collection of functions covering the following types of scientific computations:

- I/O for trajectory, coordinate, and topology files used for MD simulation
- New complex type for containing and operating molecular dynamics trajecory data
- Flexible atom selections
- Least-squares fitting of structures
- Potential mean force (PMF) or free energy profile from scattered data
- Statistical estimates (WHAM and MBAR methods) from biased data
- Dimensional reductions (Principal Component Analysis, and others)
- Other utility functions

Unique point of MDToolbox is that it is a pure Julia package; all codes are written with Julia. 
This make developers easy to maintain the package codes. 


```@contents
```

## Reduction

```@docs
pca(X::AbstractMatrix; k=nothing)
```

## Index

```@index
```