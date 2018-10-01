module MDToolbox

using NetCDF
#using FFTW

# satndard library
using Printf
using Statistics
using LinearAlgebra
using Base.Threads

# package code goes here
export AbstractTrajectory, TrjArray, select_atom
export readdcd, readnetcdf, writenetcdf, readpsf, readpdb
export centerofmass, decenter, superimpose, superimpose_serial, calcrmsd, meanstructure, calcrmsf
export calcbond, calcangle, calcdihedral
export ksdensity, ksdensity_serial, calcpmf

#TODO: introduce physical constatns, such as const K_B_in_kcalpermol

include("trjarray.jl")
include("fileio.jl")
include("structure.jl")
include("ksdensity.jl")

end # module
