module MDToolbox

using NetCDF
#using FFTW

# satndard library
using Printf
using Statistics
using Base.Threads

# package code goes here
export AbstractTrajectory, TrjArray, select_atom
export readdcd, readpsf
export centerofmass, decenter, superimpose, calcrmsd, calcbond
#fastCalcRMSDAndRotation!, innerproduct!

#TODO: introduce physical constatns, such as const K_B_in_kcalpermol

include("trjarray.jl")
include("fileio.jl")
include("structure.jl")

end # module
