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
export centerofmass, superimpose, distance

#TODO: introduce physical constatns, such as K_B

include("trjarray.jl")
include("fileio.jl")
include("structure.jl")

end # module
