module MDToolbox

using FileIO
using NetCDF
using NLsolve
#using FFTW

# satndard library
using Printf
using Statistics
using LinearAlgebra
using Base.Threads

# package code goes here
export KB_kcalpermol, KB_kjpermol
export AbstractTrajectory, TrjArray, select_atom
export readdcd, readnetcdf, writenetcdf, readpsf, writepsf, readpdb, writepdb
export centerofmass, decenter, superimpose, superimpose_serial, getrmsd, meanstructure, getrmsf
export getdistance, getangle, getdihedral
export ksdensity, ksdensity_serial, getpmf
export wham, wham_iteration

# constants
const KB_kcalpermol = 0.0019872041 #Boltzmann constant taken from wikipedia
const KB_kjpermol = 0.0083144621 #Boltzmann constant taken from wikipedia

# codes
include("trjarray.jl")
include("fileio.jl")
include("structure.jl")
include("ksdensity.jl")
include("wham.jl")

end # module
