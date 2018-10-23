module MDToolbox

using FileIO
using NetCDF
#using FFTW

# satndard library
using Printf
using Statistics
using LinearAlgebra
using Base.Threads

# package code goes here
export AbstractTrajectory, TrjArray, select_atom
export load_dcd, load_netcdf, save_netcdf, load_psf, load_pdb
export centerofmass, decenter, superimpose, superimpose_serial, get_rmsd, meanstructure, get_rmsf
export get_distance, get_angle, get_dihedral
export ksdensity, ksdensity_serial, get_pmf

#TODO: introduce physical constatns, such as const K_B_in_kcalpermol

include("trjarray.jl")
include("fileio.jl")
include("structure.jl")
include("ksdensity.jl")

end # module
