module MDToolbox

using FileIO
using NetCDF
using NLsolve
using FFTW
using Bio3DView
using CuArrays

# satndard library
using Printf
using Statistics
using LinearAlgebra
using Dates
using Base.Threads

# package code goes here
export KB_kcalpermol, KB_kjpermol
export AbstractTrajectory, TrjArray, select_atom
export readdcd, readnetcdf, writenetcdf, readpsf, writepsf, readpdb, writepdb, readcrd
export centerofmass, decenter, decenter!, superimpose, superimpose_serial, getrmsd, meanstructure, getrmsf, rotate, rotate!
export getdistance, getangle, getdihedral
export getmsd
export propagate_mcmc, propagate_md
export ksdensity, ksdensity_serial, getpmf
export wham, wham_iteration
export mbar
export sp_delta_pmf, sp_design_matrix, sp_design_matrix_atom, sp_lsquares, sp_admm, sp_descent, sp_standardize!, sp_standardize, sp_cumulate_pmf, sp_cumulate_pmf_atom
export afmize, AfmizeConfig, translateafm, getafmposterior
export Asd, readasd
export viewstruc, gpu, gpu32, cpu, cpu64

# constants
const KB_kcalpermol = 0.0019872041 #Boltzmann constant taken from wikipedia
const KB_kjpermol = 0.0083144621 #Boltzmann constant taken from wikipedia

# codes
include("trjarray.jl")
include("fileio.jl")
include("structure.jl")
include("dynamics.jl")
include("ksdensity.jl")
include("wham.jl")
include("mbar.jl")
include("sparse.jl")
include("afm.jl")
include("asd.jl")
include("utils.jl")

end # module
