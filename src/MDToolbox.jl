module MDToolbox

using FileIO
using NetCDF
using NLsolve
using FFTW
using Bio3DView
using CUDA
using ProgressMeter
using SparseArrays
using MetaGraphs, LightGraphs, GraphRecipes, Plots

# satndard library
using Printf
using Statistics
using LinearAlgebra
#using SparseArrays
using Dates
using Distributed
using Base.Threads

# package code goes here
export KB_kcalpermol, KB_kjpermol
export AbstractTrajectory, TrjArray, to3, select_atom
export mdload, mdsave, readdcd, readnetcdf, writenetcdf, readpsf, writepsf, readpdb, writepdb, readcrd
export centerofmass, decenter, decenter!, superimpose, superimpose_serial, compute_rmsd, meanstructure, compute_rmsf
export rotate, rotate!, rotate_with_matrix
export compute_distance, compute_distancemap, compute_contactmap, compute_angle, compute_dihedral
export compute_qscore, compute_drms, compute_pairlist, compute_pairlist_bruteforce
export clusterkcenters, compute_cov, rsvd, pca, tica
export propagate_mcmc, propagate_md
export ksdensity, ksdensity_serial, compute_pmf
export wham, wham_iteration
export mbar
export msmplot, msmgenerate, msmcountmatrix, msmtransitionmatrix, msmforward, msmbackward, msmbaumwelch, msmviterbi, msmimpliedtime
export sp_delta_pmf, sp_design_matrix, sp_design_matrix_atom, sp_lsquares, sp_admm, sp_descent, sp_standardize!, sp_standardize, sp_cumulate_pmf, sp_cumulate_pmf_atom
export afmize, AfmizeConfig, translateafm, getafmposterior, getposterior_parallel
export Asd, readasd
export get_residues, gpu, gpu32, cpu, cpu64, logsumexp

# constants
const KB_kcalpermol = 0.0019872041 #Boltzmann constant taken from wikipedia
const KB_kjpermol = 0.0083144621 #Boltzmann constant taken from wikipedia

# codes
include("trjarray.jl")
include("fileio.jl")
include("structure.jl")
include("reduction.jl")
include("dynamics.jl")
include("ksdensity.jl")
include("wham.jl")
include("mbar.jl")
include("msm.jl")
include("sparse.jl")
include("afm.jl")
include("asd.jl")
include("docking.jl")
include("utils.jl")

end # module
