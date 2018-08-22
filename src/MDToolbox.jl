module MDToolbox

using NetCDF
using FFTW
using Printf

# package code goes here
export AbstractTrajectory, TrjArray, to3
export readdcd
#export readdcd, writedcd, to3, calcbond

include("trjarray.jl")
include("fileio.jl")
#include("structure.jl")

end # module
