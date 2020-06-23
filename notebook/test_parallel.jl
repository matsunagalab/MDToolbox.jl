using ProgressMeter
using LinearAlgebra
#using MDToolbox

println(nprocs())

M = Matrix{Float64}[rand(1000,1000) for i = 1:50];

@time mvalues = @showprogress pmap(svdvals, M);

@everywhere f(x, k) = x .+ k

@time fvalues = progress_pmap(x -> f(x, 1.0), M);

