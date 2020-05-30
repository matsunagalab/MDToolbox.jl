using Plots, Printf, DelimitedFiles, BenchmarkTools, DelimitedFiles, Serialization, Random
include("../src/MDToolbox.jl")
using .MDToolbox

seed_num, sigma, num = map(x->parse(Int,x),split(readline()))
seed = MersenneTwister(seed_num)

q_array = readdlm("data/quaternion/QUATERNION_LIST_576_Orient")
model_array = readpdb("data/t1r/cluster.pdb");
for iatom = 1:model_array.natom
    model_array.atomname[iatom] = model_array.resname[iatom]
end
MDToolbox.decenter!(model_array)
model_num = size(model_array)[1]
q_num = size(q_array)[1]

function make_rand_tada(seed, sigmam, num)
    ret = []
    for i in 1:num
        model = model_array[rand(seed, 1:model_num), :]
        radius = rand(seed, 10:30)
        quate = q_array[rand(seed, 1:q_num), :]
        dx = rand(seed, 1:25) * 0.25
        dy = rand(seed, 1:25) * 0.25
        config = AfmizeConfig(10.0 * (pi / 180), 
                      radius, 
                      MDToolbox.Point2D(-250 + dx, -200 + dy), 
                      MDToolbox.Point2D(250 + dx, 200 + dy), 
                      MDToolbox.Point2D(6.25, 6.25), 
                      MDToolbox.defaultParameters())
        model = MDToolbox.rotate(model, quate)
        afm = MDToolbox.afmize(model, config)
        h, w = size(afm)
        afm .+= randn(seed, h, w) * sigmam
        push!(ret, afm)
    end
    return ret
end

param_array = [];
for r in [15, 20, 25, 30]
  param_array = [param_array; AfmizeConfig(10.0 * (pi / 180),
                                            r, 
                                            MDToolbox.Point2D(-250, -200), 
                                            MDToolbox.Point2D(250, 200), 
                                            MDToolbox.Point2D(6.25, 6.25), 
                                            MDToolbox.defaultParameters())]
end

frames = make_rand_tada(seed, sigma, num)
result = getafmposteriors_alpha(frames, model_array, q_array, param_array)

MDToolbox.outputResults(result, "afm_test_seed_$(seed_num)_sigma_$(sigma)_num_$(num).txt")