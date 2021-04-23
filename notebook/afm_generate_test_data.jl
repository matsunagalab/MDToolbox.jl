using Plots, Printf, DelimitedFiles, Random, Statistics, ProgressMeter
using MDToolbox

qs = readdlm("data/quaternion/QUATERNION_LIST_576_Orient")
models = readpdb("data/t1r/cluster.pdb");
for iatom = 1:models.natom
    models.atomname[iatom] = models.resname[iatom]
end
MDToolbox.decenter!(models)

probe_radius = 20.0
param = AfmizeConfig(10.0 * (pi / 180),
    probe_radius, 
    MDToolbox.Point2D(-250, -200), 
    MDToolbox.Point2D(250, 200), 
    MDToolbox.Point2D(6.25, 6.25), 
    MDToolbox.defaultParameters())

myseed = 777
sigma_noise = 3.0
nframe = 1000

seed = MersenneTwister(myseed)
afms = []
imodel_array = []
iq_array = []
dxdy_array = []
@showprogress 1 "Computing..." for i in 1:nframe
    imodel = rand(seed, 1:models.nframe)
    model = models[imodel, :]
    iq = rand(seed, 1:size(qs, 1))
    q = qs[iq, :]
    model = MDToolbox.rotate(model, q)
    dx = rand(seed)*150.0 - 50.0
    dy = rand(seed)*150.0 - 50.0
    model.x .+=  dx
    model.y .+=  dy
    afm = MDToolbox.afmize(model, param)
    h, w = size(afm)
    afm .+= randn(seed, h, w) * sigma_noise
    push!(afms, afm)
    push!(imodel_array, imodel)
    push!(iq_array, iq)
    push!(dxdy_array, (dx, dy))
end

anim = @animate for i = 1:length(afms)
    heatmap(afms[i])
end
gif(anim, "afm_generate_test_data.gif", fps = 10)

using BSON: @save, @load
@save "afm_generate_test_data.bson" models afms qs param imodel_array iq_array dxdy_array sigma_noise





