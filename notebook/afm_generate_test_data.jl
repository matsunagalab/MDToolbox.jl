using Plots, Printf, DelimitedFiles, Random, Statistics, ProgressMeter
using MDToolbox

qs = readdlm("data/quaternion/QUATERNION_LIST_576_Orient")
models_template = readpdb("data/t1r/cluster.pdb");
for iatom = 1:models_template.natom
    models_template.atomname[iatom] = models_template.resname[iatom]
end
MDToolbox.decenter!(models_template)

probe_radius = 20.0
afmize_config = AfmizeConfig(10.0 * (pi / 180),
    probe_radius, 
    MDToolbox.Point2D(-250, -200), 
    MDToolbox.Point2D(250, 200), 
    MDToolbox.Point2D(6.25, 6.25), 
    MDToolbox.defaultParameters())

myseed = 777
sigma_noise = 3.0
nframe =100

seed = MersenneTwister(myseed)
afms = []
models = []
@showprogress 1 "Computing..." for i in 1:nframe
    model = models_template[rand(seed, 1:models_template.nframe), :]
    q = qs[rand(seed, 1:size(qs, 1)), :]
    dx = rand(seed)*150.0 - 50.0
    dy = rand(seed)*150.0 - 50.0
    model.x .+=  dx
    model.y .+=  dy
    model = MDToolbox.rotate(model, q)
    afm = MDToolbox.afmize(model, afmize_config)
    #afm = zeros(80, 64)
    h, w = size(afm)
    afm .+= randn(seed, h, w) * sigma_noise
    push!(afms, afm)
    push!(models, model)
end

anim = @animate for i = 1:length(afms)
    heatmap(afms[i])
end
gif(anim, "afm_generate_test_data.gif", fps = 10)

using BSON: @save, @load
@save "afm_generate_test_data.bson" afms models models_template qs myseed sigma_noise nframe





