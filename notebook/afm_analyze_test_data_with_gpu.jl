using Plots, Printf, DelimitedFiles, Random, Statistics, ProgressMeter
using MDToolbox
using CUDA

using BSON: @save, @load
@load "afm_generate_test_data.bson" afms models models_template qs myseed sigma_noise nframe

configs = [AfmizeConfig(10.0 * (pi / 180),
    r, 
    MDToolbox.Point2D(-250, -200), 
    MDToolbox.Point2D(250, 200), 
    MDToolbox.Point2D(6.25, 6.25), 
    MDToolbox.defaultParameters())
    for r in [10.0, 15.0, 20.0, 25.0, 30.0]]

for iatom = 1:models_template.natom
    push!(models_template.mass, MDToolbox.defaultParameters()[models_template.resname[iatom]])
end

afms_d = cu.(afms)
models_template_d = gpu32(models_template)
qs_d = cu(qs)

# single frame
r = MDToolbox.getafmposterior_gpu(afms_d[1], models_template_d, qs_d, configs)

# all frames with parallel map
#r = pmap(x -> getafmposterior(x, models_template, qs, configs), afms[1:20])

#@save "afm_analyze_test_data.bson" r configs





