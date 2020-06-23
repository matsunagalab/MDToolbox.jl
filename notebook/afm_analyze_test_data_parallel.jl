using Printf, ProgressMeter
@everywhere using MDToolbox

using BSON: @save, @load
@load "afm_generate_test_data.bson" afms models models_template qs myseed sigma_noise nframe

configs = [AfmizeConfig(10.0 * (pi / 180),
    r, 
    MDToolbox.Point2D(-250, -200), 
    MDToolbox.Point2D(250, 200), 
    MDToolbox.Point2D(6.25, 6.25), 
    MDToolbox.defaultParameters())
    for r in [10.0, 15.0, 20.0, 25.0, 30.0]]

# single frame
#r = getafmposterior(afms[1], models_template, qs, configs)

# all frames with parallel map
@time r = @showprogress pmap(x -> getafmposterior(x, models_template, qs, configs), afms[1:20])

@save "afm_analyze_test_data_parallel.bson" r configs





