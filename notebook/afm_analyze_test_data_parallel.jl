using Printf, ProgressMeter
@everywhere using MDToolbox

using BSON: @save, @load
@load "afm_generate_test_data.bson" models afms qs param imodel_array iq_array dxdy_array sigma_noise

params = [AfmizeConfig(10.0 * (pi / 180),
    r, 
    MDToolbox.Point2D(-250, -200), 
    MDToolbox.Point2D(250, 200), 
    MDToolbox.Point2D(6.25, 6.25), 
    MDToolbox.defaultParameters())
    for r in [10.0, 15.0, 20.0, 25.0, 30.0]]

# single model
# p = MDToolbox.logprob_eachmodel(models[1], afms, qs, params)

# all models with parallel map
@time r = getposterior_parallel(models, afms, qs, params)

@save "afm_analyze_test_data_parallel.bson" params r





