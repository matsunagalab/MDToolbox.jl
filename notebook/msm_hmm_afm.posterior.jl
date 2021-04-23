using Printf, ProgressMeter
@everywhere using MDToolbox

using BSON: @save, @load
@load "msm_hmm_afm.bson" afms models state qs afmize_config probe_radius sigma_noise

configs = [AfmizeConfig(10.0 * (pi / 180),
    r, 
    MDToolbox.Point2D(-250, -200), 
    MDToolbox.Point2D(250, 200), 
    MDToolbox.Point2D(6.25, 6.25), 
    MDToolbox.defaultParameters())
    for r in [10.0, 15.0, 20.0, 25.0, 30.0]]

# all frames with parallel map
@time res = @showprogress pmap(x -> getafmposterior(x, models, qs, configs), afms)

@save "msm_hmm_afm.posterior.bson" res configs

