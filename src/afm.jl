"""
afmize

simulate afm image from a given structure

example::
pdb = readpbd("1ake.pdb");
decenter!(pdb);
afm_image = afmize(pdb, (150.0, 150.0), (16, 16));
heatmap(afm_image) # using Plots
"""
function afmize_alpha(ta::TrjArray, box, npixel)::Matrix{Float64}
    rpixel = box ./ npixel
    nx = floor.(Int, (ta.x[1, :] .+ 0.5*box[1]) ./ rpixel[1]);
    ny = floor.(Int, (ta.y[1, :] .+ 0.5*box[2]) ./ rpixel[2]);
    #@show nx
    #@show ny
    afm_image = zeros(Float64, npixel)
    z_min = minimum(ta.z[1, :])
    z = ta.z[1, :] .- z_min
    for iatom in 1:ta.natom
        if (1 <= nx[iatom] <= npixel[1]) && (1 <= ny[iatom] <= npixel[2])
            if z[iatom] > afm_image[nx[iatom], ny[iatom]]
                afm_image[nx[iatom], ny[iatom]] = z[iatom]
            end
        end
    end
    return afm_image
end

function defaultParameters()
    return Dict(
     "H" => 1.20,
    "HE" => 1.40,
     "B" => 1.92,
     "C" => 1.70,
     "N" => 1.55,
     "O" => 1.52,
     "F" => 1.47,
    "NE" => 1.54,
    "NA" => 2.27,
    "MG" => 1.73,
    "AL" => 1.84,
    "SI" => 2.10,
     "P" => 1.80,
     "S" => 1.80,
    "CL" => 1.75,
    "AR" => 1.88,
     "K" => 2.75,
    "CA" => 2.31,
    "CYS" => 2.75,
    "PHE" => 3.2,
    "LEU" => 3.1,
    "TRP" => 3.4,
    "VAL" => 2.95,
    "ILE" => 3.1,
    "MET" => 3.1,
    "HIS" => 3.05,
    "HSD" => 3.05,
    "TYR" => 3.25,
    "ALA" => 2.5,
    "GLY" => 2.25,
    "PRO" => 2.8,
    "ASN" => 2.85,
    "THR" => 2.8,
    "SER" => 2.6,
    "ARG" => 3.3,
    "GLN" => 3.0,
    "ASP" => 2.8,
    "LYS" => 3.2,
    "GLU" => 2.95,
    )
end

mutable struct Point2D{T<:Real}
    x::T
    y::T
end

mutable struct Point3D{T<:Real}
    x::T
    y::T
    z::T
end

struct Probe
    x::Float64
    y::Float64
    r::Float64
    angle::Float64
end

mutable struct Sphere
    x::Float64
    y::Float64
    z::Float64
    r::Float64
end

mutable struct AfmizeConfig
    probeAngle::Float64        # Radian
    probeRadius::Float64
    range_min::Point2D
    range_max::Point2D
    resolution::Point2D
    atomRadiusDict::Dict{String, Float64}
end

function AfmizeConfig(r::Float64)
    return AfmizeConfig(10.0 * (pi / 180),
                        r,
                        Point2D(-90, -90),
                        Point2D(90, 90),
                        Point2D(15, 15),
                        defaultParameters())
end

function defaultConfig()
    return AfmizeConfig(10.0 * (pi / 180),
                        20.0,
                        Point2D(-90, -90),
                        Point2D(90, 90),
                        Point2D(15, 15),
                        defaultParameters())
end

function checkConfig(tra::TrjArray ,config::AfmizeConfig)::Union{String, Nothing}
    if (config.range_max.x - config.range_min.x) % config.resolution.x != 0
        return "resolution.x must divide range"
    end
    if (config.range_max.y - config.range_min.y) % config.resolution.y != 0
        return "resolution.y must divide range"
    end
    if config.probeAngle == 0
        return "probeAngle must be positive"
    end

    for name in tra.atomname
        if !haskey(config.atomRadiusDict, name)
            return "config dosen't know atom name $(name)"
        end
    end

    return nothing
end

# 底を0にする
function moveBottom(atoms::Array{Sphere})
    bottom = atoms[1].z - atoms[1].r
    for atom in atoms
        bottom = min(bottom, atom.z - atom.r)
    end
    for atom in atoms
        atom.z -= bottom
    end
end

# 針を球として見たとき、どこで衝突するかの計算
function calcCollisionAsSphere(probe::Probe, atom::Sphere)
    distXY = sqrt((probe.x - atom.x)^2 + (probe.y - atom.y)^2)
    dr = probe.r + atom.r
    if 0 < dr^2 - distXY^2
        return atom.z + sqrt(dr^2 - distXY^2) - probe.r
    else
        return 0.0
    end
end

# 針を円形推台として見たとき、どこで衝突するかの計算
function calcCollisionAsCircularThrusters(probe::Probe, atom::Sphere)
    distXY = sqrt((probe.x - atom.x)^2 + (probe.y - atom.y)^2)
    collisionDist = probe.r + atom.r * cos(probe.angle)
    if collisionDist < distXY
        return atom.z + atom.r * sin(probe.angle) - (distXY - collisionDist) / tan(probe.angle) - probe.r
    elseif probe.r < distXY
        return atom.z + sqrt(atom.r^2 - (distXY - probe.r)^2) - probe.r
    else
        return 0.0
    end
end

# 衝突計算をする範囲を返す
function calRectangle(min_point::Point2D, max_point::Point2D, config::AfmizeConfig)
    if max_point.x <= config.range_min.x || max_point.y <= config.range_min.y
        return nothing
    end
    if config.range_max.x <= min_point.x || config.range_max.y <= min_point.y
        return nothing
    end

    min_point = Point2D(max(min_point.x, config.range_min.x),
                        max(min_point.y, config.range_min.y))
    max_point = Point2D(min(max_point.x, config.range_max.x),
                        min(max_point.y, config.range_max.y))
    width  = floor(Int, (config.range_max.x - config.range_min.x) / config.resolution.x)
    height = floor(Int, (config.range_max.y - config.range_min.y) / config.resolution.y)
    return (max(1, floor(Int, (min_point.x - config.range_min.x + config.resolution.x - 1) / config.resolution.x)),
            max(1, floor(Int, (min_point.y - config.range_min.y + config.resolution.y - 1) / config.resolution.y)),
            min(width , floor(Int, (max_point.x - config.range_min.x + config.resolution.x - 1) / config.resolution.x)),
            min(height, floor(Int, (max_point.y - config.range_min.y + config.resolution.y - 1) / config.resolution.y)))
end


function calRectangle_sphere(atom::Sphere, config::AfmizeConfig)
    min_point = Point2D(atom.x - atom.r - config.probeRadius,
                        atom.y - atom.r - config.probeRadius)
    max_point = Point2D(atom.x + atom.r + config.probeRadius,
                        atom.y + atom.r + config.probeRadius)

    return calRectangle(min_point, max_point, config)
end

function calRectangle_circularThrusters(atom::Sphere, config::AfmizeConfig)
    dist_xy = (tan(config.probeAngle) * (atom.z + atom.r * sin(config.probeAngle) - config.probeRadius)
                + config.probeRadius + atom.r * cos(config.probeAngle))
    min_point = Point2D(atom.x - dist_xy,
                        atom.y - dist_xy)
    max_point = Point2D(atom.x + dist_xy,
                        atom.y + dist_xy)

    return calRectangle(min_point, max_point, config)
end

function afmize(tra::TrjArray, config::AfmizeConfig)
    message = checkConfig(tra, config)
    if !isnothing(message)
        println(message)
        return zeros(1, 1)
    end

    width = floor(Int, (config.range_max.x - config.range_min.x) / config.resolution.x)
    height = floor(Int, (config.range_max.y - config.range_min.y) / config.resolution.y)
    atoms = [Sphere(tra.x[i], tra.y[i], tra.z[i],
            config.atomRadiusDict[tra.atomname[i]]) for i in 1:tra.natom]
    moveBottom(atoms)

    stage = zeros(height, width)
    probes = [Probe(config.range_min.x + (w-0.5) * config.resolution.x,
                    config.range_min.y + (h-0.5) * config.resolution.y,
                    config.probeRadius, config.probeAngle)
             for h in 1:height, w in 1:width]

    # 各原子事に計算する価値のある矩形を求めて、その中で探索を行う
    for atom in atoms
        rectangle = calRectangle_sphere(atom, config)

        if isnothing(rectangle) continue end

        for h in rectangle[2]:rectangle[4], w in rectangle[1]:rectangle[3]
            probe = probes[h, w]

            stage[h, w] = max(stage[h, w], calcCollisionAsSphere(probe, atom))
        end
    end

    for atom in atoms
        rectangle = calRectangle_circularThrusters(atom, config)

        if isnothing(rectangle) continue end

        for h in rectangle[2]:rectangle[4], w in rectangle[1]:rectangle[3]
            probe = probes[h, w]

            stage[h, w] = max(stage[h, w], calcCollisionAsCircularThrusters(probe, atom))
        end
    end

    return stage
end

## -------------------------------------------------- afmize -----------------------------------------------------------

struct Probe_beta{T}
    x::AbstractArray{T}
    y::AbstractArray{T}
    r::AbstractArray{T}
    angle::AbstractArray{T}
end

mutable struct Sphere_beta{T}
    x::AbstractArray{T}
    y::AbstractArray{T}
    z::AbstractArray{T}
    r::AbstractArray{T}
end

mutable struct AfmizeConfig_beta{T}
    probeAngle::T        # Radian
    probeRadius::T
    range_min::Tuple{T,T}
    range_max::Tuple{T,T}
    resolution::Tuple{T,T}
    atomRadiusDict::Dict{String, Float64}
end

function AfmizeConfig_beta(r::Float64)
    return AfmizeConfig_beta(10.0 * (pi / 180),
                            r,
                            (-90, -90),
                            (90, 90),
                            (15, 15),
                            defaultParameters())
end

function checkConfig_beta(tra::TrjArray ,config::AfmizeConfig_beta)::Union{String, Nothing}
    if (config.range_max[1] - config.range_min[1]) % config.resolution[1] != 0
        return "resolution.x must divide range"
    end
    if (config.range_max[2] - config.range_min[2]) % config.resolution[2] != 0
        return "resolution.y must divide range"
    end
    if config.probeAngle == 0
        return "probeAngle must be positive"
    end

    for name in tra.atomname
        if !haskey(config.atomRadiusDict, name)
            return "config dosen't know atom name $(name)"
        end
    end

    return nothing
end

# 底を0にする
function moveBottom_beta(atoms::Sphere_beta)
    bottom = minimum(atoms.z .- atoms.r)
    atoms.z .-= bottom
end

# 針を球として見たとき、どこで衝突するかの計算
function calcCollisionAsSphere_beta(probes::Probe_beta, atoms::Sphere_beta, id::Int, h::Int, w::Int)
    distXY = sqrt((probes.x[h, w] - atoms.x[id])^2 + (probes.y[h, w] - atoms.y[id])^2)
    dr = probes.r[h, w] + atoms.r[id]
    if 0 < dr^2 - distXY^2
        return atoms.z[id] + sqrt(dr^2 - distXY^2) - probes.r[h, w]
    else
        return 0.0
    end
end

# 針を円形推台として見たとき、どこで衝突するかの計算
function calcCollisionAsCircularThrusters_beta(probes::Probe_beta, atoms::Sphere_beta, id::Int, h::Int, w::Int)
    distXY = sqrt((probes.x[h, w] - atoms.x[id])^2 + (probes.y[h, w] - atoms.y[id])^2)
    collisionDist = probes.r[h, w] + atoms.r[id] * cos(probes.angle[h, w])
    if collisionDist < distXY
        return atoms.z[id] + atoms.r[id] * sin(probes.angle[h, w]) - (distXY - collisionDist) / tan(probes.angle[h, w]) - probes.r[h, w]
    elseif probes.r[h, w] < distXY
        return atoms.z[id] + sqrt(atoms.r[id]^2 - (distXY - probes.r[h, w])^2) - probes.r[h, w]
    else
        return 0.0
    end
end

# 衝突計算をする範囲を返す
function calRectangle_beta(min_x, min_y, max_x, max_y, config::AfmizeConfig_beta)
    if max_x <= config.range_min[1] || max_y <= config.range_min[2]
        return nothing
    end
    if config.range_max[1] <= min_x || config.range_max[2] <= min_y
        return nothing
    end

    min_x = max(min_x, config.range_min[1])
    min_y = max(min_y, config.range_min[2])
    max_x = min(max_x, config.range_max[1])
    max_y = min(max_y, config.range_max[2])
    
    width  = floor(Int, (config.range_max[1] - config.range_min[1]) / config.resolution[1])
    height = floor(Int, (config.range_max[2] - config.range_min[2]) / config.resolution[2])
    return (max(1, floor(Int, (min_x - config.range_min[1] + config.resolution[1] - 1) / config.resolution[1])),
            max(1, floor(Int, (min_y - config.range_min[2] + config.resolution[2] - 1) / config.resolution[2])),
            min(width , floor(Int, (max_x - config.range_min[1] + config.resolution[1] - 1) / config.resolution[1])),
            min(height, floor(Int, (max_y - config.range_min[2] + config.resolution[2] - 1) / config.resolution[2])))
end


function calRectangle_sphere_beta(atoms::Sphere_beta, config::AfmizeConfig_beta, id::Int)
    return calRectangle_beta(atoms.x[id] - atoms.r[id] - config.probeRadius,
                             atoms.y[id] - atoms.r[id] - config.probeRadius,
                             atoms.x[id] + atoms.r[id] + config.probeRadius,
                             atoms.y[id] + atoms.r[id] + config.probeRadius, config)
end

function calRectangle_circularThrusters_beta(atoms::Sphere_beta, config::AfmizeConfig_beta, id::Int)
    dist_xy = (tan(config.probeAngle) * (atoms.z[id] + 
               atoms.r[id] * sin(config.probeAngle) - 
               config.probeRadius) + 
               config.probeRadius + 
               atoms.r[id] * cos(config.probeAngle))

    return calRectangle_beta(atoms.x[id] - dist_xy, 
                             atoms.y[id] - dist_xy,
                             atoms.x[id] + dist_xy, 
                             atoms.y[id] + dist_xy, config)
end

function afmize_beta(tra::TrjArray, config::AfmizeConfig_beta)
    message = checkConfig_beta(tra, config)
    if !isnothing(message)
        println(message)
        return zeros(1, 1)
    end

    atom_num = tra.natom
    width = floor(Int, (config.range_max[1] - config.range_min[1]) / config.resolution[1])
    height = floor(Int, (config.range_max[2] - config.range_min[2]) / config.resolution[2])
    atoms = Sphere_beta([tra.x[i] for i in 1:atom_num], 
                        [tra.y[i] for i in 1:atom_num], 
                        [tra.z[i] for i in 1:atom_num], 
                        [config.atomRadiusDict[tra.atomname[i]] for i in 1:atom_num])
    println(typeof(atoms))
    
    moveBottom_beta(atoms)

    stage = zeros(height, width)
    probes = Probe_beta([config.range_min[1] + (w-0.5) * config.resolution[1] for h in 1:height, w in 1:width],
                        [config.range_min[2] + (h-0.5) * config.resolution[2] for h in 1:height, w in 1:width],
                        [config.probeRadius for h in 1:height, w in 1:width],
                        [config.probeAngle  for h in 1:height, w in 1:width])

    # 各原子事に計算する価値のある矩形を求めて、その中で探索を行う
    
    for atom_id in 1:atom_num
        rectangle = calRectangle_sphere_beta(atoms, config, atom_id)

        if isnothing(rectangle) continue end

        for h in rectangle[2]:rectangle[4], w in rectangle[1]:rectangle[3]
            stage[h, w] = max(stage[h, w], calcCollisionAsSphere_beta(probes, atoms, atom_id, h, w))
        end
    end

    for atom_id in 1:atom_num
        rectangle = calRectangle_circularThrusters_beta(atoms, config, atom_id)

        if isnothing(rectangle) continue end

        for h in rectangle[2]:rectangle[4], w in rectangle[1]:rectangle[3]
            stage[h, w] = max(stage[h, w], calcCollisionAsCircularThrusters_beta(probes, atoms, atom_id, h, w))
        end
    end

    return stage
end

## -------------------------------------------------- afmize -----------------------------------------------------------

function translateafm(afm, (dx, dy))
    afm_translated = zeros(eltype(afm), size(afm))
    (nx, ny) = size(afm)
    for i in maximum([1, 1-dx]):minimum([nx, nx-dx])
        for j in maximum([1, 1-dy]):minimum([ny, ny-dy])
            afm_translated[i+dx, j+dy] = afm[i, j]
        end
    end
    afm_translated
end

function translateafm_periodic(afm, (dx, dy))
    afm_translated = zeros(eltype(afm), size(afm))
    (nx, ny) = size(afm)
    for x in 0:nx-1
        for y in 0:ny-1
            xx = (x + dx + nx) % nx
            yy = (y + dy + ny) % ny
            afm_translated[xx+1, yy+1] = afm[x+1, y+1]
        end
    end
    afm_translated
end

function direct_convolution(observed, calculated)
    H, W = size(observed)
    ret = zeros(H, W)
    for y in 1:H, x in 1:W
        for dy in 1:H, dx in 1:W
            ny = y + dy - 1
            nx = x + dx - 1
            if ny > H ny -= H end
            if nx > W nx -= W end
            ret[y, x] += observed[dy, dx] * calculated[ny, nx]
        end
    end
    return ret
end

function fft_convolution(observed, calculated)
    return real.(ifft(fft(observed).*conj.(fft(calculated))))
end

function calcLogProb(observed, calculated, convolution_func)
    npix = Float64(size(observed, 1) * size(observed, 2))
    C_o  = sum(observed)
    C_c  = sum(calculated)
    C_cc = sum(calculated.^2)
    C_oo = sum(observed.^2)
    C_oc = convolution_func(observed, calculated)

    log01 = npix .* (C_cc .* C_oo .- C_oc.^2) .+ 2.0 .* C_o .* C_oc .* C_c .- C_cc .* C_o.^2 .- C_oo .* C_c.^2
    log01[log01 .<= 0.0] .= eps(Float64)
    log02 = (npix .- 2.0) .* (npix .* C_cc .- C_c.^2)
    log02 = log02 <= 0 ? eps(Float64) : log02
    logprob = 0.5 .* (3.0 .- npix) .* log.(log01) .+ (0.5 .* npix .- 2.0) .* log.(log02)

    return logprob
end

function getafmposterior(afm::Matrix{Float64}, model_array::TrjArray, q_array::Matrix{Float64}, param_array)
    println("in getafmposterior")
    imax_model = 0
    imax_q = 0
    best_param = param_array[1]
    best_translate = [0, 0]
    best_afm = similar(afm)
    best_posterior = -Inf

    observed = afm
    npix = Float64(size(observed, 1) * size(observed, 2))

    decenter!(model_array)

    Threads.@threads for imodel in 1:size(model_array, 1)
        @show imodel
        model = model_array[imodel, :]
        for iq in 1:size(q_array, 1)
            q = q_array[iq, :]
            model_rotated = MDToolbox.rotate(model, q)
            for iparam in 1:length(param_array)
                param = param_array[iparam]
                calculated = afmize(model_rotated, param)
                C_o  = sum(observed)
                C_c  = sum(calculated)
                #@btime C_oc = sum(observed_translated .* calculated)
                C_oc = sum(observed .* calculated)
                C_cc = sum(calculated.^2)
                C_oo = sum(observed.^2)
                #@btime C_oc_dxdy = real.(ifftshift(ifft(fft(observed_translated).*conj.(fft(calculated)))))
                C_oc_dxdy = real.(ifftshift(ifft(fft(observed).*conj.(fft(calculated)))))
                log01 = npix .* (C_cc .* C_oo .- C_oc_dxdy.^2) .+ 2.0 .* C_o .* C_oc_dxdy .* C_c .- C_cc .* C_o.^2 .- C_oo .* C_c.^2
                log01[log01 .<= 0.0] .= eps(Float64)
                log02 = (npix .- 2.0) .* (npix .* C_cc .- C_c.^2)
                log02 = log02 <= 0 ? eps(Float64) : log02
                logprob = 0.5 .* (3.0 .- npix) .* log.(log01) .+ (0.5 .* npix .- 2.0) .* log.(log02)
                maximum_logprob = maximum(logprob)
                if best_posterior < maximum_logprob
                    best_posterior = maximum_logprob
                    imax_model = imodel
                    imax_q = iq
                    best_param = param
                    x_center = ceil(Int32, (size(observed,1)/2)+1.0)
                    y_center = ceil(Int32, (size(observed,2)/2)+1.0)
                    dx_estimated = argmax(logprob)[1] - x_center
                    dy_estimated = argmax(logprob)[2] - y_center
                    best_translate = (dx_estimated, dy_estimated)
                    best_afm = translateafm(calculated, best_translate)
                end
            end
        end
    end

    return imax_model, imax_q, best_param, best_translate, best_afm, best_posterior
end

mutable struct posteriorResult
    posteriors
    each_best
    afm_results
    posterior_results
    best_translate
    best_posterior
    best_model
    best_quate
    best_model_rotated
    best_param
    best_afm
end

function getafmposteriors_alpha(afm_frames, model_array, quate_array, param_array, opt = "")
    frame_num = size(afm_frames)[1]
    model_num = size(model_array)[1]
    quate_num = size(quate_array)[1]
    param_num = size(param_array)[1]

    convolution_func = fft_convolution
    if opt == "direct"
        convolution_func = direct_convolution
    end

    println("func is $convolution_func")

    results = []
    for i in 1:frame_num
        # 各モデルについて、一つの角度と半径ごとのlogprobの最高値を保持しておく
        posteriors = zeros(model_num, quate_num * param_num)
        # 各モデルについて、最も良い値をだした時のafm画像
        afm_results = [zeros(size(afm_frames[1])) for i in 1:model_num]
        each_best = ones(model_num) .* -Inf

        best_posterior = -Inf
        best_model = model_array[1]
        best_quate = quate_array[1, :]
        best_model_rotated = model_array[1]
        best_param = param_array[1]
        best_afm = zeros(size(afm_frames[1]))
        posterior_results = zeros(model_num)
        push!(results, posteriorResult(posteriors, 
                                        each_best, 
                                        afm_results, 
                                        posterior_results,
                                        (0, 0), 
                                        best_posterior, 
                                        best_model, 
                                        best_quate,
                                        best_model_rotated, 
                                        best_param, 
                                        best_afm))
    end

    for (model_id, model) in zip(1:model_num, model_array)
        @show model_id
        for quate_id in 1:quate_num
            rotated_model = MDToolbox.rotate(model, quate_array[quate_id, :])
            for param_id in 1:param_num
                cal_frame = afmize(rotated_model, param_array[param_id])
                for frame_id in 1:frame_num
                    prob_mat = calcLogProb(afm_frames[frame_id], cal_frame, convolution_func)
                    max_prob = maximum(prob_mat)
                    id = (quate_id - 1) * param_num + param_id
                    results[frame_id].posteriors[model_id, id] = max_prob

                    if results[frame_id].each_best[model_id] < max_prob
                        results[frame_id].each_best[model_id] = max_prob
                        results[frame_id].afm_results[model_id] = cal_frame
                    end

                    if results[frame_id].best_posterior < max_prob
                        results[frame_id].best_posterior = max_prob
                        results[frame_id].best_model = model
                        results[frame_id].best_quate = quate_array[quate_id, :]
                        results[frame_id].best_model_rotated = rotated_model
                        results[frame_id].best_param = param_array[param_id]
                        results[frame_id].best_translate = Tuple(argmax(prob_mat))
                        results[frame_id].best_afm = translateafm_periodic(cal_frame, results[frame_id].best_translate)
                    end
                end
            end
        end
    end

    for frame_id in 1:frame_num
        results[frame_id].posteriors .-= maximum(results[frame_id].posteriors)
        results[frame_id].posteriors = exp.(results[frame_id].posteriors)
        results[frame_id].posterior_results = sum(results[frame_id].posteriors, dims = 2)
        results[frame_id].posterior_results ./= sum(results[frame_id].posterior_results)
    end

    return results
end
