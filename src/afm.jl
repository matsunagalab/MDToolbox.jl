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

function defaultConfig()
    return AfmizeConfig(10.0 * (pi / 180),
                        10.0,
                        Point2D(-100, -100),
                        Point2D(100, 100),
                        Point2D(10, 10),
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

"""
高速化前

function afmize_beta(tra::TrjArray, config::AfmizeConfig)
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

    for h in 1:height, w in 1:width
        probe = Probe(config.range_min.x + (w-0.5) * config.resolution.x,
                      config.range_min.y + (h-0.5) * config.resolution.y,
                      config.probeRadius, config.probeAngle)
        for atom in atoms
            stage[h, w] = max(stage[h, w], calcCollisionAsSphere(probe, atom))
            stage[h, w] = max(stage[h, w], calcCollisionAsCircularThrusters(probe, atom))
        end
    end

    return stage
end
"""

function getafmposterior(afm_frame, model_list = nothing, quaternion_list = nothing, radius_list = nothing)
    for model = 1:3
        for rotation = 1:3
            for probe_radius = 1:3
                for dx = 1:3
                    nothing
                end
            end
        end
    end
end
