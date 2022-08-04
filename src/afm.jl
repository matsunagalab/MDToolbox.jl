"""
Blind tip reconstruction from AFM image
"""

function compute_xc_yc(tip)
    tip_xsiz, tip_ysiz = size(tip)
    xc = round(Int, tip_xsiz/2, RoundNearestTiesUp) - 1
    yc = round(Int, tip_ysiz/2, RoundNearestTiesUp) - 1
    return xc, yc
end

function ireflect(surface)
    surface2 = - surface[end:-1:1, end:-1:1]
    return surface2
end

function icmap(image, surface, tip; thresh=0.1)
    xc, yc = compute_xc_yc(tip)
    surf_xsiz, surf_ysiz = size(surface)
    tip_xsiz, tip_ysiz = size(tip)
    r = zeros(Int64, size(surface))
    x = y = 0
    num_contact = 0
    for i = 1:surf_xsiz
        for j = 1:surf_ysiz
            num_contact = 0
            pxmin = max(xc-tip_xsiz, 1-i)
            pymin = max(yc-tip_ysiz, 1-j)
            pxmax = min(xc-1, surf_xsiz-i)
            pymax = min(yc-1, surf_ysiz-j)
            for px = pxmin:pxmax
                for py = pymin:pymax
                    temp1 = image[i, j] - tip[xc-px, yc-py]
                    temp2 = surface[i+px, j+py]
                    if abs(temp1 - temp2) < thresh
                        num_contact += 1
                        x = i+px
                        y = j+py
                    end
                end
            end
            if num_contact == 1
                r[x, y] = 1
            end
        end
    end
    return r
end

function idilation(surface, tip)
    xc, yc = compute_xc_yc(tip)
    surf_xsiz, surf_ysiz = size(surface)
    tip_xsiz, tip_ysiz = size(tip)
    r = zeros(eltype(surface), size(surface))
    for i = 1:surf_xsiz
        for j = 1:surf_ysiz
            pxmin = max(i-surf_xsiz, -xc+1)
            pymin = max(j-surf_ysiz, -yc+1)
            pxmax = min(i-1, -xc+tip_xsiz)
            pymax = min(j-1, -yc+tip_ysiz)
            dil_max = surface[i-pxmin, j-pymin] + tip[xc+pxmin, yc+pymin]
            for px = pxmin:pxmax
                for py = pymin:pymax
                    temp = surface[i-px, j-py] + tip[xc+px, yc+py]
                    dil_max = max(temp, dil_max)
                end
            end
            r[i, j] = dil_max
        end
    end
    return r
end

function ChainRulesCore.rrule(::typeof(idilation), surface::AbstractArray, tip::AbstractArray)
    tip_xsiz, tip_ysiz = size(tip)
    xc = round(Int, tip_xsiz/2, RoundNearestTiesUp) - 1
    yc = round(Int, tip_ysiz/2, RoundNearestTiesUp) - 1
    surf_xsiz, surf_ysiz = size(surface)
    image = similar(surface)
    MAX_i = zeros(Int, size(surface))
    MAX_j = zeros(Int, size(surface))
    MAX_u = zeros(Int, size(surface))
    MAX_v = zeros(Int, size(surface))
    for i = 1:surf_xsiz
        for j = 1:surf_ysiz
            pxmin = max(i-surf_xsiz, -xc+1)
            pymin = max(j-surf_ysiz, -yc+1)
            pxmax = min(i-1, -xc+tip_xsiz)
            pymax = min(j-1, -yc+tip_ysiz)
            dil_max = surface[i-pxmin, j-pymin] + tip[xc+pxmin, yc+pymin]
            i_max = i-pxmin
            j_max = j-pymin
            u_max = xc+pxmin
            v_max = yc+pymin
            for px = pxmin:pxmax
                for py = pymin:pymax
                    temp = surface[i-px, j-py] + tip[xc+px, yc+py]
                    if temp > dil_max
                        dil_max = temp
                        i_max = i-px
                        j_max = j-py
                        u_max = xc+px
                        v_max = yc+py
                    end
                end
            end
	    if (xc+pxmin) == 1 && (xc+pxmax) == tip_xsiz && (yc+pymin) == 1 && (yc+pymax) == tip_ysiz
              MAX_i[i, j] = i_max
              MAX_j[i, j] = j_max
              MAX_u[i, j] = u_max
              MAX_v[i, j] = v_max
	    end
            image[i, j] = dil_max
        end
    end

    function idilation_pullback(dI)
        dS = zeros(eltype(surface), size(surface))
        dP = zeros(eltype(tip), size(tip))
        for i = 1:surf_xsiz
            for j = 1:surf_ysiz
    	        if MAX_i[i,j] != 0
                    dS[MAX_i[i,j], MAX_j[i,j]] += dI[i,j]
                    dP[MAX_u[i,j], MAX_v[i,j]] += dI[i,j]
        		end
            end
        end
        #dP .= dP .- dP[xc, yc]
        #dP .= max.(0.0, dP)
        #dP[xc, yc] = 0.0
        return NoTangent(), dS, dP
    end
    return image, idilation_pullback
end

function ierosion(image, tip)
    xc, yc = compute_xc_yc(tip)
    im_xsiz, im_ysiz = size(image)
    tip_xsiz, tip_ysiz = size(tip)
    r = similar(image)
    for i = 1:im_xsiz
        for j = 1:im_ysiz
            pxmin = max(-i+1, -xc+1)
            pymin = max(-j+1, -yc+1)
            pxmax = min(-i+im_xsiz, -xc+tip_xsiz)
            pymax = min(-j+im_ysiz, -yc+tip_ysiz)
            eros_min = image[i+pxmin, j+pymin] - tip[xc+pxmin, yc+pymin]
            for px = pxmin:pxmax
                for py = pymin:pymax
                    temp = image[i+px, j+py] - tip[xc+px, yc+py]
                    eros_min = min(temp, eros_min)
                end
            end
            r[i, j] = eros_min
        end
    end
    return r
end

function ChainRulesCore.rrule(::typeof(ierosion), image::AbstractArray, tip::AbstractArray)
    tip_xsiz, tip_ysiz = size(tip)
    xc = round(Int, tip_xsiz/2, RoundNearestTiesUp) - 1
    yc = round(Int, tip_ysiz/2, RoundNearestTiesUp) - 1
    im_xsiz, im_ysiz = size(image)
    surface = similar(image)
    MIN_i = zeros(Int, size(image))
    MIN_j = zeros(Int, size(image))
    MIN_u = zeros(Int, size(image))
    MIN_v = zeros(Int, size(image))
    for i = 1:im_xsiz
        for j = 1:im_ysiz
            pxmin = max(-i+1, -xc+1)
            pymin = max(-j+1, -yc+1)
            pxmax = min(-i+im_xsiz, -xc+tip_xsiz)
            pymax = min(-j+im_ysiz, -yc+tip_ysiz)
            eros_min = image[i+pxmin, j+pymin] - tip[xc+pxmin, yc+pymin]
            i_min = i+pxmin
            j_min = j+pymin
            u_min = xc+pxmin
            v_min = yc+pymin
            for px = pxmin:pxmax
                for py = pymin:pymax
                    temp = image[i+px, j+py] - tip[xc+px, yc+py]
                    if temp < eros_min
                        eros_min = temp
                        i_min = i+px
                        j_min = j+py
                        u_min = xc+px
                        v_min = yc+py
                    end
                end
            end
	    if (xc+pxmin) == 1 && (xc+pxmax) == tip_xsiz && (yc+pymin) == 1 && (yc+pymax) == tip_ysiz
              MIN_i[i, j] = i_min
              MIN_j[i, j] = j_min
              MIN_u[i, j] = u_min
              MIN_v[i, j] = v_min
	    end
            surface[i, j] = eros_min
        end
    end

    function ierosion_pullback(dS)
        dI = zeros(eltype(image), size(image))
        dP = zeros(eltype(tip), size(tip))
        for i = 1:im_xsiz
            for j = 1:im_ysiz
    	        if MIN_i[i,j] != 0
                      dI[MIN_i[i,j], MIN_j[i,j]] += dS[i,j]
                      dP[MIN_u[i,j], MIN_v[i,j]] -= dS[i,j]
        		end
            end
        end
        #dP .= dP .- dP[xc, yc]
        #dP .= max.(0.0, dP)
        #dP[xc, yc] = 0.0
        return NoTangent(), dI, dP
    end
    return surface, ierosion_pullback
end

function iopen(image, tip)
    r = ierosion(image, tip)
    r = idilation(r, tip)
    return r
end

function itip_estimate_point!(tip0, image, ixp, jxp; thresh=0.0, iregularize=false)
    xc, yc = compute_xc_yc(tip0)
    im_xsiz, im_ysiz = size(image)
    tip_xsiz, tip_ysiz = size(tip0)
    interior = ((tip_xsiz-1) <= ixp) & (ixp <= (im_xsiz-tip_xsiz)) & ((tip_ysiz-1) <= jxp) & (jxp <= (im_ysiz-tip_ysiz))
    count = 0
    if interior
        for ix = 0:(tip_xsiz-1)
            for jx = 0:(tip_ysiz-1)
                if ((ix+1) == xc) & ((jx+1) == yc)
                    continue
                end
                imagep = image[ixp+1, jxp+1]
                dil = - eltype(image)(Inf)
                for id = 0:(tip_xsiz-1)
                    for jd = 0:(tip_ysiz-1)
                        #if (imagep - image[ixp+xc-id, jxp+yc-jd]) > tip0[id+1, jd+1]
                        #    continue
                        #end
                        #if ((imagep - image[ixp+xc-id, jxp+yc-jd]) > 0.0) | (xc == (id+1) & yc == (jd+1))
                        #    continue
                        #end
                        if iregularize
                            if ((imagep - image[ixp+xc-id, jxp+yc-jd] + thresh) > tip0[id+1, jd+1])
                                continue
                            end
                        else
                            #if ((imagep - image[ixp+xc-id, jxp+yc-jd]) > 0.0)
                            if ((imagep - image[ixp+xc-id, jxp+yc-jd]) > tip0[id+1, jd+1])
                                continue
                            end
                        end
                        #temp = image[ix+ixp-id+1, jx+jxp-jd+1] + tip0[id+1, jd+1] - imagep + thresh
                        temp = image[ix+ixp-id+1, jx+jxp-jd+1] - imagep + thresh
                        dil = max(dil, temp)
                    end
                end
                if isinf(-dil)
                    continue
                end
                if dil < tip0[ix+1, jx+1]
                    tip0[ix+1, jx+1] = dil
                    count += 1
                end
            end
        end
    end
    return count
end

function itip_estimate_iter!(tip0, image::Matrix{T}; thresh=0.0, iregularize=false) where {T <: Number}
    xc, yc = compute_xc_yc(tip0)
    im_xsiz, im_ysiz = size(image)
    tip_xsiz, tip_ysiz = size(tip0)
    count = 0
    open_image = iopen(image, tip0)
    for ixp = (tip_xsiz-xc):(im_xsiz-xc)
        for jxp = (tip_ysiz-yc):(im_ysiz-yc)
            if (image[ixp+1, jxp+1] - open_image[ixp+1, jxp+1]) > thresh
                c = itip_estimate_point!(tip0, image, ixp, jxp, thresh=thresh, iregularize=iregularize)
                if c > 0
                    count += 1
                end
            end
        end
    end
    return count
end

function itip_estimate_single_image!(tip0, image::Matrix{T}; thresh=0.0, iregularize=false) where {T <: Number}
    iter = 0
    count = 1
    count_cumsum = 0
    while count > 0
        iter += 1
        count = itip_estimate_iter!(tip0, image, thresh=thresh, iregularize=iregularize)
        count_cumsum += count
    end
    return count_cumsum
end

function itip_estimate!(tip0, images::AbstractVector; thresh=0.0, iregularize=false)
    nframe = length(images)
    count_cumsum2 = 0
    for iframe = 1:nframe
        count_cumsum = itip_estimate_single_image!(tip0, images[iframe], thresh=thresh, iregularize=iregularize)
        count_cumsum2 += count_cumsum
    end
    @printf "Processed %d image\n" nframe
    @printf "%d refinements \n" count_cumsum2
    return
end

function itip_estimate!(tip0, image::AbstractMatrix; thresh=0.0, iregularize=false)
    itip_estimate!(tip0, [image], thresh=thresh, iregularize=iregularize)
    return
end

function _logsumexp(x; tau=1)
    if tau > 0.0
        c = maximum(x)
    else
        c = minimum(x)
    end
    return c + tau * log(sum(exp.((x .- c) ./ tau)))
end

function _softmax(x; tau=1)
    m = _logsumexp(x; tau=tau)
    if tau > 0.0
        return exp.((x .- m) ./ tau)
    else
        return exp.((x .- m) ./ tau)
    end
end

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
    nx = floor.(Int, (ta.xyz[1, 1:3:end] .+ 0.5*box[1]) ./ rpixel[1]);
    ny = floor.(Int, (ta.xyz[1, 2:3:end] .+ 0.5*box[2]) ./ rpixel[2]);
    #@show nx
    #@show ny
    afm_image = zeros(Float64, npixel)
    z_min = minimum(ta.xyz[1, 3:3:end])
    z = ta.xyz[1, 3:3:end] .- z_min
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
     "CA" => 1.70,
     "CB" => 1.70,
     "CG" => 1.70,
     "CG1" => 1.70,
     "CG2" => 1.70,
     "CG3" => 1.70,
     "CD" => 1.70,
     "CD1" => 1.70,
     "CD2" => 1.70,
     "CD3" => 1.70,
     "CZ" => 1.70,
     "CZ1" => 1.70,
     "CZ2" => 1.70,
     "CZ3" => 1.70,
     "CD" => 1.70,
     "CD1" => 1.70,
     "CD2" => 1.70,
     "CD3" => 1.70,
     "CE" => 1.70,
     "CE1" => 1.70,
     "CE2" => 1.70,
     "CE3" => 1.70,
     "CH" => 1.70,
     "CH1" => 1.70,
     "CH2" => 1.70,
     "CH3" => 1.70,
     "N" => 1.55,
     "NE" => 1.55,
     "NZ" => 1.55,
     "ND1" => 1.55,
     "ND2" => 1.55,
     "NE1" => 1.55,
     "NE2" => 1.55,
     "NH1" => 1.55,
     "NH2" => 1.55,
      "O" => 1.52,
     "OH" => 1.52,
     "OG" => 1.52,
     "OE1" => 1.52,
     "OE2" => 1.52,
     "OG1" => 1.52,
     "OG2" => 1.52,
     "OD1" => 1.52,
     "OD2" => 1.52,
     "OXT" => 1.52,
     "F" => 1.47,
    #"NE" => 1.54,
    "MG" => 1.73,
    "AL" => 1.84,
    "SI" => 2.10,
     "P" => 1.80,
     "S" => 1.80,
     "SD" => 1.80,
     "SG" => 1.80,
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
    probeAngle::Float64  # Radian
    probeRadius::Float64
    range_min::Point2D
    range_max::Point2D
    resolution::Point2D
    atomRadiusDict::Dict{String, Float64}
end

function afmize!(tip::Matrix, config::AfmizeConfig)
    tip_xsiz, tip_ysiz = size(tip)
    xc, yc = compute_xc_yc(tip)
    #xc += 1
    #yc += 1
    @show xc, yc
    for ix = 1:tip_xsiz
        for iy = 1:tip_ysiz
            x = config.resolution.x * abs(ix - xc)
            y = config.resolution.y * abs(iy - yc)
            d = sqrt(x^2 + y^2)
            if d <= config.probeRadius
                z = sqrt(config.probeRadius^2 - d^2)
            else
                #theta = (0.5*pi) - (0.5*config.probeAngle)
                theta = (0.5*pi) - (config.probeAngle)
                z = - tan(theta) * (d - config.probeRadius)
            end
            tip[ix, iy] = z
        end
    end
    tip .= tip .- maximum(tip)
    return 0
end

function surfing(t::TrjArray, config::AfmizeConfig)
    message = checkConfig(t, config)
    if !isnothing(message)
        println(message)
        return zeros(1, 1)
    end

    width = floor(Int, (config.range_max.x - config.range_min.x) / config.resolution.x)
    height = floor(Int, (config.range_max.y - config.range_min.y) / config.resolution.y)
    stage = zeros(height, width)

    #radius_max = maximum(t.radius)
    radius = zeros(Float64, t.natom)
    for iatom = 1:t.natom
        radius[iatom] = config.atomRadiusDict[t.atomname[iatom]]
    end

    x = t.xyz[1, 1:3:end]
    y = t.xyz[1, 2:3:end]
    z = t.xyz[1, 3:3:end]
    z .= z .- minimum(z)
    for w in 1:width
        grid_x = config.range_min.x + (w-0.5) * config.resolution.x
        dx = abs.(grid_x .- x)
        index_x = dx .< radius
        for h in 1:height
            grid_y = config.range_min.y + (h-0.5) * config.resolution.y
            dy = abs.(grid_y .- y)
            index_y = dy .< radius
            index = index_x .& index_y
            if any(index)
                d = sqrt.(dx[index].^2 + dy[index].^2)
                r = radius[index]
                index2 = d .< r
                if any(index2)
                    z_surface = z[index][index2] .+ sqrt.(r[index2].^2 .- d[index2].^2)
                    stage[h, w] = maximum(z_surface)
                end
            end
        end
    end

    return stage
end

function ChainRulesCore.rrule(::typeof(surfing), xyz::AbstractArray, radius::AbstractArray, max_x::Number, min_x::Number, max_y::Number, min_y::Number, resolution_x::Number, resolution_y::Number)
    width = floor(Int, (max_x - min_x) / resolution_x)
    height = floor(Int, (max_y - min_y) / resolution_y)
    S = zeros(height, width)
    natom3 = length(xyz)
    natom = Int(natom3 / 3)

    #radius_max = maximum(t.radius)

    x = xyz[1:3:end]
    y = xyz[2:3:end]
    z = xyz[3:3:end]
    z .= z .- minimum(z)
    
    S_iatom = zeros(Int, height, width)

    for w in 1:width
        grid_x = min_x + (w-0.5) * resolution_x
        dx = abs.(grid_x .- x)
        index_x = dx .< radius
        for h in 1:height
            grid_y = min_y + (h-0.5) * resolution_y
            dy = abs.(grid_y .- y)
            index_y = dy .< radius
            index = index_x .& index_y
            if any(index)
                d = sqrt.(dx[index].^2 + dy[index].^2)
                r = radius[index]
                index2 = d .< r
                if any(index2)
                    z_surface = z[index][index2] .+ sqrt.(r[index2].^2 .- d[index2].^2)
                    max_s, max_iatom = findmax(z_surface)
                    S[h, w] = max_s
                    S_iatom[h, w] = (1:natom)[index][index2][max_iatom]
                end
            end
        end
    end

    function surfing_pullback(dS)
        dxyz = zeros(eltype(xyz), natom3)
        dradius = zeros(eltype(radius), natom)
        for w in 1:width
            grid_x = min_x + (w-0.5) * resolution_x
            for h in 1:height
                grid_y = min_y + (h-0.5) * resolution_y
                iatom = S_iatom[h, w]
                if iatom != 0
                    ix = 3 * (iatom - 1) + 1
                    dxyz[ix] += - 2.0 * (x[iatom] - grid_x) / sqrt(radius[iatom]^2 - (x[iatom] - grid_x)^2 - (y[iatom] - grid_y)^2) * dS[h, w]
                    iy = 3 * (iatom - 1) + 2
                    dxyz[iy] += - 2.0 * (y[iatom] - grid_y) / sqrt(radius[iatom]^2 - (x[iatom] - grid_x)^2 - (y[iatom] - grid_y)^2) * dS[h, w]
                    iz = 3 * (iatom - 1) + 3
                    dxyz[iz] += dS[h, w]
                    dradius[iatom] += 2.0 * radius[iatom] / sqrt(radius[iatom]^2 - (x[iatom] - grid_x)^2 - (y[iatom] - grid_y)^2) * dS[h, w]
                end
            end
        end
        return NoTangent(), dxyz, dradius, ZeroTangent(), ZeroTangent(), ZeroTangent(), ZeroTangent(), ZeroTangent(), ZeroTangent()
    end
    
    return S, surfing_pullback
end

###################################################
"""
Afmize: calculate pseudo AFM image from a given molecule and a tip shape

the following codes are licensed under the MIT license (Copyright (c) 2018-2021 Tohru Niina), see LICENSE.md
"""

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
            return "config doesn't know atom name $(name)"
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

function afmize(tra::TrjArray, config::AfmizeConfig; removeBottom=true)
    message = checkConfig(tra, config)
    if !isnothing(message)
        println(message)
        return zeros(1, 1)
    end

    width = floor(Int, (config.range_max.x - config.range_min.x) / config.resolution.x)
    height = floor(Int, (config.range_max.y - config.range_min.y) / config.resolution.y)
    atoms = [Sphere(tra.xyz[1, 3*(i-1)+1], tra.xyz[1, 3*(i-1)+2], tra.xyz[1, 3*(i-1)+3],
            config.atomRadiusDict[tra.atomname[i]]) for i in 1:tra.natom]
    if removeBottom
        moveBottom(atoms)
    end

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
            #stage[h, w] = calcCollisionAsSphere(probe, atom)
        end
    end

    for atom in atoms
        rectangle = calRectangle_circularThrusters(atom, config)

        if isnothing(rectangle) continue end

        for h in rectangle[2]:rectangle[4], w in rectangle[1]:rectangle[3]
            probe = probes[h, w]

            stage[h, w] = max(stage[h, w], calcCollisionAsCircularThrusters(probe, atom))
            #stage[h, w] = calcCollisionAsCircularThrusters(probe, atom)
        end
    end

    return stage
end

# 高速化前
function afmize_beta(tra::TrjArray, config::AfmizeConfig)
    message = checkConfig(tra, config)
    if !isnothing(message)
        println(message)
        return zeros(1, 1)
    end

    width = floor(Int, (config.range_max.x - config.range_min.x) / config.resolution.x)
    height = floor(Int, (config.range_max.y - config.range_min.y) / config.resolution.y)
    atoms = [Sphere(tra.xyz[1, 3*(i-1)+1], tra.xyz[1, 3*(i-1)+2], tra.xyz[1, 3*(i-1)+3],
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

function afmize_gpu(tra::TrjArray, config::AfmizeConfig)
    message = checkConfig(tra, config)
    if !isnothing(message)
        println(message)
        return zeros(1, 1)
    end

    width = floor(Int, (config.range_max.x - config.range_min.x) / config.resolution.x)
    height = floor(Int, (config.range_max.y - config.range_min.y) / config.resolution.y)
    atom_r = tra.mass
    atom_z = tra.z[:]
    atom_z = atom_z .- (minimum(atom_z) - atom_r[argmin(atom_z)])
    #atom_z .= atom_z .- (minimum(atom_z))
    atom_z_max = maximum(atom_z)
    stage = similar(tra.x, 1, height*width)
    stage .= 0.0

    probe_x = similar(tra.x, 1, width)
    probe_y = similar(tra.x, 1, height)
    probe_x[:] .= range(config.range_min.x, step=config.resolution.x, stop=config.range_min.x + (width-1)*config.resolution.x) .+ (0.5*config.resolution.x)
    probe_y[:] .= range(config.range_min.y, step=config.resolution.y, stop=config.range_min.y + (height-1)*config.resolution.y) .+ (0.5*config.resolution.y)
    probe_r = config.probeRadius
    probe_angle = config.probeAngle

    # collision with sphere
    dist_x2 = (tra.x[:] .- probe_x).^2
    dist_y2 = (tra.y[:] .- probe_y).^2
    index_x_to_xy = [j for j=1:width for i=1:height]
    index_y_to_xy = [i for j=1:width for i=1:height]

    dist_xy2 = dist_x2[:, index_x_to_xy] .+ dist_y2[:, index_y_to_xy]
    #stage_xy = similar(dist_xy2)
    #stage_xy .= 0.0
    dr = probe_r .+ atom_r
    dr2 = dr.^2
    s = dr2 .- dist_xy2
    @show typeof(s)
    index_collide = s .> 0.0
    index_not_collide = .!index_collide
    @show typeof(index_collide)
    if any(index_collide)
        s[Array(index_collide)] .= sqrt.(s[index_collide])
    end
    s .= atom_z .+ s .- probe_r
    if any(index_not_collide)
        s[Array(index_not_collide)] .= 0.0
    end
    s_max = maximum(s, dims=1)
    stage .= max.(stage, s_max)

    # collision with circular thruster
    dist_xy2 .= sqrt.(dist_xy2)
    dist_collision = probe_r .+ atom_r .* cos(probe_angle)
    index_side   = dist_collision .< dist_xy2
    index_corner = (probe_r .< dist_xy2) .& .!index_side
    index_not_collide = .!(index_corner .| index_side)
    s[Array(index_side)] .= (atom_r .* sin.(probe_angle) .- (dist_xy2 .- dist_collision) ./ tan.(probe_angle))[index_side]
    s[Array(index_corner)] .= sqrt.( (atom_r.^2 .- (dist_xy2 .- probe_r).^2)[index_corner]  )
    s .= s .+ atom_z .- probe_r
    s[Array(index_not_collide)] .= 0.0
    s_max = maximum(s, dims=1)
    stage .= max.(stage, s_max)

    stage_reshape = reshape(stage, (height, width))
    return stage_reshape
end

###################################################
"""
AFM image analysis: calculate likelihood of many parameters, such as the orientation of molecule from AFM images. 
"""

function translateafm(afm, (dx, dy))
    afm_translated = zeros(eltype(afm), size(afm))
    (ny, nx) = size(afm)
    for i in maximum([1, 1-dx]):minimum([nx, nx-dx])
        for j in maximum([1, 1-dy]):minimum([ny, ny-dy])
            afm_translated[j+dy, i+dx] = afm[j, i]
        end
    end
    afm_translated
end

function translateafm_periodic(afm, (dx, dy))
    afm_translated = zeros(eltype(afm), size(afm))
    (nx, ny) = size(afm)
    for x in 1:nx
        for y in 1:ny
            xx = x + dx
            yy = y + dy
            if xx > nx
                xx -= nx
            end
            if yy > ny
                yy -= ny
            end
            afm_translated[xx, yy] = afm[x, y]
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

function direct_convolution_non_periodic(observed, calculated)
    H, W = size(observed)
    ret = zeros(H, W)
    for y in 1:H, x in 1:W
        for dy in 1:(H-y+1), dx in 1:(W-x+1)
            ny = y + dy - 1
            nx = x + dx - 1
            ret[y, x] += observed[dy, dx] * calculated[ny, nx]
        end
    end
    return ret
end

function fft_convolution(observed, calculated)
    #return real.(ifft(fft(observed).*conj.(fft(calculated))))
    return real.(ifftshift(ifft(fft(observed).*conj.(fft(calculated)))))
end

function calcLogProb(observed::AbstractMatrix{T}, calculated::AbstractMatrix{T}, convolution_func=fft_convolution) where {T}
    npix = eltype(observed)(size(observed, 1) * size(observed, 2))
    C_o  = sum(observed)
    C_c  = sum(calculated)
    C_cc = sum(calculated.^2)
    C_oo = sum(observed.^2)
    C_oc = convolution_func(observed, calculated)

    log01 = npix .* (C_cc .* C_oo .- C_oc.^2) .+ eltype(observed)(2) .* C_o .* C_oc .* C_c .- C_cc .* C_o.^2 .- C_oo .* C_c.^2
    #log01[Array(log01 .< log(eps(eltype(observed))))] .= eps(eltype(observed))
    log02 = (npix .- eltype(observed)(2)) .* (npix .* C_cc .- C_c.^2)
    #log02 = log02 <= log(eps(eltype(observed))) ? eps(eltype(observed)) : log02
    logprob = eltype(observed)(0.5) .* (eltype(observed)(3.0) .- npix) .* log.(log01) .+ (eltype(observed)(0.5) .* npix .- eltype(observed)(2.0)) .* log(log02)

    return logprob
end

function getafmposterior(afm::Matrix{Float64}, model_array::TrjArray, q_array::Matrix{Float64}, param_array)
    imax_model = 0
    imax_q = 0
    best_model = model_array[1, :]
    best_param = param_array[1]
    best_translate = [0, 0]
    best_afm = similar(afm)
    best_posterior = -Inf

    observed = afm
    npix = size(observed, 1) * size(observed, 2)
    decenter!(model_array)

    posterior_all = zeros(Float64, size(q_array, 1)*length(param_array)*npix, size(model_array, 1))
    posterior = zeros(Float64, size(model_array, 1))

    ### loop over models(structures)
    Threads.@threads for imodel in 1:size(model_array, 1)
        @show imodel
        icount = 0
        model = model_array[imodel, :]
        ### loop over rotations
        for iq in 1:size(q_array, 1)
            q = q_array[iq, :]
            model_rotated = MDToolbox.rotate(model, q)
            ### loop over afmize parameters
            for iparam in 1:length(param_array)
                param = param_array[iparam]
                calculated = afmize(model_rotated, param)
                logprob = calcLogProb(observed, calculated)
                posterior_all[(icount+1):(icount+npix), imodel] .= logprob[:]
                icount += npix
                maximum_logprob = maximum(logprob)
                if best_posterior < maximum_logprob
                    best_posterior = maximum_logprob
                    imax_model = imodel
                    best_model = model_rotated
                    imax_q = iq
                    best_param = param
                    x_center = ceil(Int64, (size(observed,1)/2)+1.0)
                    y_center = ceil(Int64, (size(observed,2)/2)+1.0)
                    dx_estimated = argmax(logprob)[1] - x_center
                    dy_estimated = argmax(logprob)[2] - y_center
                    best_translate = (dx_estimated, dy_estimated)
                    best_afm = translateafm(calculated, best_translate)
                end
            end
        end
    end

    for imodel in 1:size(model_array, 1)
        posterior[imodel] = logsumexp(posterior_all[:, imodel])
    end
    posterior .= posterior .- maximum(posterior)
    posterior .= exp.(posterior)
    posterior .= posterior ./ sum(posterior)

    return imax_model, imax_q, best_param, best_translate, best_afm, best_posterior, posterior
end


function getafmposterior_gpu(afm::AbstractMatrix{T}, model_array::TrjArray{T, U}, q_array::AbstractMatrix{T}, param_array) where {T, U}
    imax_model = 0
    imax_q = 0
    best_param = param_array[1]
    best_translate = [0, 0]
    best_afm = similar(afm)
    best_posterior = -Inf

    observed = afm
    decenter!(model_array)

    for imodel in 1:size(model_array, 1)
        @show imodel
        model = model_array[imodel, :]
        models_rotated = MDToolbox.rotate(model, q_array)
        for iq in 1:size(q_array, 1)
            model_rotated = models_rotated[iq, :]
            for iparam in 1:length(param_array)
                param = param_array[iparam]
                #@show typeof(model_rotated.x)
                calculated = MDToolbox.afmize_gpu(model_rotated, param)
                #@show typeof(calculated)
                logprob = MDToolbox.calcLogProb(observed, calculated, MDToolbox.fft_convolution)
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


function logprob_eachmodel(model::TrjArray, afm_array, q_array::Matrix{Float64}, param_array)
    nq = size(q_array, 1)
    nparam = length(param_array)
    nafm = length(afm_array)
    npix = size(afm_array[1], 1) * size(afm_array[1], 2)
    logprob_array = zeros(eltype(afm_array[1]), nafm, nq, nparam)
    dxdy_array = Array{Tuple{Int,Int}}(undef, nafm, nq, nparam)
    dxdy = Array{Tuple{Int,Int}}(undef, nafm)
    x_center = ceil(Int64, (size(afm_array[1],2)/2)+1.0)
    y_center = ceil(Int64, (size(afm_array[2],1)/2)+1.0)

    decenter!(model)    
    ### loop over rotations
    Threads.@threads for iq in 1:nq
        @show iq
        q = q_array[iq, :]
        model_rotated = MDToolbox.rotate(model, q)
        ### loop over afmize parameters
        for iparam in 1:nparam
            param = param_array[iparam]
            calculated = afmize(model_rotated, param)
            for iafm in 1:nafm
                observed = afm_array[iafm]
                logprob = calcLogProb(observed, calculated)
                logprob_array[iafm, iq, iparam] = logsumexp(logprob)
                dx = argmax(logprob)[2] - x_center
                dy = argmax(logprob)[1] - y_center
                dxdy_array[iafm, iq, iparam] = (dx, dy)
            end
        end
    end

    for iafm in 1:nafm
        ind = argmax(logprob_array[iafm, :, :])
        dxdy[iafm] = dxdy_array[iafm, :, :][ind]
    end

    return (logprob_array=logprob_array, dxdy=dxdy)
end


function getposterior_parallel(models::TrjArray, afm_array, q_array::Matrix{Float64}, param_array)
    nmodel = models.nframe
    nafm = length(afm_array)
    nq = size(q_array, 1)
    nparam = length(param_array)

    p = @showprogress pmap(x -> logprob_eachmodel(x, afm_array, q_array, param_array), models)

    logprob_all = []
    logprob_model = []
    logprob_q = []
    logprob_param = []
    dxdy_best = []
    for iafm = 1:nafm
        pp = zeros(Float64, nmodel, nq, nparam)
        for imodel = 1:nmodel
            for iq = 1:nq
                for iparam = 1:nparam
                    pp[imodel, iq, iparam] = p[imodel].logprob_array[iafm, iq, iparam]
                end
            end
        end

        push!(logprob_all, pp)

        t = zeros(Float64, nmodel)
        for imodel = 1:nmodel
            t[imodel] = logsumexp(pp[imodel, :, :][:])
        end
        push!(logprob_model, t)

        imax = argmax(t)
        push!(dxdy_best, p[imax].dxdy[iafm])

        t = zeros(Float64, nq)
        for iq = 1:nq
            t[iq] = logsumexp(pp[:, iq, :][:])
        end
        push!(logprob_q, t)

        t = zeros(Float64, nparam)
        for iparam = 1:nparam
            t[iparam] = logsumexp(pp[:, :, iparam][:])
        end
        push!(logprob_param, t)
    end

    return (all=logprob_all, 
            model=logprob_model, 
            q=logprob_q, 
            param=logprob_param,
            dxdy=dxdy_best)
end


mutable struct posteriorResult
    each_quate_id
    each_param_id
    posterior_results
    best_translate
    best_posterior
    best_model_id
    best_quate_id
    best_param_id
    best_afm
end

mutable struct posteriorTmpData
    posteriors
    each_best
end

function outputResults(posteriorResults, fileName)
    open(fileName, "w") do out
        N = size(posteriorResults)[1]
        model_size = size(posteriorResults[1].each_quate_id)[1]
        frame_h, frame_w = size(posteriorResults[1].best_afm)
        println(out, "$(N) $(model_size)")
        for result in posteriorResults
            for i in 1:model_size
                if i != 1 print(out, " ") end
                print(out, "$(result.each_quate_id[i])")
            end
            println(out, "")
            for i in 1:model_size
                if i != 1 print(out, " ") end
                print(out, "$(result.each_param_id[i])")
            end
            println(out, "")
            for i in 1:model_size
                if i != 1 print(out, " ") end
                print(out, "$(result.posterior_results[i])")
            end
            println(out, "")
            println(out, "$(result.best_translate[1]) $(result.best_translate[2])")
            println(out, "$(result.best_posterior)")
            println(out, "$(result.best_model_id)")
            println(out, "$(result.best_quate_id)")
            println(out, "$(result.best_param_id)")
            println(out, "$(frame_h) $(frame_w)")
            for h in 1:frame_h
                for w in 1:frame_w
                    if w != 1 print(out, " ") end
                    print(out, "$(result.best_afm[h, w])")
                end
                println(out, "")
            end
        end
    end
end

function inputResults(fileName)
    ret = []
    open(fileName, "r") do io
        N, model_size = map(x->parse(Int,x),split(readline(io)))
        for result in 1:N
            each_quate_id = map(x->parse(Int,x),split(readline(io)))
            each_param_id = map(x->parse(Int,x),split(readline(io)))
            posterior_results = map(x->parse(Float64,x),split(readline(io)))
            best_translate = Tuple(map(x->parse(Int,x),split(readline(io))))
            best_posterior = parse(Float64, readline(io))
            best_model_id = parse(Int, readline(io))
            best_quate_id = parse(Int, readline(io))
            best_param_id = parse(Int, readline(io))
            frame_h, frame_w = map(x->parse(Int,x),split(readline(io)))
            frame = zeros(Float64, frame_h, frame_w)
            for h in 1:frame_h
                frame[h, :] = map(x->parse(Float64,x),split(readline(io)))
            end
            push!(ret, posteriorResult(each_quate_id,
                                        each_param_id,
                                        posterior_results,
                                        best_translate,
                                        best_posterior,
                                        best_model_id,
                                        best_quate_id,
                                        best_param_id,
                                        frame))
        end
    end
    return ret
end
