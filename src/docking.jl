function golden_section_spiral(n)
    points = zeros(Float64, n, 3)
    inc = pi * (3.0 - sqrt(5.0))
    offset = 2.0 / Float64(n)
    for k = 1:n
        y = (k-1) * offset - 1.0 + (offset / 2.0)
        r = sqrt(1.0 - y*y)
        phi = (k-1) * inc
        points[k, 1] = cos(phi) * r
        points[k, 2] = y
        points[k, 3] = sin(phi) * r
    end
    return points
end

function set_radius(ta::TrjArray{T, U}) where {T, U}
    r = 1.0
end

function compute_sasa(ta::TrjArray{T, U}; probe_radius=8.0::T, npoints=960::Int, frame=1::Int) where {T, U}
    # construct neighbor rist
    #radius_max = maximum(maximum(ta.radius), probe_radius)
    radius_max = 2.0
    pairlist = compute_pairlist(ta[frame, :], radius_max*4.0)
    neighbor_list = []
    dist_list = []
    for iatom in 1:ta.natom
        push!(neighbor_list, Array{U}(undef, 0))
        push!(dist_list, Array{T}(undef, 0))
    end
    for ipair in 1:size(pairlist.pair, 1)
        i = pairlist.pair[ipair, 1]
        j = pairlist.pair[ipair, 2]
        d = pairlist.dist[ipair]
        push!(neighbor_list[i], j)
        push!(neighbor_list[j], i)
        push!(dist_list[i], d)
        push!(dist_list[j], d)
    end

    points = golden_section_spiral(npoints)
    is_accessible = true
    sasa = Array{T}(undef, ta.natom)

    for iatom = 1:ta.natom
        n_accessible_point = 0
        for ipoint in 1:npoints
            is_accessible = true
            point = points[ipoint, :] .* (ta.radius[iatom] + probe_radius) ################
            point[1] =  point[1] + ta.x[frame, iatom]
            point[2] =  point[2] + ta.y[frame, iatom]
            point[3] =  point[3] + ta.z[frame, iatom]
            for j in 1:length(neighbor_list[i])
                jatom = neighbor_list[j]
                d = dist_list[j]
                if d < (ta.radius[jatom] + probe_radius)
                    is_accessible = false
                    break
                end
            end
        end
        if is_accessible
            n_accessible_point += 1
        end
        sasa[iatom] = 4.0 * pi * probe_radius * probe_radius * n_accessible_point / npoints
    end

    return sasa
end
