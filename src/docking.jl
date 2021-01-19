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
    radius = Array{T}(undef, ta.natom)
    element = Array{String}(undef, ta.natom)
    for iatom = 1:ta.natom
        if !isnothing(match(r"H.*", ta.atomname[iatom]))
            element[iatom] = "H"
        elseif !isnothing(match(r"C.*", ta.atomname[iatom]))
            element[iatom] = "C"
        elseif !isnothing(match(r"N.*", ta.atomname[iatom]))
            element[iatom] = "N"
        elseif !isnothing(match(r"O.*", ta.atomname[iatom]))
            element[iatom] = "O"
        elseif !isnothing(match(r"S.*", ta.atomname[iatom]))
            element[iatom] = "S"
        else
            error("failed to assign element: " * ta.atomname[iatom])
        end
    end

    radius_dict = Dict("H" => 1.20, 
                       "C" => 1.70, 
                       "N" => 1.55, 
                       "O" => 1.52, 
                       "S" => 1.80)

    for iatom = 1:ta.natom
        radius[iatom] = radius_dict[element[iatom]]
    end

    return TrjArray(ta, radius=radius)
end

function compute_sasa(ta::TrjArray{T, U}, probe_radius=1.4::T; npoint=960::Int, iframe=1::Int, candicate = 10) where {T, U}
    # construct neighbor rist
    max_radius = 2.0 * maximum(ta.radius) + 2.0 * probe_radius ############# TODO
    pairlist = compute_pairlist(ta[iframe, :], max_radius)
    neighbor_list = []
    for iatom in 1:ta.natom
        push!(neighbor_list, Array{U}(undef, 0))
    end
    for ipair in 1:size(pairlist.pair, 1)
        i = pairlist.pair[ipair, 1]
        j = pairlist.pair[ipair, 2]
        push!(neighbor_list[i], j)
        push!(neighbor_list[j], i)
    end

    # generate uniform points on a unit sphere
    points = golden_section_spiral(npoint)

    # compute the ratio of exposed area for each sphere
    sasa = Array{T}(undef, ta.natom)
    Threads.@threads for iatom = 1:ta.natom
        n_accessible_point = 0
        neighbor_list_iatom = neighbor_list[iatom]
        for ipoint in 1:npoint
            is_accessible = true
            point = points[ipoint, :] .* (ta.radius[iatom] + probe_radius)
            point[1] += ta.x[iframe, iatom]
            point[2] += ta.y[iframe, iatom]
            point[3] += ta.z[iframe, iatom]
            for j in 1:length(neighbor_list_iatom)
                jatom = neighbor_list_iatom[j]
                d = 0.0
                d += (point[1] - ta.x[iframe, jatom])^2
                d += (point[2] - ta.y[iframe, jatom])^2
                d += (point[3] - ta.z[iframe, jatom])^2
                d = sqrt(d)
                if d < (ta.radius[jatom] + probe_radius)
                    is_accessible = false
                    break
                end
            end
            if is_accessible
                n_accessible_point += 1
            end
        end
        sasa[iatom] = 4.0 * pi * (ta.radius[iatom] + probe_radius)^2 * n_accessible_point / npoint
    end

    return sasa
end

function assign_shape_complementarity!(grid, ta::TrjArray{T, U}, grid_space, 
                                       rcut1, rcut2, x_grid, y_grid, z_grid, iframe) where {T,U}
    grid .= 0.0
    nx, ny, nz = size(grid)

    for iatom = 1:ta.natom
        x = ta.x[iframe, iatom]
        y = ta.y[iframe, iatom]
        z = ta.z[iframe, iatom]

        rcut = rcut1[iatom]

        dx = x - x_grid[1]
        ix_min = floor(U, (dx - rcut)/grid_space) + 1
        ix_min = ix_min >= 1 ? ix_min : 1
        ix_max = floor(U, (dx + rcut)/grid_space) + 2
        ix_max = ix_max <= nx ? ix_max : nx

        dy = y - y_grid[1]
        iy_min = floor(U, (dy - rcut)/grid_space) + 1
        iy_min = iy_min >= 1 ? iy_min : 1
        iy_max = floor(U, (dy + rcut)/grid_space) + 2
        iy_max = iy_max <= ny ? iy_max : ny

        dz = z - z_grid[1]
        iz_min = floor(U, (dz - rcut)/grid_space) + 1
        iz_min = iz_min >= 1 ? iz_min : 1
        iz_max = floor(U, (dz + rcut)/grid_space) + 2
        iz_max = iz_max <= nz ? iz_max : nz

        for ix = ix_min:ix_max
            for iy = iy_min:iy_max
                for iz = iz_min:iz_max
                    if imag(grid[ix, iy, iz]) < 0.0001
                        dist = sqrt((x - x_grid[ix])^2 + (y - y_grid[ix])^2 + (z - z_grid[ix])^2)
                        if dist < rcut
                            grid[ix, iy, iz] = 9.0im
                        end
                    end
                end
            end
        end

        rcut = rcut2[iatom]

        dx = x - x_grid[1]
        ix_min = floor(U, (dx - rcut)/grid_space) + 1
        ix_min = ix_min >= 1 ? ix_min : 1
        ix_max = floor(U, (dx + rcut)/grid_space) + 2
        ix_max = ix_max <= nx ? ix_max : nx

        dy = y - y_grid[1]
        iy_min = floor(U, (dy - rcut)/grid_space) + 1
        iy_min = iy_min >= 1 ? iy_min : 1
        iy_max = floor(U, (dy + rcut)/grid_space) + 2
        iy_max = iy_max <= ny ? iy_max : ny

        dz = z - z_grid[1]
        iz_min = floor(U, (dz - rcut)/grid_space) + 1
        iz_min = iz_min >= 1 ? iz_min : 1
        iz_max = floor(U, (dz + rcut)/grid_space) + 2
        iz_max = iz_max <= nz ? iz_max : nz

        for ix = ix_min:ix_max
            for iy = iy_min:iy_max
                for iz = iz_min:iz_max
                    if imag(grid[ix, iy, iz]) < 0.0001
                        dist = sqrt((x - x_grid[ix])^2 + (y - y_grid[ix])^2 + (z - z_grid[ix])^2)
                        if dist < rcut
                            grid[ix, iy, iz] = 1.0
                        end
                    end
                end
            end
        end
    end
end

function compute_docking_score_with_fft(quaternion, grid_RSC, grid_LSC, ligand2, grid_space, rcut1, rcut2, x_grid, y_grid, z_grid, iframe, nx, ny, nz, tops, iq)
    ligand2_rotated = rotate(ligand2, quaternion)
    assign_shape_complementarity!(grid_LSC, ligand2_rotated, grid_space, rcut1, rcut2, x_grid, y_grid, z_grid, iframe)
    if CUDA.functional()
        grid_RSC_gpu = cu(grid_RSC)
        grid_LSC_gpu = cu(grid_LSC)
        t_gpu = ifftshift(ifft(ifft(grid_RSC_gpu) .* fft(grid_LSC_gpu)))
        score_gpu = real(t_gpu) .- imag(t_gpu)
        score = Array(score_gpu)
    else
        t = ifftshift(ifft(ifft(grid_RSC) .* fft(grid_LSC)))
        score = real(t) .- imag(t)
    end
    
    ret = []
    x_center = ceil(Int64, (nx/2)+1.0)
    y_center = ceil(Int64, (ny/2)+1.0)
    z_center = ceil(Int64, (nz/2)+1.0)
    for t in 1:tops
        id = argmax(score)
        dx_estimated = id[1] - x_center
        dy_estimated = id[2] - y_center  
        dz_estimated = id[3] - z_center 
        
        push!(ret, (score[id], dx_estimated, dy_estimated, dz_estimated, iq))        
        score[id] = -Inf
    end
    
    return ret
end

function dock_fft(receptor::TrjArray{T, U}, ligand::TrjArray{T, U}, quaternions; grid_space=1.2, iframe=1, tops=10) where {T, U}
    # generate grid coordinates for receptor
    receptor2, _dummy = decenter(receptor)
    ligand2, _dummy = decenter(ligand)

    x_min, x_max = minimum(ligand2.x), maximum(ligand2.x)
    y_min, y_max = minimum(ligand2.y), maximum(ligand2.y)
    z_min, z_max = minimum(ligand2.z), maximum(ligand2.z)
    size_ligand = sqrt((x_max - x_min)^2 + (y_max - y_min)^2 + (z_max - z_min)^2)

    x_min = minimum(receptor2.x) - size_ligand - grid_space
    y_min = minimum(receptor2.y) - size_ligand - grid_space
    z_min = minimum(receptor2.z) - size_ligand - grid_space

    x_max = maximum(receptor2.x) + size_ligand + grid_space
    y_max = maximum(receptor2.y) + size_ligand + grid_space
    z_max = maximum(receptor2.z) + size_ligand + grid_space

    x_grid = collect(x_min:grid_space:x_max)
    y_grid = collect(y_min:grid_space:y_max)
    z_grid = collect(z_min:grid_space:z_max)

    nx, ny, nz = length(x_grid), length(y_grid), length(z_grid)

    sasa_receptor = receptor2.mass
    sasa_ligand = ligand2.mass

    # spape complementarity of receptor
    iatom_surface = sasa_receptor .> 1.0
    iatom_core = .!iatom_surface
    rcut1 = zeros(T, receptor2.natom)
    rcut2 = zeros(T, receptor2.natom)
    
    rcut1[iatom_core] .= receptor2.radius[iatom_core] * sqrt(1.5)
    rcut2[iatom_core] .= 0.0
    rcut1[iatom_surface] .= receptor2.radius[iatom_surface] * sqrt(0.8)
    rcut2[iatom_surface] .= receptor2.radius[iatom_surface] .+ 3.4

    grid_RSC = zeros(complex(T), nx, ny, nz)
    assign_shape_complementarity!(grid_RSC, receptor2, grid_space, rcut1, rcut2, x_grid, y_grid, z_grid, iframe)

    # spape complementarity of ligand
    iatom_surface = sasa_ligand .> 1.0
    iatom_core = .!iatom_surface
    rcut1 = zeros(T, ligand2.natom)
    rcut2 = zeros(T, ligand2.natom)
        
    rcut1[iatom_core] .= ligand2.radius[iatom_core] * sqrt(1.5)
    rcut2[iatom_core] .= 0.0
    rcut1[iatom_surface] .= 0.0
    rcut2[iatom_surface] .= ligand2.radius[iatom_surface]
    
    grid_LSC = zeros(complex(T), nx, ny, nz)

    # compute score with FFT
    nq = size(quaternions, 1)
    score = zeros(T, nx, ny, nz, nq)
    score .= -Inf
    #for iq in 1:size(quaternions, 1)
    #    ligand2_rotated = rotate(ligand2, quaternions[iq, :])
    #    assign_shape_complementarity!(grid_LSC, ligand2_rotated, grid_space, rcut1, rcut2, x_grid, y_grid, z_grid, iframe)
    # 
    #    t = ifftshift(ifft(ifft(grid_RSC) .* fft(grid_LSC)))
    #    score[:, :, :, iq] .= real(t) .- imag(t)
    #end
    #for iq in 1:size(quaternions, 1)
    #    s = compute_docking_score_with_fft(quaternions[iq, :], grid_RSC, grid_LSC, ligand2, grid_space, rcut1, rcut2, x_grid, y_grid, z_grid, iframe)
    #    score[:, :, :, iq] .= s
    #end
    s = @showprogress pmap(x -> compute_docking_score_with_fft(quaternions[x, :], grid_RSC, grid_LSC, ligand2, grid_space, rcut1, rcut2, x_grid, y_grid, z_grid, iframe, nx, ny, nz, tops, x), 1:nq)
    
    result_tops = []
    for iq = 1:nq
        result_tops = [result_tops; s[iq]]
    end
    sort!(result_tops, rev=true)
    resize!(result_tops, tops)

    score_tops = zeros(T, tops)
    trans_tops = zeros(Int64, tops, 3)
    quate_tops = zeros(Float64, tops, 4)

    for t in 1:tops
        score_tops[t] = result_tops[t][1]
        trans_tops[t, :] = [result_tops[t][2], result_tops[t][3], result_tops[t][4]]
        quate_tops[t, :] = quaternions[result_tops[t][5], :]
    end

    itop = 1
    ligand_return = MDToolbox.rotate(ligand2, quate_tops[itop, :])
    ligand_return.x .+= (trans_tops[itop, 1] * grid_space)
    ligand_return.y .+= (trans_tops[itop, 2] * grid_space)
    ligand_return.z .+= (trans_tops[itop, 3] * grid_space)

    for itop = 2:tops
        ligand_tmp = MDToolbox.rotate(ligand2, quate_tops[itop, :])
        ligand_tmp.x .+= (trans_tops[itop, 1] * grid_space)
        ligand_tmp.y .+= (trans_tops[itop, 2] * grid_space)
        ligand_tmp.z .+= (trans_tops[itop, 3] * grid_space)    
        ligand_return = [ligand_return; ligand_tmp]
    end

    return (receptor=receptor2, ligand=ligand_return, score=score_tops)
end