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

function compute_sasa(ta::TrjArray{T, U}, probe_radius=8.0::T; npoint=960::Int, iframe=1::Int, candicate = 10) where {T, U}
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

function assign_shape_complementarity!(thread_id, grid_space, rcut1, rcut2, 
    x, y, z, x_grid, y_grid, z_grid, nx, ny, nz, grid_private)

    rcut = rcut1 > rcut2 ? rcut1 : rcut2
    dx = x - x_grid[1]
    ix_min = floor(Int, (dx - rcut)/grid_space) + 1
    ix_min = ix_min >= 1 ? ix_min : 1
    ix_max = floor(Int, (dx + rcut)/grid_space) + 2
    ix_max = ix_max <= nx ? ix_max : nx

    dy = y - y_grid[1]
    iy_min = floor(Int, (dy - rcut)/grid_space) + 1
    iy_min = iy_min >= 1 ? iy_min : 1
    iy_max = floor(Int, (dy + rcut)/grid_space) + 2
    iy_max = iy_max <= ny ? iy_max : ny

    dz = z - z_grid[1]
    iz_min = floor(Int, (dz - rcut)/grid_space) + 1
    iz_min = iz_min >= 1 ? iz_min : 1
    iz_max = floor(Int, (dz + rcut)/grid_space) + 2
    iz_max = iz_max <= nz ? iz_max : nz

    for ix = ix_min:ix_max
        for iy = iy_min:iy_max
            for iz = iz_min:iz_max
                dist = sqrt((x - x_grid[ix])^2 + (y - y_grid[ix])^2 + (z - z_grid[ix])^2)
                if dist < rcut1
                    grid_private[ix, iy, iz, thread_id] = 9.0im
                elseif dist < rcut2
                    grid_private[ix, iy, iz, thread_id] = 1.0
                end
            end
        end
    end
end

function assign_desolvation_free_energy!(thread_id, grid_space, rcut1, rcut2, 
    x, y, z, x_grid, y_grid, z_grid, nx, ny, nz, grid_private)

    rcut = rcut1 > rcut2 ? rcut1 : rcut2
    dx = x - x_grid[1]
    ix_min = floor(Int, (dx - rcut)/grid_space) + 1
    ix_min = ix_min >= 1 ? ix_min : 1
    ix_max = floor(Int, (dx + rcut)/grid_space) + 2
    ix_max = ix_max <= nx ? ix_max : nx

    dy = y - y_grid[1]
    iy_min = floor(Int, (dy - rcut)/grid_space) + 1
    iy_min = iy_min >= 1 ? iy_min : 1
    iy_max = floor(Int, (dy + rcut)/grid_space) + 2
    iy_max = iy_max <= ny ? iy_max : ny

    dz = z - z_grid[1]
    iz_min = floor(Int, (dz - rcut)/grid_space) + 1
    iz_min = iz_min >= 1 ? iz_min : 1
    iz_max = floor(Int, (dz + rcut)/grid_space) + 2
    iz_max = iz_max <= nz ? iz_max : nz

    for ix = ix_min:ix_max
        for iy = iy_min:iy_max
            for iz = iz_min:iz_max
                dist = sqrt((x - x_grid[ix])^2 + (y - y_grid[ix])^2 + (z - z_grid[ix])^2)
                if dist < rcut1
                    grid_private[ix, iy, iz, thread_id] = 9.0im
                elseif dist < rcut2
                    grid_private[ix, iy, iz, thread_id] = 1.0
                end
            end
        end
    end


    # buffer
    for ix = ix_min:ix_max
        for iy = iy_min:iy_max
            for iz = iz_min:iz_max
                if grid_private[ix, iy, iz, thread_id] == 9.0im continue end
                
                for dx in -1:1
                    for dy in -1:1
                        for dz in -1:1
                            if dx == 0 && dy == 0 && dz == 0 continue end
                            now_x = ix + dx
                            now_y = iy + dy
                            now_z = iz + dz
                            if now_x < 1 || now_y < 1 || now_z < 1 continue end
                            if now_x > size(grid_private, 1) || now_y > size(grid_private, 2) || now_z > size(grid_private, 3) continue end
                            
                            if grid_private[now_x, now_y, now_z, thread_id] == 9.0im
                                # println("in $(now_x) $(now_y) $(now_z)")
                                grid_private[ix, iy, iz, thread_id] += 1
                            end
                        end
                    end
                end
            end
        end
    end

end

function dock_fft(receptor::TrjArray{T, U}, ligand::TrjArray{T, U}, quaternions; grid_space=1.2, iframe=1, tops=10) where {T, U}
    # generate grid coordinates for receptor
    receptor2, = decenter(receptor)
    ligand2, = decenter(ligand)

    x_min, x_max = minimum(ligand2.x), maximum(ligand2.x)
    y_min, y_max = minimum(ligand2.y), maximum(ligand2.y)
    z_min, z_max = minimum(ligand2.z), maximum(ligand2.z)
    size_ligand = sqrt((x_max - x_min)^2 + (y_max - y_min)^2 + (z_max - z_min)^2)

    x_min, x_max = minimum(receptor2.x) - size_ligand, maximum(receptor2.x) + size_ligand
    y_min, y_max = minimum(receptor2.y) - size_ligand, maximum(receptor2.y) + size_ligand
    z_min, z_max = minimum(receptor2.z) - size_ligand, maximum(receptor2.z) + size_ligand

    x_grid = collect((x_min - grid_space):grid_space:(x_max + grid_space))
    y_grid = collect((y_min - grid_space):grid_space:(y_max + grid_space))
    z_grid = collect((z_min - grid_space):grid_space:(z_max + grid_space))

    nx, ny, nz = length(x_grid), length(y_grid), length(z_grid)

    sasa_receptor = receptor2.mass
    sasa_ligand = ligand2.mass

    # spape complementarity of receptor
    iatom_surface = sasa_receptor .> 1.0
    iatom_core = .!iatom_surface
    rcut1 = zeros(T, receptor2.natom)
    rcut2 = zeros(T, receptor2.natom)
    
    if size(iatom_core, 1) > 0
        rcut1[iatom_core] .= receptor2.radius[iatom_core] * sqrt(1.5)
        rcut2[iatom_core] .= 0.0
    end
    if size(iatom_surface, 1) > 0
        rcut1[iatom_surface] .= receptor2.radius[iatom_surface] * sqrt(0.8)
        rcut2[iatom_surface] .= receptor2.radius[iatom_surface] .+ 3.4
    end

    
    nthread = Threads.nthreads()
    grid_private = zeros(complex(T), nx, ny, nz, nthread)
    Threads.@threads for iatom = 1:receptor2.natom
        tid = Threads.threadid()
        assign_shape_complementarity!(tid, grid_space, rcut1[iatom], rcut2[iatom], 
                                      receptor2.x[iframe, iatom], receptor2.y[iframe, iatom], receptor2.z[iframe, iatom], 
                                      x_grid, y_grid, z_grid, nx, ny, nz, grid_private)
    end
    grid_RSC = dropdims(sum(grid_private, dims=4), dims=4)

    # spape complementarity of ligand
    iatom_surface = sasa_ligand .> 1.0
    iatom_core = .!iatom_surface
    rcut1 = zeros(T, ligand2.natom)
    rcut2 = zeros(T, ligand2.natom)
    
    if size(iatom_core, 1) > 0        
        rcut1[iatom_core] .= ligand2.radius[iatom_core] * sqrt(1.5)
        rcut2[iatom_core] .= 0.0
    end
    if size(iatom_surface, 1) > 0
        rcut1[iatom_surface] .= 0.0
        rcut2[iatom_surface] .= ligand2.radius[iatom_surface]
    end

    grid_LSC = zeros(complex(T), nx, ny, nz)
    score = zeros(T, nx, ny, nz, size(quaternions, 1))
    for iq in 1:size(quaternions, 1)
        ligand2_rotated = rotate(ligand2, quaternions[iq, :])
        grid_private .= 0.0
        Threads.@threads for iatom = 1:ligand2.natom
            tid = Threads.threadid()
            assign_shape_complementarity!(tid, grid_space, rcut1[iatom], rcut2[iatom], 
                                          ligand2_rotated.x[iframe, iatom], ligand2_rotated.y[iframe, iatom], 
                                          ligand2_rotated.z[iframe, iatom], 
                                          x_grid, y_grid, z_grid, nx, ny, nz, grid_private)
        end
        grid_LSC .= dropdims(sum(grid_private, dims=4), dims=4)    
        #@show iq sum(grid_LSC)

        t = ifftshift(ifft(fft(grid_RSC) .* conj.(fft(grid_LSC))))
        score[:, :, :, iq] .= real(t) .- imag(t)
    end

    # desolvation_free_energy of receptor
    iatom_surface = sasa_receptor .> 1.0
    iatom_core    = .!iatom_surface
    rcut1                 = zeros(T, receptor2.natom)
    rcut2                 = zeros(T, receptor2.natom)
    if size(iatom_core, 1) > 0
        rcut1[iatom_core]    .= receptor2.radius[iatom_core] * sqrt(1.5)
        rcut2[iatom_core]    .= 0.0
    end
    if size(iatom_surface, 1) > 0
        rcut1[iatom_surface] .= receptor2.radius[iatom_surface] * sqrt(0.8)
        rcut2[iatom_surface] .= receptor2.radius[iatom_surface] .+ 3.4
    end


    nthread = Threads.nthreads()
    grid_private = zeros(complex(T), nx, ny, nz, nthread)
    Threads.@threads for iatom = 1:receptor2.natom
        tid = Threads.threadid()
        assign_desolvation_free_energy!(tid, grid_space, rcut1[iatom], rcut2[iatom], 
                                      receptor2.x[iframe, iatom], receptor2.y[iframe, iatom], receptor2.z[iframe, iatom], 
                                      x_grid, y_grid, z_grid, nx, ny, nz, grid_private)
    end
    grid_RDS = dropdims(sum(grid_private, dims=4), dims=4) 

    # desolvation_free_energy of ligand
    iatom_surface = sasa_ligand .> 1.0
    iatom_core = .!iatom_surface
    rcut1 = zeros(T, ligand2.natom)
    rcut2 = zeros(T, ligand2.natom)
    
    if size(iatom_core, 1) > 0
        rcut1[iatom_core] .= ligand2.radius[iatom_core] * sqrt(1.5)
        rcut2[iatom_core] .= 0.0
    end

    if size(iatom_surface, 1) > 0
        rcut1[iatom_surface] .= 0.0
        rcut2[iatom_surface] .= ligand2.radius[iatom_surface]
    end


    grid_LDS = zeros(complex(T), nx, ny, nz)
    score = zeros(T, nx, ny, nz, size(quaternions, 1))
    for iq in 1:size(quaternions, 1)
        ligand2_rotated = rotate(ligand2, quaternions[iq, :])
        grid_private .= 0.0
        Threads.@threads for iatom = 1:ligand2.natom
            tid = Threads.threadid()
            assign_desolvation_free_energy!(tid, grid_space, rcut1[iatom], rcut2[iatom], 
                                          ligand2_rotated.x[iframe, iatom], ligand2_rotated.y[iframe, iatom],
                                          ligand2_rotated.z[iframe, iatom], 
                                          x_grid, y_grid, z_grid, nx, ny, nz, grid_private)
        end
        grid_LDS .= dropdims(sum(grid_private, dims=4), dims=4)    
        #@show iq sum(grid_LDS)

        t = ifftshift(ifft(fft(grid_RDS) .* conj.(fft(grid_LDS))))
        score[:, :, :, iq] .= score[:, :, :, iq] .+ (imag(t) ./ 2)
    end
    
    score_tops = zeros(T, tops)
    trans_tops = zeros(Int64, tops, 3)
    quate_tops = zeros(Float64, tops, 4)
    x_center = ceil(Int64, (nx/2)+1.0)
    y_center = ceil(Int64, (ny/2)+1.0)
    z_center = ceil(Int64, (nz/2)+1.0)
     
    for t in 1:tops
        id = argmax(score)
        dx_estimated = id[1] - x_center
        dy_estimated = id[2] - y_center  
        dz_estimated = id[3] - z_center 
        
        score_tops[t] = score[id]
        trans_tops[t, :] = [dx_estimated, dy_estimated, dz_estimated]
        quate_tops[t, :] = quaternions[id[4], :]
        
        score[id] = -Inf
    end
        

    return (receptor = receptor2, ligand = ligand2, score = score_tops, trans = trans_tops, quate = quate_tops)
end