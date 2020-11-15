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

function compute_sasa(ta::TrjArray{T, U}, probe_radius=8.0::T; npoint=960::Int, iframe=1::Int) where {T, U}
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
                    grid_private[ix, iy, iz, thread_id] = 9.0 im
                elseif dist < rcut2
                    grid_private[ix, iy, iz, thread_id] = 1.0
                end
            end
        end
    end
end

function dock_fft(receptor::TrjArray{T, U}, ligand::TrjArray{T, U}, quaternions; grid_space=1.2, iframe=1) where {T, U}
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

    nx = length(x_grid)
    ny = length(y_grid)
    nz = length(z_grid)
    ngrid = nx*ny*nz
    grid = zeros(T, nx, ny, nz)

    x_gpoints = zeros(T, ngrid, 1)
    y_gpoints = zeros(T, ngrid, 1)
    z_gpoints = zeros(T, ngrid, 1)
    for iz = 1:nz, iy = 1:ny, ix = 1:nx
        ipoint = ix + (iy-1)*nx + (iz-1)*nx*ny
        x_gpoints[ipoint] = x_grid[ix]
        y_gpoints[ipoint] = y_grid[iy]
        z_gpoints[ipoint] = z_grid[iz]
    end

    x_receptor = zeros(T, 1, receptor2.natom)
    x_receptor[1, :] .= receptor2.x[iframe, :]

    y_receptor = zeros(T, 1, receptor2.natom)
    y_receptor[1, :] .= receptor2.y[iframe, :]

    z_receptor = zeros(T, 1, receptor2.natom)
    z_receptor[1, :] .= receptor2.z[iframe, :]

    radius_receptor = zeros(T, 1, receptor2.natom)
    radius_receptor[1, :] .= receptor2.radius

    sasa_receptor = zeros(T, 1, receptor2.natom)
    sasa_receptor[1, :] .= receptor2.mass

    # spape complementarity
    d = zeros(T, ngrid, receptor2.natom)
    d .= (x_gpoints .- x_receptor).^2 .+ (y_gpoints .- y_receptor).^2 .+ (z_gpoints .- z_receptor).^2

    iatom_surface = sasa_receptor .> 1.0
    iatom_core = .!iatom_surface
    
    grid_RSC = zeros(complex(T), nx, ny, nz)
    id = iatom_core .& (d .< (radius_receptor .* sqrt(1.5) ))
    #grid_RSC[ ] .= 1.0


    # generate grid coordinates for ligand with quaternions

        # convolution
    
    return id
end
