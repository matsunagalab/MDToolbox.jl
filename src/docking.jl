################ docking score parameters
function get_acescore()
    ace_score = Array{Float64}(undef, 18)
    ace_score[1] = -0.495 # ATOM TYPE "N"
    ace_score[2] = -0.553 # ATOM TYPE "CA"
    ace_score[3] = -0.464 # ATOM TYPE "C"
    ace_score[4] = -0.079 # ATOM TYPE "O"
    ace_score[5] = 0.008 # ATOM TYPE "GCA"
    ace_score[6] = -0.353 # ATOM TYPE "CB"
    ace_score[7] = 1.334 # ATOM TYPE "KNZ"
    ace_score[8] = 1.046 # ATOM TYPE "KCD"
    ace_score[9] = 0.933 # ATOM TYPE "DOD"
    ace_score[10] = 0.726 # ATOM TYPE "RNH"
    ace_score[11] = 0.693 # ATOM TYPE "NND"
    ace_score[12] = 0.606 # ATOM TYPE "RNE"
    ace_score[13] = 0.232 # ATOM TYPE "SOG"
    ace_score[14] = 0.061 # ATOM TYPE "HNE"
    ace_score[15] = -0.289 # ATOM TYPE "YCZ"
    ace_score[16] = -0.432 # ATOM TYPE "FCZ"
    ace_score[17] = -0.987 # ATOM TYPE "LCD"
    ace_score[18] = -1.827 # ATOM TYPE "CSG"

    return ace_score
end

function set_atomtype_id(ta::TrjArray{T,U}) where {T,U}
    atomtype_id = Array{Int64}(undef, ta.natom)

    for iatom = 1:ta.natom
        # ATOM TYPE "N"
        if ta.atomname[iatom] == "N"
            atomtype_id[iatom] = 1

            # ATOM TYPE "C"
        elseif ta.atomname[iatom] == "C"
            atomtype_id[iatom] = 3

            # ATOM TYPE "O"
        elseif ta.atomname[iatom] == "O" || ta.atomname[iatom] == "OXT"
            atomtype_id[iatom] = 4

            # ATOM TYPE "GCA"
        elseif ta.resname[iatom] == "GLY" && ta.atomname[iatom] == "CA"
            atomtype_id[iatom] = 5

            # ATOM TYPE "CA"
        elseif ta.atomname[iatom] == "CA"
            atomtype_id[iatom] = 2

            # ATOM TYPE "CB"
        elseif ta.resname[iatom] == "ALA" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "ARG" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "ASN" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "ASP" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "CYS" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "GLN" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "GLU" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "HIS" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "ILE" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "LEU" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "LYS" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "MET" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "PHE" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "PRO" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "PRO" && ta.atomname[iatom] == "CG"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "PRO" && ta.atomname[iatom] == "CD"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "THR" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "TRP" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "TYR" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6
        elseif ta.resname[iatom] == "VAL" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 6

            # ATOM TYPE "KNZ"
        elseif ta.resname[iatom] == "LYS" && ta.atomname[iatom] == "CE"
            atomtype_id[iatom] = 7
        elseif ta.resname[iatom] == "LYS" && ta.atomname[iatom] == "NZ"
            atomtype_id[iatom] = 7

            # ATOM TYPE "KCD"
        elseif ta.resname[iatom] == "LYS" && ta.atomname[iatom] == "CD"
            atomtype_id[iatom] = 8

            # ATOM TYPE "DOD"
        elseif ta.resname[iatom] == "ASP" && ta.atomname[iatom] == "CG"
            atomtype_id[iatom] = 9
        elseif ta.resname[iatom] == "ASP" && ta.atomname[iatom] == "OD1"
            atomtype_id[iatom] = 9
        elseif ta.resname[iatom] == "ASP" && ta.atomname[iatom] == "OD2"
            atomtype_id[iatom] = 9
        elseif ta.resname[iatom] == "GLU" && ta.atomname[iatom] == "CD"
            atomtype_id[iatom] = 9
        elseif ta.resname[iatom] == "GLU" && ta.atomname[iatom] == "OE1"
            atomtype_id[iatom] = 9
        elseif ta.resname[iatom] == "GLU" && ta.atomname[iatom] == "OE2"
            atomtype_id[iatom] = 9

            # ATOM TYPE "RNH"
        elseif ta.resname[iatom] == "ARG" && ta.atomname[iatom] == "CZ"
            atomtype_id[iatom] = 10
        elseif ta.resname[iatom] == "ARG" && ta.atomname[iatom] == "NH1"
            atomtype_id[iatom] = 10
        elseif ta.resname[iatom] == "ARG" && ta.atomname[iatom] == "NH2"
            atomtype_id[iatom] = 10

            # ATOM TYPE "NND"
        elseif ta.resname[iatom] == "ASN" && ta.atomname[iatom] == "CG"
            atomtype_id[iatom] = 11
        elseif ta.resname[iatom] == "ASN" && ta.atomname[iatom] == "OD1"
            atomtype_id[iatom] = 11
        elseif ta.resname[iatom] == "ASN" && ta.atomname[iatom] == "ND2"
            atomtype_id[iatom] = 11
        elseif ta.resname[iatom] == "GLN" && ta.atomname[iatom] == "CD"
            atomtype_id[iatom] = 11
        elseif ta.resname[iatom] == "GLN" && ta.atomname[iatom] == "OE1"
            atomtype_id[iatom] = 11
        elseif ta.resname[iatom] == "GLN" && ta.atomname[iatom] == "NE2"
            atomtype_id[iatom] = 11

            # ATOM TYPE "RNE"
        elseif ta.resname[iatom] == "ARG" && ta.atomname[iatom] == "CD"
            atomtype_id[iatom] = 12
        elseif ta.resname[iatom] == "ARG" && ta.atomname[iatom] == "NE"
            atomtype_id[iatom] = 12

            # ATOM TYPE "SOG"
        elseif ta.resname[iatom] == "SER" && ta.atomname[iatom] == "CB"
            atomtype_id[iatom] = 13
        elseif ta.resname[iatom] == "SER" && ta.atomname[iatom] == "OG"
            atomtype_id[iatom] = 13
        elseif ta.resname[iatom] == "THR" && ta.atomname[iatom] == "OG1"
            atomtype_id[iatom] = 13
        elseif ta.resname[iatom] == "TYR" && ta.atomname[iatom] == "OH"
            atomtype_id[iatom] = 13

            # ATOM TYPE "HNE"
        elseif ta.resname[iatom] == "HIS" && ta.atomname[iatom] == "CG"
            atomtype_id[iatom] = 14
        elseif ta.resname[iatom] == "HIS" && ta.atomname[iatom] == "ND1"
            atomtype_id[iatom] = 14
        elseif ta.resname[iatom] == "HIS" && ta.atomname[iatom] == "CD2"
            atomtype_id[iatom] = 14
        elseif ta.resname[iatom] == "HIS" && ta.atomname[iatom] == "CE1"
            atomtype_id[iatom] = 14
        elseif ta.resname[iatom] == "HIS" && ta.atomname[iatom] == "NE2"
            atomtype_id[iatom] = 14
        elseif ta.resname[iatom] == "TRP" && ta.atomname[iatom] == "NE1"
            atomtype_id[iatom] = 14

            # ATOM TYPE "YCZ"
        elseif ta.resname[iatom] == "TYR" && ta.atomname[iatom] == "CE1"
            atomtype_id[iatom] = 15
        elseif ta.resname[iatom] == "TYR" && ta.atomname[iatom] == "CE2"
            atomtype_id[iatom] = 15
        elseif ta.resname[iatom] == "TYR" && ta.atomname[iatom] == "CZ"
            atomtype_id[iatom] = 15

            # ATOM TYPE "FCZ"
        elseif ta.resname[iatom] == "ARG" && ta.atomname[iatom] == "CG"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "GLN" && ta.atomname[iatom] == "CG"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "GLU" && ta.atomname[iatom] == "CG"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "ILE" && ta.atomname[iatom] == "CG1"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "LEU" && ta.atomname[iatom] == "CG"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "LYS" && ta.atomname[iatom] == "CG"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "MET" && ta.atomname[iatom] == "CG"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "MET" && ta.atomname[iatom] == "SD"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "PHE" && ta.atomname[iatom] == "CG"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "PHE" && ta.atomname[iatom] == "CD1"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "PHE" && ta.atomname[iatom] == "CD2"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "PHE" && ta.atomname[iatom] == "CE1"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "PHE" && ta.atomname[iatom] == "CE2"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "PHE" && ta.atomname[iatom] == "CZ"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "THR" && ta.atomname[iatom] == "CG2"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "TRP" && ta.atomname[iatom] == "CG"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "TRP" && ta.atomname[iatom] == "CD1"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "TRP" && ta.atomname[iatom] == "CD2"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "TRP" && ta.atomname[iatom] == "CE2"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "TRP" && ta.atomname[iatom] == "CE3"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "TRP" && ta.atomname[iatom] == "CZ2"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "TRP" && ta.atomname[iatom] == "CZ3"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "TRP" && ta.atomname[iatom] == "CH2"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "TYR" && ta.atomname[iatom] == "CG"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "TYR" && ta.atomname[iatom] == "CD1"
            atomtype_id[iatom] = 16
        elseif ta.resname[iatom] == "TYR" && ta.atomname[iatom] == "CD2"
            atomtype_id[iatom] = 16

            # ATOM TYPE "LCD"
        elseif ta.resname[iatom] == "ILE" && ta.atomname[iatom] == "CG2"
            atomtype_id[iatom] = 17
        elseif ta.resname[iatom] == "ILE" && ta.atomname[iatom] == "CD"
            atomtype_id[iatom] = 17
        elseif ta.resname[iatom] == "ILE" && ta.atomname[iatom] == "CD1"
            atomtype_id[iatom] = 17
        elseif ta.resname[iatom] == "LEU" && ta.atomname[iatom] == "CD1"
            atomtype_id[iatom] = 17
        elseif ta.resname[iatom] == "LEU" && ta.atomname[iatom] == "CD2"
            atomtype_id[iatom] = 17
        elseif ta.resname[iatom] == "MET" && ta.atomname[iatom] == "CE"
            atomtype_id[iatom] = 17
        elseif ta.resname[iatom] == "VAL" && ta.atomname[iatom] == "CG1"
            atomtype_id[iatom] = 17
        elseif ta.resname[iatom] == "VAL" && ta.atomname[iatom] == "CG2"
            atomtype_id[iatom] = 17

            # ATOM TYPE "CSG"
        elseif ta.resname[iatom] == "CYS" && ta.atomname[iatom] == "SG"
            atomtype_id[iatom] = 18
        else
            error("error: faled to assign atom type " * ta.resname[iatom] * "-" * ta.atomname[iatom])
        end
    end

    return TrjArray(ta, atomtype_id=atomtype_id)
end

function set_radius(ta::TrjArray{T,U}) where {T,U}
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

function set_charge(ta)
    charge = Array{Float64}(undef, ta.natom)

    is_first = true
    for iatom = 1:ta.natom
        # ATOM TYPE "N"
        if ta.atomname[iatom] == "N"
            if is_first
                charge[iatom] = 1.0
                is_first = false
            else
                charge[iatom] = 0.5
            end

            # ATOM TYPE "O"
        elseif ta.atomname[iatom] == "O"
            charge[iatom] = -0.5
        elseif ta.atomname[iatom] == "OXT"
            charge[iatom] = -1.0

        elseif ta.resname[iatom] == "ARG" && ta.atomname[iatom] == "NH1"
            charge[iatom] = 0.5
        elseif ta.resname[iatom] == "ARG" && ta.atomname[iatom] == "NH2"
            charge[iatom] = 0.5
        elseif ta.resname[iatom] == "GLU" && ta.atomname[iatom] == "OE1"
            charge[iatom] = -0.5
        elseif ta.resname[iatom] == "GLU" && ta.atomname[iatom] == "OE2"
            charge[iatom] = -0.5
        elseif ta.resname[iatom] == "ASP" && ta.atomname[iatom] == "OD1"
            charge[iatom] = -0.5
        elseif ta.resname[iatom] == "ASP" && ta.atomname[iatom] == "OD2"
            charge[iatom] = -0.5
        elseif ta.resname[iatom] == "LYS" && ta.atomname[iatom] == "NZ"
            charge[iatom] = 1.0
        elseif ta.resname[iatom] == "PRO" && ta.atomname[iatom] == "N"
            charge[iatom] = -0.1

        else
            charge[iatom] = 0.0
        end
    end

    return TrjArray(ta, charge=charge)
end

################ solvent accessible surface area
function golden_section_spiral(n)
    points = zeros(Float64, n, 3)
    inc = pi * (3.0 - sqrt(5.0))
    offset = 2.0 / Float64(n)
    for k = 1:n
        y = (k - 1) * offset - 1.0 + (offset / 2.0)
        r = sqrt(1.0 - y * y)
        phi = (k - 1) * inc
        points[k, 1] = cos(phi) * r
        points[k, 2] = y
        points[k, 3] = sin(phi) * r
    end
    return points
end

function compute_sasa(ta::TrjArray{T,U}, probe_radius=1.4::T; npoint=960::Int, iframe=1::Int, candicate=10) where {T,U}
    # construct pair rist
    if isempty(ta.radius)
        error("radius is empty.")
    end
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
    for iatom = 1:ta.natom
        n_accessible_point = 0
        neighbor_list_iatom = neighbor_list[iatom]
        for ipoint in 1:npoint
            is_accessible = true
            point = points[ipoint, :] .* (ta.radius[iatom] + probe_radius)
            point[1] += ta.xyz[iframe, 3*(iatom-1)+1]
            point[2] += ta.xyz[iframe, 3*(iatom-1)+2]
            point[3] += ta.xyz[iframe, 3*(iatom-1)+3]
            for j in 1:length(neighbor_list_iatom)
                jatom = neighbor_list_iatom[j]
                dx = point[1] - ta.xyz[iframe, 3*(jatom-1)+1]
                dy = point[2] - ta.xyz[iframe, 3*(jatom-1)+2]
                dz = point[3] - ta.xyz[iframe, 3*(jatom-1)+3]
                if !isempty(ta.boxsize)
                    dx = dx - round(dx / ta.boxsize[iframe, 1]) * ta.boxsize[iframe, 1]
                    dy = dy - round(dy / ta.boxsize[iframe, 2]) * ta.boxsize[iframe, 2]
                    dz = dz - round(dz / ta.boxsize[iframe, 3]) * ta.boxsize[iframe, 3]
                end
                d = sqrt(dx^2 + dy^2 + dz^2)
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

    return TrjArray(ta, sasa=sasa)
end

################ grid
function generate_grid(receptor_org::TrjArray{T,U}, ligand_org::TrjArray{T,U}; iframe=1, spacing=1.2) where {T,U}
    receptor = deepcopy(receptor_org)
    ligand = deepcopy(ligand_org)
    decenter!(receptor)
    orient!(ligand)
    xmin_ligand = minimum(ligand.xyz[iframe, 1:3:end])
    xmax_ligand = maximum(ligand.xyz[iframe, 1:3:end])
    size_ligand = xmax_ligand - xmin_ligand

    xmin_receptor = minimum(receptor.xyz[iframe, 1:3:end])
    xmax_receptor = maximum(receptor.xyz[iframe, 1:3:end])
    xmin_grid = xmin_receptor - size_ligand - spacing
    xmax_grid = xmax_receptor + size_ligand + spacing

    ymin_receptor = minimum(receptor.xyz[iframe, 2:3:end])
    ymax_receptor = maximum(receptor.xyz[iframe, 2:3:end])
    ymin_grid = ymin_receptor - size_ligand - spacing
    ymax_grid = ymax_receptor + size_ligand + spacing

    zmin_receptor = minimum(receptor.xyz[iframe, 3:3:end])
    zmax_receptor = maximum(receptor.xyz[iframe, 3:3:end])
    zmin_grid = zmin_receptor - size_ligand - spacing
    zmax_grid = zmax_receptor + size_ligand + spacing

    x_grid = Array{T,1}(range(xmin_grid, xmax_grid, step=spacing))
    if typeof(receptor_org.xyz) <: CuArray
        x_grid = CuArray(x_grid)
    end

    y_grid = Array{T,1}(range(ymin_grid, ymax_grid, step=spacing))
    if typeof(receptor_org.xyz) <: CuArray
        y_grid = CuArray(y_grid)
    end

    z_grid = Array{T,1}(range(zmin_grid, zmax_grid, step=spacing))
    if typeof(receptor_org.xyz) <: CuArray
        z_grid = CuArray(z_grid)
    end

    nx = length(x_grid)
    ny = length(y_grid)
    nz = length(z_grid)

    grid_real = similar(receptor_org.xyz, (nx, ny, nz))
    grid_real .= zero(T)
    grid_imag = similar(receptor_org.xyz, (nx, ny, nz))
    grid_imag .= zero(T)

    return grid_real, grid_imag, x_grid, y_grid, z_grid
end

function spread_nearest_add!(grid::AbstractArray{T},
    x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T},
    x_grid::AbstractVector{T}, y_grid::AbstractVector{T}, z_grid::AbstractVector{T},
    weight::AbstractVector{T}) where {T}

    natom = length(x)
    x_grid_delta = x_grid[2] - x_grid[1]
    y_grid_delta = y_grid[2] - y_grid[1]
    z_grid_delta = z_grid[2] - z_grid[1]

    x_grid_min = x_grid[1]
    y_grid_min = y_grid[1]
    z_grid_min = z_grid[1]

    for iatom = 1:natom
        ix = ceil(Int, (x[iatom] - x_grid_min) / x_grid_delta)
        iy = ceil(Int, (y[iatom] - y_grid_min) / y_grid_delta)
        iz = ceil(Int, (z[iatom] - z_grid_min) / z_grid_delta)
        grid[ix, iy, iz] += weight[iatom]
    end

    return nothing
end

function spread_nearest_add!(grid::CuArray{T},
    x::CuArray{T}, y::CuArray{T}, z::CuArray{T},
    x_grid::CuArray{T}, y_grid::CuArray{T}, z_grid::CuArray{T},
    weight::CuArray{T}) where {T}

    natom = length(x)
    nthreads = 256
    @cuda blocks = ceil(Int, natom / nthreads) threads = nthreads spread_nearest_add_kernel!(grid, x, y, z, x_grid, y_grid, z_grid, weight)

    return nothing
end

function spread_nearest_add_kernel!(grid::CuDeviceArray{T},
    x::CuDeviceVector{T}, y::CuDeviceVector{T}, z::CuDeviceVector{T},
    x_grid::CuDeviceVector{T}, y_grid::CuDeviceVector{T}, z_grid::CuDeviceVector{T},
    weight::CuDeviceVector{T}) where {T}

    natom = length(x)
    tid = threadIdx().x
    gtid = (blockIdx().x - 1) * blockDim().x + tid  # global thread id

    x_grid_delta = x_grid[2] - x_grid[1]
    y_grid_delta = y_grid[2] - y_grid[1]
    z_grid_delta = z_grid[2] - z_grid[1]

    x_grid_min = x_grid[1]
    y_grid_min = y_grid[1]
    z_grid_min = z_grid[1]

    iatom = gtid
    if iatom <= natom
        ix = ceil(Int, (x[iatom] - x_grid_min) / x_grid_delta)
        iy = ceil(Int, (y[iatom] - y_grid_min) / y_grid_delta)
        iz = ceil(Int, (z[iatom] - z_grid_min) / z_grid_delta)
        CUDA.@atomic grid[ix, iy, iz] += weight[iatom]
    end

    return nothing
end

function spread_nearest_substitute!(grid::AbstractArray{T},
    x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T},
    x_grid::AbstractVector{T}, y_grid::AbstractVector{T}, z_grid::AbstractVector{T},
    weight::AbstractVector{T}) where {T}

    natom = length(x)
    x_grid_delta = x_grid[2] - x_grid[1]
    y_grid_delta = y_grid[2] - y_grid[1]
    z_grid_delta = z_grid[2] - z_grid[1]

    x_grid_min = x_grid[1]
    y_grid_min = y_grid[1]
    z_grid_min = z_grid[1]

    for iatom = 1:natom
        ix = ceil(Int, (x[iatom] - x_grid_min) / x_grid_delta)
        iy = ceil(Int, (y[iatom] - y_grid_min) / y_grid_delta)
        iz = ceil(Int, (z[iatom] - z_grid_min) / z_grid_delta)
        grid[ix, iy, iz] = weight[iatom]
    end

    return nothing
end

function spread_nearest_substitute!(grid::CuArray{T},
    x::CuArray{T}, y::CuArray{T}, z::CuArray{T},
    x_grid::CuArray{T}, y_grid::CuArray{T}, z_grid::CuArray{T},
    weight::CuArray{T}) where {T}

    natom = length(x)
    nthreads = 256
    @cuda blocks = ceil(Int, natom / nthreads) threads = nthreads spread_nearest_substitute_kernel!(grid, x, y, z, x_grid, y_grid, z_grid, weight)

    return nothing
end

function spread_nearest_substitute_kernel!(grid::CuDeviceArray{T},
    x::CuDeviceVector{T}, y::CuDeviceVector{T}, z::CuDeviceVector{T},
    x_grid::CuDeviceVector{T}, y_grid::CuDeviceVector{T}, z_grid::CuDeviceVector{T},
    weight::CuDeviceVector{T}) where {T}

    natom = length(x)
    tid = threadIdx().x
    gtid = (blockIdx().x - 1) * blockDim().x + tid  # global thread id

    x_grid_delta = x_grid[2] - x_grid[1]
    y_grid_delta = y_grid[2] - y_grid[1]
    z_grid_delta = z_grid[2] - z_grid[1]

    x_grid_min = x_grid[1]
    y_grid_min = y_grid[1]
    z_grid_min = z_grid[1]

    iatom = gtid
    if iatom <= natom
        ix = ceil(Int, (x[iatom] - x_grid_min) / x_grid_delta)
        iy = ceil(Int, (y[iatom] - y_grid_min) / y_grid_delta)
        iz = ceil(Int, (z[iatom] - z_grid_min) / z_grid_delta)
        grid[ix, iy, iz] = weight[iatom]
    end

    return nothing
end

function spread_neighbors_add!(grid::AbstractArray{T},
    x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T},
    x_grid::AbstractVector{T}, y_grid::AbstractVector{T}, z_grid::AbstractVector{T},
    weight::AbstractVector{T}, rcut::AbstractVector{T}) where {T}

    natom = length(x)
    nx, ny, nz = size(grid)
    for iatom = 1:natom
        for ix = 1:nx
            dx = x[iatom] - x_grid[ix]
            if abs(dx) > rcut[iatom]
                continue
            end
            for iy = 1:ny
                dy = y[iatom] - y_grid[iy]
                if abs(dy) > rcut[iatom]
                    continue
                end
                for iz = 1:nz
                    dz = z[iatom] - z_grid[iz]
                    if abs(dz) > rcut[iatom]
                        continue
                    end
                    d = dx * dx + dy * dy + dz * dz
                    if d < rcut[iatom] * rcut[iatom]
                        grid[ix, iy, iz] += weight[iatom]
                    end
                end
            end
        end
    end

    return nothing
end

function spread_neighbors_add!(grid::CuArray{T},
    x::CuArray{T}, y::CuArray{T}, z::CuArray{T},
    x_grid::CuArray{T}, y_grid::CuArray{T}, z_grid::CuArray{T},
    weight::CuArray{T}, rcut::CuArray{T}) where {T}

    natom = length(x)
    nthreads = 256
    @cuda blocks = ceil(Int, natom / nthreads) threads = nthreads spread_neighbors_add_kernel!(grid, x, y, z, x_grid, y_grid, z_grid, weight, rcut)

    return nothing
end

function spread_neighbors_add_kernel!(grid::CuDeviceArray{T},
    x::CuDeviceVector{T}, y::CuDeviceVector{T}, z::CuDeviceVector{T},
    x_grid::CuDeviceVector{T}, y_grid::CuDeviceVector{T}, z_grid::CuDeviceVector{T},
    weight::CuDeviceVector{T}, rcut::CuDeviceVector{T}) where {T}

    natom = length(x)
    tid = threadIdx().x
    gtid = (blockIdx().x - 1) * blockDim().x + tid  # global thread id
    nx, ny, nz = size(grid)

    iatom = gtid
    if iatom <= natom
        for ix = 1:nx
            dx = x[iatom] - x_grid[ix]
            if abs(dx) > rcut[iatom]
                continue
            end
            for iy = 1:ny
                dy = y[iatom] - y_grid[iy]
                if abs(dy) > rcut[iatom]
                    continue
                end
                for iz = 1:nz
                    dz = z[iatom] - z_grid[iz]
                    if abs(dz) > rcut[iatom]
                        continue
                    end
                    d = dx * dx + dy * dy + dz * dz
                    if d < rcut[iatom] * rcut[iatom]
                        CUDA.@atomic grid[ix, iy, iz] += weight[iatom]
                    end
                end
            end
        end
    end

    return nothing
end

function spread_neighbors_substitute!(grid::AbstractArray{T2},
    x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T},
    x_grid::AbstractVector{T}, y_grid::AbstractVector{T}, z_grid::AbstractVector{T},
    weight::AbstractVector{T}, rcut::AbstractVector{T}) where {T2,T}

    natom = length(x)
    nx, ny, nz = size(grid)

    for iatom = 1:natom
        for ix = 1:nx
            dx = x[iatom] - x_grid[ix]
            if abs(dx) > rcut[iatom]
                continue
            end
            for iy = 1:ny
                dy = y[iatom] - y_grid[iy]
                if abs(dy) > rcut[iatom]
                    continue
                end
                for iz = 1:nz
                    dz = z[iatom] - z_grid[iz]
                    if abs(dz) > rcut[iatom]
                        continue
                    end
                    d = dx * dx + dy * dy + dz * dz
                    if d < rcut[iatom] * rcut[iatom]
                        grid[ix, iy, iz] = weight[iatom]
                    end
                end
            end
        end
    end

    return nothing
end

function spread_neighbors_substitute!(grid::CuArray{T},
    x::CuArray{T}, y::CuArray{T}, z::CuArray{T},
    x_grid::CuArray{T}, y_grid::CuArray{T}, z_grid::CuArray{T},
    weight::CuArray{T}, rcut::CuArray{T}) where {T}

    natom = length(x)
    nthreads = 256
    @cuda blocks = ceil(Int, natom / nthreads) threads = nthreads spread_neighbors_substitute_kernel!(grid, x, y, z, x_grid, y_grid, z_grid, weight, rcut)

    return nothing
end

function spread_neighbors_substitute_kernel!(grid::CuDeviceArray{T},
    x::CuDeviceVector{T}, y::CuDeviceVector{T}, z::CuDeviceVector{T},
    x_grid::CuDeviceVector{T}, y_grid::CuDeviceVector{T}, z_grid::CuDeviceVector{T},
    weight::CuDeviceVector{T}, rcut::CuDeviceVector{T}) where {T}
    natom = length(x)
    tid = threadIdx().x
    gtid = (blockIdx().x - 1) * blockDim().x + tid  # global thread id
    nx, ny, nz = size(grid)

    iatom = gtid
    if iatom <= natom
        for ix = 1:nx
            dx = x[iatom] - x_grid[ix]
            if abs(dx) > rcut[iatom]
                continue
            end
            for iy = 1:ny
                dy = y[iatom] - y_grid[iy]
                if abs(dy) > rcut[iatom]
                    continue
                end
                for iz = 1:nz
                    dz = z[iatom] - z_grid[iz]
                    if abs(dz) > rcut[iatom]
                        continue
                    end
                    d = dx * dx + dy * dy + dz * dz
                    if d < rcut[iatom] * rcut[iatom]
                        grid[ix, iy, iz] = weight[iatom]
                    end
                end
            end
        end
    end

    return nothing
end

function assign_sc_receptor!(grid_real::AbstractArray{T}, grid_imag::AbstractArray{T},
    x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T},
    x_grid::AbstractVector{T}, y_grid::AbstractVector{T}, z_grid::AbstractVector{T},
    radius::AbstractVector{T}, id_surface::AbstractVector) where {T}

    grid_real .= zero(T)
    grid_imag .= zero(T)

    x_s = x[id_surface]
    y_s = y[id_surface]
    z_s = z[id_surface]
    radius_s = radius[id_surface]
    weight_s = similar(radius_s)

    x_c = x[.!id_surface]
    y_c = y[.!id_surface]
    z_c = z[.!id_surface]
    radius_c = radius[.!id_surface]
    weight_c = similar(radius_c)

    weight_s .= T(1.0)
    spread_neighbors_substitute!(grid_real, x_s, y_s, z_s, x_grid, y_grid, z_grid, weight_s, radius_s .+ T(3.4))
    weight_s .= T(0.0)
    spread_neighbors_substitute!(grid_real, x_s, y_s, z_s, x_grid, y_grid, z_grid, weight_s, radius_s .* T(sqrt(0.8)))
    weight_c .= T(0.0)
    spread_neighbors_substitute!(grid_real, x_c, y_c, z_c, x_grid, y_grid, z_grid, weight_c, radius_c .* T(sqrt(1.5)))

    weight_s .= T(9.0)
    spread_neighbors_substitute!(grid_imag, x_s, y_s, z_s, x_grid, y_grid, z_grid, weight_s, radius_s .* T(sqrt(0.8)))
    weight_c .= T(9.0)
    spread_neighbors_substitute!(grid_imag, x_c, y_c, z_c, x_grid, y_grid, z_grid, weight_c, radius_c .* T(sqrt(1.5)))

    return nothing
end

function assign_sc_ligand!(grid_real::AbstractArray{T}, grid_imag::AbstractArray{T},
    x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T},
    x_grid::AbstractVector{T}, y_grid::AbstractVector{T}, z_grid::AbstractVector{T},
    radius::AbstractVector{T}, id_surface::AbstractVector) where {T}

    grid_real .= zero(T)
    grid_imag .= zero(T)

    x_s = x[id_surface]
    y_s = y[id_surface]
    z_s = z[id_surface]
    radius_s = radius[id_surface]
    weight_s = similar(radius_s)

    x_c = x[.!id_surface]
    y_c = y[.!id_surface]
    z_c = z[.!id_surface]
    radius_c = radius[.!id_surface]
    weight_c = similar(radius_c)

    weight_s .= T(1.0)
    spread_neighbors_substitute!(grid_real, x_s, y_s, z_s, x_grid, y_grid, z_grid, weight_s, radius_s)

    weight_c .= T(9.0)
    spread_neighbors_substitute!(grid_imag, x_c, y_c, z_c, x_grid, y_grid, z_grid, weight_c, radius_c .* T(sqrt(1.5)))

    return nothing
end

function assign_ds!(grid_real::AbstractArray{T}, grid_imag::AbstractArray{T},
    x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T},
    x_grid::AbstractVector{T}, y_grid::AbstractVector{T}, z_grid::AbstractVector{T},
    ace_score::AbstractVector{T}) where {T}

    grid_real .= zero(T)
    grid_imag .= zero(T)

    radius = similar(ace_score)
    radius .= T(6.0)
    spread_neighbors_add!(grid_real, x, y, z, x_grid, y_grid, z_grid, ace_score, radius)

    radius .= T(1.0)
    spread_nearest_substitute!(grid_imag, x, y, z, x_grid, y_grid, z_grid, radius)

    return nothing
end

################ rotate
function rotate!(x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T},
    q::AbstractVector{T}) where {T}
    natom = length(x)
    r1 = 1.0 - 2.0 * q[2] * q[2] - 2.0 * q[3] * q[3]
    r2 = 2.0 * (q[1] * q[2] + q[3] * q[4])
    r3 = 2.0 * (q[1] * q[3] - q[2] * q[4])
    r4 = 2.0 * (q[1] * q[2] - q[3] * q[4])
    r5 = 1.0 - 2.0 * q[1] * q[1] - 2.0 * q[3] * q[3]
    r6 = 2.0 * (q[2] * q[3] + q[1] * q[4])
    r7 = 2.0 * (q[1] * q[3] + q[2] * q[4])
    r8 = 2.0 * (q[2] * q[3] - q[1] * q[4])
    r9 = 1.0 - 2.0 * q[1] * q[1] - 2.0 * q[2] * q[2]
    for iatom = 1:natom
        x_new = r1 * x[iatom] + r2 * y[iatom] + r3 * z[iatom]
        y_new = r4 * x[iatom] + r5 * y[iatom] + r6 * z[iatom]
        z_new = r7 * x[iatom] + r8 * y[iatom] + r9 * z[iatom]
        x[iatom] = x_new
        y[iatom] = y_new
        z[iatom] = z_new
    end
    return nothing
end

function rotate!(x::CuVector{T}, y::CuVector{T}, z::CuVector{T},
                 q::CuVector{T}) where {T}
    natom = length(x)
    nthreads = 256
    @cuda blocks = ceil(Int, natom / nthreads) threads = nthreads rotate_kernel!(x, y, z, q)
    return nothing
end

function rotate_kernel!(x::CuDeviceVector{T}, y::CuDeviceVector{T}, z::CuDeviceVector{T},
                        q::CuDeviceVector{T}) where {T}
    natom = length(x)
    tid = threadIdx().x
    gtid = (blockIdx().x - 1) * blockDim().x + tid  # global thread id

    r1 = 1.0 - 2.0 * q[2] * q[2] - 2.0 * q[3] * q[3]
    r2 = 2.0 * (q[1] * q[2] + q[3] * q[4])
    r3 = 2.0 * (q[1] * q[3] - q[2] * q[4])
    r4 = 2.0 * (q[1] * q[2] - q[3] * q[4])
    r5 = 1.0 - 2.0 * q[1] * q[1] - 2.0 * q[3] * q[3]
    r6 = 2.0 * (q[2] * q[3] + q[1] * q[4])
    r7 = 2.0 * (q[1] * q[3] + q[2] * q[4])
    r8 = 2.0 * (q[2] * q[3] - q[1] * q[4])
    r9 = 1.0 - 2.0 * q[1] * q[1] - 2.0 * q[2] * q[2]

    iatom = gtid
    if iatom <= natom
        x_new = r1 * x[iatom] + r2 * y[iatom] + r3 * z[iatom]
        y_new = r4 * x[iatom] + r5 * y[iatom] + r6 * z[iatom]
        z_new = r7 * x[iatom] + r8 * y[iatom] + r9 * z[iatom]
        x[iatom] = x_new
        y[iatom] = y_new
        z[iatom] = z_new
    end

    return nothing
end

################ docking
function compute_docking_score_with_fft(quaternion, grid_RSC, grid_LSC, ligand2, grid_space, rcut1, rcut2, x_grid, y_grid, z_grid, iframe, tops, iq)
    ligand2_rotated = rotate(ligand2, quaternion)
    assign_shape_complementarity!(grid_LSC, ligand2_rotated, grid_space, rcut1, rcut2, x_grid, y_grid, z_grid, iframe)
    #grid_LSC .= grid_LSC[end:-1:1, end:-1:1, end:-1:1]

    if CUDA.functional()
        grid_RSC_gpu = cu(grid_RSC)
        grid_LSC_gpu = cu(grid_LSC)
        t_gpu = ifft(fft(grid_RSC_gpu) .* conj.(fft(conj.(grid_LSC_gpu))))
        #score_gpu = real(t_gpu) .- imag(t_gpu)
        score_gpu = real(t_gpu)
        score = Array(score_gpu)
    else
        t = ifft(fft(grid_RSC) .* conj.(fft(conj.(grid_LSC))))
        #score = real(t) .- imag(t)
        score = real(t)
    end

    ret = []
    for t in 1:tops
        id = argmax(score)
        dx_estimated = id[1]
        dy_estimated = id[2]
        dz_estimated = id[3]
        push!(ret, (score[id], dx_estimated, dy_estimated, dz_estimated, iq))
        score[id] = -Inf
    end

    return ret
end

function filter_tops!(score_tops, cartesian_tops, iq_tops, score, iq, tops)
    if any(score .< score_tops[tops+1])
        c = findall(score .< score_tops[tops+1])
        s = score[c]
        if isempty(s)
            return
        end
        nrows = length(s)
        nrows = min(nrows, tops)
        score_tops[(tops+1):(tops+nrows)] .= Array(s[1:nrows])
        cartesian_tops[(tops+1):(tops+nrows)] .= Array(c[1:nrows])
        iq_tops[(tops+1):(tops+nrows)] .= iq
        id = sortperm(score_tops)
        score_tops .= score_tops[id]
        cartesian_tops .= cartesian_tops[id]
        iq_tops .= iq_tops[id]
    end
    return
end


function filter_tops!(score_tops, cartesian_tops, rotz_tops, rotx_tops, score, rotz, rotx, tops)
    nz = size(score, 3)
    iz_center = ceil(Int, (nz / 2.0) + 1.0)
    score_xy = score[:, :, iz_center]
    if any(score_xy .< score_tops[tops+1])
        c = findall(score_xy .< score_tops[tops+1])
        s = score_xy[c]
        if isempty(s)
            return
        end
        nrows = length(s)
        nrows = min(nrows, tops)
        score_tops[(tops+1):(tops+nrows)] .= Array(s[1:nrows])
        cartesian_tops[(tops+1):(tops+nrows)] .= Array(c[1:nrows])
        rotz_tops[(tops+1):(tops+nrows)] .= rotz
        rotx_tops[(tops+1):(tops+nrows)] .= rotx
        id = sortperm(score_tops)
        score_tops .= score_tops[id]
        cartesian_tops .= cartesian_tops[id]
        rotz_tops .= rotz_tops[id]
        rotx_tops .= rotx_tops[id]
    end
    return
end


function dock!(receptor::TrjArray{T,U}, ligand::TrjArray{T,U}, quaternions::Matrix{T}; deg=15.0, grid_space=1.2, iframe=1, tops=100, alpha=0.01, beta=0.06) where {T,U}
    decenter!(receptor)
    decenter!(ligand)
    orient!(receptor)
    orient!(ligand)

    # Assign atom radius
    println("step1: assigning atom radius")
    receptor = set_radius(receptor)
    ligand = set_radius(ligand)

    # SASA
    println("step2: computing SASA")
    receptor = compute_sasa(receptor, 1.4)
    ligand = compute_sasa(ligand, 1.4)

    # ACE scores
    println("step3: assigning ACE scores")
    receptor = set_acescore(receptor)
    ligand = set_acescore(ligand)

    # Determine grid size and coordinates
    println("step4: computing grid size and coordinates")
    size_ligand = maximum(ligand.xyz[iframe, 1:3:end]) - minimum(ligand.xyz[iframe, 1:3:end])
    #size_ligand = size_ligand*2

    x_min = minimum(receptor.xyz[iframe, 1:3:end]) - size_ligand - grid_space
    y_min = minimum(receptor.xyz[iframe, 2:3:end]) - size_ligand - grid_space
    z_min = minimum(receptor.xyz[iframe, 3:3:end]) - size_ligand - grid_space

    x_max = maximum(receptor.xyz[iframe, 1:3:end]) + size_ligand + grid_space
    y_max = maximum(receptor.xyz[iframe, 2:3:end]) + size_ligand + grid_space
    z_max = maximum(receptor.xyz[iframe, 3:3:end]) + size_ligand + grid_space

    x_grid = collect(x_min:grid_space:x_max)
    y_grid = collect(y_min:grid_space:y_max)
    z_grid = collect(z_min:grid_space:z_max)

    nx, ny, nz = length(x_grid), length(y_grid), length(z_grid)
    nxyz = nx * ny * nz

    println("step5: assigning values to grid points")
    # receptor grid: shape complementarity
    grid_RSC = zeros(Complex{T}, (nx, ny, nz))

    id_surface = receptor.sasa .> 1.0
    id_core = .!id_surface

    x = receptor.xyz[iframe, 1:3:end]
    y = receptor.xyz[iframe, 2:3:end]
    z = receptor.xyz[iframe, 3:3:end]
    rcut = zeros(T, receptor.natom)
    rcut[id_core] .= receptor.radius[id_core] .* sqrt(1.5)
    rcut[id_surface] .= receptor.radius[id_surface] .* sqrt(0.8)
    value_core = fill(Complex{T}(9im), receptor.natom)

    x_surface = x[id_surface]
    y_surface = y[id_surface]
    z_surface = z[id_surface]
    rcut_surface = receptor.radius[id_surface] .+ 3.4
    value_surface = fill(Complex{T}(1.0), length(rcut_surface))

    spread_neighbors_substitute!(grid_RSC, x_surface, y_surface, z_surface, x_grid, y_grid, z_grid, rcut_surface, value_surface)
    spread_neighbors_substitute!(grid_RSC, x, y, z, x_grid, y_grid, z_grid, rcut, value_core)

    # receptor grid: desolvation free energy
    grid_RDS = zeros(Complex{T}, (nx, ny, nz))
    rcut_ds = fill(T(6.0), receptor.natom)
    spread_nearest!(grid_RDS, x, y, z, x_grid, y_grid, z_grid)
    spread_neighbors_add!(grid_RDS, x, y, z, x_grid, y_grid, z_grid, rcut_ds, receptor.mass)

    # ligand grid: shape complementarity
    grid_LSC = zeros(Complex{T}, (nx, ny, nz))

    id_surface = ligand.sasa .> 1.0
    id_core = .!id_surface

    x = ligand.xyz[iframe, 1:3:end]
    y = ligand.xyz[iframe, 2:3:end]
    z = ligand.xyz[iframe, 3:3:end]
    rcut = zeros(T, ligand.natom)
    rcut[id_core] .= ligand.radius[id_core] .* sqrt(1.5)
    rcut[id_surface] .= ligand.radius[id_surface] .* sqrt(0.8)
    value_core = fill(Complex{T}(9im), ligand.natom)
    rcut_surface = ligand.radius[id_surface]
    value_surface = fill(Complex{T}(1.0), length(rcut_surface))

    x_surface = x[id_surface]
    y_surface = y[id_surface]
    z_surface = z[id_surface]
    spread_neighbors_substitute!(grid_LSC, x, y, z, x_grid, y_grid, z_grid, rcut, value_core)
    spread_neighbors_substitute!(grid_LSC, x_surface, y_surface, z_surface, x_grid, y_grid, z_grid, rcut_surface, value_surface)

    # ligand grid: desolvation free energy
    grid_LDS = zeros(Complex{T}, (nx, ny, nz))
    rcut_ds = fill(T(6.0), ligand.natom)
    spread_nearest!(grid_LDS, x, y, z, x_grid, y_grid, z_grid)
    spread_neighbors_add!(grid_LDS, x, y, z, x_grid, y_grid, z_grid, rcut_ds, ligand.mass)

    println("step7: docking")
    score_sc = similar(grid_RSC, T)
    score_ds = similar(grid_RSC, T)
    score_tops = fill(typemax(eltype(Float32)), 2 * tops)
    cartesian_tops = Vector{CartesianIndex{3}}(undef, 2 * tops)
    iq_tops = fill(-1, 2 * tops)
    if CUDA.functional()
        grid_RSC_d = cu(grid_RSC)
        grid_LSC_d = cu(grid_LSC)
        grid_RDS_d = cu(grid_RDS)
        grid_LDS_d = cu(grid_LDS)
        grid_LDS_real = real.(grid_LDS_d)
        grid_LDS_imag = real.(grid_LDS_d)
        x_grid_d = cu(x_grid)
        y_grid_d = cu(y_grid)
        z_grid_d = cu(z_grid)
        x_org = cu(x)
        y_org = cu(y)
        z_org = cu(z)
        x_d = cu(x)
        y_d = cu(y)
        z_d = cu(z)

        id_surface_d = cu(id_surface)
        rcut_d = cu(rcut)
        rcut_surface_d = cu(rcut_surface)
        value_core_d = cu(value_core)
        value_surface_d = cu(value_surface)

        rcut_ds_d = cu(rcut_ds)
        ligand_mass_d = cu(ligand.mass)

        t_d = similar(grid_RSC_d)
        score_sc_d = real(t_d)
        score_ds_d = real(t_d)
        score_d = real(t_d)

        quaternions_d = cu(quaternions)
        nblocks = ceil(Int, length(x) / 256)
        @time @showprogress for iq = 1:size(quaternions, 1)
            # rotate ligand
            x_d .= x_org
            y_d .= y_org
            z_d .= z_org
            #@time CUDA.@sync @cuda threads=256 blocks=nblocks rotate_with_matrix_gpu!(x_d, y_d, z_d, R_d)
            @cuda threads = 256 blocks = nblocks rotate_gpu!(x_d, y_d, z_d, quaternions_d[iq, :])

            # shape complementarity
            grid_LSC_d .= zero(eltype(grid_LSC_d))
            @cuda threads = 256 blocks = nblocks spread_neighbors_substitute_gpu!(grid_LSC_d, x_d, y_d, z_d, x_grid_d, y_grid_d, z_grid_d, rcut_d, value_core_d)
            @cuda threads = 256 blocks = nblocks spread_neighbors_substitute_gpu!(grid_LSC_d, x_d[id_surface_d], y_d[id_surface_d], z_d[id_surface_d], x_grid_d, y_grid_d, z_grid_d, rcut_surface_d, value_surface_d)
            #t_d .= ifftshift(ifft(ifft(grid_RSC_d) .* fft(grid_LSC_d))) .* nxyz
            #score_sc_d .= - real(t_d) .+ imag(t_d)
            t_d .= ifftshift(ifft(ifft(grid_RSC_d) .* fft(grid_LSC_d))) .* nxyz
            score_sc_d .= -real(t_d) .+ imag(t_d)

            # desolvation free energy
            grid_LDS_d .= zero(eltype(grid_LDS_d))
            grid_LDS_real .= zero(eltype(grid_LDS_real))
            grid_LDS_imag .= zero(eltype(grid_LDS_imag))
            @cuda threads = 256 blocks = nblocks spread_nearest_gpu!(grid_LDS_imag, x_d, y_d, z_d, x_grid_d, y_grid_d, z_grid_d)
            #@show typeof(grid_LDS_real)
            #@show typeof(ligand_mass_d)
            @cuda threads = 256 blocks = nblocks spread_neighbors_add_gpu!(grid_LDS_real, x_d, y_d, z_d, x_grid_d, y_grid_d, z_grid_d, rcut_ds_d, ligand_mass_d)
            grid_LDS_d .= grid_LDS_real .+ grid_LDS_imag .* im
            t_d .= 0.5 .* ifftshift(ifft(ifft(grid_RDS_d) .* fft(grid_LDS_d))) .* nxyz
            score_ds_d .= -imag(t_d)

            # filter top scores
            #score_d .= alpha .* score_sc_d
            score_d .= alpha .* score_sc_d .+ score_ds_d
            filter_tops!(score_tops, cartesian_tops, iq_tops, score_d, iq, tops)
        end
        grid_RSC .= Array(grid_RSC_d)
        grid_LSC .= Array(grid_LSC_d)
        grid_RDS .= Array(grid_RDS_d)
        grid_LDS .= Array(grid_LDS_d)
        score_sc .= Array(score_sc_d)
        score_ds .= Array(score_ds_d)
    else
        x_org = deepcopy(x)
        y_org = deepcopy(y)
        z_org = deepcopy(z)
        t = similar(grid_LSC)
        score = real(t)

        @time @showprogress for iq = 1:size(quaternions, 1)
            # rotate ligand
            x .= x_org
            y .= y_org
            z .= z_org
            rotate!(x, y, z, quaternions[iq, :])

            # shape complementarity
            grid_LSC .= zero(eltype(grid_LSC))
            spread_neighbors_substitute!(grid_LSC, x, y, z, x_grid, y_grid, z_grid, rcut, value_core)
            spread_neighbors_substitute!(grid_LSC, x[id_surface], y[id_surface], z[id_surface], x_grid, y_grid, z_grid, rcut_surface, value_surface)
            #t .= ifft(fft(grid_RSC) .* conj.(fft(conj.(grid_LSC))))
            #score_sc .= real(t)
            t .= ifftshift(ifft(ifft(grid_RSC) .* fft(grid_LSC))) .* nxyz
            score_sc .= -real(t) .+ imag(t)

            # desolvation free energy
            grid_LDS .= zero(eltype(grid_LDS))
            spread_nearest!(grid_LDS, x, y, z, x_grid, y_grid, z_grid)
            spread_neighbors_add!(grid_LDS, x, y, z, x_grid, y_grid, z_grid, rcut_ds, ligand.mass)
            t .= 0.5 .* ifftshift(ifft(ifft(grid_RDS) .* fft(grid_LDS))) .* nxyz
            score_ds .= -imag(t)

            # filter top scores
            #score .= score_sc
            score .= alpha .* score_sc .+ score_ds
            filter_tops!(score_tops, cartesian_tops, iq_tops, score, iq, tops)
        end
    end

    ligand_init = deepcopy(ligand[iframe, :])
    ligand_return = deepcopy(ligand[0, :])
    ix_center = ceil(Int, (nx / 2.0) + 1.0)
    iy_center = ceil(Int, (ny / 2.0) + 1.0)
    iz_center = ceil(Int, (nz / 2.0) + 1.0)
    for itop = 1:tops
        iq = iq_tops[itop]
        ligand_tmp = rotate(ligand_init, quaternions[iq, :])
        dx = (cartesian_tops[itop][1] - ix_center) * grid_space
        dy = (cartesian_tops[itop][2] - iy_center) * grid_space
        dz = (cartesian_tops[itop][3] - iz_center) * grid_space
        ligand_tmp.xyz[1, 1:3:end] .-= dx
        ligand_tmp.xyz[1, 2:3:end] .-= dy
        ligand_tmp.xyz[1, 3:3:end] .-= dz
        ligand_return = [ligand_return; ligand_tmp]
    end

    return (receptor=receptor, ligand=ligand_return, score=score_tops, iq=iq_tops, cartesian=cartesian_tops,
        grid_RSC=grid_RSC, grid_LSC=grid_LSC,
        grid_RDS=grid_RDS, grid_LDS=grid_LDS,
        score_sc, score_ds)
end

