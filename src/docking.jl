function set_acescore(ta)
    atom_type = Array{String}(undef, ta.natom)
    ace_score = Array{Float64}(undef, ta.natom)

    for iatom = 1:ta.natom
        # ATOM TYPE "N"
        if ta.atomname[iatom] == "N"
            atom_type[iatom] = "N"
            ace_score[iatom] = -0.495

            # ATOM TYPE "C"
        elseif ta.atomname[iatom] == "C"
            atom_type[iatom] = "C"
            ace_score[iatom] = -0.464

            # ATOM TYPE "O"
        elseif ta.atomname[iatom] == "O" || ta.atomname[iatom] == "OXT"
            atom_type[iatom] = "O"
            ace_score[iatom] = -0.079

            # ATOM TYPE "GCA"
        elseif ta.resname[iatom] == "GLY" && ta.atomname[iatom] == "CA"
            atom_type[iatom] = "GCA"
            ace_score[iatom] = -0.008

            # ATOM TYPE "CA"
        elseif ta.atomname[iatom] == "CA"
            atom_type[iatom] = "CA"
            ace_score[iatom] = -0.553

            # ATOM TYPE "CB"
        elseif ta.resname[iatom] ==  "ALA" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "ARG" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "ASN" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "ASP" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "CYS" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "GLN" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "GLU" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "HIS" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "ILE" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "LEU" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "LYS" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "MET" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "PHE" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "PRO" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "PRO" && ta.atomname[iatom] == "CG"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "PRO" && ta.atomname[iatom] == "CD"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "THR" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "TRP" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "TYR" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353
        elseif ta.resname[iatom] ==  "VAL" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "CB"
            ace_score[iatom] = -0.353

            # ATOM TYPE "KNZ"
        elseif ta.resname[iatom] ==  "LYS" && ta.atomname[iatom] == "CE"
            atom_type[iatom] = "KNZ"
            ace_score[iatom] = 1.334
        elseif ta.resname[iatom] ==  "LYS" && ta.atomname[iatom] == "NZ"
            atom_type[iatom] = "KNZ"
            ace_score[iatom] = 1.334

            # ATOM TYPE "KCD"
        elseif ta.resname[iatom] ==  "LYS" && ta.atomname[iatom] == "CD"
            atom_type[iatom] = "KCD"
            ace_score[iatom] = 1.046

            # ATOM TYPE "DOD"
        elseif ta.resname[iatom] ==  "ASP" && ta.atomname[iatom] == "CG"
            atom_type[iatom] = "DOD"
            ace_score[iatom] = 0.933
        elseif ta.resname[iatom] ==  "ASP" && ta.atomname[iatom] == "OD1"
            atom_type[iatom] = "DOD"
            ace_score[iatom] = 0.933
        elseif ta.resname[iatom] ==  "ASP" && ta.atomname[iatom] == "OD2"
            atom_type[iatom] = "DOD"
            ace_score[iatom] = 0.933
        elseif ta.resname[iatom] ==  "GLU" && ta.atomname[iatom] == "CD"
            atom_type[iatom] = "DOD"
            ace_score[iatom] = 0.933
        elseif ta.resname[iatom] ==  "GLU" && ta.atomname[iatom] == "OE1"
            atom_type[iatom] = "DOD"
            ace_score[iatom] = 0.933
        elseif ta.resname[iatom] ==  "GLU" && ta.atomname[iatom] == "OE2"
            atom_type[iatom] = "DOD"
            ace_score[iatom] = 0.933

            # ATOM TYPE "RNH"
        elseif ta.resname[iatom] ==  "ARG" && ta.atomname[iatom] == "CZ"
            atom_type[iatom] = "RNH"
            ace_score[iatom] = 0.726
        elseif ta.resname[iatom] ==  "ARG" && ta.atomname[iatom] == "NH1"
            atom_type[iatom] = "RNH"
            ace_score[iatom] = 0.726
        elseif ta.resname[iatom] ==  "ARG" && ta.atomname[iatom] == "NH2"
            atom_type[iatom] = "RNH"
            ace_score[iatom] = 0.726

            # ATOM TYPE "NND"
        elseif ta.resname[iatom] ==  "ASN" && ta.atomname[iatom] == "CG"
            atom_type[iatom] = "NND"
            ace_score[iatom] = 0.693
        elseif ta.resname[iatom] ==  "ASN" && ta.atomname[iatom] == "OD1"
            atom_type[iatom] = "NND"
            ace_score[iatom] = 0.693
        elseif ta.resname[iatom] ==  "ASN" && ta.atomname[iatom] == "ND2"
            atom_type[iatom] = "NND"
            ace_score[iatom] = 0.693
        elseif ta.resname[iatom] ==  "GLN" && ta.atomname[iatom] == "CD"
            atom_type[iatom] = "NND"
            ace_score[iatom] = 0.693
        elseif ta.resname[iatom] ==  "GLN" && ta.atomname[iatom] == "OE1"
            atom_type[iatom] = "NND"
            ace_score[iatom] = 0.693
        elseif ta.resname[iatom] ==  "GLN" && ta.atomname[iatom] == "NE2"
            atom_type[iatom] = "NND"
            ace_score[iatom] = 0.693

            # ATOM TYPE "RNE"
        elseif ta.resname[iatom] ==  "ARG" && ta.atomname[iatom] == "CD"
            atom_type[iatom] = "RNE"
            ace_score[iatom] = 0.606
        elseif ta.resname[iatom] ==  "ARG" && ta.atomname[iatom] == "NE"
            atom_type[iatom] = "RNE"
            ace_score[iatom] = 0.606

            # ATOM TYPE "SOG"
        elseif ta.resname[iatom] ==  "SER" && ta.atomname[iatom] == "CB"
            atom_type[iatom] = "SOG"
            ace_score[iatom] = 0.232
        elseif ta.resname[iatom] ==  "SER" && ta.atomname[iatom] == "OG"
            atom_type[iatom] = "SOG"
            ace_score[iatom] = 0.232
        elseif ta.resname[iatom] ==  "THR" && ta.atomname[iatom] == "OG1"
            atom_type[iatom] = "SOG"
            ace_score[iatom] = 0.232
        elseif ta.resname[iatom] ==  "TYR" && ta.atomname[iatom] == "OH"
            atom_type[iatom] = "SOG"
            ace_score[iatom] = 0.232

            # ATOM TYPE "HNE"
        elseif ta.resname[iatom] ==  "HIS" && ta.atomname[iatom] == "CG"
            atom_type[iatom] = "HNE"
            ace_score[iatom] = 0.061
        elseif ta.resname[iatom] ==  "HIS" && ta.atomname[iatom] == "ND1"
            atom_type[iatom] = "HNE"
            ace_score[iatom] = 0.061
        elseif ta.resname[iatom] ==  "HIS" && ta.atomname[iatom] == "CD2"
            atom_type[iatom] = "HNE"
            ace_score[iatom] = 0.061
        elseif ta.resname[iatom] ==  "HIS" && ta.atomname[iatom] == "CE1"
            atom_type[iatom] = "HNE"
            ace_score[iatom] = 0.061
        elseif ta.resname[iatom] ==  "HIS" && ta.atomname[iatom] == "NE2"
            atom_type[iatom] = "HNE"
            ace_score[iatom] = 0.061
        elseif ta.resname[iatom] ==  "TRP" && ta.atomname[iatom] == "NE1"
            atom_type[iatom] = "HNE"
            ace_score[iatom] = 0.061

            # ATOM TYPE "YCZ"
        elseif ta.resname[iatom] ==  "TYR" && ta.atomname[iatom] == "CE1"
            atom_type[iatom] = "YCZ"
            ace_score[iatom] = -0.289
        elseif ta.resname[iatom] ==  "TYR" && ta.atomname[iatom] == "CE2"
            atom_type[iatom] = "YCZ"
            ace_score[iatom] = -0.289
        elseif ta.resname[iatom] ==  "TYR" && ta.atomname[iatom] == "CZ"
            atom_type[iatom] = "YCZ"
            ace_score[iatom] = -0.289

            # ATOM TYPE "FCZ"
        elseif ta.resname[iatom] ==  "ARG" && ta.atomname[iatom] == "CG"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "GLN" && ta.atomname[iatom] == "CG"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "GLU" && ta.atomname[iatom] == "CG"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "ILE" && ta.atomname[iatom] == "CG1"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "LEU" && ta.atomname[iatom] == "CG"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "LYS" && ta.atomname[iatom] == "CG"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "MET" && ta.atomname[iatom] == "CG"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "MET" && ta.atomname[iatom] == "SD"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "PHE" && ta.atomname[iatom] == "CG"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "PHE" && ta.atomname[iatom] == "CD1"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "PHE" && ta.atomname[iatom] == "CD2"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "PHE" && ta.atomname[iatom] == "CE1"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "PHE" && ta.atomname[iatom] == "CE2"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "PHE" && ta.atomname[iatom] == "CZ"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "THR" && ta.atomname[iatom] == "CG2"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "TRP" && ta.atomname[iatom] == "CG"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "TRP" && ta.atomname[iatom] == "CD1"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "TRP" && ta.atomname[iatom] == "CD2"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "TRP" && ta.atomname[iatom] == "CE2"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "TRP" && ta.atomname[iatom] == "CE3"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "TRP" && ta.atomname[iatom] == "CZ2"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "TRP" && ta.atomname[iatom] == "CZ3"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "TRP" && ta.atomname[iatom] == "CH2"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "TYR" && ta.atomname[iatom] == "CG"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "TYR" && ta.atomname[iatom] == "CD1"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432
        elseif ta.resname[iatom] ==  "TYR" && ta.atomname[iatom] == "CD2"
            atom_type[iatom] = "FCZ"
            ace_score[iatom] = -0.432

            # ATOM TYPE "LCD"
        elseif ta.resname[iatom] ==  "ILE" && ta.atomname[iatom] == "CG2"
            atom_type[iatom] = "LCD"
            ace_score[iatom] = -0.987
        elseif ta.resname[iatom] ==  "ILE" && ta.atomname[iatom] == "CD"
            atom_type[iatom] = "LCD"
            ace_score[iatom] = -0.987
        elseif ta.resname[iatom] ==  "ILE" && ta.atomname[iatom] == "CD1"
            atom_type[iatom] = "LCD"
            ace_score[iatom] = -0.987
        elseif ta.resname[iatom] ==  "LEU" && ta.atomname[iatom] == "CD1"
            atom_type[iatom] = "LCD"
            ace_score[iatom] = -0.987
        elseif ta.resname[iatom] ==  "LEU" && ta.atomname[iatom] == "CD2"
            atom_type[iatom] = "LCD"
            ace_score[iatom] = -0.987
        elseif ta.resname[iatom] ==  "MET" && ta.atomname[iatom] == "CE"
            atom_type[iatom] = "LCD"
            ace_score[iatom] = -0.987
        elseif ta.resname[iatom] ==  "VAL" && ta.atomname[iatom] == "CG1"
            atom_type[iatom] = "LCD"
            ace_score[iatom] = -0.987
        elseif ta.resname[iatom] ==  "VAL" && ta.atomname[iatom] == "CG2"
            atom_type[iatom] = "LCD"
            ace_score[iatom] = -0.987

            # ATOM TYPE "CSG"
        elseif ta.resname[iatom] ==  "CYS" && ta.atomname[iatom] == "SG"
            atom_type[iatom] = "CSG"
            ace_score[iatom] = -1.827
        else
            println("error: faled to assign atom type " * ta.atomname[iatom] * "-" * ta.resname[iatom])
        end
    end

    return TrjArray(ta, mass=ace_score)
end

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
        elseif  ta.atomname[iatom] == "OXT"
            charge[iatom] = -1.0
            
        elseif  ta.resname[iatom] ==  "ARG" && ta.atomname[iatom] == "NH1"
            charge[iatom] = 0.5
        elseif ta.resname[iatom] ==  "ARG" && ta.atomname[iatom] == "NH2"
            charge[iatom] = 0.5
        elseif ta.resname[iatom] ==  "GLU" && ta.atomname[iatom] == "OE1"
            charge[iatom] = -0.5
        elseif ta.resname[iatom] ==  "GLU" && ta.atomname[iatom] == "OE2"
            charge[iatom] = -0.5
        elseif ta.resname[iatom] ==  "ASP" && ta.atomname[iatom] == "OD1"
            charge[iatom] = -0.5
        elseif ta.resname[iatom] ==  "ASP" && ta.atomname[iatom] == "OD2"
            charge[iatom] = -0.5
        elseif ta.resname[iatom] ==  "LYS" && ta.atomname[iatom] == "NZ"
            charge[iatom] = 1.0
        elseif ta.resname[iatom] ==  "PRO" && ta.atomname[iatom] == "N"
            charge[iatom] = -0.1
            
        else
            charge[iatom] = 0.0
        end
    end
    
    return TrjArray(ta, charge=charge)
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
                d = 0.0
                d += (point[1] - ta.xyz[iframe, 3*(jatom-1)+1])^2
                d += (point[2] - ta.xyz[iframe, 3*(jatom-1)+2])^2
                d += (point[3] - ta.xyz[iframe, 3*(jatom-1)+3])^2
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
    
    return TrjArray(ta, sasa=sasa)
end

function assign_shape_complementarity!(grid, ta::TrjArray{T, U}, grid_space, 
    rcut1, rcut2, x_grid, y_grid, z_grid, iframe) where {T,U}
    grid .= 0.0 + 0.0im
    nx, ny, nz = size(grid)
    
    for iatom = 1:ta.natom
        rcut = rcut1[iatom]
        if rcut < 0.0
            continue
        end
        
        x = ta.xyz[iframe, 3*(iatom-1)+1]
        y = ta.xyz[iframe, 3*(iatom-1)+2]
        z = ta.xyz[iframe, 3*(iatom-1)+3]
        
        dx = x - x_grid[1]
        ix_min = floor(U, (dx - rcut)/grid_space) + 1
        ix_min = max(ix_min, 1)
        ix_max = floor(U, (dx + rcut)/grid_space) + 2
        ix_max = min(ix_max, nx)
        
        dy = y - y_grid[1]
        iy_min = floor(U, (dy - rcut)/grid_space) + 1
        iy_min = max(iy_min, 1)
        iy_max = floor(U, (dy + rcut)/grid_space) + 2
        iy_max = min(iy_max, ny)
        
        dz = z - z_grid[1]
        iz_min = floor(U, (dz - rcut)/grid_space) + 1
        iz_min = max(iz_min, 1)
        iz_max = floor(U, (dz + rcut)/grid_space) + 2
        iz_max = min(iz_max, nz)
        
        for ix = ix_min:ix_max
            for iy = iy_min:iy_max
                for iz = iz_min:iz_max
                    if imag(grid[ix, iy, iz]) < 0.0001
                        dist = sqrt((x - x_grid[ix])^2 + (y - y_grid[iy])^2 + (z - z_grid[iz])^2)
                        if dist < rcut
                            grid[ix, iy, iz] = 0.0 + 9.0im
                        end
                    end
                end
            end
        end
    end
    
    for iatom = 1:ta.natom
        rcut = rcut2[iatom]
        if rcut < 0.0
            continue
        end
        
        x = ta.xyz[iframe, 3*(iatom-1)+1]
        y = ta.xyz[iframe, 3*(iatom-1)+2]
        z = ta.xyz[iframe, 3*(iatom-1)+3]
        
        dx = x - x_grid[1]
        ix_min = floor(U, (dx - rcut)/grid_space) + 1
        ix_min = max(ix_min, 1)
        ix_max = floor(U, (dx + rcut)/grid_space) + 2
        ix_max = min(ix_max, nx)
        
        dy = y - y_grid[1]
        iy_min = floor(U, (dy - rcut)/grid_space) + 1
        iy_min = max(iy_min, 1)
        iy_max = floor(U, (dy + rcut)/grid_space) + 2
        iy_max = min(iy_max, ny)
        
        dz = z - z_grid[1]
        iz_min = floor(U, (dz - rcut)/grid_space) + 1
        iz_min = max(iz_min, 1)
        iz_max = floor(U, (dz + rcut)/grid_space) + 2
        iz_max = min(iz_max, nz)
        
        for ix = ix_min:ix_max
            for iy = iy_min:iy_max
                for iz = iz_min:iz_max
                    if imag(grid[ix, iy, iz]) < 0.0001
                        dist = sqrt((x - x_grid[ix])^2 + (y - y_grid[iy])^2 + (z - z_grid[iz])^2)
                        if dist < rcut
                            grid[ix, iy, iz] = 1.0 + 0.0im
                        end
                    end
                end
            end
        end
    end
end

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

function compute_docking_score_on_xyplane(grid_RSC, grid_LSC)
    if CUDA.functional()
        grid_RSC_gpu = cu(grid_RSC)
        grid_LSC_gpu = cu(grid_LSC)
        t_gpu = ifft(fft(grid_RSC_gpu) .* conj.(fft(conj.(grid_LSC_gpu))))
        score_gpu = real(t_gpu)
        score = Array(score_gpu)
    else
        t = ifft(fft(grid_RSC) .* conj.(fft(conj.(grid_LSC))))
        score = real(t)
    end
    
    return score[:, :, 1]
end

function spread_nearest!(grid, x, y, z, grid_x, grid_y, grid_z)
    T = eltype(grid)
    natom = length(x)
    grid .= zero(T)

    for iatom = 1:natom
        dx = x[iatom] .- grid_x
        dy = y[iatom] .- grid_y
        dz = z[iatom] .- grid_z
        ix = argmin(abs.(dx))
        iy = argmin(abs.(dy))
        iz = argmin(abs.(dz))
        grid[ix, iy, iz] = T(1.0*im)
    end
            
    return
end

function spread_nearest_gpu!(grid, x, y, z, grid_x, grid_y, grid_z)
    T = eltype(grid)
    natom = length(x)
    ngid = gridDim().x*blockDim().x
    gid = (blockIdx().x - 1)*blockDim().x + threadIdx().x - 1

    grid_x_min = grid_x[1]
    grid_y_min = grid_y[1]
    grid_z_min = grid_z[1]

    grid_x_delta = grid_x[2] - grid_x[1]
    grid_y_delta = grid_y[2] - grid_y[1]
    grid_z_delta = grid_z[2] - grid_z[1]

    grid_x_delta_half = grid_x_delta*0.5
    grid_y_delta_half = grid_y_delta*0.5
    grid_z_delta_half = grid_z_delta*0.5


    ntile = ceil(Int, natom/ngid)

    for itile = 1:ntile
        iatom = ngid*(itile-1) + gid + 1
        if iatom > natom
            return nothing
        end
        atom_x = x[iatom]
        atom_y = y[iatom]
        atom_z = z[iatom]
        
        ix = ceil(Int, (atom_x - grid_x_min)/grid_x_delta + grid_x_delta_half) + 1
        iy = ceil(Int, (atom_y - grid_y_min)/grid_y_delta + grid_y_delta_half) + 1
        iz = ceil(Int, (atom_z - grid_z_min)/grid_z_delta + grid_z_delta_half) + 1
        grid[ix,iy,iz] = T(1.0)
    end

    return nothing
end

function spread_neighbors_add!(grid, x, y, z, grid_x, grid_y, grid_z, rcut, charge)
    T = eltype(grid)
    natom = length(x)
    nx, ny, nz = size(grid)

    for iatom = 1:natom
        for ix = 1:nx
            dx = x[iatom] - grid_x[ix]
            if abs(dx) > rcut[iatom]
                continue
            end
            for iy = 1:ny
                dy = y[iatom] - grid_y[iy]
                if abs(dy) > rcut[iatom]
                    continue
                end
                for iz = 1:nz
                    dz = z[iatom] - grid_z[iz]
                    if abs(dz) > rcut[iatom]
                        continue
                    end
                    d = dx*dx + dy*dy + dz*dz
                    if d < rcut[iatom]*rcut[iatom]
                        grid[ix, iy, iz] += charge[iatom]
                    end
                end
            end
        end
    end
            
    return
end

function spread_neighbors_add_gpu!(grid, x, y, z, grid_x, grid_y, grid_z, rcut, charge)
    T = eltype(grid)
    natom = length(x)
    ngid = gridDim().x*blockDim().x
    gid = (blockIdx().x - 1)*blockDim().x + threadIdx().x - 1

    grid_x_min = grid_x[1]
    grid_y_min = grid_y[1]
    grid_z_min = grid_z[1]

    grid_x_delta = grid_x[2] - grid_x[1]
    grid_y_delta = grid_y[2] - grid_y[1]
    grid_z_delta = grid_z[2] - grid_z[1]

    nx, ny, nz = size(grid)

    ntile = ceil(Int, natom/ngid)

    for itile = 1:ntile
        iatom = ngid*(itile-1) + gid + 1
        if iatom > natom
            return
        end
        atom_x = x[iatom]
        atom_y = y[iatom]
        atom_z = z[iatom]
        r = rcut[iatom]
        r2 = r*r
        c = charge[iatom]

        ix_min = max(floor(Int, (atom_x - r - grid_x_min)/grid_x_delta) - 1, 1)
        iy_min = max(floor(Int, (atom_y - r - grid_y_min)/grid_y_delta) - 1, 1)
        iz_min = max(floor(Int, (atom_z - r - grid_z_min)/grid_z_delta) - 1, 1)

        ix_max = min(floor(Int, (atom_x + r - grid_x_min)/grid_x_delta) + 2, nx)
        iy_max = min(floor(Int, (atom_y + r - grid_y_min)/grid_y_delta) + 2, ny)
        iz_max = min(floor(Int, (atom_z + r - grid_z_min)/grid_z_delta) + 2, nz)        

        for ix = ix_min:ix_max
            for iy = iy_min:iy_max
                for iz = iz_min:iz_max
                    xx = grid_x_min + grid_x_delta*(ix-1)
                    yy = grid_y_min + grid_y_delta*(iy-1)
                    zz = grid_z_min + grid_z_delta*(iz-1)
                    d = (atom_x - xx)^2 + (atom_y - yy)^2 + (atom_z - zz)^2
                    if d < r2
                        #@cuprintln("called1")
                        #@cuprintln("c = $(c)")
                        #@cuprintln("called2")
                        #@cuprintln("grid[ix,iy,iz] = $(grid[ix,iy,iz])")
                        @atomic grid[ix,iy,iz] += c
                    end
                end
            end
        end
    end

    return
end

function spread_neighbors_substitute!(grid, x, y, z, grid_x, grid_y, grid_z, rcut, charge)
    T = eltype(grid)
    natom = length(x)
    nx, ny, nz = size(grid)

    for iatom = 1:natom
        for ix = 1:nx
            dx = x[iatom] - grid_x[ix]
            if abs(dx) > rcut[iatom]
                continue
            end
            for iy = 1:ny
                dy = y[iatom] - grid_y[iy]
                if abs(dy) > rcut[iatom]
                    continue
                end
                for iz = 1:nz
                    dz = z[iatom] - grid_z[iz]
                    if abs(dz) > rcut[iatom]
                        continue
                    end
                    d = dx*dx + dy*dy + dz*dz
                    if d < rcut[iatom]*rcut[iatom]
                        grid[ix, iy, iz] = charge[iatom]
                    end
                end
            end
        end
    end
            
    return
end

function spread_neighbors_substitute_gpu!(grid, x, y, z, grid_x, grid_y, grid_z, rcut, charge)
    T = eltype(grid)
    natom = length(x)
    ngid = gridDim().x*blockDim().x
    gid = (blockIdx().x - 1)*blockDim().x + threadIdx().x - 1

    grid_x_min = grid_x[1]
    grid_y_min = grid_y[1]
    grid_z_min = grid_z[1]

    grid_x_delta = grid_x[2] - grid_x[1]
    grid_y_delta = grid_y[2] - grid_y[1]
    grid_z_delta = grid_z[2] - grid_z[1]

    nx, ny, nz = size(grid)

    ntile = ceil(Int, natom/ngid)

    for itile = 1:ntile
        iatom = ngid*(itile-1) + gid + 1
        if iatom > natom
            return
        end
        atom_x = x[iatom]
        atom_y = y[iatom]
        atom_z = z[iatom]
        r = rcut[iatom]
        r2 = r*r
        c = charge[iatom]

        ix_min = max(floor(Int, (atom_x - r - grid_x_min)/grid_x_delta) - 1, 1)
        iy_min = max(floor(Int, (atom_y - r - grid_y_min)/grid_y_delta) - 1, 1)
        iz_min = max(floor(Int, (atom_z - r - grid_z_min)/grid_z_delta) - 1, 1)

        ix_max = min(floor(Int, (atom_x + r - grid_x_min)/grid_x_delta) + 2, nx)
        iy_max = min(floor(Int, (atom_y + r - grid_y_min)/grid_y_delta) + 2, ny)
        iz_max = min(floor(Int, (atom_z + r - grid_z_min)/grid_z_delta) + 2, nz)        

        for ix = ix_min:ix_max
            for iy = iy_min:iy_max
                for iz = iz_min:iz_max
                    xx = grid_x_min + grid_x_delta*(ix-1)
                    yy = grid_y_min + grid_y_delta*(iy-1)
                    zz = grid_z_min + grid_z_delta*(iz-1)
                    d = (atom_x - xx)^2 + (atom_y - yy)^2 + (atom_z - zz)^2
                    if d < r2
                        #@cuprintln("called")
                        #@cuprintln("grid[ix,iy,iz] = $(grid[ix,iy,iz])")
                        #@cuprintln("c = $(c)")
                        grid[ix,iy,iz] = c
                    end
                end
            end
        end
    end

    return
end

function rotate_with_matrix!(x, y, z, R)
    natom = length(x)
    for iatom = 1:natom
        x_curr = x[iatom]
        y_curr = y[iatom]
        z_curr = z[iatom]
        x_new = R[1, 1]*x_curr + R[1, 2]*y_curr + R[1, 3]*z_curr
        y_new = R[2, 1]*x_curr + R[2, 2]*y_curr + R[2, 3]*z_curr
        z_new = R[3, 1]*x_curr + R[3, 2]*y_curr + R[3, 3]*z_curr
        x[iatom] = x_new
        y[iatom] = y_new
        z[iatom] = z_new
    end
    return
end

#function rotate_with_matrix_gpu!(x, y, z, R)
#    natom = length(x)
#    ngid = gridDim().x*blockDim().x
#    gid = (blockIdx().x - 1)*blockDim().x + threadIdx().x - 1
#
#    ntile = ceil(Int, natom/ngid)
#    for itile = 1:ntile
#        iatom = ngid*(itile-1) + gid + 1
#        if iatom > natom
#            return
#        end
#        x_curr = x[iatom]
#        y_curr = y[iatom]
#        z_curr = z[iatom]
#        x_new = R[1, 1]*x_curr + R[1, 2]*y_curr + R[1, 3]*z_curr
#        y_new = R[2, 1]*x_curr + R[2, 2]*y_curr + R[2, 3]*z_curr
#        z_new = R[3, 1]*x_curr + R[3, 2]*y_curr + R[3, 3]*z_curr
#        x[iatom] = x_new
#        y[iatom] = y_new
#        z[iatom] = z_new
#    end
#    return
#end

function rotate!(x, y, z, q::AbstractVector{T}) where {T, U}
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

function rotate_with_matrix!(x, y, z, R::AbstractMatrix{T}) where {T, U}
    natom = length(x)
    r1 = R[1, 1]
    r2 = R[1, 2]
    r3 = R[1, 3]
    r4 = R[2, 1]
    r5 = R[2, 2]
    r6 = R[2, 3]
    r7 = R[3, 1]
    r8 = R[3, 2]
    r9 = R[3, 3]
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

function rotate_gpu!(x, y, z, q)
    natom = length(x)
    ngid = gridDim().x*blockDim().x
    gid = (blockIdx().x - 1)*blockDim().x + threadIdx().x - 1

    r1 = 1.0 - 2.0 * q[2] * q[2] - 2.0 * q[3] * q[3]
    r2 = 2.0 * (q[1] * q[2] + q[3] * q[4])
    r3 = 2.0 * (q[1] * q[3] - q[2] * q[4])
    r4 = 2.0 * (q[1] * q[2] - q[3] * q[4])
    r5 = 1.0 - 2.0 * q[1] * q[1] - 2.0 * q[3] * q[3]
    r6 = 2.0 * (q[2] * q[3] + q[1] * q[4])
    r7 = 2.0 * (q[1] * q[3] + q[2] * q[4])
    r8 = 2.0 * (q[2] * q[3] - q[1] * q[4])
    r9 = 1.0 - 2.0 * q[1] * q[1] - 2.0 * q[2] * q[2]

    ntile = ceil(Int, natom/ngid)
    for itile = 1:ntile
        iatom = ngid*(itile-1) + gid + 1
        if iatom > natom
            return
        end
        x_new = r1 * x[iatom] + r2 * y[iatom] + r3 * z[iatom]
        y_new = r4 * x[iatom] + r5 * y[iatom] + r6 * z[iatom]
        z_new = r7 * x[iatom] + r8 * y[iatom] + r9 * z[iatom]
        x[iatom] = x_new
        y[iatom] = y_new
        z[iatom] = z_new
    end
    return nothing
end

function rotate_with_matrix_gpu!(x, y, z, q)
    natom = length(x)
    ngid = gridDim().x*blockDim().x
    gid = (blockIdx().x - 1)*blockDim().x + threadIdx().x - 1

    r1 = R[1, 1]
    r2 = R[1, 2]
    r3 = R[1, 3]
    r4 = R[2, 1]
    r5 = R[2, 2]
    r6 = R[2, 3]
    r7 = R[3, 1]
    r8 = R[3, 2]
    r9 = R[3, 3]

    ntile = ceil(Int, natom/ngid)
    for itile = 1:ntile
        iatom = ngid*(itile-1) + gid + 1
        if iatom > natom
            return
        end
        x_new = r1 * x[iatom] + r2 * y[iatom] + r3 * z[iatom]
        y_new = r4 * x[iatom] + r5 * y[iatom] + r6 * z[iatom]
        z_new = r7 * x[iatom] + r8 * y[iatom] + r9 * z[iatom]
        x[iatom] = x_new
        y[iatom] = y_new
        z[iatom] = z_new
    end
    return nothing
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
    iz_center = ceil(Int, (nz/2.0)+1.0)
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

function dock!(receptor::TrjArray{T, U}, ligand::TrjArray{T, U}, quaternions::Matrix{T}; deg=15.0, grid_space=1.2, iframe=1, tops=100, alpha=0.01, beta=0.06) where {T, U}
    decenter!(receptor)
    decenter!(ligand)
    orient!(receptor)
    orient!(ligand)

    # Assign atom radius
    println("step1: assigning atom radius")
    receptor = set_radius(receptor)
    ligand = set_radius(ligand)    

    # Solvent accessible surface
    println("step2: computing SASA")
    receptor = compute_sasa(receptor, 1.4)
    ligand = compute_sasa(ligand, 1.4)

    # ACE scores
    println("step3: assigning ACE scores")
    receptor = set_acescore(receptor)
    ligand = set_acescore(ligand)

    # Determine grid size and coordinates
    println("step5: computing grid size and coordinates")
    size_ligand = maximum(ligand.xyz[iframe, 1:3:end]) - minimum(ligand.xyz[iframe, 1:3:end])
    
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
    nxyz = nx*ny*nz

    println("step6: assigning values to grid points")
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
    score_tops = fill(typemax(eltype(Float32)), 2*tops)
    cartesian_tops = Vector{CartesianIndex{3}}(undef, 2*tops)
    iq_tops = fill(-1, 2*tops);
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
        nblocks = ceil(Int, length(x)/256)
        @time @showprogress for iq = 1:size(quaternions, 1)
            # rotate ligand
            x_d .= x_org
            y_d .= y_org
            z_d .= z_org
            #@time CUDA.@sync @cuda threads=256 blocks=nblocks rotate_with_matrix_gpu!(x_d, y_d, z_d, R_d)
            @cuda threads=256 blocks=nblocks rotate_gpu!(x_d, y_d, z_d, quaternions_d[iq, :])

            # shape complementarity
            grid_LSC_d .= zero(eltype(grid_LSC_d))
            @cuda threads=256 blocks=nblocks spread_neighbors_substitute_gpu!(grid_LSC_d, x_d, y_d, z_d, x_grid_d, y_grid_d, z_grid_d, rcut_d, value_core_d)
            @cuda threads=256 blocks=nblocks spread_neighbors_substitute_gpu!(grid_LSC_d, x_d[id_surface_d], y_d[id_surface_d], z_d[id_surface_d], x_grid_d, y_grid_d, z_grid_d, rcut_surface_d, value_surface_d)
            #t_d .= ifftshift(ifft(ifft(grid_RSC_d) .* fft(grid_LSC_d))) .* nxyz
            #score_sc_d .= - real(t_d) .+ imag(t_d)
            t_d .= ifftshift(ifft(ifft(grid_RSC_d) .* fft(grid_LSC_d))) .* nxyz
            score_sc_d .= - real(t_d) .+ imag(t_d)

            # desolvation free energy
            grid_LDS_d .= zero(eltype(grid_LDS_d))
            grid_LDS_real .= zero(eltype(grid_LDS_real))
            grid_LDS_imag .= zero(eltype(grid_LDS_imag))
            @cuda threads=256 blocks=nblocks spread_nearest_gpu!(grid_LDS_imag, x_d, y_d, z_d, x_grid_d, y_grid_d, z_grid_d)
            #@show typeof(grid_LDS_real)
            #@show typeof(ligand_mass_d)
            @cuda threads=256 blocks=nblocks spread_neighbors_add_gpu!(grid_LDS_real, x_d, y_d, z_d, x_grid_d, y_grid_d, z_grid_d, rcut_ds_d, ligand_mass_d)
            grid_LDS_d .= grid_LDS_real .+ grid_LDS_imag .* im
            t_d .= 0.5 .* ifftshift(ifft(ifft(grid_RDS_d) .* fft(grid_LDS_d))) .* nxyz
            score_ds_d .= - imag(t_d)
        
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
            score_sc .= - real(t) .+ imag(t)

            # desolvation free energy
            grid_LDS .= zero(eltype(grid_LDS))
            spread_nearest!(grid_LDS, x, y, z, x_grid, y_grid, z_grid)
            spread_neighbors_add!(grid_LDS, x, y, z, x_grid, y_grid, z_grid, rcut_ds, ligand.mass)
            t .= 0.5 .* ifftshift(ifft(ifft(grid_RDS) .* fft(grid_LDS))) .* nxyz
            score_ds .= - imag(t)

            # filter top scores
            #score .= score_sc
            score .= alpha .* score_sc .+ score_ds
            filter_tops!(score_tops, cartesian_tops, iq_tops, score, iq, tops)
        end
    end

    ligand_init = deepcopy(ligand[iframe, :])
    ligand_return = deepcopy(ligand[0, :])
    ix_center = ceil(Int, (nx/2.0)+1.0)
    iy_center = ceil(Int, (ny/2.0)+1.0)
    iz_center = ceil(Int, (nz/2.0)+1.0)
    for itop = 1:tops
        iq = iq_tops[itop]
        ligand_tmp = rotate(ligand_init, quaternions[iq, :])
        dx = (cartesian_tops[itop][1]-ix_center) * grid_space
        dy = (cartesian_tops[itop][2]-iy_center) * grid_space
        dz = (cartesian_tops[itop][3]-iz_center) * grid_space
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

function dock_multimer!(receptor::TrjArray{T, U}; deg=15.0, grid_space=1.2, iframe=1, tops=100, alpha=0.01, beta=0.06, nfold=3, rot_space=10.0) where {T, U}
    decenter!(receptor)
    orient!(receptor)

    # Assign atom radius
    println("step1: assigning atom radius")
    receptor = set_radius(receptor)

    # Solvent accessible surface
    println("step2: computing SASA")
    receptor = compute_sasa(receptor, 1.4)

    # ACE scores
    println("step3: assigning ACE scores")
    receptor = set_acescore(receptor)

    # Determine grid size and coordinates
    println("step5: computing grid size and coordinates")
    size_ligand = maximum(receptor.xyz[iframe, 1:3:end]) - minimum(receptor.xyz[iframe, 1:3:end])
    
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
    nxyz = nx*ny*nz

    println("step6: assigning values to grid points")
    # receptor grid: shape complementarity
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
    
    # receptor grid: desolvation free energy
    grid_RSC = zeros(Complex{T}, (nx, ny, nz))
    grid_LSC = zeros(Complex{T}, (nx, ny, nz))
    grid_RDS = zeros(Complex{T}, (nx, ny, nz))
    grid_LDS = zeros(Complex{T}, (nx, ny, nz))
    rcut_ds = fill(T(6.0), receptor.natom)
    #spread_nearest!(grid_RDS, x, y, z, x_grid, y_grid, z_grid)
    #spread_neighbors_add!(grid_RDS, x, y, z, x_grid, y_grid, z_grid, rcut_ds, receptor.mass)

    println("step7: docking")
    score_sc = similar(grid_RSC, T)
    score_ds = similar(grid_RSC, T)
    score_tops = fill(typemax(eltype(Float32)), 2*tops)
    cartesian_tops = Vector{CartesianIndex{2}}(undef, 2*tops)
    rotz_tops = fill(-1.0, 2*tops);
    rotx_tops = fill(-1.0, 2*tops);
    rot_z_pi = collect((0.0:rot_space:360.0)./ 360.0 .* (2.0*pi))
    rot_x_pi = collect((0.0:rot_space:90.0)./ 360.0 .* (2.0*pi))
    #if CUDA.functional()
    if 1 == 0
    else
        x_org = deepcopy(x)
        y_org = deepcopy(y)
        z_org = deepcopy(z)
        t = similar(grid_LSC)
        score = real(t)
        @showprogress for rotz in rot_z_pi
            R = [cos(rotz) -sin(rotz) 0.0; sin(rotz) cos(rotz) 0.0; 0.0 0.0 1.0]
            receptor_rot = rotate_with_matrix(receptor, R)
            for rotx in rot_x_pi
                R = [1.0 0.0 0.0; 0.0 cos(rotx) -sin(rotx); 0.0 sin(rotx) cos(rotx)]
                receptor_rot = rotate_with_matrix(receptor_rot, R)
                R = [cos(2.0*pi/nfold) -sin(2.0*pi/nfold) 0.0; sin(2.0*pi/nfold) cos(2.0*pi/nfold) 0.0; 0.0 0.0 1.0]
                ligand_rot = rotate_with_matrix(receptor_rot, R)

                x = receptor_rot.xyz[iframe, 1:3:end]
                y = receptor_rot.xyz[iframe, 2:3:end]
                z = receptor_rot.xyz[iframe, 3:3:end]            
                grid_RSC .= zero(eltype(grid_RSC))
                spread_neighbors_substitute!(grid_RSC, x[id_surface], y[id_surface], z[id_surface], x_grid, y_grid, z_grid, rcut_surface, value_surface)
                spread_neighbors_substitute!(grid_RSC, x, y, z, x_grid, y_grid, z_grid, rcut, value_core)
    
                x = ligand_rot.xyz[iframe, 1:3:end]
                y = ligand_rot.xyz[iframe, 2:3:end]
                z = ligand_rot.xyz[iframe, 3:3:end]            
                grid_LSC .= zero(eltype(grid_LSC))
                spread_neighbors_substitute!(grid_LSC, x, y, z, x_grid, y_grid, z_grid, rcut, value_core)
                spread_neighbors_substitute!(grid_LSC, x[id_surface], y[id_surface], z[id_surface], x_grid, y_grid, z_grid, rcut_surface, value_surface)

                t .= ifftshift(ifft(ifft(grid_RSC) .* fft(grid_LSC))) .* nxyz
                score_sc .= - real(t) .+ imag(t)

                x = receptor_rot.xyz[iframe, 1:3:end]
                y = receptor_rot.xyz[iframe, 2:3:end]
                z = receptor_rot.xyz[iframe, 3:3:end]            
                grid_RDS .= zero(eltype(grid_RDS))
                spread_nearest!(grid_RDS, x, y, z, x_grid, y_grid, z_grid)
                spread_neighbors_add!(grid_RDS, x, y, z, x_grid, y_grid, z_grid, rcut_ds, receptor.mass)
                    
                x = ligand_rot.xyz[iframe, 1:3:end]
                y = ligand_rot.xyz[iframe, 2:3:end]
                z = ligand_rot.xyz[iframe, 3:3:end]            
                grid_LDS .= zero(eltype(grid_LDS))
                spread_nearest!(grid_LDS, x, y, z, x_grid, y_grid, z_grid)
                spread_neighbors_add!(grid_LDS, x, y, z, x_grid, y_grid, z_grid, rcut_ds, receptor.mass)

                t .= 0.5 .* ifftshift(ifft(ifft(grid_RDS) .* fft(grid_LDS))) .* nxyz
                score_ds .= - imag(t)

                score .= alpha .* score_sc .+ score_ds
                filter_tops!(score_tops, cartesian_tops, rotz_tops, rotx_tops, score, rotz, rotx, tops)
            end
        end
    end

    receptor_init = deepcopy(receptor[iframe, :])
    receptor_return = deepcopy(receptor[0, :])
    ligand_return = deepcopy(receptor[0, :])
    ix_center = ceil(Int, (nx/2.0)+1.0)
    iy_center = ceil(Int, (ny/2.0)+1.0)
    iz_center = ceil(Int, (nz/2.0)+1.0)
    for itop = 1:tops
        rotz = rotz_tops[itop]
        R = [cos(rotz) -sin(rotz) 0.0; sin(rotz) cos(rotz) 0.0; 0.0 0.0 1.0]
        receptor_rot = rotate_with_matrix(receptor_init, R)
        rotx = rotx_tops[itop]
        R = [1.0 0.0 0.0; 0.0 cos(rotx) -sin(rotx); 0.0 sin(rotx) cos(rotx)]
        receptor_rot = rotate_with_matrix(receptor_rot, R)
        receptor_return = [receptor_return; receptor_rot]
        R = [cos(2.0*pi/nfold) -sin(2.0*pi/nfold) 0.0; sin(2.0*pi/nfold) cos(2.0*pi/nfold) 0.0; 0.0 0.0 1.0]
        ligand_rot = rotate_with_matrix(receptor_rot, R)
        dx = (cartesian_tops[itop][1]-ix_center) * grid_space
        dy = (cartesian_tops[itop][2]-iy_center) * grid_space
        #dz = (cartesian_tops[itop][3]-iz_center) * grid_space
        ligand_rot.xyz[1, 1:3:end] .-= dx
        ligand_rot.xyz[1, 2:3:end] .-= dy
        #ligand_rot.xyz[1, 3:3:end] .-= dz
        ligand_return = [ligand_return; ligand_rot]
    end

    natom = receptor_return.natom
    cc = TrjArray(receptor_return[1, :], chainname=fill(string(1), natom), chainid=fill(1, natom))
    for ifold = 2:nfold
        cc = [cc TrjArray(receptor_return[1, :], chainname=fill(string(ifold), natom), chainid=fill(ifold, natom))]
    end
    @show cc
    complex_return = deepcopy(cc[0, :])
    for itop = 1:tops
        com1 = centerofmass(receptor_return[itop, :])
        com2 = centerofmass(ligand_return[itop, :])
        L = [com2.xyz[1, 1] - com1.xyz[1, 1], com2.xyz[1, 2] - com1.xyz[1, 2]]
        L_norm = sqrt(L[1]^2 + L[2]^2)
        alpha = (180.0 - (360.0/nfold))/2.0
        alpha = (alpha / 360.0) * 2.0 * pi
        d_norm = L_norm / (2.0*cos(alpha))

        cc = TrjArray(receptor_return[itop, :], chainname=fill(string(1), natom), chainid=fill(1, natom))
        cc.xyz[:, 1:3:end] .= cc.xyz[:, 1:3:end] .+ d_norm
        tmp = deepcopy(cc)
        for ifold = 2:nfold
            beta = (2.0 * pi / nfold) * (ifold-1)
            R = [cos(beta) -sin(beta) 0.0; sin(beta) cos(beta) 0.0; 0.0 0.0 1.0]    
            r = rotate_with_matrix(tmp, R)
            cc = [cc TrjArray(r, chainname=fill(string(ifold), natom), chainid=fill(ifold, natom))]
        end
        complex_return = [complex_return; cc]
    end

    return (receptor=receptor, ligand=ligand_return, multimer=complex_return, score=score_tops, rotz=rotz_tops, rotx=rotx_tops, cartesian=cartesian_tops, 
            grid_RSC=grid_RSC, grid_LSC=grid_LSC, 
            grid_RDS=grid_RDS, grid_LDS=grid_LDS,
            score_sc, score_ds)
end

function dock_multimer(receptor::TrjArray{T, U}; rot_space=10.0, grid_space=1.2, iframe=1, tops=10, nfold=3) where {T, U}
    # generate grid coordinates for receptor
    receptor2, dummy = decenter(receptor)
    
    @show "step 1"
    x_min, x_max = minimum(receptor2.xyz[iframe, 1:3:end]), maximum(receptor2.xyz[iframe, 1:3:end])
    y_min, y_max = minimum(receptor2.xyz[iframe, 2:3:end]), maximum(receptor2.xyz[iframe, 2:3:end])
    z_min, z_max = minimum(receptor2.xyz[iframe, 3:3:end]), maximum(receptor2.xyz[iframe, 3:3:end])
    size_receptor = sqrt((x_max - x_min)^2 + (y_max - y_min)^2 + (z_max - z_min)^2)
    
    x_min = minimum(receptor2.xyz[iframe, 1:3:end]) - size_receptor - grid_space
    y_min = minimum(receptor2.xyz[iframe, 2:3:end]) - size_receptor - grid_space
    z_min = minimum(receptor2.xyz[iframe, 3:3:end]) - size_receptor - grid_space
    
    x_max = maximum(receptor2.xyz[iframe, 1:3:end]) + size_receptor + grid_space
    y_max = maximum(receptor2.xyz[iframe, 2:3:end]) + size_receptor + grid_space
    z_max = maximum(receptor2.xyz[iframe, 3:3:end]) + size_receptor + grid_space
    
    x_grid = collect(x_min:grid_space:x_max)
    y_grid = collect(y_min:grid_space:y_max)
    z_grid = collect(z_min:grid_space:z_max)
    
    @show "step 2"
    nx, ny, nz = length(x_grid), length(y_grid), length(z_grid)
    
    # shape complementarity of receptor
    @show "step 3"
    iatom_surface = receptor2.sasa .> 1.0
    iatom_core = .!iatom_surface
    rcut1 = zeros(T, receptor2.natom)
    rcut2 = zeros(T, receptor2.natom)
    
    rcut1[iatom_core] .= receptor2.radius[iatom_core] * sqrt(1.5)
    rcut2[iatom_core] .= -1.0
    rcut1[iatom_surface] .= receptor2.radius[iatom_surface] * sqrt(0.8)
    rcut2[iatom_surface] .= receptor2.radius[iatom_surface] .+ 3.4
    
    @show nx ny nz
    @show "step 4"
    grid_RSC = zeros(complex(T), nx, ny, nz)
    grid_LSC = zeros(complex(T), nx, ny, nz)
    #assign_shape_complementarity!(grid_RSC, receptor2, grid_space, rcut1, rcut2, x_grid, y_grid, z_grid, iframe)
    
    @show "step 4.1"
    # compute score with FFT
    #nq = size(quaternions, 1)
    #s = @showprogress pmap(q -> compute_docking_score_with_fft(quaternions[q, :], grid_RSC, grid_LSC, ligand2, grid_space, rcut1, rcut2, x_grid, y_grid, z_grid, iframe, tops, q), 1:nq)
    rot_z_pi = collect((0.0:rot_space:360.0)./ 360.0 .* (2.0*pi))
    rot_x_pi = collect((0.0:rot_space:360.0)./ 360.0 .* (2.0*pi))
    score_best = -Inf
    @show "step 4.2"
    receptor_best = deepcopy(receptor2)
    ligand_best = deepcopy(receptor2)
    dx_estimated = 0
    dy_estimated = 0
    @show "step 5"
    for z in rot_z_pi
        @printf "%f \n" z
        R = [cos(z) -sin(z) 0.0; sin(z) cos(z) 0.0; 0.0 0.0 1.0]
        receptor_rot = rotate_with_matrix(receptor2, R)
        for x in rot_x_pi
            R = [1.0 0.0 0.0; 0.0 cos(x) -sin(x); 0.0 sin(x) cos(x)]
            receptor_rot = rotate_with_matrix(receptor_rot, R)
            R = [cos(2.0*pi/nfold) -sin(2.0*pi/nfold) 0.0; sin(2.0*pi/nfold) cos(2.0*pi/nfold) 0.0; 0.0 0.0 1.0]
            ligand_rot = rotate_with_matrix(receptor_rot, R)
            assign_shape_complementarity!(grid_RSC, receptor_rot, grid_space, rcut1, rcut2, x_grid, y_grid, z_grid, iframe)
            assign_shape_complementarity!(grid_LSC, ligand_rot, grid_space, rcut1, rcut2, x_grid, y_grid, z_grid, iframe)
            score = compute_docking_score_on_xyplane(grid_RSC, grid_LSC)
            score_max = maximum(score)
            if score_best < score_max
                id = argmax(score)
                dx_estimated = id[1]
                dy_estimated = id[2]
                receptor_best = deepcopy(receptor_rot)
                ligand_best = deepcopy(ligand_rot)
                score_best = score_max
            end
        end    
    end
    
    ligand_return = deepcopy(ligand_best)
    dx = (dx_estimated-1) * grid_space
    if dx > (nx*grid_space / 2.0)
        dx = dx - (nx*grid_space)
    end
    dy = (dy_estimated-1) * grid_space
    if dy > (ny*grid_space / 2.0)
        dy = dy - (ny*grid_space)
    end
    ligand_return.xyz[iframe, 1:3:end] .+= dx
    ligand_return.xyz[iframe, 2:3:end] .+= dy
    
    return (receptor=receptor_best, ligand=ligand_return, score=score_best)
end

function docking_by_electrostatic_energy(receptor::TrjArray{T, U}, ligand::TrjArray{T, U}, quaternions; grid_space=1.2, iframe=1, tops=10) where {T, U}
    
    # generate grid coordinates of receptor
    decenter!(receptor)
    decenter!(ligand)
    
    # determine size of ligand
    x_min = minimum(ligand.xyz[iframe, 1:3:end])
    y_min = minimum(ligand.xyz[iframe, 2:3:end])
    z_min = minimum(ligand.xyz[iframe, 3:3:end])
    x_max = maximum(ligand.xyz[iframe, 1:3:end])
    y_max = maximum(ligand.xyz[iframe, 2:3:end])
    z_max = maximum(ligand.xyz[iframe, 3:3:end])
    size_ligand = sqrt((x_max-x_min)^2 + (y_max-y_min)^2 + (z_max-z_min)^2)
    
    # extension grid of receptor using size of ligand
    x_min = minimum(receptor.xyz[iframe, 1:3:end]) - size_ligand - grid_space
    y_min = minimum(receptor.xyz[iframe, 2:3:end]) - size_ligand - grid_space
    z_min = minimum(receptor.xyz[iframe, 3:3:end]) - size_ligand - grid_space
    x_max = maximum(receptor.xyz[iframe, 1:3:end]) + size_ligand + grid_space
    y_max = maximum(receptor.xyz[iframe, 2:3:end]) + size_ligand + grid_space
    z_max = maximum(receptor.xyz[iframe, 3:3:end]) + size_ligand + grid_space
    
    x_grid = collect(x_min:grid_space:x_max)
    y_grid = collect(y_min:grid_space:y_max)
    z_grid = collect(z_min:grid_space:z_max)
    
    nx = length(x_grid)
    ny = length(y_grid)
    nz = length(z_grid)
    
    grid_elec = zeros(T, nx, ny, nz)      
    # assign ace score for RDS
    for iatom = 1:receptor.natom
        # create atom coordinates of receptor
        x_atom = receptor.xyz[iframe, 3*(iatom-1)+1]
        y_atom = receptor.xyz[iframe, 3*(iatom-1)+2]
        z_atom = receptor.xyz[iframe, 3*(iatom-1)+3]
        
        # determin Real part of RDS
        rcut = 20.0
        ix_min = findfirst(abs.(x_atom .- x_grid) .< rcut)
        iy_min = findfirst(abs.(y_atom .- y_grid) .< rcut)
        iz_min = findfirst(abs.(z_atom .- z_grid) .< rcut)
        ix_max = findlast(abs.(x_atom .- x_grid) .< rcut)
        iy_max = findlast(abs.(y_atom .- y_grid) .< rcut)
        iz_max = findlast(abs.(z_atom .- z_grid) .< rcut)
        
        # ６Åの判別（まだ条件分岐をつくってないので注意）    
        for ix = ix_min:ix_max
            for iy = iy_min:iy_max
                for iz = iz_min:iz_max
                    dist = sqrt((x_atom - x_grid[ix])^2 + (y_atom - y_grid[iy])^2 + (z_atom - z_grid[iz])^2)
                    if dist < 2.0
                        grid_elec[ix, iy, iz] += 0.0
                    elseif dist < 6.0
                        dielec = 4.0
                        grid_elec[ix, iy, iz] += receptor.charge[iatom] / (dielec * dist^2)
                    elseif dist < 8.0
                        dielec = 38.0*dist - 224.0
                        grid_elec[ix, iy, iz] += receptor.charge[iatom] / (dielec * dist^2)
                    else
                        dielec = 80.0
                        grid_elec[ix, iy, iz] += receptor.charge[iatom] / (dielec * dist^2)
                    end
                end
            end
        end
    end
    
    for iatom = 1:receptor.natom
        if receptor.radius[iatom] > 1.0
            rcut = receptor.radius[iatom] * sqrt(0.8)
        else
            rcut = receptor.radius[iatom] * sqrt(1.5)
        end
        
        x_atom = receptor.xyz[iframe, 3*(iatom-1)+1]
        y_atom = receptor.xyz[iframe, 3*(iatom-1)+2]
        z_atom = receptor.xyz[iframe, 3*(iatom-1)+3]
        
        ix_min = findfirst(abs.(x_atom .- x_grid) .< rcut)
        iy_min = findfirst(abs.(y_atom .- y_grid) .< rcut)
        iz_min = findfirst(abs.(z_atom .- z_grid) .< rcut)
        ix_max = findlast(abs.(x_atom .- x_grid) .< rcut)
        iy_max = findlast(abs.(y_atom .- y_grid) .< rcut)
        iz_max = findlast(abs.(z_atom .- z_grid) .< rcut)
        
        for ix = ix_min:ix_max
            for iy = iy_min:iy_max
                for iz = iz_min:iz_max
                    grid_elec[ix, iy, iz] = 0.0
                end
            end
        end
    end  
    
    # generate grid coordinates of ligand at every rotation by quaternions
    grid_charge = zeros(T, nx, ny, nz)
    score = similar(grid_charge)
    nq = size(quaternions, 1)
    @showprogress for iq = 1:nq
        grid_charge .= 0.0
        
        # rotate ligand by quaternions
        ligand_rotated = rotate(ligand, quaternions[iq, :])
        
        # assign ace score for LDS
        for iatom = 1:ligand_rotated.natom
            
            # create atom coordinates of ligand
            x_atom = ligand_rotated.xyz[iframe, 3*(iatom-1)+1]
            y_atom = ligand_rotated.xyz[iframe, 3*(iatom-1)+2]
            z_atom = ligand_rotated.xyz[iframe, 3*(iatom-1)+3]
            
            # invert atoms coordinates into grid coordinates
            ix_min = floor(Int, (x_atom - x_grid[1]) / grid_space)
            iy_min = floor(Int, (y_atom - y_grid[1]) / grid_space)
            iz_min = floor(Int, (z_atom - z_grid[1]) / grid_space)
            
            a = x_atom - x_grid[ix_min] - (0.5*grid_space)
            b = y_atom - y_grid[iy_min] - (0.5*grid_space)
            c = z_atom - z_grid[iz_min] - (0.5*grid_space)
            
            for ix = 1:2
                for iy = 1:2
                    for iz = 1:2
                        X = (-1)^ix * (0.5*grid_space)
                        Y = (-1)^iy * (0.5*grid_space)
                        Z = (-1)^iz * (0.5*grid_space)
                        weight = (X + a)*(Y + b)*(Z + c) / (X*Y*Z) / 8.0
                        grid_charge[ix_min+ix-1, iy_min+iy-1, iz_min+iz-1] = weight * ligand.charge[iatom]
                    end
                end
            end
        end
        
        # compute elec socre with FFT
        t = ifft(fft(grid_elec) .* conj.(fft(conj.(grid_charge))))
        score = real(t)
    end
    
    return score
end