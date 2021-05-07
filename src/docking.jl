function get_atom_type(ta::TrjArray{T,U}) where {T,U}
    ace_score = Array{T}(undef, ta.natom)
    atom_type = Array{String}(undef, ta.natom)
    
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
    return TrjArray(ta, ace_score=ace_score)
end

function get_ace_score(ta::TrjArray{T,U}, iatom) where {T,U}
    atom_type = Array{String}(undef, ta.natom)
    ace_score = Array{Float64}(undef, ta.natom)
    
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

    return ace_score[iatom]
end

function docking_by_desolvation_energy(receptor::TrjArray{T, U}, ligand::TrjArray{T, U}, quaternions; grid_space=1.2, iframe=1, tops=10) where {T, U}
    
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
    #size_ligand = size_ligand*2 

    LDS_x_grid = collect(x_min:grid_space:x_max)
    LDS_y_grid = collect(y_min:grid_space:y_max)
    LDS_z_grid = collect(z_min:grid_space:z_max)

    LDS_nx = length(LDS_x_grid)
    LDS_ny = length(LDS_y_grid)
    LDS_nz = length(LDS_z_grid)

    println("grid size of ligand is ", LDS_nx, "," LDS_ny, "," LDS_nz)

  # extension grid of receptor using size of ligand
    x_min = minimum(receptor.xyz[iframe, 1:3:end]) - size_ligand - grid_space
    y_min = minimum(receptor.xyz[iframe, 2:3:end]) - size_ligand - grid_space
    z_min = minimum(receptor.xyz[iframe, 3:3:end]) - size_ligand - grid_space
    x_max = maximum(receptor.xyz[iframe, 1:3:end]) + size_ligand + grid_space
    y_max = maximum(receptor.xyz[iframe, 2:3:end]) + size_ligand + grid_space
    z_max = maximum(receptor.xyz[iframe, 3:3:end]) + size_ligand + grid_space

    RDS_x_grid = collect(x_min:grid_space:x_max)
    RDS_y_grid = collect(y_min:grid_space:y_max)
    RDS_z_grid = collect(z_min:grid_space:z_max)
    
    RDS_nx = length(RDS_x_grid)
    RDS_ny = length(RDS_y_grid)
    RDS_nz = length(RDS_z_grid)

    println("RDS x grid size is " , RDS_nx)
    println("RDS y grid size is " , RDS_ny)
    println("RDS z grid size is " , RDS_nz)

    grid_RDS = zeros(complex(T), RDS_nx, RDS_ny, RDS_nz)
    grid_LDS = zeros(complex(T), RDS_nx, RDS_ny, RDS_nz)
    
    
  # generate grid coordinates of ligand at every rotation by quaternions
    nq = size(quaternions, 1)
    score_DS_max = Vector{T}(undef, nq)
    score_DS_min = Vector{T}(undef, nq)
    for iq = 1:nq
        grid_RDS .= 0.0 + 0.0im
        grid_LDS .= 0.0 + 0.0im
        
      # rotate ligand by quaternions
        ligand = rotate(ligand, quaternions[iq, :])

      # assign ace score for RDS
        for iatom = 1:receptor.natom

          # create atom coordinates of receptor
            x_atom = receptor.xyz[iframe, 3*(iatom-1)+1]
            y_atom = receptor.xyz[iframe, 3*(iatom-1)+2]
            z_atom = receptor.xyz[iframe, 3*(iatom-1)+3]

          # invert atom coordinates into grid coordinates
            x = round(Int, x_atom / grid_space + RDS_nx/2, RoundNearestTiesAway)
            y = round(Int, y_atom / grid_space + RDS_ny/2, RoundNearestTiesAway)
            z = round(Int, z_atom / grid_space + RDS_nz/2, RoundNearestTiesAway)
        
          # determin Real part of RDS
            ix_min = round(Int, (x_atom - 6) / grid_space + RDS_nx/2, RoundToZero)
            iy_min = round(Int, (y_atom - 6) / grid_space + RDS_ny/2, RoundToZero)
            iz_min = round(Int, (z_atom - 6) / grid_space + RDS_nz/2, RoundToZero)
            ix_max = round(Int, (x_atom + 6) / grid_space + RDS_nx/2, RoundToZero)
            iy_max = round(Int, (y_atom + 6) / grid_space + RDS_ny/2, RoundToZero)
            iz_max = round(Int, (z_atom + 6) / grid_space + RDS_nz/2, RoundToZero)
 

            for ix = ix_min:ix_max
                for iy = iy_min:iy_max
                    for iz = iz_min:iz_max
                        #println(grid_RDS[Int(ix), Int(iy), Int(iz)])
                        dist = sqrt((x_atom - RDS_x_grid[ix])^2 + (y_atom - RDS_y_grid[iy])^2 + (z_atom - RDS_z_grid[iz])^2)
                        if dist <= 6
                            grid_RDS[ix, iy, iz] += receptor.ace_score[iatom]
                        end
                    end
                end
            end


          # determine Imag part of RDS
            grid_RDS[x, y, z] += 1.0im


        end

        @show sum(imag(grid_RDS))

      # assign ace score for LDS
        for iatom = 1:ligand.natom

          # create atom coordinates of ligand
            x_atom = ligand.xyz[iframe, 3*(iatom-1)+1]
            y_atom = ligand.xyz[iframe, 3*(iatom-1)+2]
            z_atom = ligand.xyz[iframe, 3*(iatom-1)+3]
            
          # invert atoms coordinates into grid coordinates
            x = round(Int, x_atom / grid_space + RDS_nx/2, RoundNearestTiesAway)
            y = round(Int, y_atom / grid_space + RDS_ny/2, RoundNearestTiesAway)
            z = round(Int, z_atom / grid_space + RDS_nz/2, RoundNearestTiesAway)
      
          # determin Real part of LDS
            ix_min = round(Int, (x_atom - 6) / grid_space + RDS_nx/2, RoundToZero)
            iy_min = round(Int, (y_atom - 6) / grid_space + RDS_ny/2, RoundToZero)
            iz_min = round(Int, (z_atom - 6) / grid_space + RDS_nz/2, RoundToZero)
            ix_max = round(Int, (x_atom + 6) / grid_space + RDS_nx/2, RoundToZero)
            iy_max = round(Int, (y_atom + 6) / grid_space + RDS_ny/2, RoundToZero)
            iz_max = round(Int, (z_atom + 6) / grid_space + RDS_nz/2, RoundToZero)

            for ix = ix_min:ix_max
                for iy = iy_min:iy_max
                    for iz = iz_min:iz_max
                        dist = sqrt((x_atom - RDS_x_grid[ix])^2 + (y_atom - RDS_y_grid[iy])^2 + (z_atom - RDS_z_grid[iz])^2)
                        if dist <= 6
                            grid_LDS[ix, iy, iz] += ligand.ace_score[iatom]
                        end
                    end
                end
            end

            # determine Imag part of LDS
            grid_LDS[x, y, z] += 1.0im
        end

      # compute DS socre with FFT
        t = ifft(ifft(grid_RDS) .* fft(grid_LDS))
        score_DS_max[iq] = maximum(imag(t)) / 2 * RDS_nx * RDS_ny * RDS_nz
        score_DS_min[iq] = minimum(imag(t)) / 2 * LDS_nx * LDS_ny * LDS_nz

    end

    return score_DS_max, score_DS_min
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
            point[iframe] += ta.xyz[iframe, 3*(iatom-1)+1]
            point[2] += ta.xyz[iframe, 3*(iatom-1)+2]
            point[3] += ta.xyz[iframe, 3*(iatom-1)+3]
            for j in 1:length(neighbor_list_iatom)
                jatom = neighbor_list_iatom[j]
                d = 0.0
                d += (point[iframe] - ta.xyz[iframe, 3*(jatom-1)+1])^2
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

        dx = x - x_grid[iframe]
        ix_min = floor(U, (dx - rcut)/grid_space) + 1
        ix_min = max(ix_min, 1)
        ix_max = floor(U, (dx + rcut)/grid_space) + 2
        ix_max = min(ix_max, nx)

        dy = y - y_grid[iframe]
        iy_min = floor(U, (dy - rcut)/grid_space) + 1
        iy_min = max(iy_min, 1)
        iy_max = floor(U, (dy + rcut)/grid_space) + 2
        iy_max = min(iy_max, ny)

        dz = z - z_grid[iframe]
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

        dx = x - x_grid[iframe]
        ix_min = floor(U, (dx - rcut)/grid_space) + 1
        ix_min = max(ix_min, 1)
        ix_max = floor(U, (dx + rcut)/grid_space) + 2
        ix_max = min(ix_max, nx)

        dy = y - y_grid[iframe]
        iy_min = floor(U, (dy - rcut)/grid_space) + 1
        iy_min = max(iy_min, 1)
        iy_max = floor(U, (dy + rcut)/grid_space) + 2
        iy_max = min(iy_max, ny)

        dz = z - z_grid[iframe]
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
        dx_estimated = id[iframe]
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

function dock_fft(receptor::TrjArray{T, U}, ligand::TrjArray{T, U}, quaternions; grid_space=1.2, iframe=1, tops=10) where {T, U}
    # generate grid coordinates for receptor
    receptor2, _dummy = decenter(receptor)
    ligand2, _dummy = decenter(ligand)

    x_min, x_max = minimum(ligand2.xyz[iframe, 1:3:end]), maximum(ligand2.xyz[iframe, 1:3:end])
    y_min, y_max = minimum(ligand2.xyz[iframe, 2:3:end]), maximum(ligand2.xyz[iframe, 2:3:end])
    z_min, z_max = minimum(ligand2.xyz[iframe, 3:3:end]), maximum(ligand2.xyz[iframe, 3:3:end])
    size_ligand = sqrt((x_max - x_min)^2 + (y_max - y_min)^2 + (z_max - z_min)^2)
    size_ligand = size_ligand*2

    x_min = minimum(receptor2.xyz[iframe, 1:3:end]) - size_ligand - grid_space
    y_min = minimum(receptor2.xyz[iframe, 2:3:end]) - size_ligand - grid_space
    z_min = minimum(receptor2.xyz[iframe, 3:3:end]) - size_ligand - grid_space

    x_max = maximum(receptor2.xyz[iframe, 1:3:end]) + size_ligand + grid_space
    y_max = maximum(receptor2.xyz[iframe, 2:3:end]) + size_ligand + grid_space
    z_max = maximum(receptor2.xyz[iframe, 3:3:end]) + size_ligand + grid_space

    x_grid = collect(x_min:grid_space:x_max)
    y_grid = collect(y_min:grid_space:y_max)
    z_grid = collect(z_min:grid_space:z_max)

    nx, ny, nz = length(x_grid), length(y_grid), length(z_grid)

    # spape complementarity of receptor
    iatom_surface = receptor2.sasa .> 1.0
    iatom_core = .!iatom_surface
    rcut1 = zeros(T, receptor2.natom)
    rcut2 = zeros(T, receptor2.natom)
    
    rcut1[iatom_core] .= receptor2.radius[iatom_core] * sqrt(1.5)
    rcut2[iatom_core] .= -1.0
    rcut1[iatom_surface] .= receptor2.radius[iatom_surface] * sqrt(0.8)
    rcut2[iatom_surface] .= receptor2.radius[iatom_surface] .+ 3.4

    grid_RSC = zeros(complex(T), nx, ny, nz)
    assign_shape_complementarity!(grid_RSC, receptor2, grid_space, rcut1, rcut2, x_grid, y_grid, z_grid, iframe)

    # spape complementarity of ligand
    iatom_surface = ligand2.sasa .> 1.0
    iatom_core = .!iatom_surface
    rcut1 = zeros(T, ligand2.natom)
    rcut2 = zeros(T, ligand2.natom)
    rcut1[iatom_core] .= ligand2.radius[iatom_core] * sqrt(1.5)
    rcut2[iatom_core] .= -1.0
    rcut1[iatom_surface] .= -1.0
    rcut2[iatom_surface] .= ligand2.radius[iatom_surface]
    
    grid_LSC = zeros(complex(T), nx, ny, nz)

    # compute score with FFT
    nq = size(quaternions, 1)
    s = @showprogress pmap(q -> compute_docking_score_with_fft(quaternions[q, :], grid_RSC, grid_LSC, ligand2, grid_space, rcut1, rcut2, x_grid, y_grid, z_grid, iframe, tops, q), 1:nq)
    
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
        score_tops[t] = result_tops[t][iframe]
        trans_tops[t, 1] = result_tops[t][2]
        trans_tops[t, 2] = result_tops[t][3]
        trans_tops[t, 3] = result_tops[t][4]
        quate_tops[t, :] .= quaternions[result_tops[t][5], :]
    end

    itop = 1
    ligand_return = rotate(ligand2, quate_tops[itop, :])
    dx = (trans_tops[itop, 1]-1) * grid_space
    if dx > (nx*grid_space / 2.0)
        dx = dx - (nx*grid_space)
    end
    dy = (trans_tops[itop, 2]-1) * grid_space
    if dy > (ny*grid_space / 2.0)
        dy = dy - (ny*grid_space)
    end
    dz = (trans_tops[itop, 3]-1) * grid_space
    if dz > (nz*grid_space / 2.0)
        dz = dz - (nz*grid_space)
    end
    ligand_return.xyz[iframe, 1:3:end] .+= dx
    ligand_return.xyz[iframe, 2:3:end] .+= dy
    ligand_return.xyz[iframe, 3:3:end] .+= dz

    for itop = 2:tops
        ligand_tmp = rotate(ligand2, quate_tops[itop, :])
        dx = (trans_tops[itop, 1]-1) * grid_space
        if dx > (nx*grid_space / 2.0)
            dx = dx - (nx*grid_space)
        end
        dy = (trans_tops[itop, 2]-1) * grid_space
        if dy > (ny*grid_space / 2.0)
            dy = dy - (ny*grid_space)
        end
        dz = (trans_tops[itop, 3]-1) * grid_space
        if dz > (nz*grid_space / 2.0)
            dz = dz - (nz*grid_space)
        end
        ligand_tmp.xyz[iframe, 1:3:end] .+= dx
        ligand_tmp.xyz[iframe, 2:3:end] .+= dy
        ligand_tmp.xyz[iframe, 3:3:end] .+= dz
        ligand_return = [ligand_return; ligand_tmp]
    end

    return (receptor=receptor2, ligand=ligand_return, score=score_tops, grid_RSC=grid_RSC, grid_LSC=grid_LSC)
end


function dock_multimer(receptor::TrjArray{T, U}; rot_space=10.0, radius=3.0:1.2:20.0, grid_space=1.2, iframe=1, tops=10, nfold=3) where {T, U}
    # generate grid coordinates for receptor
    receptor2, _dummy = decenter(receptor)

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

    # spape complementarity of receptor
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
            R = [cos(2*pi/nfold) -sin(2*pi/nfold) 0.0; sin(2*pi/nfold) cos(2*pi/nfold) 0.0; 0.0 0.0 1.0]
            ligand_rot = rotate_with_matrix(receptor_rot, R)
            assign_shape_complementarity!(grid_RSC, receptor_rot, grid_space, rcut1, rcut2, x_grid, y_grid, z_grid, iframe)
            assign_shape_complementarity!(grid_LSC, ligand_rot, grid_space, rcut1, rcut2, x_grid, y_grid, z_grid, iframe)
            score = compute_docking_score_on_xyplane(grid_RSC, grid_LSC)
            score_max = maximum(score)
            if score_best < score_max
                id = argmax(score)
                dx_estimated = id[iframe]
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