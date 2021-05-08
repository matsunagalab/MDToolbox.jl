import Base: convert, copy, show, getindex, lastindex, isempty,
             size, length, eachindex, ==, isequal, hash, vcat, hcat, merge, map

abstract type AbstractTrajectory end

struct TrjArray{T, U} <: AbstractTrajectory
    #DONE: support for x, y, z = nothing, nothing, nothing
    #DONE: natom, nframe
    #TODO: revert to parametric type, fix slowness of constructor
    #TODO: simple and easy constructions for other functions which outputs a new trjarray (such as centerofmass)
    #TODO: should allow one-dimensional array (vector) for x, y, z constructor input
    xyz::AbstractArray{T}
    boxsize::AbstractArray{T}
    chainname::AbstractArray{String}
    chainid::AbstractArray{U}
    resname::AbstractArray{String}
    resid::AbstractArray{U}
    atomname::AbstractArray{String}
    atomid::AbstractArray{U}
    mass::AbstractArray{T}
    radius::AbstractArray{T}
    charge::AbstractArray{T}
    sasa::AbstractArray{T}
    list_bond::AbstractArray{U}
    list_angle::AbstractArray{U}
    list_dihedral::AbstractArray{U}
    list_improper::AbstractArray{U}
    list_cmap::AbstractArray{U}
    natom::Int64
    nframe::Int64

    function TrjArray{T, U}(xyz, boxsize, chainname, chainid, resname, resid, atomname, atomid, mass, radius, charge, sasa,
                      list_bond, list_angle, list_dihedral, list_improper, list_cmap) where {T, U}
        # nrow, ncol = (size(trj, 1), size(trj, 2))
        # natom = Int64(ncol/3)
        # ischecked && return new(trj, atomname, atomid, meta)
        # natom != length(chainname) && throw(DimensionMismatch("chainname must match width of trajectory"))
        # natom != length(chainid) && throw(DimensionMismatch("chainid must match width of trajectory"))

        # check size and define nframe
        nframe = 0
        if !isempty(xyz)
            nframe = size(xyz, 1)
            if !isempty(boxsize)
                @assert nframe == size(boxsize, 1)
                @assert 3 == size(boxsize, 2)
            end
        end
        nframe = U(nframe)

        # check size and define natom
        natom = 0
        vec_collection = [chainname, chainid, resname, resid, atomname, atomid, mass, radius, charge, sasa]
        if !isempty(xyz)
            natom = U(size(xyz, 2)/3)
        else
            for vec in vec_collection
                if !isempty(vec)
                    natom = length(vec)
                    break
                end
            end
        end
        if !isempty(xyz)
            @assert (natom*3) == size(xyz, 2)
        end
        for vec in vec_collection
            if !isempty(vec)
                @assert natom == length(vec)
            end
        end
        natom = U(natom)

        return new(xyz, boxsize, chainname, chainid, resname, resid, atomname, atomid, mass, radius, charge, sasa,
                   list_bond, list_angle, list_dihedral, list_improper, list_cmap, natom, nframe)
    end
end

###### outer constructors ########
function TrjArray{T, U}(;xyz = Matrix{T}(undef, 0, 0), 
                  boxsize = Matrix{T}(undef, 0, 0),
                  chainname = Vector{String}(undef, 0), chainid = Vector{U}(undef, 0),
                  resname = Vector{String}(undef, 0), resid = Vector{U}(undef, 0),
                  atomname = Vector{String}(undef, 0), atomid = Vector{U}(undef, 0),
                  mass = Vector{T}(undef, 0), radius = Vector{T}(undef, 0), charge = Vector{T}(undef, 0), sasa = Vector{T}(undef, 0),
                  list_bond = Matrix{U}(undef, 0, 0), list_angle = Matrix{U}(undef, 0, 0),
                  list_dihedral = Matrix{U}(undef, 0, 0), list_improper = Matrix{U}(undef, 0, 0),
                  list_cmap = Matrix{U}(undef, 0, 0)) where {T, U}
    TrjArray{T, U}(xyz, boxsize, chainname, chainid, resname, resid, atomname, atomid, mass, radius, charge, sasa,
             list_bond, list_angle, list_dihedral, list_improper, list_cmap)
end

TrjArray(ta::TrjArray{T, U}; xyz = nothing, 
    boxsize = nothing,
    chainname = nothing, chainid = nothing,
    resname = nothing, resid = nothing,
    atomname = nothing, atomid = nothing,
    mass = nothing, radius = nothing, charge = nothing, sasa = nothing,
    list_bond = nothing, list_angle = nothing,
    list_dihedral = nothing, list_improper = nothing,
    list_cmap = nothing) where {T, U} = TrjArray{T, U}(
        xyz = isnothing(xyz) ? ta.xyz : xyz,
        boxsize = isnothing(boxsize) ? ta.boxsize : boxsize,
        chainname = isnothing(chainname) ? ta.chainname : chainname,
        chainid = isnothing(chainid) ? ta.chainid : chainid,
        resname = isnothing(resname) ? ta.resname : resname,
        resid = isnothing(resid) ? ta.resid : resid,
        atomname = isnothing(atomname) ? ta.atomname : atomname,
        atomid = isnothing(atomid) ? ta.atomid : atomid,
        mass = isnothing(mass) ? ta.mass : mass, 
        radius = isnothing(radius) ? ta.radius : radius, 
        charge = isnothing(charge) ? ta.charge : charge,
        sasa = isnothing(sasa) ? ta.sasa : sasa,
        list_bond = isnothing(list_bond) ? ta.list_bond : list_bond, 
        list_angle = isnothing(list_angle) ? ta.list_angle : list_angle, 
        list_dihedral = isnothing(list_dihedral) ? ta.list_dihedral : list_dihedral, 
        list_improper = isnothing(list_improper) ? ta.list_improper : list_improper, 
        list_cmap = isnothing(list_cmap) ? ta.list_cmap : list_cmap)

TrjArray(xyz::AbstractArray{T}, boxsize::AbstractArray{T}, ta::TrjArray{T, U}) where {T, U} =
       TrjArray{T, U}(xyz = xyz, boxsize = boxsize,
                      chainname = ta.chainname, chainid = ta.chainid,
                      resname = ta.resname, resid = ta.resid,
                      atomname = ta.atomname, atomid = ta.atomid,
                      mass = ta.mass, radius = ta.radius, charge = ta.charge, sasa = ta.sasa,
                      list_bond = ta.list_bond, list_angle = ta.list_angle,
                      list_dihedral = ta.list_dihedral, list_improper = ta.list_improper,
                      list_cmap = ta.list_cmap)

TrjArray(xyz::AbstractArray{T}, ta::TrjArray{T, U}) where {T, U} =
       TrjArray{T, U}(xyz = xyz, boxsize = ta.boxsize,
                      chainname = ta.chainname, chainid = ta.chainid,
                      resname = ta.resname, resid = ta.resid,
                      atomname = ta.atomname, atomid = ta.atomid,
                      mass = ta.mass, radius = ta.radius, charge = ta.charge, sasa = ta.sasa,
                      list_bond = ta.list_bond, list_angle = ta.list_angle,
                      list_dihedral = ta.list_dihedral, list_improper = ta.list_improper,
                      list_cmap = ta.list_cmap)

TrjArray(x...; y...) = TrjArray{Float64,Int64}(x...; y...)
TrjArray() = TrjArray{Float64, Int64}()

###### size, length #################
size(ta::TrjArray) = (ta.nframe, ta.natom)

size(ta::TrjArray, dim::Int) = (dim == 1) ? ta.nframe : (dim == 2) ? ta.natom : 1

length(ta::TrjArray) = ta.nframe

###### getindex #################
# all element indexing
getindex(ta::TrjArray, ::Colon) = ta
getindex(ta::TrjArray, ::Colon, ::Colon) = ta

# single row
function getindex(ta::TrjArray{T, U}, n::Int) where {T, U}
    if iszero(n)
        return TrjArray(ta, xyz = Array{T, 2}(undef, 0, 0), boxsize = Array{T, 2}(undef, 0, 0)) ## TODO need to fix here
    elseif !isempty(ta.boxsize)
        return TrjArray(ta, xyz = ta.xyz[n:n, :], boxsize = ta.boxsize[n:n, :])
    else
        return TrjArray(ta, xyz = ta.xyz[n:n, :])
    end
end
getindex(ta::TrjArray{T, U}, n::Int, ::Colon) where {T, U} = getindex(ta, n)

# range of rows
function getindex(ta::TrjArray{T, U}, r::UnitRange{Int}) where {T, U}
    if !isempty(ta.boxsize)
        TrjArray(ta, xyz = ta.xyz[r, :], boxsize = ta.boxsize[r, :])
    else
        TrjArray(ta, xyz = ta.xyz[r, :])
    end
end
getindex(ta::TrjArray{T, U}, r::UnitRange{Int}, ::Colon) where {T, U} = getindex(ta, r)

# array of rows (integer)
function getindex(ta::TrjArray{T, U}, a::AbstractVector{S}) where {T, U, S <: Integer}
    if !isempty(ta.boxsize)
        TrjArray(ta, xyz = ta.xyz[a, :], boxsize = ta.boxsize[a, :])
    else
        TrjArray(ta, xyz = ta.xyz[a, :])
    end
end
getindex(ta::TrjArray{T, U}, a::AbstractVector{S}, ::Colon) where {T, U, S <: Integer} = getindex(ta, a)

# array of rows (bool)
#getindex(ta::TrjArray, a::AbstractVector{S}) where {S <: Bool} = TrjArray(ta.x[a, :], ta.y[a, :], ta.z[a, :], ta)
#getindex(ta::TrjArray, a::AbstractVector{S}, ::Colon) where {S <: Bool} = getindex(ta, a)
function getindex(ta::TrjArray{T, U}, a::AbstractVector{Bool}) where {T, U}
    if !isempty(ta.boxsize)
        TrjArray(ta, xyz = ta.xyz[a, :], boxsize = ta.boxsize[a, :])
    else
        TrjArray(ta, xyz = ta.xyz[a, :])
    end
end
getindex(ta::TrjArray{T, U}, a::AbstractVector{Bool}, ::Colon) where {T, U} = getindex(ta, a)

# define function for updating list_bond, list_angle, list_dihedral, list_improper, list_cmap
function reindex_list(natom, list_some, index)
  if isempty(list_some) || isempty(index)
    return Matrix{eltype(list_some)}(undef, 0, 0)
  end
  if typeof(index) <: AbstractVector{Bool}
    index_logical = index
    reindex = zeros(eltype(list_some), natom)
    reindex[index] = 1:sum(index)
  else
    index_logical = falses(natom)
    index_logical[index] .= true
    reindex = zeros(eltype(list_some), natom)
    reindex[index] = 1:length(index)
  end
  nlist = size(list_some, 1)
  list_some_new = similar(list_some)
  icount = 0
  for i = 1:nlist
    if all(index_logical[list_some[i, :]])
      icount += 1
      list_some_new[icount, :] .= reindex[list_some[i, :]]
    end
  end
  if icount == 0
    return Matrix{eltype(list_some)}(undef, 0, 0)
  else
    return list_some_new[1:icount, :]
  end
end

# single column
to3(r::Int) = [3*(r-1)+1, 3*(r-1)+2, 3*(r-1)+3]

getindex(ta::TrjArray{T, U}, ::Colon, r::Int) where {T, U} = TrjArray{T, U}(
             xyz = isempty(ta.xyz) ? ta.xyz : ta.xyz[:, to3(r)],
             boxsize = ta.boxsize,
             chainname = isempty(ta.chainname) ? ta.chainname : ta.chainname[r:r],
             chainid = isempty(ta.chainid) ? ta.chainid : ta.chainid[r:r],
             resname = isempty(ta.resname) ? ta.resname : ta.resname[r:r],
             resid = isempty(ta.resid) ? ta.resid : ta.resid[r:r],
             atomname = isempty(ta.atomname) ? ta.atomname : ta.atomname[r:r],
             #atomid = isempty(ta.atomid) ? ta.atomid : ta.atomid[r:r],
             atomid = isempty(ta.atomid) ? ta.atomid : collect(Int64, 1:1),
             mass = isempty(ta.mass) ? ta.mass : ta.mass[r:r],
             radius = isempty(ta.radius) ? ta.radius : ta.radius[r:r],
             charge = isempty(ta.charge) ? ta.charge : ta.charge[r:r],
             sasa = isempty(ta.sasa) ? ta.sasa : ta.sasa[r:r],
             list_bond = isempty(ta.list_bond) ? ta.list_bond : reindex_list(ta.natom, ta.list_bond, r:r),
             list_angle = isempty(ta.list_angle) ? ta.list_angle : reindex_list(ta.natom, ta.list_angle, r:r),
             list_dihedral = isempty(ta.list_dihedral) ? ta.list_dihedral : reindex_list(ta.natom, ta.list_dihedral, r:r),
             list_improper = isempty(ta.list_improper) ? ta.list_improper : reindex_list(ta.natom, ta.list_improper, r:r),
             list_cmap = isempty(ta.list_cmap) ? ta.list_cmap : reindex_list(ta.natom, ta.list_cmap, r:r))

# range of columns
function to3(r::UnitRange{Int})
    n = length(r)
    r3 = similar(r, 3*n)
    r3[1:3:end] .= 3 .* (r .- 1) .+ 1
    r3[2:3:end] .= 3 .* (r .- 1) .+ 2
    r3[3:3:end] .= 3 .* (r .- 1) .+ 3
    return r3
end

getindex(ta::TrjArray{T, U}, ::Colon, r::UnitRange{Int}) where {T, U} = TrjArray{T, U}(
             xyz = isempty(ta.xyz) ? ta.xyz : ta.xyz[:, to3(r)],
             boxsize = ta.boxsize,
             chainname = isempty(ta.chainname) ? ta.chainname : ta.chainname[r],
             chainid = isempty(ta.chainid) ? ta.chainid : ta.chainid[r],
             resname = isempty(ta.resname) ? ta.resname : ta.resname[r],
             resid = isempty(ta.resid) ? ta.resid : ta.resid[r],
             atomname = isempty(ta.atomname) ? ta.atomname : ta.atomname[r],
             #atomid = isempty(ta.atomid) ? ta.atomid : ta.atomid[r],
             atomid = isempty(ta.atomid) ? ta.atomid : collect(Int64, 1:length(r)),
             mass = isempty(ta.mass) ? ta.mass : ta.mass[r],
             radius = isempty(ta.radius) ? ta.radius : ta.radius[r],
             charge = isempty(ta.charge) ? ta.charge : ta.charge[r],
             sasa = isempty(ta.sasa) ? ta.sasa : ta.sasa[r],
             list_bond = isempty(ta.list_bond) ? ta.list_bond : reindex_list(ta.natom, ta.list_bond, r),
             list_angle = isempty(ta.list_angle) ? ta.list_angle : reindex_list(ta.natom, ta.list_angle, r),
             list_dihedral = isempty(ta.list_dihedral) ? ta.list_dihedral : reindex_list(ta.natom, ta.list_dihedral, r),
             list_improper = isempty(ta.list_improper) ? ta.list_improper : reindex_list(ta.natom, ta.list_improper, r),
             list_cmap = isempty(ta.list_cmap) ? ta.list_cmap : reindex_list(ta.natom, ta.list_cmap, r))

# array of columns (integer)
# range of columns
function to3(r::AbstractVector{S}) where {S <: Integer}
    n = length(r)
    r3 = similar(r, 3*n)
    r3[1:3:end] .= 3 .* (r .- 1) .+ 1
    r3[2:3:end] .= 3 .* (r .- 1) .+ 2
    r3[3:3:end] .= 3 .* (r .- 1) .+ 3
    return r3
end

getindex(ta::TrjArray{T, U}, ::Colon, r::AbstractVector{S}) where {T, U, S <: Integer} = TrjArray{T, U}(
             xyz = isempty(ta.xyz) ? ta.xyz : ta.xyz[:, to3(r)],
             boxsize = ta.boxsize,
             chainname = isempty(ta.chainname) ? ta.chainname : ta.chainname[r],
             chainid = isempty(ta.chainid) ? ta.chainid : ta.chainid[r],
             resname = isempty(ta.resname) ? ta.resname : ta.resname[r],
             resid = isempty(ta.resid) ? ta.resid : ta.resid[r],
             atomname = isempty(ta.atomname) ? ta.atomname : ta.atomname[r],
             #atomid = isempty(ta.atomid) ? ta.atomid : ta.atomid[r],
             atomid = isempty(ta.atomid) ? ta.atomid : collect(Int64, 1:length(r)),
             mass = isempty(ta.mass) ? ta.mass : ta.mass[r],
             radius = isempty(ta.radius) ? ta.radius : ta.radius[r],
             charge = isempty(ta.charge) ? ta.charge : ta.charge[r],
             sasa = isempty(ta.sasa) ? ta.sasa : ta.sasa[r],
             list_bond = isempty(ta.list_bond) ? ta.list_bond : reindex_list(ta.natom, ta.list_bond, r),
             list_angle = isempty(ta.list_angle) ? ta.list_angle : reindex_list(ta.natom, ta.list_angle, r),
             list_dihedral = isempty(ta.list_dihedral) ? ta.list_dihedral : reindex_list(ta.natom, ta.list_dihedral, r),
             list_improper = isempty(ta.list_improper) ? ta.list_improper : reindex_list(ta.natom, ta.list_improper, r),
             list_cmap = isempty(ta.list_cmap) ? ta.list_cmap : reindex_list(ta.natom, ta.list_cmap, r))

# array of rows (bool)
function to3(r::AbstractVector{Bool})
    n = length(r)
    r3 = similar(r, 3*n)
    r3[1:3:end] .= r
    r3[2:3:end] .= r
    r3[3:3:end] .= r
    return r3
end

getindex(ta::TrjArray{T, U}, ::Colon, r::AbstractVector{Bool}) where {T, U} = TrjArray{T, U}(
             xyz = isempty(ta.xyz) ? ta.xyz : ta.xyz[:, to3(r)],
             boxsize = ta.boxsize,
             chainname = isempty(ta.chainname) ? ta.chainname : ta.chainname[r],
             chainid = isempty(ta.chainid) ? ta.chainid : ta.chainid[r],
             resname = isempty(ta.resname) ? ta.resname : ta.resname[r],
             resid = isempty(ta.resid) ? ta.resid : ta.resid[r],
             atomname = isempty(ta.atomname) ? ta.atomname : ta.atomname[r],
             #atomid = isempty(ta.atomid) ? ta.atomid : ta.atomid[r],
             atomid = isempty(ta.atomid) ? ta.atomid : collect(Int64, 1:sum(r)),
             mass = isempty(ta.mass) ? ta.mass : ta.mass[r],
             radius = isempty(ta.radius) ? ta.radius : ta.radius[r],
             charge = isempty(ta.charge) ? ta.charge : ta.charge[r],
             sasa = isempty(ta.sasa) ? ta.sasa : ta.sasa[r],
             list_bond = isempty(ta.list_bond) ? ta.list_bond : reindex_list(ta.natom, ta.list_bond, r),
             list_angle = isempty(ta.list_angle) ? ta.list_angle : reindex_list(ta.natom, ta.list_angle, r),
             list_dihedral = isempty(ta.list_dihedral) ? ta.list_dihedral : reindex_list(ta.natom, ta.list_dihedral, r),
             list_improper = isempty(ta.list_improper) ? ta.list_improper : reindex_list(ta.natom, ta.list_improper, r),
             list_cmap = isempty(ta.list_cmap) ? ta.list_cmap : reindex_list(ta.natom, ta.list_cmap, r))

# combinations
getindex(ta::TrjArray{T, U}, rows, cols) where {T, U} = ta[rows, :][:, cols]

# end index
lastindex(ta::TrjArray, dims::Int = 1) = (dims == 1) ? ta.nframe : (dims == 2) ? ta.natom : 1

###### atom selection #################
function unfold(A)
    V = []
    if typeof(A) <: AbstractString
        push!(V, A)
    else
        for x in A
            if x === A
                push!(V, x)
            else
                append!(V, unfold(x))
            end
        end
    end
    V
end

function match_query(some_array, query)
    natom = length(some_array)
    #index = fill(false, natom)
    index = falses(natom)
    if isempty(some_array)
        return index
    elseif typeof(some_array[1]) == String
        query = split(query)
    else
        query = Meta.eval(Meta.parse(replace("[" * query * "]", r"(\d)\s" => s"\1, ")))
    end
    query = unique(unfold(query))
    for q in query
        index = index .| (some_array .== q)
    end
    index
end

function replace_ex!(ex::Expr, ta::TrjArray)
    chainid = ta.chainid
    chainname = ta.chainname
    resid = ta.resid
    resname = ta.resname
    atomid = ta.atomid
    atomname = ta.atomname
    for (i, arg) in enumerate(ex.args)
        if arg == :(match_query)
            ex.args[i] = :($match_query)
        elseif arg == :(chainid)
            ex.args[i] = :($chainid)
        elseif arg == :(chainname)
            ex.args[i] = :($chainname)
        elseif arg == :(resid)
            ex.args[i] = :($resid)
        elseif arg == :(resname)
            ex.args[i] = :($resname)
        elseif arg == :(atomid)
            ex.args[i] = :($atomid)
        elseif arg == :(atomname)
            ex.args[i] = :($atomname)
        end
        if isa(arg, Expr)
            replace_ex!(arg, ta)
        end
    end
    ex
end

function select_atom(ta::TrjArray, s::AbstractString)
    #s = strip(s)

    # alias
    s = replace(s, "waters" => "water")
    s = replace(s, "solvents" => "solvent")
    s = replace(s, "proteins" => "protein")

    """
    this code is licensed under the GNU LGPL version 2.1, or at your option a later version of the license.  (Copyright (c) R.T McGibbon et al.), see LICENSE.md
    """
    s = replace(s, "water" => "resname H2O SOL WAT HHO OHH HOH OH2 TIP TIP2 TIP3 TIP4")

    """
    this code is licensed under the GNU LGPL version 2.1, or at your option a later version of the license.  (Copyright (c) R.T McGibbon et al.), see LICENSE.md
    """
    s = replace(s, "hydrogen" => "atomname 1H 2H 3H H1 H2 H3 H4 H5 H HA 1HA 2HA HB 1HB " *
    "2HB 3HB HG 1HG 2HG 1HG1 2HG1 3HG1 1HG2 2HG2 3HG2 1HD 2HD 1HD1 2HD1 3HD1 1HD2 2HD2 " *
    "3HD2 1HH1 2HH1 1HH2 2HH2 HE 1HE 2HE 3HE 1HE1 2HE1 1HE2 2HE2 HD1 HD2 HE1 HNZ1 HNZ2 " *
    "HNZ3 HG1 HH HH1 HH2 ")

    """
    this code is licensed under the GNU LGPL version 2.1, or at your option a later version of the license.  (Copyright (c) R.T McGibbon et al.), see LICENSE.md
    """
    s = replace(s, "solvent" => "resname 118 119 1AL 1CU 2FK 2HP 2OF 3CO 3MT 3NI 3OF 4MO " * 
    "543 6MO ACT AG AL ALF ATH AU AU3 AUC AZI Ag BA BAR BCT BEF BF4 BO4 BR BS3 BSY Be CA " *
    "CA+2 Ca+2 CAC CAD CAL CD CD1 CD3 CD5 CE CES CHT CL CL- CLA Cl- CO CO3 CO5 CON CR CS " *
    "CSB CU CU1 CU3 CUA CUZ CYN Cl- Cr DME DMI DSC DTI DY E4N EDR EMC ER3 EU EU3 F FE FE2 " *
    "FPO GA GD3 GEP HAI HG HGC HOH IN IOD ION IR IR3 IRI IUM K K+ KO4 LA LCO LCP LI LIT LU " *
    "MAC MG MH2 MH3 MLI MMC MN MN3 MN5 MN6 MO1 MO2 MO3 MO4 MO5 MO6 MOO MOS MOW MW1 MW2 MW3 " *
    "NA NA+2 NA2 NA5 NA6 NAO NAW Na+2 NET NH4 NI NI1 NI2 NI3 NO2 NO3 NRU Na+ O4M OAA OC1 OC2 " *
    "OC3 OC4 OC5 OC6 OC7 OC8 OCL OCM OCN OCO OF1 OF2 OF3 OH OS OS4 OXL PB PBM PD PER PI PO3 " *
    "PO4 POT PR PT PT4 PTN RB RH3 RHD RU RUB Ra SB SCN SE4 SEK SM SMO SO3 SO4 SOD SR Sm Sn " *
    "T1A TB TBA TCN TEA THE TL TMA TRA UNX V V2+ VN3 VO4 W WO5 Y1 YB YB2 YH YT3 ZN ZN2 ZN3 ZNA ZNO ZO3")

    """
    this code is licensed under the GNU LGPL version 2.1, or at your option a later version of the license.  (Copyright (c) R.T McGibbon et al.), see LICENSE.md
    """
    s = replace(s, "protein" => "resname ACE NME 00C 01W 02K 02L 03Y 07O 08P 0A0 0A1 0A2 " *
    "0A8 0AA 0AB 0AC 0AF 0AG 0AH 0AK 0BN 0CS 0E5 0EA 0FL 0NC 0WZ 0Y8 143 193 1OP 1PA 1PI " *
    "1TQ 1TY 1X6 200 23F 23S 26B 2AD 2AG 2AO 2AS 2CO 2DO 2FM 2HF 2KK 2KP 2LU 2ML 2MR 2MT " *
    "2OR 2PI 2QZ 2R3 2SI 2TL 2TY 2VA 2XA 32S 32T 33X 3AH 3AR 3CF 3GA 3MD 3NF 3QN 3TY 3XH " *
    "4BF 4CF 4CY 4DP 4FB 4FW 4HT 4IN 4MM 4PH 4U7 56A 5AB 5CS 5CW 5HP 6CL 6CW 6GL 6HN 7JA " *
    "9NE 9NF 9NR 9NV A5N A66 AA3 AA4 AAR AB7 ABA ACB ACL ADD AEA AEI AFA AGM AGT AHB AHH " *
    "AHO AHP AHS AHT AIB AKL AKZ ALA ALC ALM ALN ALO ALS ALT ALV ALY AN8 APE APH API APK " *
    "APM APP AR2 AR4 AR7 ARG ARM ARO ARV AS2 AS9 ASA ASB ASI ASK ASL ASM ASN ASP ASQ ASX " *
    "AVN AYA AZK AZS AZY B1F B2A B2F B2I B2V B3A B3D B3E B3K B3L B3M B3Q B3S B3T B3U B3X " *
    "B3Y BB6 BB7 BB8 BB9 BBC BCS BE2 BFD BG1 BH2 BHD BIF BIL BIU BJH BL2 BLE BLY BMT BNN " *
    "BNO BOR BPE BSE BTA BTC BTR BUC BUG C1X C22 C3Y C4R C5C C66 C6C CAF CAL CAS CAV CAY " *
    "CCL CCS CDE CDV CEA CGA CGU CHF CHG CHP CHS CIR CLE CLG CLH CME CMH CML CMT CPC CPI " *
    "CR5 CS0 CS1 CS3 CS4 CSA CSB CSD CSE CSJ CSO CSP CSR CSS CSU CSW CSX CSZ CTE CTH CUC " *
    "CWR CXM CY0 CY1 CY3 CY4 CYA CYD CYF CYG CYJ CYM CYQ CYR CYS CZ2 CZZ D11 D3P D4P DA2 " *
    "DAB DAH DAL DAR DAS DBB DBS DBU DBY DBZ DC2 DCL DCY DDE DFI DFO DGH DGL DGN DHA DHI " *
    "DHL DHN DHP DHV DI7 DIL DIR DIV DLE DLS DLY DM0 DMH DMK DMT DNE DNL DNP DNS DOA DOH " *
    "DON DPL DPN DPP DPQ DPR DSE DSG DSN DSP DTH DTR DTY DVA DYS ECC EFC EHP ESB ESC EXY " *
    "EYS F2F FAK FB5 FB6 FCL FGA FGL FGP FH7 FHL FHO FLA FLE FLT FME FOE FP9 FRD FT6 FTR " *
    "FTY FVA FZN GAU GCM GFT GGL GHG GHP GL3 GLH GLJ GLK GLM GLN GLQ GLU GLX GLY GLZ GMA " *
    "GND GPL GSC GSU GT9 GVL H14 H5M HAC HAR HBN HCS HFA HGL HHI HIA HIC HIP HIQ HIS HL2 " *
    "HLU HMR HPC HPE HPH HPQ HQA HRG HRP HS8 HS9 HSE HSL HSO HTI HTN HTR HV5 HVA HY3 HYP " *
    "HZP I2M I58 IAM IAR IAS IEL IGL IIL ILE ILG ILX IML IOY IPG IT1 IYR IYT IZO JJJ JJK " *
    "JJL K1R KCX KGC KNB KOR KPI KST KYN KYQ L2A LA2 LAA LAL LBY LCK LCX LCZ LDH LED LEF " *
    "LEH LEI LEM LEN LET LEU LEX LHC LLP LLY LME LMF LMQ LP6 LPD LPG LPL LPS LSO LTA LTR " *
    "LVG LVN LYF LYK LYM LYN LYR LYS LYX LYZ M0H M2L M2S M30 M3L MA MAA MAI MBQ MC1 MCG " *
    "MCL MCS MD3 MD6 MDF MDH MEA MED MEG MEN MEQ MET MEU MF3 MGG MGN MGY MHL MHO MHS MIS " *
    "MK8 ML3 MLE MLL MLY MLZ MME MMO MND MNL MNV MOD MP8 MPH MPJ MPQ MSA MSE MSL MSO MSP " *
    "MT2 MTY MVA N10 N2C N7P N80 N8P NA8 NAL NAM NB8 NBQ NC1 NCB NCY NDF NEM NEP NFA NHL " *
    "NIY NLE NLN NLO NLP NLQ NMC NMM NNH NPH NPI NSK NTR NTY NVA NYS NZH O12 OAR OAS OBF " *
    "OBS OCS OCY OHI OHS OIC OLE OLT OLZ OMT ONH ONL OPR ORN ORQ OSE OTB OTH OXX P1L P2Y " *
    "PAQ PAS PAT PAU PBB PBF PCA PCC PCE PCS PDL PEC PF5 PFF PFX PG1 PG9 PGL PGY PH6 PHA " *
    "PHD PHE PHI PHL PHM PIV PLE PM3 POM PPN PR3 PR9 PRO PRS PSA PSH PTA PTH PTM PTR PVH " *
    "PVL PYA PYL PYX QCS QMM QPA QPH R1A R4K RE0 RE3 RON RVX RZ4 S1H S2C S2D S2P SAC SAH " *
    "SAR SBL SCH SCS SCY SD2 SDP SE7 SEB SEC SEG SEL SEM SEN SEP SER SET SGB SHC SHP SHR " *
    "SIB SLR SLZ SMC SME SMF SNC SNN SOC SOY SRZ STY SUB SUN SVA SVV SVW SVX SVY SVZ SYS " *
    "T11 T66 TA4 TAV TBG TBM TCQ TCR TDD TFQ TH6 THC THO THR THZ TIH TMB TMD TNB TNR TOQ " *
    "TPH TPL TPO TPQ TQI TQQ TRF TRG TRN TRO TRP TRQ TRW TRX TRY TST TTQ TTS TXY TY1 TY2 " *
    "TY3 TY5 TYB TYI TYJ TYN TYO TYQ TYR TYS TYT TYW TYX TYY TZB TZO UMA UN1 UN2 UNK VAD " *
    "VAF VAL VB1 VDL VLL VLM VMS VOL WLU WPA WRP WVL X2W XCN XCP XDT XPL XPR XSN XX1 YCM " *
    "YOF YTH Z01 ZAL ZCL ZFB ZU0 ZZJ")

    # attributes for selection
    s = replace(s, "chainid" => "match_query(chainid, \" ")
    s = replace(s, "chainname" => "match_query(chainname, \" ")
    s = replace(s, "resid" => "match_query(resid, \" ")
    s = replace(s, "resname" => "match_query(resname, \" ")
    s = replace(s, "atomid" => "match_query(atomid, \" ")
    s = replace(s, "atomname" => "match_query(atomname, \" ")
    s = replace(s, r"(^|(\s+))name(\s+)" => "match_query(atomname, \" ")
    s = replace(s, r"(^|(\s+))chain(\s+)" => "match_query(chainname, \" ")

    # parentheses (only 1-depth) and logical operators (and, or , not)
    s = replace(s, r"\)(\s+)and" => "\" ) ) .& ")
    s = replace(s, r"([^\)])(\s+)and" => s"\1 \" ) .& ")
    s = replace(s, r"\)(\s+)or" => "\" ) ) .| ")
    s = replace(s, r"([^\)])(\s+)or" => s"\1 \" ) .| ")
    s = replace(s, r"(^|(\s+))not(\s+)" => " .! ")

    s = strip(s)
    if s[end] == ')'
        s = s[1:end-1] * " \" ) )"
    elseif s[end-length("all")+1:end] != "all" && s[end-length("backbone")+1:end] != "backbone"
        s = s * " \" )"
    end

    #resname = ta.resname
    #println(s)
    ex = Meta.parse(s)
    #Meta.dump(ex)
    replace_ex!(ex, ta)
    #Meta.dump(ex)
    #println(s)
    index = Meta.eval(ex)
    #println(index')
    findall(index)
end

function getindex(ta::TrjArray, s::AbstractString)
    index = select_atom(ta, s)
    ta[:, index]
end

# combinations
getindex(ta::TrjArray, ::Colon, s::AbstractString) = getindex(ta, s)
getindex(ta::TrjArray, rows, s::AbstractString) = ta[rows, :][:, s]

###### equality #################
==(x::TrjArray, y::TrjArray) = all(f -> getfield(x, f) == getfield(y, f), fieldnames(TrjArray))
isequal(x::TrjArray, y::TrjArray) = all(f -> isequal(getfield(x, f), getfield(y, f)), fieldnames(TrjArray))

###### copy ###############
copy(ta::TrjArray)::TrjArray =
    TrjArray(ta.xyz, ta.boxsize,
             ta.chainname, ta.chainid,
             ta.resname, ta.resid,
             ta.atomname, ta.atomid,
             ta.mass, ta.radius, ta.charge, ta.sasa,
             ta.list_bond, ta.list_angle, ta.list_dihedral, ta.list_improper, ta.list_cmap)

###### accessors to field values #################

###### vcat, hcat, merge #################
function vcat(ta_collection::TrjArray...)
    natom = ta_collection[1].natom
    isbox_empty = false

    T = Float64
    if isempty(ta_collection[1].boxsize)
        isbox_empty = true
    else
        T = eltype(ta_collection[1].boxsize)
    end

    for i = 2:length(ta_collection)
        if natom != ta_collection[i].natom
            throw(ArgumentError("number of atoms doesn't match"))
        end
    end
    xyz, boxsize = ta_collection[1].xyz, ta_collection[1].boxsize
    for i = 2:length(ta_collection)
        if !isempty(ta_collection[i].xyz)
            if isempty(xyz)
                xyz = ta_collection[i].xyz
            else
                xyz = [xyz; ta_collection[i].xyz]
            end
        end
        if !isempty(ta_collection[i].boxsize) & !isbox_empty
            if isempty(boxsize)
                boxsize = ta_collection[i].boxsize
            else
                boxsize = [boxsize; ta_collection[i].boxsize]
            end
        else
            @printf "Warning: boxsize information discarded in some trajectories\n"
            isbox_empty = true
            boxsize = Matrix{T}(undef, 0, 0)
        end
    end
    TrjArray(ta_collection[1], xyz=xyz, boxsize=boxsize)
end

###### end keyword #################
# endof(ta::TrjArray) = isempty(ta.x) ? nothing : size(ta.x, 1)
# eachindex(ta::TrjArray) = Base.OneTo(size(ta.x, 1))

###### iterator #################
Base.iterate(ta::TrjArray, state=1) = state > ta.nframe ? nothing : (ta[state], state + 1)

###### conversion #################
function convert(::Type{T}, ta::TrjArray) where {T<:AbstractArray}
    c = Matrix{eltype(T)}(undef, ta.nframe, 3*ta.natom)
    c .= ta.xyz
    #for iatom = 1:ta.natom
    #    c[1:ta.nframe, (3*(iatom-1)+1):(3*(iatom-1)+1)] .= eltype(T).(ta.x[1:ta.nframe, iatom:iatom])
    #    c[1:ta.nframe, (3*(iatom-1)+2):(3*(iatom-1)+2)] .= eltype(T).(ta.y[1:ta.nframe, iatom:iatom])
    #    c[1:ta.nframe, (3*(iatom-1)+3):(3*(iatom-1)+3)] .= eltype(T).(ta.z[1:ta.nframe, iatom:iatom])
    #end
    return c
end

###### broadcast #################

###### show #####################
function alignment_xyz(io::IO, X::AbstractVecOrMat,
        rows::AbstractVector, cols::AbstractVector,
        cols_if_complete::Integer, cols_otherwise::Integer, sep::Integer)
    a = Tuple{Int, Int}[]
    for j in cols
        l = 14
        r = 13
        push!(a, (l, r))
        if length(a) > 1 && sum(map(sum,a)) + sep*length(a) >= cols_if_complete
            pop!(a)
            break
        end
    end
    if 1 < length(a) < length(axes(X,2))
        while sum(map(sum,a)) + sep*length(a) >= cols_otherwise
            pop!(a)
        end
    end
    return a
end


function print_matrix_xyz(io::IO,
        X::AbstractVecOrMat, Y::AbstractVecOrMat, Z::AbstractVecOrMat,
        columnwidth::Vector,
        i::Integer, cols::AbstractVector, sep::AbstractString)
    for (k, j) = enumerate(cols)
        k > length(columnwidth) && break
        # if X[i,j] < 10^6 && X[i,j] > -10^5
        #     print_x = Printf.@sprintf " %8.2f" X[i,j]
        # else
        #     print_x = Printf.@sprintf " %8.2e" X[i,j]
        # end
        # TODO: treatment for very large or small values
        print_x = Printf.@sprintf " %8.2f" X[i,j]
        print_y = Printf.@sprintf " %8.2f" Y[i,j]
        print_z = Printf.@sprintf " %8.2f" Z[i,j]
        print(io, print_x, print_y, print_z, sep) # 27 characters + sep
        #if k < length(columnwidth); print(io, sep); end
    end
end

function string_column(x, y, j)
    if isempty(x)
        string_x = ""
    else
        string_x = string(x[j])
    end
    if isempty(y)
        string_y = ""
    else
        string_y = string(y[j])
    end
    s = " " * string_x * string_y
    if length(s) > 24
        s = " " * s[1:24] * "~" * " "
    else
        s = rpad(s, 27)
    end
    s
end

function print_column(io::IO,
                      x, y, columnwidth::Vector, cols::AbstractVector,
                      sep::AbstractString)
    for (k, j) = enumerate(cols)
        k > length(columnwidth) && break
        s = string_column(x, y, j)
        print(io, s, sep) # 27 characters + sep
    end
end


function print_matrix_vdots(io::IO, vdots::AbstractString,
        A::Vector, sep::AbstractString)
    for k = 1:length(A)
        w = A[k][1] + A[k][2]
        l = repeat(" ", max(0, A[k][1]-length(vdots)))
        r = repeat(" ", max(0, w-length(vdots)-length(l)))
        print(io, l, vdots, r)
        if k <= length(A); print(io, sep); end
    end
end


function show(io::IO, ta::TrjArray)
    # for jupyter notebook
    #ENV["COLUMNS"] = 150
    #ENV["LINES"] = 30

    #TODO: fix double empty lines when nframe=0
    pre = "|"  # pre-matrix string
    sep = " |" # separator between elements
    hdots = "  \u2026  "
    vdots = "\u22ee"
    ddots = "  \u22f1  "
    hmod = 5

    # summary line
    @printf(io, "%dx%d %s\n", ta.nframe, ta.natom, typeof(ta))

    X = ta.xyz[:, 1:3:end]
    Y = ta.xyz[:, 2:3:end]
    Z = ta.xyz[:, 3:3:end]

    sz = displaysize(io)
    screenheight, screenwidth = sz[1] - 4, sz[2]
    screenwidth -= length(pre)

    if !isempty(ta.chainid) || !isempty(ta.chainname)
        screenheight -= 1
    end

    if !isempty(ta.resid) || !isempty(ta.resname)
        screenheight -= 1
    end

    if !isempty(ta.atomid) || !isempty(ta.atomname)
        screenheight -= 1
    end

    @assert textwidth(hdots) == textwidth(ddots)

    #rowsA, colsA = UnitRange(axes(X,1)), UnitRange(axes(X,2))
    rowsA, colsA = 1:ta.nframe, 1:ta.natom
    m, n = length(rowsA), length(colsA)

    # To figure out alignments, only need to look at as many rows as could
    # fit down screen. If screen has at least as many rows as A, look at A.
    # If not, then we only need to look at the first and last chunks of A,
    # each half a screen height in size.
    halfheight = div(screenheight,2)
    if m > screenheight
        rowsA = [rowsA[(0:halfheight-1) .+ firstindex(rowsA)]; rowsA[(end-div(screenheight-1,2)+1):end]]
    end
    # Similarly for columns, only necessary to get alignments for as many
    # columns as could conceivably fit across the screen
    maxpossiblecols = div(screenwidth, 3+length(sep))
    if n > maxpossiblecols
        colsA = [colsA[(1:maxpossiblecols)]; colsA[(end-maxpossiblecols+1):end]]
    end
    columnwidth = alignment_xyz(io, X, rowsA, colsA, screenwidth, screenwidth, length(sep))

    # Nine-slicing is accomplished using print_matrix_row repeatedly
    if n <= length(columnwidth) # rows and cols fit so just print whole matrix in one piece
        if !isempty(ta.chainid) || !isempty(ta.chainname)
            print(io, pre)
            print_column(io, ta.chainid, ta.chainname, columnwidth, colsA, sep)
            println(io)
        end
        if !isempty(ta.resid) || !isempty(ta.resname)
            print(io, pre)
            print_column(io, ta.resid, ta.resname, columnwidth, colsA, sep)
            println(io)
        end
        if !isempty(ta.atomid) || !isempty(ta.atomname)
            print(io, pre)
            print_column(io, ta.atomid, ta.atomname, columnwidth, colsA, sep)
            println(io)
        end
    else # rows fit down screen but cols don't, so need horizontal ellipsis
        c = div(screenwidth-length(hdots)+1,2)+1  # what goes to right of ellipsis
        Rcolumnwidth = reverse(alignment_xyz(io, X, rowsA, reverse(colsA), c, c, length(sep))) # alignments for right
        c = screenwidth - sum(map(sum,Rcolumnwidth)) - (length(Rcolumnwidth)-1)*length(sep) - length(hdots)
        Lcolumwidth = alignment_xyz(io, X, rowsA, colsA, c, c, length(sep)) # alignments for left of ellipsis
        if !isempty(ta.chainid) || !isempty(ta.chainname)
            print(io, pre)
            #print_matrix_xyz(io, X, Y, Z, Lcolumwidth, i, colsA[1:length(Lcolumwidth)], sep)
            print_column(io, ta.chainid, ta.chainname, Lcolumwidth, colsA[1:length(Lcolumwidth)], sep)
            print(io, hdots)
            #print_matrix_xyz(io, X, Y, Z, Rcolumnwidth, i, (n - length(Rcolumnwidth)) .+ colsA, sep)
            print_column(io, ta.chainid, ta.chainname, Rcolumnwidth, (n - length(Rcolumnwidth)) .+ colsA, sep)
            println(io)
        end
        if !isempty(ta.resid) || !isempty(ta.resname)
            print(io, pre)
            #print_matrix_xyz(io, X, Y, Z, Lcolumwidth, i, colsA[1:length(Lcolumwidth)], sep)
            print_column(io, ta.resid, ta.resname, Lcolumwidth, colsA[1:length(Lcolumwidth)], sep)
            print(io, hdots)
            #print_matrix_xyz(io, X, Y, Z, Rcolumnwidth, i, (n - length(Rcolumnwidth)) .+ colsA, sep)
            print_column(io, ta.resid, ta.resname, Rcolumnwidth, (n - length(Rcolumnwidth)) .+ colsA, sep)
            println(io)
        end
        if !isempty(ta.atomid) || !isempty(ta.atomname)
            print(io, pre)
            #print_matrix_xyz(io, X, Y, Z, Lcolumwidth, i, colsA[1:length(Lcolumwidth)], sep)
            print_column(io, ta.atomid, ta.atomname, Lcolumwidth, colsA[1:length(Lcolumwidth)], sep)
            print(io, hdots)
            #print_matrix_xyz(io, X, Y, Z, Rcolumnwidth, i, (n - length(Rcolumnwidth)) .+ colsA, sep)
            print_column(io, ta.atomid, ta.atomname, Rcolumnwidth, (n - length(Rcolumnwidth)) .+ colsA, sep)
            println(io)
        end
    end

    if m <= screenheight # rows fit vertically on screen
        if n <= length(columnwidth) # rows and cols fit so just print whole matrix in one piece
            for i in rowsA
                print(io, pre)
                print_matrix_xyz(io, X, Y, Z, columnwidth, i, colsA, sep)
                if i != last(rowsA); println(io); end
            end
        else # rows fit down screen but cols don't, so need horizontal ellipsis
            c = div(screenwidth-length(hdots)+1,2)+1  # what goes to right of ellipsis
            Rcolumnwidth = reverse(alignment_xyz(io, X, rowsA, reverse(colsA), c, c, length(sep))) # alignments for right
            c = screenwidth - sum(map(sum,Rcolumnwidth)) - (length(Rcolumnwidth)-1)*length(sep) - length(hdots)
            Lcolumwidth = alignment_xyz(io, X, rowsA, colsA, c, c, length(sep)) # alignments for left of ellipsis
            for i in rowsA
                print(io, pre)
                print_matrix_xyz(io, X, Y, Z, Lcolumwidth, i, colsA[1:length(Lcolumwidth)], sep)
                print(io, (i - first(rowsA)) % hmod == 0 ? hdots : repeat(" ", length(hdots)))
                print_matrix_xyz(io, X, Y, Z, Rcolumnwidth, i, (n - length(Rcolumnwidth)) .+ colsA, sep)
                if i != last(rowsA); println(io); end
            end
        end
    else # rows don't fit so will need vertical ellipsis
        if n <= length(columnwidth) # rows don't fit, cols do, so only vertical ellipsis
            for i in rowsA
                print(io, pre)
                print_matrix_xyz(io, X, Y, Z, columnwidth, i, colsA, sep)
                if i != rowsA[end] || i == rowsA[halfheight]; println(io); end
                if i == rowsA[halfheight]
                    print(io, pre)
                    print_matrix_vdots(io, vdots, columnwidth, sep)
                    println(io)
                end
            end
        else # neither rows nor cols fit, so use all 3 kinds of dots
            c = div(screenwidth-length(hdots)+1,2)+1
            Rcolumnwidth = reverse(alignment_xyz(io, X, rowsA, reverse(colsA), c, c, length(sep)))
            c = screenwidth - sum(map(sum,Rcolumnwidth)) - (length(Rcolumnwidth)-1)*length(sep) - length(hdots)
            Lcolumwidth = alignment_xyz(io, X, rowsA, colsA, c, c, length(sep))
            for i in rowsA
                print(io, pre)
                print_matrix_xyz(io, X, Y, Z, Lcolumwidth, i, colsA[1:length(Lcolumwidth)], sep)
                print(io, (i - first(rowsA)) % hmod == 0 ? hdots : repeat(" ", length(hdots)))
                print_matrix_xyz(io, X, Y, Z, Rcolumnwidth, i, (n-length(Rcolumnwidth)).+colsA, sep)
                if i != rowsA[end] || i == rowsA[halfheight]; println(io); end
                if i == rowsA[halfheight]
                    print(io, pre)
                    print_matrix_vdots(io, vdots, Lcolumwidth, sep)
                    print(io, ddots)
                    print_matrix_vdots(io, vdots, Rcolumnwidth, sep)
                    println(io)
                end
            end
        end
    end
end
