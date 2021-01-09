"""
load xplor or charmm (namd) format dcd file
"""
function readdcd(filename::String; index=nothing, stride=1, isbox=true)
    #TODO: endian
    header_ischarmm_4dims = false
    header_ischarmm_extrablock = false
    x = y = z = Matrix{Float64}(undef, 0, 0)
    boxsize = Matrix{Float64}(undef, 0, 0)

    open(filename, "r") do io
        seekend(io)
        file_size = position(io) # get file size in byte
        seekstart(io)
        blocksize = read(io, Int32) # should be 84
        # header (4 chars) either "CORD" or "VELD"; position(io) == 4 bytes
        header_hdr = Array{Char, 1}(undef, 4)
        for i in 1:4
            header_hdr[i] = read(io, Char)
        end
        # the total # of frames (snapshots); position(io) == 8 bytes
        header_nset = read(io, Int32)
        # starting time-step; position(io) == 12 bytes
        header_istrt = read(io, Int32)
        # frequency to save trajectory; position(io) == 16 bytes
        header_nsavc = read(io, Int32)
        # the total # of simulation step; position(io) == 20 bytes
        header_nstep = read(io, Int32)
        # null4 (int*4); position(io) == 24 bytes
        header_null4 = Array{Int32, 1}(undef, 4)
        read!(io, header_null4)
        # of degrees of freedom; position(io) == 40 bytes
        header_nfreat = read(io, Int32)
        # step-size of simulation; position(io) == 44 bytes
        header_delta = read(io, Float32)
        # null9 (int*9); position(io) == 48 bytes
        header_null9 = Array{Int32, 1}(undef, 9)
        read!(io, header_null9)
        # version; position(io) == 84 bytes
        header_version = read(io, Int32)
        # charmm extension format
        if header_version > 0
            header_ischarmm = true
        end
        # delta is double precision in xplor format
        if header_ischarmm == false
            cof = position(io)
            # reread delta in double precision
            fseek(fid, 44)
            header_delta = read(io, Float64)
            seek(io, cof)
        end
        # check charmm extensions
        if header_ischarmm == true
            cof = position(io)
            # charmm extrablock extension
            seek(io, 48)
            n = read(io, Int32)
            if n == 1
                header_ischarmm_extrablock = true
            end
            # charmm 4dims extension
            seek(io, 52);
            n = read(io, Int32)
            if n == 1
                header_ischarmm_4dims = true
            end
            seek(io, cof)
        end
        # blocksize1; position(io) == 88 bytes
        blocksize1 = read(io, Int32)
        #### read block 2 (title)
        # blocksize2; position(io) == 92 bytes
        blocksize2 = read(io, Int32)
        # of title lines; position(io) == 96 bytes
        header_ntitle = read(io, Int32)
        # title
        header_title = []
        for i in 1:header_ntitle
            t = Array{Char, 1}(undef, 80)
            for i in 1:80
                t[i] = read(io, Char)
            end
            push!(header_title, String(t))
        end
        # blocksize2
        blocksize2 = read(io, Int32)
        #### read block 3 (natom)
        # blocksize3
        blocksize3 = read(io, Int32)
        # # of atoms
        header_natom = read(io, Int32)
        # blocksize3
        blocksize3 = read(io, Int32)
        ## read coordinates
        header_size = position(io)
        if header_ischarmm_extrablock == true
            extrablocksize = 4*2 + 8*6
        else
            extrablocksize = 0
        end
        coordblocksize = (4*2 + 4*header_natom)*3
        nframe = (file_size - header_size) / (extrablocksize + coordblocksize)
        nframe = Int64(floor(nframe))
        nframe_stride = length(1:stride:nframe)
        if typeof(index) <: AbstractVector{Bool}
            index2 = findall(index)
        elseif index == nothing
            index2 = collect(1:header_natom)
        else
            index2 = index
        end
        x = zeros(Float64, (nframe_stride, length(index2)))
        y = zeros(Float64, (nframe_stride, length(index2)))
        z = zeros(Float64, (nframe_stride, length(index2)))
        if header_ischarmm_extrablock == true
            boxsize = zeros(Float64, (nframe_stride, 3))
        end
        # read next frames
        dummy = Array{Float64, 1}(undef, 6)
        crd_x = Array{Float32, 1}(undef, header_natom)
        crd_y = Array{Float32, 1}(undef, header_natom)
        crd_z = Array{Float32, 1}(undef, header_natom)
        for iframe in 1:nframe_stride
            # read charmm extrablock (unitcell info)
            if header_ischarmm_extrablock == true
                blocksize_box = read(io, Int32)
                read!(io, dummy)
                blocksize_box = read(io, Int32)
            end
            # read x coordinates
            blocksize_x = read(io, Int32)
            read!(io, crd_x)
            blocksize_x = read(io, Int32)
            # read y coordinates
            blocksize_y = read(io, Int32)
            read!(io, crd_y)
            blocksize_y = read(io, Int32)
            # read z coordinates
            blocksize_z = read(io, Int32)
            read!(io, crd_z)
            blocksize_z = read(io, Int32)
            # skip charmm 4dims extension
            if header_ischarmm_4dims == true
                blocksize_4dims = read(io, Int32)
                skip(io, blocksize_4dims)
                blocksize_4dims = read(io, Int32)
            end
            # get data
            x[iframe, :] = crd_x[index2];
            y[iframe, :] = crd_y[index2];
            z[iframe, :] = crd_z[index2];
            if header_ischarmm_extrablock == true
                boxsize[iframe, :] = dummy[[1 3 6]];
            end

            # skip stride-1 frames
            if iframe < nframe_stride
                for ii = 1:(stride-1)
                    if header_ischarmm_extrablock == true
                        blocksize_box = read(io, Int32)
                        skip(io, blocksize_box)
                        blocksize_box = read(io, Int32)
                    end
                    # read x coordinates
                    blocksize_x = read(io, Int32)
                    skip(io, blocksize_x)
                    blocksize_x = read(io, Int32)
                    # read y coordinates
                    blocksize_y = read(io, Int32)
                    skip(io, blocksize_y)
                    blocksize_y = read(io, Int32)
                    # read z coordinates
                    blocksize_z = read(io, Int32)
                    skip(io, blocksize_z)
                    blocksize_z = read(io, Int32)
                    # skip charmm 4dims extension
                    if header_ischarmm_4dims == true
                        blocksize_4dims = read(io, Int32)
                        skip(io, blocksize_4dims)
                        blocksize_4dims = read(io, Int32)
                    end
                end
            end
        end
    end

    if isempty(boxsize) | !isbox
        TrjArray{Float64, Int64}(x=x, y=y, z=z)
    else
        TrjArray{Float64, Int64}(x=x, y=y, z=z, boxsize=boxsize)
    end
end


############################################################################
"""
read netcdf file
"""
function readnetcdf(filename::String; index=nothing)
    #finfo = ncinfo(filename)
    #attributes = finfo.gatts
    #dimensions = finfo.dims
    nc = NetCDF.open(filename, mode=NetCDF.NC_64BIT_OFFSET);
    nframe = Int64(nc.dim["frame"].dimlen)
    natom = Int64(nc.dim["atom"].dimlen)

    is_trj = haskey(nc.vars, "coordinates") ? true : false
    is_box = haskey(nc.vars, "cell_lengths") ? true : false
    is_vel = haskey(nc.vars, "velocities") ? true : false
    is_temp = haskey(nc.vars, "temp0") ? true : false

    NetCDF.close(nc)

    if typeof(index) <: AbstractVector{Bool}
        index2 = findall(index)
    elseif index == nothing
        index2 = collect(1:natom)
    else
        index2 = index
    end

    start_atom = minimum(index2)
    count_atom = maximum(index2) - start_atom + 1
    index3 = index2 .- start_atom .+ 1
    start_time = 1
    count_time = nframe

    if is_trj
        d = ncread(filename, "coordinates", start=[1, start_atom, start_time], count=[3, count_atom, count_time])
        d = map(Float64, d)
        x = d[1, index3, :]'
        y = d[2, index3, :]'
        z = d[3, index3, :]'
    else
        x = y = z = Matrix{Float64}(undef, 0, 0)
    end

    if is_box
        d = ncread(filename, "cell_lengths", start=[1, start_time], count=[3, count_time]);
        d = map(Float64, d)
        boxsize = d'
    else
        boxsize = Matrix{Float64}(undef, 0, 0)
    end

    if is_vel
        d = ncread(filename, "velocities", start=[1, start_atom, start_time], count=[3, count_atom, count_time])
        d = map(Float64, d)
        vx = d[1, index3, :]'
        vy = d[2, index3, :]'
        vz = d[3, index3, :]'
    else
        vx = vy = vz = Matrix{Float64}(undef, 0, 0)
    end

    if is_temp
        d = ncread(filename, "temp0", start=[start_time], count=[count_time]);
        temp = map(Float64, d)
    else
        temp = Vector{Float64}(undef, 0)
    end

    TrjArray{Float64, Int64}(x=x, y=y, z=z, boxsize=boxsize)
end

function writenetcdf2(filename::String, ta::TrjArray; velocity = nothing, force = nothing)
    scale_factor = 20.455
    natom = ta.natom
    nframe = ta.nframe

    isfile(filename) && rm(filename)

end

############################################################################
"""
write netcdf file
"""
function writenetcdf(filename::String, ta::TrjArray; velocity = nothing, force = nothing)
    scale_factor = 20.455
    natom = ta.natom
    nframe = ta.nframe

    # global attributes
    gatts = Dict(
        "ConventionVersion" => "1.0",
        "Conventions" => "AMBER",
        "application" => "Julia",
        "program" => "MDToolbox.jl",
        "programVersion" => "0.1",
        "title" => "Created by MDToolbox.jl"
    )

    # dimensions
    ncdim_frame = NcDim("frame", nframe, unlimited=true)
    ncdim_atom = NcDim("atom", natom)
    #ncdim_dimvalue = NcDim("dimvalue", 1)
    ncdim_spatial = NcDim("spatial", 3)
    ncdim_cell_spatial = NcDim("cell_spatial", 3)
    ncdim_cell_angular = NcDim("cell_angular", 3)
    ncdim_label = NcDim("label", 5)

    # variables
    ncvar_spatial = NcVar("spatial", ncdim_spatial, t=NetCDF.NC_CHAR)
    ncvar_time = NcVar("time", ncdim_frame, t=NetCDF.NC_FLOAT, atts=Dict("units" => "picosecond"))
    ncvar_coordinates = NcVar("coordinates", [ncdim_spatial, ncdim_atom, ncdim_frame], t=NetCDF.NC_FLOAT, atts=Dict("units" => "angstrom"))
    if !isempty(ta.boxsize)
        ncvar_cell_spatial = NcVar("cell_spatial", ncdim_cell_spatial, t=NetCDF.NC_CHAR)
        ncvar_cell_angular = NcVar("cell_angular", [ncdim_label, ncdim_cell_angular], t=NetCDF.NC_CHAR)
        ncvar_cell_lengths = NcVar("cell_lengths", [ncdim_cell_spatial, ncdim_frame], t=NetCDF.NC_FLOAT, atts=Dict("units" => "angstrom"))
        ncvar_cell_angles = NcVar("cell_angles", [ncdim_cell_angular, ncdim_frame], t=NetCDF.NC_FLOAT, atts=Dict("units" => "degree"))
    end
    if velocity != nothing
        ncvar_velocities = NcVar("velocities", [ncdim_spatial, ncdim_atom, ncdim_frame], t=NetCDF.NC_FLOAT, atts=Dict("units" => "angstrom/picosecond"))
    end
    if force != nothing
        ncvar_forces = NcVar("forces", [ncdim_spatial, ncdim_atom, ncdim_frame], t=NetCDF.NC_FLOAT, atts=Dict("units" => "amu*angstrom/picosecond^2"))
    end
    #ncvar_temp0 = NcVar("temp0", ncdim_frame, t=NetCDF.NC_FLOAT, atts=Dict("units" => "kelvin"))

    # create the NetCDF file
    varlist = [ncvar_spatial, ncvar_time, ncvar_coordinates]
    if !isempty(ta.boxsize)
        push!(varlist, ncvar_cell_spatial)
        push!(varlist, ncvar_cell_angular)
        push!(varlist, ncvar_cell_lengths)
        push!(varlist, ncvar_cell_angles)
    end
    if velocity != nothing
        push!(varlist, ncvar_velocities)
    end
    if force != nothing
        push!(varlist, ncvar_forces)
    end

    NetCDF.create(filename, varlist, gatts=gatts, mode=NetCDF.NC_64BIT_OFFSET) do nc
        # write data
        # you can check by ncread("tmp.nc", "spatial")
        #c = zeros(UInt8, 3); copyto!(c, 1, "xyz", 1); #@show c
        NetCDF.putvar(nc, "spatial", nc_string2char(["xyz"])[1:3] )
        NetCDF.putvar(nc, "time", zeros(Float32, nframe))
        trj = zeros(Float32, 3, natom, nframe)
        for iframe in 1:nframe
            for iatom in 1:natom
                trj[1, iatom, iframe] = ta.x[iframe, iatom]
                trj[2, iatom, iframe] = ta.y[iframe, iatom]
                trj[3, iatom, iframe] = ta.z[iframe, iatom]
            end
        end
        NetCDF.putvar(nc, "coordinates", trj)

        if !isempty(ta.boxsize)
            #c = zeros(UInt8, 3); copyto!(c, 1, "abc", 1); #@show c
            NetCDF.putvar(nc, "cell_spatial", nc_string2char(["abc"])[1:3])
            #c = zeros(UInt8, 5, 3); copyto!(c, 1, "alpha", 1); copyto!(c, 6, "beta ", 1); copyto!(c, 11, "gamma", 1); #@show c
            NetCDF.putvar(nc, "cell_angular", nc_string2char(["alpha", "beta ", "gamma"])[1:5, :])
            d = zeros(Float32, 3, nframe)
            d[:, :] .= ta.boxsize'
            NetCDF.putvar(nc, "cell_lengths", d)
            d[:, :] .= 90.0
            NetCDF.putvar(nc, "cell_angles", d)
        end

        if velocity != nothing
            for iframe in 1:nframe
                for iatom in 1:natom
                    trj[1, iatom, iframe] = velocity.x[iframe, iatom]/scale_factor
                    trj[2, iatom, iframe] = velocity.y[iframe, iatom]/scale_factor
                    trj[3, iatom, iframe] = velocity.z[iframe, iatom]/scale_factor
                end
            end
            NetCDF.putvar(nc, "velocities", trj)
        end
        if force != nothing
            for iframe in 1:nframe
                for iatom in 1:natom
                    trj[1, iatom, iframe] = force.x[iframe, iatom]; #TODO: scale_factor?
                    trj[2, iatom, iframe] = force.y[iframe, iatom]
                    trj[3, iatom, iframe] = force.z[iframe, iatom]
                end
            end
            NetCDF.putvar(nc, "forces", trj)
        end
    end

    return nothing
end


############################################################################
function parse_line(line::String, index, mytype::DataType, default_value)
    try
        if mytype == String
            return strip(string(line[index]))
        else
            return parse(mytype, line[index])
        end
    catch
        return default_value
    end
end


############################################################################
"""
read charmm or xplor type psf file
"""
function readpsf(filename::String)
    isPSF = false
    isEXT = false
    isCMAP = false
    isCHEQ = false

    lines = open(filename, "r" ) do fp
        readlines(fp)
    end

    # first line
    line = lines[1]
    line_size = length(line)
    if line_size >= 3
        if line[1:3] == "PSF"
            isPSF = true
        end
    end
    for i = 1:4:line_size
        i_max = minimum([i+3, line_size])
        isPSF = strip(line[i:i_max]) == "PSF" ? true : isPSF
        isEXT = strip(line[i:i_max]) == "EXT" ? true : isEXT
        isCMAP = strip(line[i:i_max]) == "CMAP" ? true : isCMAP
        isCHEQ = strip(line[i:i_max]) == "CHEQ" ? true : isCHEQ
    end
    if !isPSF
        print("Sorry, this seems not be a PSF file")
        return false
    end
    if isEXT
        #fmt_atom = "%10d| %8s| %8d| %8s| %8s| %4s| %14f|%14f|%8d %*f %*f" %octave
        #fmt_atom = [1:10, 12:19, 21:28, 30:37, 39:46, 48:51, 53:66, 67:80, 81:88]
        fmt_atom = [1:10, 12:19, 21:28, 30:37, 39:46, 48:51, 53:68, 69:80, 81:88]
        #fmt_list = "%10d"
    else
        #fmt_atom = "%8d %4s %4d %4s| %4s| %4s| %14f%14f%8d %*f %*f" #octave
        fmt_atom = [1:8, 10:13, 15:18, 20:23, 25:28, 30:33, 34:47, 48:61, 62:69]
        #fmt_list = "%8d"
        fmt_list = "%8d"
    end

    psf_title = Vector{String}(undef, 0)
    psf_segment_name = Vector{String}(undef, 0)
    psf_residue_name = Vector{String}(undef, 0)
    psf_residue_id = Vector{Int64}(undef, 0)
    psf_atom_name = Vector{String}(undef, 0)
    psf_atom_type_letter = Vector{String}(undef, 0)
    psf_atom_type_number = Vector{Int64}(undef, 0)
    psf_atom_id = Vector{Int64}(undef, 0)
    psf_charge = Vector{Float64}(undef, 0)
    psf_mass = Vector{Float64}(undef, 0)
    list_bond = Array{Int64}(undef, 0, 0)
    list_angle = Array{Int64}(undef, 0, 0)
    list_dihedral = Array{Int64}(undef, 0, 0)
    list_improper = Array{Int64}(undef, 0, 0)
    list_cmap = Array{Int64}(undef, 0, 0)

    iline = 2
    while iline <= length(lines)
        line = lines[iline]; iline += 1

        if occursin(r".*\d.*!\w.*", line) && !occursin("NGRP", line)
            line_splitted = split(line, "!")
            #num = parse(Int64, line_splitted[1])
            num = parse_line(String(line_splitted[1]), 1:length(line_splitted[1]), Int64, 0)
            key = strip(line_splitted[2])
            if occursin(r"^NATOM", key)
                for ii = iline:(iline+num-1)
                    line = lines[ii]
                    push!(psf_atom_id, parse_line(line, fmt_atom[1], Int64, 0))
                    push!(psf_segment_name, parse_line(line, fmt_atom[2], String, "None"))
                    push!(psf_residue_id, parse_line(line, fmt_atom[3], Int64, 0))
                    push!(psf_residue_name, parse_line(line, fmt_atom[4], String, "None"))
                    push!(psf_atom_name, parse_line(line, fmt_atom[5], String, "None"))
                    if isletter(strip(line[fmt_atom[6]])[1])
                        push!(psf_atom_type_letter, parse_line(line, fmt_atom[6], String, "None")) # CHARMM format
                    else
                        push!(psf_atom_type_number, parse_line(line, fmt_atom[6], Int64, 0)) # XPLOR format
                    end
                    push!(psf_charge, parse_line(line, fmt_atom[7], Float64, 0.0))
                    push!(psf_mass, parse_line(line, fmt_atom[8], Float64, 0.0))
                end
                iline += num

            elseif occursin(r"^NBOND", key)
                a = Vector{Int64}(undef, 0)
                while length(a) < (num*2)
                    line = lines[iline]; iline += 1
                    append!(a, [parse(Int64, y) for y in split(line)])
                end
                list_bond = collect(reshape(a, 2, :)')

            elseif occursin(r"^NTHETA", key)
                a = Vector{Int64}(undef, 0)
                while length(a) < (num*3)
                    line = lines[iline]; iline += 1
                    append!(a, [parse(Int64, y) for y in split(line)])
                end
                list_angle = collect(reshape(a, 3, :)')

            elseif occursin(r"^NPHI", key)
                a = Vector{Int64}(undef, 0)
                while length(a) < (num*4)
                    line = lines[iline]; iline += 1
                    append!(a, [parse(Int64, y) for y in split(line)])
                end
                list_dihedral = collect(reshape(a, 4, :)')

            elseif occursin(r"^NIMPHI", key)
                a = Vector{Int64}(undef, 0)
                while length(a) < (num*4)
                    line = lines[iline]; iline += 1
                    append!(a, [parse(Int64, y) for y in split(line)])
                end
                list_improper = collect(reshape(a, 4, :)')

            elseif occursin(r"^NCRTERM", key)
                a = Vector{Int64}(undef, 0)
                while length(a) < (num*8)
                    line = lines[iline]; iline += 1
                    append!(a, [parse(Int64, y) for y in split(line)])
                end
                list_cmap = collect(reshape(a, 8, :)')

            end
        end
    end

    TrjArray{Float64, Int64}(chainname=psf_segment_name,
             resname=psf_residue_name, resid=psf_residue_id,
             atomname=psf_atom_name, atomid=psf_atom_id,
             mass=psf_mass, charge=psf_charge,
             list_bond=list_bond, list_angle=list_angle,
             list_dihedral=list_dihedral, list_improper=list_improper,
             list_cmap=list_cmap)
end


############################################################################
"""
write psf file
"""
function writepsf(filename::String, ta::TrjArray)
    natom = ta.natom

    open(filename, "w") do io
        Printf.@printf(io, "%4s", "PSF ")
        if !isempty(ta.list_cmap)
            Printf.@printf(io, "%4s", "CMAP")
        end
        Printf.@printf(io, "\n\n")

        Printf.@printf(io, "%8d !NTITLE\n", 1)
        Printf.@printf(io, "CREATED by MDToolbox.jl\n")
        Printf.@printf(io, "\n")

        Printf.@printf(io, "%8d !NATOM\n", natom)
        for iatom = 1:natom
            Printf.@printf(io, "%8d", iatom)
            Printf.@printf(io, " %4s", isempty(ta.chainname) ? rpad("NONE", 4) : rpad(ta.chainname[iatom], 4))
            Printf.@printf(io, " %-4d", isempty(ta.resid) ? 0 : ta.resid[iatom])
            Printf.@printf(io, " %4s", isempty(ta.resname) ? rpad("NONE", 4) : rpad(ta.resname[iatom], 4))
            Printf.@printf(io, " %4s", isempty(ta.atomname) ? rpad("NONE", 4) : rpad(ta.atomname[iatom], 4))
            Printf.@printf(io, " %4s", isempty(ta.atomname) ? rpad("NONE", 4) : rpad(ta.atomname[iatom], 4))
            Printf.@printf(io, "%14.9f", isempty(ta.charge) ? 0.0 : ta.charge[iatom])
            Printf.@printf(io, "%14.7f", isempty(ta.mass) ? 0.0 : ta.mass[iatom])
            Printf.@printf(io, "%9d", 0)
            Printf.@printf(io, "\n")
        end

        if !isempty(ta.list_bond)
            Printf.@printf(io, "%8d !NBOND\n", size(ta.list_bond, 1))
            for i = 1:size(ta.list_bond, 1)
                for j = 1:size(ta.list_bond, 2)
                    Printf.@printf(io, "%8d", ta.list_bond[i, j])
                end
                if (i % 4) == 0
                    Printf.@printf(io, "\n")
                end
            end
            Printf.@printf(io, "\n")
        end

        if !isempty(ta.list_angle)
            Printf.@printf(io, "%8d !NTHETA\n", size(ta.list_angle, 1))
            for i = 1:size(ta.list_angle, 1)
                for j = 1:size(ta.list_angle, 2)
                    Printf.@printf(io, "%8d", ta.list_angle[i, j])
                end
                if (i % 3) == 0
                    Printf.@printf(io, "\n")
                end
            end
            Printf.@printf(io, "\n")
        end

        if !isempty(ta.list_dihedral)
            Printf.@printf(io, "%8d !NPHI\n", size(ta.list_dihedral, 1))
            for i = 1:size(ta.list_dihedral, 1)
                for j = 1:size(ta.list_dihedral, 2)
                    Printf.@printf(io, "%8d", ta.list_dihedral[i, j])
                end
                if (i % 2) == 0
                    Printf.@printf(io, "\n")
                end
            end
            Printf.@printf(io, "\n")
        end

        if !isempty(ta.list_improper)
            Printf.@printf(io, "%8d !NIMPHI\n", size(ta.list_improper, 1))
            for i = 1:size(ta.list_improper, 1)
                for j = 1:size(ta.list_improper, 2)
                    Printf.@printf(io, "%8d", ta.list_improper[i, j])
                end
                if (i % 2) == 0
                    Printf.@printf(io, "\n")
                end
            end
            Printf.@printf(io, "\n")
        end

        if !isempty(ta.list_cmap)
            Printf.@printf(io, "%8d !NCRTERM\n", size(ta.list_cmap, 1))
            for i = 1:size(ta.list_cmap, 1)
                for j = 1:size(ta.list_cmap, 2)
                    Printf.@printf(io, "%8d", ta.list_cmap[i, j])
                end
                Printf.@printf(io, "\n")
            end
            Printf.@printf(io, "\n")
        end
    end
end


############################################################################
"""
read protein data bank (PDB) file
"""
function readpdb(filename::String)
    lines = open(filename, "r" ) do fp
        readlines(fp)
    end

    model = []
    lines_clean = []
    for i = 1:length(lines)
        line = lines[i]
        if match(r"^ATOM", line) != nothing
            push!(lines_clean, line)
        end
        if match(r"^HETATM", line) != nothing
            push!(lines_clean, line)
        end
        if match(r"^ENDMDL", line) != nothing || (i == length(lines) && !isempty(lines_clean))
            push!(model, lines_clean)
            lines_clean = []
        end
    end

    natom = length(model[1])
    pdb_record = Vector{String}(undef, natom)
    pdb_serial = Vector{Int64}(undef, natom)
    pdb_name = Vector{String}(undef, natom)
    pdb_altloc = Vector{String}(undef, natom)
    pdb_resname = Vector{String}(undef, natom)
    pdb_chainid = Vector{String}(undef, natom)
    pdb_resseq = Vector{Int64}(undef, natom)
    pdb_x = zeros(Float64, length(model), natom)
    pdb_y = zeros(Float64, length(model), natom)
    pdb_z = zeros(Float64, length(model), natom)
    pdb_occupancy = Vector{Float64}(undef, natom)
    pdb_tempfactor = Vector{Float64}(undef, natom)
    pdb_element = Vector{String}(undef, natom)
    pdb_charge = Vector{Float64}(undef, natom)

    lines = model[1]
    for iatom = 1:natom
        line = lines[iatom]
        pdb_record[iatom] = parse_line(line, 1:6, String, "None")
        num = parse_line(line, 7:12, String, "0")
        if iatom == 1 || pdb_serial[iatom-1] < 99999
            pdb_serial[iatom] = parse(Int64, num)
        elseif match(r"\*", num) != nothing
            pdb_serial[iatom] = iatom == 1 ? 1 : pdb_serial[iatom-1]
        else
            pdb_serial[iatom] = parse(Int64, num)
        end
        pdb_name[iatom] = parse_line(line, 13:16, String, "None")
        pdb_altloc[iatom] = parse_line(line, 17:17, String, "None")
        pdb_resname[iatom] = parse_line(line, 18:21, String, "None")
        pdb_chainid[iatom] = parse_line(line, 22:22, String, "None")
        pdb_resseq[iatom] = parse_line(line, 23:28, Int64, 0)
        pdb_x[1, iatom] = parse_line(line, 31:38, Float64, 0.0)
        pdb_y[1, iatom] = parse_line(line, 39:46, Float64, 0.0)
        pdb_z[1, iatom] = parse_line(line, 47:54, Float64, 0.0)
        pdb_occupancy[iatom] = parse_line(line, 55:60, Float64, 0.0)
        pdb_tempfactor[iatom] = parse_line(line, 61:66, Float64, 0.0)
        pdb_element[iatom] = parse_line(line, 77:78, String, "None")
        pdb_charge[iatom] = parse_line(line, 79:80, Float64, 0.0)
    end

    for imodel = 2:length(model)
        lines = model[imodel]
        for iatom = 1:natom
            line = lines[iatom]
            pdb_x[imodel, iatom] = parse_line(line, 31:38, Float64, 0.0)
            pdb_y[imodel, iatom] = parse_line(line, 39:46, Float64, 0.0)
            pdb_z[imodel, iatom] = parse_line(line, 47:54, Float64, 0.0)
        end
    end

    TrjArray{Float64, Int64}(x=pdb_x, y=pdb_y, z=pdb_z,
             chainname=pdb_chainid,
             resid=pdb_resseq, resname=pdb_resname,
             atomid=pdb_serial, atomname=pdb_name)
end


############################################################################
"""
write pdb file
"""
function writepdb(io::IO, ta::TrjArray; format_type="vmd")
    natom = ta.natom
    nframe = ta.nframe

    if format_type == "namd"
        Printf.@printf io "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n"
    end

    for iframe = 1:nframe
        for iatom = 1:natom
            Printf.@printf(io, "%6s", "ATOM  ")
            if isempty(ta.atomid)
                Printf.@printf(io, "%5d", mod(iatom, 100000))
            else
                Printf.@printf(io, "%5d", mod(ta.atomid[iatom], 100000))
            end
            Printf.@printf(io, "%1s", " ")
            if isempty(ta.atomname)
                Printf.@printf(io, "%4s", " CA ")
            else
                Printf.@printf(io, "%4s", rpad(ta.atomname[iatom], 3))
            end
            Printf.@printf(io, "%1s", " ")
            #Printf.@printf(io, "%3s", ta.resname(iatom, :))
            #Printf.@printf(io, "%1s", " ")
            if isempty(ta.resname)
                Printf.@printf(io, "%4s", "ALA ")
            else
                Printf.@printf(io, "%4s", rpad(ta.resname[iatom], 4))
            end
            if isempty(ta.chainname)
                Printf.@printf(io, "%1s", " ")
            else
                Printf.@printf(io, "%1s", ta.chainname[iatom])
            end
            if isempty(ta.chainname)
                Printf.@printf(io, "%4d", mod(iatom, 10000))
            else
                Printf.@printf(io, "%4d", mod(ta.resid[iatom], 10000))
            end
            #Printf.@printf(io, "%1s", ta.icode(iatom, :))
            Printf.@printf(io, "%1s", " ")
            Printf.@printf(io, "%1s", " ")
            Printf.@printf(io, "%1s", " ")
            Printf.@printf(io, "%1s", " ")
            #Printf.@printf(io, "%8.3f", ta.xyz(iatom, 1))
            #Printf.@printf(io, "%8.3f", ta.xyz(iatom, 2))
            #Printf.@printf(io, "%8.3f", ta.xyz(iatom, 3))
            Printf.@printf(io, "%8.3f", ta.x[iframe, iatom])
            Printf.@printf(io, "%8.3f", ta.y[iframe, iatom])
            Printf.@printf(io, "%8.3f", ta.z[iframe, iatom])
            Printf.@printf(io, "%6.2f", 0.0)
            if isempty(ta.charge)
                Printf.@printf(io, "%6s", "      ")
            else
                Printf.@printf(io, "%6.2f", ta.charge[iatom])
            end
            Printf.@printf(io, "%1s", " ")
            Printf.@printf(io, "%1s", " ")
            Printf.@printf(io, "%1s", " ")
            Printf.@printf(io, "%1s", " ")
            Printf.@printf(io, "%1s", " ")
            Printf.@printf(io, "%1s", " ")
            Printf.@printf(io, "%1s", " ")
            Printf.@printf(io, "%1s", " ")
            Printf.@printf(io, "%1s", " ")
            Printf.@printf(io, "%1s", " ")
            Printf.@printf(io, "%2s", "  "); #element
            Printf.@printf(io, "%2s", "  "); #charge
            Printf.@printf(io, "\n")

            if !isempty(ta.chainname)
                if (iatom < ta.natom) && (ta.chainname[iatom] != ta.chainname[iatom+1])
                    Printf.@printf(io, "TER\n");
                end
            end
        end
        if nframe > 1
            Printf.@printf(io, "TER\n");
            Printf.@printf(io, "ENDMDL\n");
        end
        #Printf.@printf " %8.2f" X[i,j]
    end
end

"""
write pdb file
"""
function writepdb(filename::String, ta::TrjArray; format_type="vmd")
    open(filename, "w") do io
        writepdb(io, ta, format_type=format_type)
    end
end


############################################################################
"""
read amber coordinates file
"""
function readcrd(filename::String)
    lines = open(filename, "r" ) do fp
        readlines(fp)
    end

    title = lines[1]
    natom = parse_line(lines[2], 1:length(lines[2]), Int64, 0)
    #pdb_name[iatom] = parse_line(line, 13:16, String, "None")

    xyz = Vector{Float64}(undef, natom*3)

    iatom3 = 1
    icount = 3
    while true
        n = length(lines[icount])
        for i = 1:12:n
            xyz[iatom3] = parse_line(lines[icount], i:(i+11), Float64, 0.0)
            iatom3 += 1
        end
        icount += 1
        if iatom3 > (natom*3)
            break
        end
    end

    x = zeros(Float64, 1, natom)
    y = zeros(Float64, 1, natom)
    z = zeros(Float64, 1, natom)
    x[:] .= xyz[1:3:end]
    y[:] .= xyz[2:3:end]
    z[:] .= xyz[3:3:end]

    TrjArray{Float64, Int64}(x=x, y=y, z=z)
end
