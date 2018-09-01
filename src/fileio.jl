"""
read xplor or charmm (namd) format dcd file
"""
function readdcd(filename::String; index=nothing)
    header_ischarmm_4dims = false
    x = []
    y = []
    z = []
    boxsize = []

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
                header_ischarmm_extrablock = true
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

        if typeof(index) <: AbstractVector{Bool}
            index2 = findall(index)
        elseif index == nothing
            index2 = collect(1:header_natom)
        else
            index2 = index
        end

        x = zeros(Float64, (nframe, length(index2)))
        y = zeros(Float64, (nframe, length(index2)))
        z = zeros(Float64, (nframe, length(index2)))
        if header_ischarmm_extrablock == true
            boxsize = zeros(Float64, (nframe, 3))
        end

        # read next frames
        dummy = Array{Float64, 1}(undef, 6)
        crd_x = Array{Float32, 1}(undef, header_natom)
        crd_y = Array{Float32, 1}(undef, header_natom)
        crd_z = Array{Float32, 1}(undef, header_natom)
        for iframe in 1:nframe
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
        end
    end

    #x, y, z, boxsize
    TrjArray(x=x, y=y, z=z, boxsize=boxsize)
end


"""
read netcdf file
"""
function readnetcdf(filename::String; index=nothing)
    finfo = ncinfo(filename)
    attributes = finfo.gatts
    #dimensions = finfo.dims
    nframe = Int64(finfo.dim["frame"].dimlen)
    natom = Int64(finfo.dim["atom"].dimlen)

    is_trj = haskey(finfo.vars, "coordinates") ? true : false
    is_box = haskey(finfo.vars, "cell_lengths") ? true : false
    is_vel = haskey(finfo.vars, "velocities") ? true : false
    is_temp = haskey(finfo.vars, "temp0") ? true : false

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
        boxsize = Vector{Float64}(undef, 0)
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

    TrjArray(x=x, y=y, z=z, boxsize=boxsize)
end


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
    for i = 1:4:(line_size-3)
        isPSF = line[i:i+3] == "PSF " ? true : isPSF
        isEXT = line[i:i+3] == "EXT " ? true : isEXT
        isCMAP = line[i:i+3] == "CMAP" ? true : isCMAP
        line[i:i+3] == "CHEQ" ? true : isCHEQ
    end
    if !isPSF
        print("Sorry, this seems not be a PSF file")
        return false
    end
    if isEXT
        #fmt_atom = "%10d| %8s| %8d| %8s| %8s| %4s| %14f|%14f|%8d %*f %*f" %octave
        fmt_atom = [1:10, 12:19, 21:28, 30:37, 39:46, 48:51, 53:66, 67:80, 81:88]
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
    psf_atom_type = Vector{String}(undef, 0)
    psf_atom_id = Vector{Int64}(undef, 0)
    psf_charge = Vector{Float64}(undef, 0)
    psf_mass = Vector{Float64}(undef, 0)

    i = 2
    while i <= length(lines)
        line = lines[i]; i += 1

        if occursin(r".*\d.*!\w.*", line) && !occursin("NGRP", line)
            line_splitted = split(line, "!")
            num = parse(Int64, line_splitted[1])
            if strip(line_splitted[end]) == "NATOM"
                for ii = i:(i+num-1)
                    line = lines[ii]
                    push!(psf_atom_id, parse_line(line, fmt_atom[1], Int64, 0))
                    push!(psf_segment_name, parse_line(line, fmt_atom[2], String, "None"))
                    push!(psf_residue_id, parse_line(line, fmt_atom[3], Int64, 0))
                    push!(psf_residue_name, parse_line(line, fmt_atom[4], String, "None"))
                    push!(psf_atom_name, parse_line(line, fmt_atom[5], String, "None"))
                    push!(psf_atom_type, parse_line(line, fmt_atom[6], String, "None")) # TODO: XPLOR or CHARMM format
                    push!(psf_charge, parse_line(line, fmt_atom[7], Float64, 0.0))
                    push!(psf_mass, parse_line(line, fmt_atom[8], Float64, 0.0))
                end
                i += num
            end
        end
    end
    TrjArray(chainname=psf_segment_name,
             resname=psf_residue_name, resid=psf_residue_id,
             atomname=psf_atom_name, atomid=psf_atom_id,
             mass=psf_mass, charge=psf_charge)
end
