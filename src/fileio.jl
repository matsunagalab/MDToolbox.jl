"""
read xplor or charmm (namd) format dcd file
"""
function readdcd(filename::String; index_atom=nothing)

    header_ischarmm_4dims = false
    x = nothing
    y = nothing
    z = nothing
    boxsize = nothing
    
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

        if typeof(index_atom) == Array{Bool, 1}
            index_atom2 = findall(index_atom)
        elseif index_atom == nothing
            index_atom2 = collect(1:header_natom)
        end

        # if ~exist('index', 'var') || isempty(index)
        #     index = 1:header_natom;
        # end

        x = zeros(Float64, (nframe, length(index_atom2)))
        y = zeros(Float64, (nframe, length(index_atom2)))
        z = zeros(Float64, (nframe, length(index_atom2)))
        boxsize = zeros(Float64, (nframe, 3))

        # read next frames
        dummy = Array{Float64, 1}(undef, 6)
        crd_x = Array{Float32, 1}(undef, header_natom)
        crd_y = Array{Float32, 1}(undef, header_natom)
        crd_z = Array{Float32, 1}(undef, header_natom)
        for iframe in 1:nframe
            #@show iframe
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

            if header_ischarmm_extrablock == true
                boxsize[iframe, :] = dummy[[1 3 6]];
            end
            x[iframe, :] = convert(Array{Float64, 1}, crd_x[index_atom2]);
            y[iframe, :] = convert(Array{Float64, 1}, crd_y[index_atom2]);
            z[iframe, :] = convert(Array{Float64, 1}, crd_z[index_atom2]);
        end

    end

    #x, y, z, boxsize
    TrjArray(x, y, z, boxsize=boxsize)

end


"""
write xplor or charmm (namd) format dcd file

function writedcd(filename::String, trj::Array{Float64}, boxsize::Array{Float64})

    # ## check existing file
    # if exist(filename, 'file')
    #     filename_old = sprintf('%s.old', filename);
    #     display(sprintf('existing file %s is moved to %s', filename, filename_old));
    #     movefile(filename, filename_old);
    # end

    ## initialization
    nframe, natom3 = size(trj);
    natom = Int64(natom3/3);

    if ~exist('header', 'var') || isempty(header)
        % default header in xplor format
        header.is_charmm = false;
        header.is_charmm_extrablock = false;
        header.is_charmm_4dims = false;
        header.blocksize1 = 84;
        header.hdr = 'CORD';
        header.nset = size(trj, 1);
        header.istrt = 0;
        header.nsavc = 1;
        header.nstep = 0;
        header.null4 = zeros(4, 1);
        header.nfreat = 0;
        header.delta = 1.0;
        header.null9 = zeros(9, 1);
        header.version = 0;
        header.blocksize2 = 164;
        header.ntitle = 2;
        title1 = sprintf('REMARKS FILENAME=%s CREATED BY MATLAB', filename);
        for i = (numel(title1)+1):80
            title1 = [title1 ' '];
        end
        title1 = title1(1:80);
        title2 = sprintf('REMARKS DATE: %s CREATED BY USER: %s', datestr(now, 'mm/dd/yy'), getenv('USER'));
        for i = (numel(title2)+1):80
            title2 = [title2 ' '];
        end
        title2 = title2(1:80);
        header.title = [title1; title2];
        header.blocksize3 = 4;
        header.natom = size(trj, 2) / 3;
    end

    if header.nset ~= size(trj, 1)
        header.nset = size(trj, 1);
    end

    if exist('box', 'var') && ~isempty(box)
        # charmm format
        header.is_charmm = true;
        header.is_charmm_extrablock = true;
        if header.version == 0
            header.version = 1; % is_charmm -> true
        end
        header.null9(1) = 1; % is_charmm_extrablock -> true
    else
        # xplor format
        header.is_charmm = false;
        header.is_charmm_extrablock = false;
        header.is_charmm_4dims = false;
        header.version = 0; % is_charmm -> false
        header.null9(1) = 0; % is_charmm_extrablock -> false
        header.null9(2) = 0; % is_charmm_4dims -> false
    end

    ## open file
    assert(ischar(filename), 'Please specify valid filename for the first argument')
    fid = fopen(filename, 'w');
    assert(fid > 0, 'Could not open file.');
    cleaner = onCleanup(@() fclose(fid));

    ## write block 1 (header)
    fwrite(fid, header.blocksize1, 'int32');             
    fwrite(fid, header.hdr, 'uchar');
    fwrite(fid, header.nset, 'int32');
    fwrite(fid, header.istrt, 'int32');
    fwrite(fid, header.nsavc, 'int32');
    fwrite(fid, header.nstep, 'int32');
    fwrite(fid, header.null4, 'int32');
    fwrite(fid, header.nfreat, 'int32');

    if header.is_charmm
        # charmm format
        fwrite(fid, header.delta, 'float32');
        fwrite(fid, header.null9, 'int32');
    else
        # xplor format
        fwrite(fid, header.delta, 'float64');
        fwrite(fid, header.null9(2:end), 'int32');
    end

    fwrite(fid, header.version, 'int32');
    fwrite(fid, header.blocksize1, 'int32');

    ## write block 2 (title)
    fwrite(fid, header.blocksize2, 'int32'); 
    fwrite(fid, header.ntitle, 'int32');
    fwrite(fid, header.title(1, :), 'uchar');
    fwrite(fid, header.title(2, :), 'uchar');
    fwrite(fid, header.blocksize2, 'int32');

    ## write block 3 (natom)
    fwrite(fid, header.blocksize3, 'int32');
    fwrite(fid, header.natom, 'int32');
    fwrite(fid, header.blocksize3, 'int32');

    ## write coordinates
    dummy = zeros(1, 6);
    for iframe = 1:nframe
        if header.is_charmm_extrablock
            fwrite(fid, 48, 'int32');
            dummy(1, [1 3 6]) = box(iframe, :);
            fwrite(fid, dummy, 'float64');
            fwrite(fid, 48, 'int32');
        end
  
        fwrite(fid, natom*4, 'int32');
        fwrite(fid, trj(iframe, 1:3:end), 'float32');
        fwrite(fid, natom*4, 'int32');

        fwrite(fid, natom*4, 'int32');
        fwrite(fid, trj(iframe, 2:3:end), 'float32');
        fwrite(fid, natom*4, 'int32');

        fwrite(fid, natom*4, 'int32');
        fwrite(fid, trj(iframe, 3:3:end), 'float32');
        fwrite(fid, natom*4, 'int32');
    end

end
"""
function readnetcdf(filename::String)
    
end


