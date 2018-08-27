import Base: convert, copy, length, show, getindex, start, next, done, isempty,
             endof, size, eachindex, ==, isequal, hash, vcat, hcat, merge, map

abstract type AbstractTrajectory end

struct TrjArray <: AbstractTrajectory

    #DONE: support for x, y, z = nothing, nothing, nothing
    #DONE: natom, nframe
    #TODO: revert to parametric type, fix slowness of constructor
    #TODO: simple and easy constructions for other functions which outputs a new trjarray (such as centerofmass)
    #TODO: should allow one-dimensional array (vector) for x, y, z constructor input
    x::Matrix{Float64}
    y::Matrix{Float64}
    z::Matrix{Float64}
    boxsize::Matrix{Float64}
    chainname::Vector{String}
    chainid::Vector{Int64}
    resname::Vector{String}
    resid::Vector{Int64}
    atomname::Vector{String}
    atomid::Vector{Int64}
    mass::Vector{Float64}
    charge::Vector{Float64}
    # bond::Matrix{Int64}
    # angle::Matrix{Int64}
    # dihedral::Matrix{Int64}
    meta::Any
    natom::Int64
    nframe::Int64

    function TrjArray(x, y, z, boxsize, chainname, chainid, resname, resid, atomname, atomid, mass, charge, meta)
        # nrow, ncol = (size(trj, 1), size(trj, 2))
        # natom = Int64(ncol/3)
        # ischecked && return new(trj, atomname, atomid, meta)
        # natom != length(chainname) && throw(DimensionMismatch("chainname must match width of trajectory"))
        # natom != length(chainid) && throw(DimensionMismatch("chainid must match width of trajectory"))

        # check size and define nframe
        nframe = 0
        if !isempty(x)
            nframe = size(x, 1)
            @assert nframe == size(y, 1) == size(z, 1)
            if !isempty(boxsize)
                @assert nframe == size(boxsize, 1)
                @assert 3 == size(boxsize, 2)
            end
        end
        nframe = Int64(nframe)

        # check size and define natom
        natom = 0
        mat_collection = [x, y, z]
        vec_collection = [chainname, chainid, resname, resid, atomname, atomid, mass, charge]
        if !isempty(x)
            natom = size(x, 2)
        else
            for vec in vec_collection
                if !isempty(vec)
                    natom = length(vec)
                    break
                end
            end
        end
        for mat in mat_collection
            if !isempty(mat)
                @assert natom == size(mat, 2)
            end
        end
        for vec in vec_collection
            if !isempty(vec)
                @assert natom == length(vec)
            end
        end
        natom = Int64(natom)

        # check type
        x2 = typeof(x) == Matrix{Float64} ? x : map(Float64, x)
        y2 = typeof(x) == Matrix{Float64} ? y : map(Float64, y)
        z2 = typeof(x) == Matrix{Float64} ? z : map(Float64, z)
        boxsize2 = typeof(boxsize) == Matrix{Float64} ? boxsize : map(Float64, boxsize)
        chainname2 = typeof(chainname) == Vector{String} ? chainname : map(strip, map(string, chainname))
        chainid2 = typeof(chainid) == Vector{Int64} ? chainid : map(Int64, chainid)
        resname2 = typeof(resname) == Vector{String} ? resname : map(strip, map(string, resname))
        resid2 = typeof(resid) == Vector{Int64} ? resid : map(Int64, resid)
        atomname2 = typeof(atomname) == Vector{String} ? atomname : map(strip, map(string, atomname))
        atomid2 = typeof(atomid) == Vector{Int64} ? atomid : map(Int64, atomid)
        mass2 = typeof(mass) == Vector{Float64} ? mass : map(Float64, mass)
        charge2 = typeof(charge) == Vector{Float64} ? charge : map(Float64, charge)

        return new(x2, y2, z2, boxsize2, chainname2, chainid2, resname2, resid2, atomname2, atomid2, mass2, charge2, meta, natom, nframe)
    end
end

###### outer constructors ########

function TrjArray(;x = Matrix{Float64}(undef, 0, 0), y = Matrix{Float64}(undef, 0, 0), z = Matrix{Float64}(undef, 0, 0),
                  boxsize = Matrix{Float64}(undef, 0, 0), 
                  chainname = Vector{String}(undef, 0), chainid = Vector{Int64}(undef, 0), 
                  resname = Vector{String}(undef, 0), resid = Vector{Int64}(undef, 0), 
                  atomname = Vector{String}(undef, 0), atomid = Vector{Int64}(undef, 0), 
                  mass = Vector{Float64}(undef, 0), charge = Vector{Float64}(undef, 0), 
                  meta = Any[])
    if typeof(x) <: AbstractVector
        x2 = reshape(x, length(x), 1)
    else
        x2 = x
    end
    if typeof(y) <: AbstractVector
        y2 = reshape(y, length(y), 1)
    else
        y2 = y
    end
    if typeof(z) <: AbstractVector
        z2 = reshape(z, length(z), 1)
    else
        z2 = z
    end
    TrjArray(x2, y2, z2, boxsize, chainname, chainid, resname, resid, atomname, atomid, mass, charge, meta)
end

TrjArray(x::Matrix{T}, y::Matrix{T}, z::Matrix{T}, boxsize::Matrix{T}, ta::TrjArray) where {T <: Real} =
             TrjArray(x = x, y = y, z = z, boxsize = boxsize, 
                      chainname = ta.chainname, chainid = ta.chainid,
                      resname = ta.resname, resid = ta.resid,
                      atomname = ta.atomname, atomid = ta.atomid,
                      mass = ta.mass, charge = ta.charge,
                      meta = ta.meta)

TrjArray(x::Matrix{T}, y::Matrix{T}, z::Matrix{T}, ta::TrjArray) where {T <: Real} =
             TrjArray(x = x, y = y, z = z, boxsize = ta.boxsize, 
                      chainname = ta.chainname, chainid = ta.chainid,
                      resname = ta.resname, resid = ta.resid,
                      atomname = ta.atomname, atomid = ta.atomid,
                      mass = ta.mass, charge = ta.charge,
                      meta = ta.meta)

###### getindex #################

# all element indexing
getindex(ta::TrjArray, ::Colon) = ta
getindex(ta::TrjArray, ::Colon, ::Colon) = ta

# single row
function getindex(ta::TrjArray, n::Int)
    if iszero(n)
        return TrjArray(Array{Float64, 2}(undef, 0, 0), Array{Float64, 2}(undef, 0, 0), Array{Float64, 2}(undef, 0, 0), ta)
    elseif !isempty(ta.boxsize)
        return TrjArray(ta.x[n:n, :], ta.y[n:n, :], ta.z[n:n, :], ta.boxsize[n:n, :], ta)
    else
        return TrjArray(ta.x[n:n, :], ta.y[n:n, :], ta.z[n:n, :], ta)
    end
end
getindex(ta::TrjArray, n::Int, ::Colon) = getindex(ta, n)

# range of rows
function getindex(ta::TrjArray, r::UnitRange{Int})
    if !isempty(ta.boxsize)
        TrjArray(ta.x[r, :], ta.y[r, :], ta.z[r, :], ta.boxsize[r, :], ta)
    else
        TrjArray(ta.x[r, :], ta.y[r, :], ta.z[r, :], ta)
    end
end
getindex(ta::TrjArray, r::UnitRange{Int}, ::Colon) = getindex(ta, r)

# array of rows (integer)
function getindex(ta::TrjArray, a::AbstractVector{S}) where {S <: Integer}
    if !isempty(ta.boxsize)
        TrjArray(ta.x[a, :], ta.y[a, :], ta.z[a, :], ta.boxsize[a, :], ta)
    else
        TrjArray(ta.x[a, :], ta.y[a, :], ta.z[a, :], ta)
    end
end
getindex(ta::TrjArray, a::AbstractVector{S}, ::Colon) where {S <: Integer} = getindex(ta, a)

# array of rows (bool)
#getindex(ta::TrjArray, a::AbstractVector{S}) where {S <: Bool} = TrjArray(ta.x[a, :], ta.y[a, :], ta.z[a, :], ta)
#getindex(ta::TrjArray, a::AbstractVector{S}, ::Colon) where {S <: Bool} = getindex(ta, a)
function getindex(ta::TrjArray, a::AbstractVector{Bool})
    if !isempty(ta.boxsize)
        TrjArray(ta.x[a, :], ta.y[a, :], ta.z[a, :], ta.boxsize[a, :], ta)
    else
        TrjArray(ta.x[a, :], ta.y[a, :], ta.z[a, :], ta)
    end
end
getindex(ta::TrjArray, a::AbstractVector{Bool}, ::Colon) = getindex(ta, a)

# single column
getindex(ta::TrjArray, ::Colon, n::Int) = TrjArray( 
             x = isempty(ta.x) ? [] : ta.x[:, n:n], 
             y = isempty(ta.y) ? [] : ta.y[:, n:n], 
             z = isempty(ta.z) ? [] : ta.z[:, n:n], 
             boxsize = ta.boxsize, 
             chainname = isempty(ta.chainname) ? [] : ta.chainname[n:n], 
             chainid = isempty(ta.chainid) ? [] : ta.chainid[n:n], 
             resname = isempty(ta.resname) ? [] : ta.resname[n:n], 
             resid = isempty(ta.resid) ? [] : ta.resid[n:n], 
             atomname = isempty(ta.atomname) ? [] : ta.atomname[n:n], 
             atomid = isempty(ta.atomid) ? [] : ta.atomid[n:n], 
             mass = isempty(ta.mass) ? [] : ta.mass[n:n], 
             charge = isempty(ta.charge) ? [] : ta.charge[n:n], 
             meta = ta.meta)

# range of columns
getindex(ta::TrjArray, ::Colon, r::UnitRange{Int}) = TrjArray( 
             x = isempty(ta.x) ? [] : ta.x[:, r], 
             y = isempty(ta.y) ? [] : ta.y[:, r], 
             z = isempty(ta.z) ? [] : ta.z[:, r], 
             boxsize = ta.boxsize, 
             chainname = isempty(ta.chainname) ? [] : ta.chainname[r], 
             chainid = isempty(ta.chainid) ? [] : ta.chainid[r], 
             resname = isempty(ta.resname) ? [] : ta.resname[r], 
             resid = isempty(ta.resid) ? [] : ta.resid[r], 
             atomname = isempty(ta.atomname) ? [] : ta.atomname[r], 
             atomid = isempty(ta.atomid) ? [] : ta.atomid[r], 
             mass = isempty(ta.mass) ? [] : ta.mass[r], 
             charge = isempty(ta.charge) ? [] : ta.charge[r], 
             meta = ta.meta)

# array of columns (integer)
getindex(ta::TrjArray, ::Colon, r::AbstractVector{S}) where {S <: Integer} = TrjArray( 
             x = isempty(ta.x) ? [] : ta.x[:, r], 
             y = isempty(ta.y) ? [] : ta.y[:, r], 
             z = isempty(ta.z) ? [] : ta.z[:, r], 
             boxsize = ta.boxsize, 
             chainname = isempty(ta.chainname) ? [] : ta.chainname[r], 
             chainid = isempty(ta.chainid) ? [] : ta.chainid[r], 
             resname = isempty(ta.resname) ? [] : ta.resname[r], 
             resid = isempty(ta.resid) ? [] : ta.resid[r], 
             atomname = isempty(ta.atomname) ? [] : ta.atomname[r], 
             atomid = isempty(ta.atomid) ? [] : ta.atomid[r], 
             mass = isempty(ta.mass) ? [] : ta.mass[r], 
             charge = isempty(ta.charge) ? [] : ta.charge[r], 
             meta = ta.meta)

# array of rows (bool)
getindex(ta::TrjArray, ::Colon, r::AbstractVector{Bool}) = TrjArray( 
             x = isempty(ta.x) ? [] : ta.x[:, r], 
             y = isempty(ta.y) ? [] : ta.y[:, r], 
             z = isempty(ta.z) ? [] : ta.z[:, r], 
             boxsize = ta.boxsize, 
             chainname = isempty(ta.chainname) ? [] : ta.chainname[r], 
             chainid = isempty(ta.chainid) ? [] : ta.chainid[r], 
             resname = isempty(ta.resname) ? [] : ta.resname[r], 
             resid = isempty(ta.resid) ? [] : ta.resid[r], 
             atomname = isempty(ta.atomname) ? [] : ta.atomname[r], 
             atomid = isempty(ta.atomid) ? [] : ta.atomid[r], 
             mass = isempty(ta.mass) ? [] : ta.mass[r], 
             charge = isempty(ta.charge) ? [] : ta.charge[r], 
             meta = ta.meta)

# combinations
getindex(ta::TrjArray, rows, cols) = ta[rows, :][:, cols]

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
    index = fill(false, natom)
    if isempty(some_array)
        return index
    elseif typeof(some_array[1]) == String
        query = split(query)
    else
        query = Meta.eval(Meta.parse(replace("[" * query * "]", r"(\d)\s" => s"\1, ")))
    end
    query = unique(unfold(query))
    for q in query
        index = index .| (some_array .== q) # conversion occurs because Vector{Bool} .| BitArray{1}
    end
    index
end

function replace_ex!(ex::Expr)
    for (i, arg) in enumerate(ex.args)
        if arg == :(match_query)
            ex.args[i] = :($match_query)
        end
        if isa(arg, Expr)
            replace_ex!(arg)
        end
    end
    ex
end

function select_atom(ta::TrjArray, s::AbstractString)
    s = strip(s)
    s = replace(s, "chainid" => "match_query($(ta.chainid), \" ")
    s = replace(s, "chainname" => "match_query($(ta.chainname), \" ")
    s = replace(s, "resid" => "match_query($(ta.resid), \" ")
    s = replace(s, "resname" => "match_query($(ta.resname), \" ")
    s = replace(s, "atomid" => "match_query($(ta.atomid), \" ")
    s = replace(s, "atomname" => "match_query($(ta.atomname), \" ")

    s = replace(s, r"\)(\s+)and" => "\" ) ) .& ")
    s = replace(s, r"([^\)])(\s+)and" => s"\1 \" ) .& ")
    s = replace(s, r"\)(\s+)or" => "\" ) ) .| ")
    s = replace(s, r"([^\)])(\s+)or" => s"\1 \" ) .| ")

    if s[end] == ')'
        s = s[1:end-1] * " \" ) )"
    elseif s[end-length("all")+1:end] != "all" && s[end-length("backbone")+1:end] != "backbone"
        s = s * " \" )"
    end

    ex = Meta.parse(s)
    replace_ex!(ex)
    #println(s)
    #println(index')
    index = Meta.eval(ex)
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

###### accessors to field values #################

###### vcat, hcat, merge #################
function vcat(ta_collection::TrjArray...)
    natom = ta_collection[1].natom
    for i = 2:length(ta_collection)
        if natom != ta_collection[i].natom
            throw(ArgumentError("number of atoms doesn't match"))
        end
    end
    x, y, z, boxsize = ta_collection[1].x, ta_collection[1].y, ta_collection[1].z, ta_collection[1].boxsize
    for i = 2:length(ta_collection)
        if !isempty(ta_collection[i].x)
            if isempty(x)
                x = ta_collection[i].x
            else
                x = [x; ta_collection[i].x]
            end
        end
        if !isempty(ta_collection[i].y)
            if isempty(y)
                y = ta_collection[i].y
            else
                y = [y; ta_collection[i].y]
            end
        end
        if !isempty(ta_collection[i].z)
            if isempty(z)
                z = ta_collection[i].z
            else
                z = [z; ta_collection[i].z]
            end
        end
        if !isempty(ta_collection[i].boxsize)
            if isempty(boxsize)
                boxsize = ta_collection[i].boxsize
            else
                boxsize = [boxsize; ta_collection[i].boxsize]
            end
        end
    end
    TrjArray(x, y, z, boxsize, ta_collection[1])
end

###### size #################

###### end keyword #################

endof(ta::TrjArray) = isempty(ta.x) ? nothing : size(ta.x, 1)
eachindex(ta::TrjArray) = Base.OneTo(size(ta.x, 1))

###### iterator #################

###### conversion #################

# conversion(Float64, ta:TrjArray)

###### mathematical operators, broadcast #################

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
    #TODO: fix double empty lines when nframe=0
    pre = "|"  # pre-matrix string
    sep = " |" # separator between elements
    hdots = "  \u2026  "
    vdots = "\u22ee"
    ddots = "  \u22f1  "
    hmod = 5

    # summary line
    @printf(io, "%dx%d %s\n", ta.nframe, ta.natom, typeof(ta))

    X = ta.x
    Y = ta.y
    Z = ta.z

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


