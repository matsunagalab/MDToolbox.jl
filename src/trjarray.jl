import Base: convert, copy, length, show, getindex, start, next, done, isempty,
             endof, size, eachindex, ==, isequal, hash

abstract type AbstractTrajectory end

struct TrjArray <: AbstractTrajectory

    x::Union{Matrix{Float64}, Nothing}
    y::Union{Matrix{Float64}, Nothing}
    z::Union{Matrix{Float64}, Nothing}

    chainname::Union{Vector{String}, Nothing}
    chainid::Union{Vector{Int64}, Nothing}

    resname::Union{Vector{String}, Nothing}
    resid::Union{Vector{Int64}, Nothing}

    atomname::Union{Vector{String}, Nothing}
    atomid::Union{Vector{Int64}, Nothing}

    boxsize::Union{Matrix{Float64}, Nothing}
    # mass::Union{Vector{Float64}, Nothing}
    # charge::Union{Vector{Float64}, Nothing}

    # bond::Union{Matrix{Int64}, Nothing}
    # angle::Union{Matrix{Int64}, Nothing}
    # dihedral::Union{Matrix{Int64}, Nothing}

    meta::Any

    function TrjArray(
            x::Union{Matrix{Float64}, Nothing}, 
            y::Union{Matrix{Float64}, Nothing}, 
            z::Union{Matrix{Float64}, Nothing}, 
            chainname::Union{Vector{String}, Nothing}, 
            chainid::Union{Vector{Int64}, Nothing}, 
            resname::Union{Vector{String}, Nothing}, 
            resid::Union{Vector{Int64}, Nothing}, 
            atomname::Union{Vector{String}, Nothing}, 
            atomid::Union{Vector{Int64}, Nothing}, 
            boxsize::Union{Matrix{Float64}, Nothing}, 
            meta::Any) 

        # nrow, ncol = (size(trj, 1), size(trj, 2))
        # natom = Int64(ncol/3)
        #@show natom

        # ischecked && return new(trj, atomname, atomid, meta)

        # natom != length(chainname) && throw(DimensionMismatch("chainname must match width of trajectory"))
        # natom != length(chainid) && throw(DimensionMismatch("chainid must match width of trajectory"))

        if x != nothing 
            x2 = map(Float64, x)
        else
            x2 = nothing
        end

        if y != nothing 
            y2 = map(Float64, y)
        else
            y2 = nothing
        end

        if z != nothing 
            z2 = map(Float64, z)
        else
            z2 = nothing
        end

        if chainname != nothing 
            chainname2 = map(string, chainname)
        else
            chainname2 = nothing
        end

        if chainid != nothing 
            chainid2 = map(Int64, chainid)
        else
            chainid2 = nothing
        end

        if resname != nothing 
            resname2 = map(string, resname)
        else
            resname2 = nothing
        end

        if resid != nothing 
            resid2 = map(Int64, resid)
        else
            resid2 = nothing
        end

        if atomname != nothing 
            atomname2 = map(string, atomname)
        else
            atomname2 = nothing
        end

        if atomid != nothing 
            atomid2 = map(Int64, atomid)
        else
            atomid2 = nothing
        end
        
        if boxsize != nothing 
            boxsize2 = map(Float64, boxsize)
        else
            boxsize2 = nothing
        end

        return new(x2, y2, z2, chainname2, chainid2, resname2, resid2, atomname2, atomid2, boxsize2, meta)
    end
end

###### outer constructor ########

TrjArray(x::Matrix{T},
         y::Matrix{T},
         z::Matrix{T}; 
         chainname = nothing, 
         chainid = nothing, 
         resname = nothing, 
         resid = nothing, 
         atomname = nothing, 
         atomid = nothing, 
         boxsize = nothing, 
         meta = nothing) where {T <: Real} =
             TrjArray(x, y, z, chainname, chainid, resname, resid, atomname, atomid, boxsize, meta)

TrjArray(x::Matrix{T}, y::Matrix{T}, z::Matrix{T}, ta::TrjArray) where {T <: Real} =
             TrjArray(x, y, z, 
                      ta.chainname,
                      ta.chainid,
                      ta.resname,
                      ta.resid,
                      ta.atomname,
                      ta.atomid,
                      ta.boxsize,
                      ta.meta)

###### getindex #################

# all element indexing
getindex(ta::TrjArray, ::Colon) = ta
getindex(ta::TrjArray, ::Colon, ::Colon) = ta

# single row
getindex(ta::TrjArray, n::Int) = TrjArray(ta.x[n:n, :], ta.y[n:n, :], ta.z[n:n, :], ta)
getindex(ta::TrjArray, n::Int, ::Colon) = getindex(ta, n)

# range of rows
getindex(ta::TrjArray, r::UnitRange{Int}) = TrjArray(ta.x[r, :], ta.y[r, :], ta.z[r, :], ta)
getindex(ta::TrjArray, r::UnitRange{Int}, ::Colon) = getindex(ta, r)

# array of rows
getindex(ta::TrjArray, a::AbstractVector{S}) where {S <: Integer} = TrjArray(ta.x[a, :], ta.y[a, :], ta.z[a, :], ta)
getindex(ta::TrjArray, a::AbstractVector{S}, ::Colon) where {S <: Integer} = getindex(ta, a)

# array of rows
getindex(ta::TrjArray, a::AbstractVector{S}) where {S <: Bool} = TrjArray(ta.x[a, :], ta.y[a, :], ta.z[a, :], ta)
getindex(ta::TrjArray, a::AbstractVector{S}, ::Colon) where {S <: Bool} = getindex(ta, a)

# indexing to xyz indexing
function to3(index::AbstractVector{S}) where {S <: Integer}
    index3 = zeros(eltype(index), length(index)*3)
    index3[1:3:end] = 3 .* (index .- 1) .+ 1
    index3[2:3:end] = 3 .* (index .- 1) .+ 2
    index3[3:3:end] = 3 .* (index .- 1) .+ 3
    index3
end

function to3(index::AbstractVector{S}) where {S <: Bool}
    index3 = zeros(Bool, length(index)*3)
    index3[1:3:end] = index
    index3[2:3:end] = index
    index3[3:3:end] = index
    index3
end

function to3(index::Integer)
    index3 = zeros(eltype(index), 3)
    index3[1] = 3 * (index - 1) + 1
    index3[2] = 3 * (index - 1) + 2
    index3[3] = 3 * (index - 1) + 3
    index3
end

function to3(index::Bool)
    index3 = zeros(Bool, 3)
    index3[1] = index
    index3[2] = index
    index3[3] = index
    index3
end

# 3 columns
# integer
# bool

# multiple columns
# array
# unitrange

# end keyword
endof(ta::TrjArray) = ta.x == nothing ? nothing : size(ta.x, 1)
eachindex(ta::TrjArray) = Base.OneTo(size(ta.x, 1))

###### show #####################

const undef_ref_alignment = (3,3)

"""
`alignment(X, rows, cols, cols_if_complete, cols_otherwise, sep)` returns the
alignment for specified parts of array `X`, returning the (left,right) info.
It will look in X's `rows`, `cols` (both lists of indices)
and figure out what's needed to be fully aligned, for example looking all
the way down a column and finding out the maximum size of each element.
Parameter `sep::Integer` is number of spaces to put between elements.
`cols_if_complete` and `cols_otherwise` indicate screen width to use.
Alignment is reported as a vector of (left,right) tuples, one for each
column going across the screen.
"""
function alignment(io::IO, X::AbstractVecOrMat,
        rows::AbstractVector, cols::AbstractVector,
        cols_if_complete::Integer, cols_otherwise::Integer, sep::Integer)
    a = Tuple{Int, Int}[]
    for j in cols # need to go down each column one at a time
        l = r = 0
        for i in rows # plumb down and see what largest element sizes are
            if isassigned(X,i,j)
                aij = Base.alignment(io, X[i,j])
            else
                aij = undef_ref_alignment
            end
            l = max(l, aij[1]) # left characters
            r = max(r, aij[2]) # right characters
        end
        push!(a, (l, r)) # one tuple per column of X, pruned to screen width
        if length(a) > 1 && sum(map(sum,a)) + sep*length(a) >= cols_if_complete
            pop!(a) # remove this latest tuple if we're already beyond screen width
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


"""
`print_matrix_row(io, X, A, i, cols, sep)` produces the aligned output for
a single matrix row X[i, cols] where the desired list of columns is given.
The corresponding alignment A is used, and the separation between elements
is specified as string sep.
`print_matrix_row` will also respect compact output for elements.
"""
function print_matrix_row(io::IO,
        X::AbstractVecOrMat, A::Vector,
        i::Integer, cols::AbstractVector, sep::AbstractString)
    for (k, j) = enumerate(cols)
        k > length(A) && break
        if isassigned(X,Int(i),Int(j)) # isassigned accepts only `Int` indices
            x = X[i,j]
            a = Base.alignment(io, x)
            sx = sprint(show, x, context=:compact => true, sizehint=0)
        else
            a = undef_ref_alignment
            sx = undef_ref_str
        end
        l = repeat(" ", A[k][1]-a[1]) # pad on left and right as needed
        r = repeat(" ", A[k][2]-a[2])
        prettysx = Base.replace_in_print_matrix(X, i, j, sx)
        print(io, l, prettysx, r)
        if k < length(A); print(io, sep); end
    end
end

function print_matrix_xyz(io::IO,
        X::AbstractVecOrMat, Y::AbstractVecOrMat, Z::AbstractVecOrMat,
        columnwidth::Vector,
        i::Integer, cols::AbstractVector, sep::AbstractString)
    icount = 0
    for (k, j) = enumerate(cols)
        k > length(columnwidth) && break
        # if X[i,j] < 10^6 && X[i,j] > -10^5
        #     print_x = Printf.@sprintf " %8.2f" X[i,j]
        # else
        #     print_x = Printf.@sprintf " %8.2e" X[i,j]
        # end
        print_x = Printf.@sprintf " %8.2f" X[i,j]
        print_y = Printf.@sprintf " %8.2f" Y[i,j]
        print_z = Printf.@sprintf " %8.2f" Z[i,j]
        icount += 1
        print(io, print_x, print_y, print_z, sep)
        #if k < length(columnwidth); print(io, sep); end
    end
end


"""
`print_matrix_vdots` is used to show a series of vertical ellipsis instead
of a bunch of rows for long matrices. Not only is the string vdots shown
but it also repeated every M elements if desired.
"""
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
    pre = "|"  # pre-matrix string
    sep = " |" # separator between elements
    post = ""  # post-matrix string
    hdots = "  \u2026  "
    vdots = "\u22ee"
    ddots = "  \u22f1  "
    hmod = 5

    # summary line
    nrow, ncol = (size(ta.x, 1), size(ta.x, 2))
    @printf(io, "%dx%d %s\n", nrow, ncol, typeof(ta))

    X = ta.x
    Y = ta.y
    Z = ta.z

    sz = displaysize(io)
    screenheight, screenwidth = sz[1] - 4, sz[2]
    screenwidth -= length(pre)

    @assert textwidth(hdots) == textwidth(ddots)
    
    rowsA, colsA = UnitRange(axes(X,1)), UnitRange(axes(X,2))
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
                #print_matrix_row(io, X, A, i, colsA, sep)
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
            #Ralign = reverse(alignment(io, X, rowsA, reverse(colsA), c, c, length(sep)))
            Rcolumnwidth = reverse(alignment_xyz(io, X, rowsA, reverse(colsA), c, c, length(sep)))
            #c = screenwidth - sum(map(sum,Ralign)) - (length(Ralign)-1)*length(sep) - length(hdots)
            c = screenwidth - sum(map(sum,Rcolumnwidth)) - (length(Rcolumnwidth)-1)*length(sep) - length(hdots)
            #Lalign = alignment(io, X, rowsA, colsA, c, c, length(sep))
            Lcolumwidth = alignment_xyz(io, X, rowsA, colsA, c, c, length(sep))
            #r = mod((length(Rcolumnwidth)-n+1),vmod) # where to put dots on right half
            for i in rowsA
                print(io, pre)
                #print_matrix_row(io, X,Lalign,i,colsA[1:length(Lalign)],sep)
                print_matrix_xyz(io, X, Y, Z, Lcolumwidth, i, colsA[1:length(Lcolumwidth)], sep)
                print(io, (i - first(rowsA)) % hmod == 0 ? hdots : repeat(" ", length(hdots)))
                #print_matrix_row(io, X,Ralign,i,(n-length(Ralign)).+colsA,sep)
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
        if isempty(rowsA)
            print(io, vdots)
            length(colsA) > 1 && print(io, "    ", ddots)
            print(io, post)
        end
    end
end


MISSING = "NaN"

@inline _showval(v::Any) = repr(v)
@inline _showval(v::Real) = string(v)
@inline _showval(v::AbstractFloat) =
    #ifelse(isnan(v), MISSING, string(round(v, digits=DECIMALS)))
    ifelse(isnan(v), MISSING, @sprintf "%8.3g" v)

"""
calculate the paging

```
> using MarketData
> AAPL  # this function will return `UnitRange{Int64}[1:9, 10:12]`
```
"""
@inline function _showpages(dcol::Int, colwidth::Array{Int})
    ret = UnitRange{Int}[]
    c = dcol - 4
    last_i = 1
    for i in eachindex(colwidth)
        w = colwidth[i] + 3
        if c - w < 0
            push!(ret, last_i:i-1)
            # next page
            c = dcol - 4 - w
            last_i = i
        elseif i == length(colwidth)
            push!(ret, last_i:i)
        else
            c -= w
        end
    end
    ret
end


function show_old(io::IO, ta::TrjArray)

    # summary line
    nrow, ncol = (size(ta.x, 1), size(ta.x, 2))
    @printf(io, "%dx%d %s\n", nrow, ncol, typeof(ta))

    # calculate column withs
    drow, dcol = displaysize(io)
    res_row    = 9  # number of reserved rows: summary line, lable line ... etc
    half_row   = floor(Int, (drow - res_row) / 2)
    add_row    = (drow - res_row) % 2

    if nrow > (drow - res_row)
        tophalf = 1:(half_row + add_row)
        bothalf = (nrow - half_row + 1):nrow
        strs = _showval.(@view ta.x[[tophalf; bothalf], :])
    else
        strs = _showval.(ta.x)
    end
    strs = strs[:, 1:3:end] .* " " .* strs[:, 2:3:end] .* " " .* strs[:, 3:3:end]

    natom = Int(ncol/3)

    colnames1 = Vector{String}(undef, natom)
    for i in 1:natom
        colnames1[i] = string(ta.chainid[i]) * " " * ta.chainname[i]
    end

    colnames2 = Vector{String}(undef, natom)
    for i in 1:natom
        colnames2[i] = string(ta.resid[i]) * " " * ta.resname[i]
    end

    colnames3 = Vector{String}(undef, natom)
    for i in 1:natom
        colnames3[i] = string(ta.atomid[i]) * " " * ta.atomname[i]
    end

    # NOTE: reshaping is a workaround in julia 0.6
    #       in 0.7, it can be:
    #         [textwidth.(colnames)'; textwidth.(strs); fill(5, ncol)']
    colwidth = maximum([
        reshape(textwidth.(colnames1), 1, :);
        reshape(textwidth.(colnames2), 1, :);
        reshape(textwidth.(colnames3), 1, :);
        textwidth.(strs);
        reshape(fill(5, natom), 1, :)], dims=1)

    # paging
    pages = _showpages(dcol, colwidth)

    for p in pages
        # row label line
        ## e.g. | Open  | High  | Low   | Close  |
        #print(io, "│", " "^(spacetime + 2))
        for (name, w) in zip(colnames1[p], colwidth[p])
            print(io, "│ ", rpad(name, w + 1))
        end
        println(io, "│")

        for (name, w) in zip(colnames2[p], colwidth[p])
            print(io, "│ ", rpad(name, w + 1))
        end
        println(io, "│")

        for (name, w) in zip(colnames3[p], colwidth[p])
            print(io, "│ ", rpad(name, w + 1))
        end
        println(io, "│")
        ## e.g. ├───────┼───────┼───────┼────────┤
        #print(io, "├", "─"^(spacetime + 2))
        for w in colwidth[p]
            print(io, "┼", "─"^(w + 2))
        end
        print(io, "┤")

        # timestamp and values line
        if nrow > (drow - res_row)
            for i in tophalf
                println(io)
                #print(io, "│ ", ts[i], " ")
                for j in p
                    print(io, "│ ", rpad(strs[i, j], colwidth[j] + 1))
                end
                print(io, "│")
            end

            print(io, "\n   \u22EE")

            for i in (length(bothalf) - 1):-1:0
                i = size(strs, 1) - i
                println(io)
                #print(io, "│ ", ts[i], " ")
                for j in p
                    print(io, "│ ", rpad(strs[i, j], colwidth[j] + 1))
                end
                print(io, "│")
            end

        else
            for i in 1:nrow
                println(io)
                #print(io, "│ ", ts[i], " ")
                for j in p
                    print(io, "│ ", rpad(strs[i, j], colwidth[j] + 1))
                end
                print(io, "│")
            end
        end

        if length(pages) > 1 && p != pages[end]
            print(io, "\n\n")
        end
    end  # for p in pages
end

