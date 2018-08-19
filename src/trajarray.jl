###### type definition ##########

import Base: convert, copy, length, show, getindex, start, next, done, isempty,
             endof, size, eachindex, ==, isequal, hash

abstract type AbstractTrajectory end

struct TrajArray{T, N, A <: AbstractArray{T, N}} <: AbstractTrajectory

    traj::A

    chainname::Vector{String}
    chainid::Vector{Int}

    resname::Vector{String}
    resid::Vector{Int}

    atomname::Vector{String}
    atomid::Vector{Int}

    # boxsize::Matrix{Float, 2}
    # mass::Vector{Float}
    # charge::Vector{Float}

    # bond::Matrix{Int, 2}
    # angle::Matrix{Int, 2}
    # dihedral::Matrix{Int, 2}

    meta::Any

    function TrajArray{T, N, A}(
            traj::A, 
            chainname::Vector{String}, 
            chainid::Vector{Int}, 
            resname::Vector{String}, 
            resid::Vector{Int}, 
            atomname::Vector{String}, 
            atomid::Vector{Int}, 
            meta::Any; 
            ischecked = false) where {T, N, A <: AbstractArray{T, N}}

        nrow, ncol = (size(traj, 1), size(traj, 2))
        natom = Int(ncol/3)

        ischecked && return new(traj, atomname, atomid, meta)

        natom != length(chainname) && throw(DimensionMismatch("chainname must match width of trajectory"))
        natom != length(chainid) && throw(DimensionMismatch("chainid must match width of trajectory"))
        natom != length(resname) && throw(DimensionMismatch("resname must match width of trajectory"))
        natom != length(resid) && throw(DimensionMismatch("resid must match width of trajectory"))
        natom != length(atomname) && throw(DimensionMismatch("atomname must match width of trajectory"))
        natom != length(atomid) && throw(DimensionMismatch("atomid must match width of trajectory"))

        return new(traj, chainname, chainid, resname, resid, atomname, atomid, meta)
    end
end

###### outer constructor ########

TrajArray(traj::AbstractArray{T, N},
          chainname::Vector{S}=fill("ChainName", Int(size(traj,2)/3)), chainid::Vector{I}=collect(1:Int(size(traj,2)/3)),
          resname::Vector{S}=fill("ResName", Int(size(traj,2)/3)), resid::Vector{I}=collect(1:Int(size(traj,2)/3)),
          atomname::Vector{S}=fill("AtomName", Int(size(traj,2)/3)), atomid::Vector{I}=collect(1:Int(size(traj,2)/3)),
          m::Any=nothing; args...) where {T, N, S <: AbstractString, I <: Int} =
              TrajArray{T, N, typeof(traj)}(traj,
                                            map(String, chainname), map(Int64, chainid),
                                            map(String, resname), map(Int64, resid),
                                            map(String, atomname), map(Int64, atomid),
                                            m; args...)

###### show #####################

DECIMALS = 2
MISSING = "NaN"

@inline _showval(v::Any) = repr(v)
@inline _showval(v::Number) = string(v)
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

function show(io::IO, ta::TrajArray{T}) where {T}

    # summary line
    nrow, ncol = (size(ta.traj, 1), size(ta.traj, 2))
    @printf(io, "%dx%d %s\n", nrow, ncol, typeof(ta))

    # calculate column withs
    drow, dcol = displaysize(io)
    res_row    = 9  # number of reserved rows: summary line, lable line ... etc
    half_row   = floor(Int, (drow - res_row) / 2)
    add_row    = (drow - res_row) % 2

    if nrow > (drow - res_row)
        tophalf = 1:(half_row + add_row)
        bothalf = (nrow - half_row + 1):nrow
        strs = _showval.(@view ta.traj[[tophalf; bothalf], :])
    else
        strs = _showval.(ta.traj)
    end
    strs = strs[:, 1:3:end] .* " " .* strs[:, 2:3:end] .* " " .* strs[:, 3:3:end]

    natom = Int(ncol/3)

    colnames1 = Vector{String}(undef, natom)
    for i in 1:natom
        colnames1[i] = string(ta.chainid[i]) * ta.chainname[i]
    end

    colnames2 = Vector{String}(undef, natom)
    for i in 1:natom
        colnames2[i] = string(ta.resid[i]) * ta.resname[i]
    end

    colnames3 = Vector{String}(undef, natom)
    for i in 1:natom
        colnames3[i] = string(ta.atomid[i]) * ta.atomname[i]
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

###### getindex #################

# range of rows
function getindex(ta::TrajArray, r::UnitRange{Int})
    TrajArray(ta.traj[r, :],
              ta.chainname, ta.chainid,
              ta.resname, ta.resid,
              ta.atomname, ta.atomid,
              ta.meta)
end

