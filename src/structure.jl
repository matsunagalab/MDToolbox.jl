function to3(index::Int64)
    index3 = collect((3*(index-1) + 1):(3*index))
    index3
end

"""
calcbond

calculate distances between pair atom from their Cartesian coordinates
"""
function calcbond(trj::Array{Float64}, pair::Array{Int64})

    ## initialization
    @show "called"
    nframe = size(trj, 1)
    npair = size(pair, 1)

    ## calculation
    bond = zeros(Float64, (nframe, npair))
    for ipair = 1:npair
        index1 = to3(pair[ipair, 1])
        index2 = to3(pair[ipair, 2])
        for iframe = 1:nframe
            d = norm(trj[iframe, index1] - trj[iframe, index2])
            bond[iframe, ipair] = d
        end
    end

    bond
end
