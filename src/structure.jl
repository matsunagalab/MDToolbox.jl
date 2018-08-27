###### center of mass #################
"""
centerofmass

calculate center of mass of given trajectory
"""
function centerofmass(ta::TrjArray; isweight::Bool=true, index=nothing)::TrjArray
    nframe = ta.nframe
    natom = ta.natom
    if isweight
        @assert isempty(ta.mass) == false
        weight = ta.mass
    else
        weight = ones(Float64, natom)
    end
    if index == nothing
        index = 1:natom
    end
    com_x = Vector{Float64}(undef, nframe)
    com_y = Vector{Float64}(undef, nframe)
    com_z = Vector{Float64}(undef, nframe)
    weight_sum = sum(weight)
    @threads for iframe in 1:nframe
        com_x[iframe] = sum(ta.x[iframe, index] .* weight) / weight_sum
        com_y[iframe] = sum(ta.y[iframe, index] .* weight) / weight_sum
        com_z[iframe] = sum(ta.z[iframe, index] .* weight) / weight_sum
    end
    TrjArray(x=com_x, y=com_y, z=com_z)
end

###### superimpose #################
"""
centerofmass

calculate center of mass of given trajectory
"""
function superimpose(ref::TrjArray, ta::TrjArray; isweight::Bool=true, index=nothing)
    nframe = ta.nframe
    natom = ta.natom
    if isweight
        @assert isempty(ta.mass) == false
        weight = ta.mass
    else
        weight = ones(Float64, natom)
    end
    if index == nothing
        index = 1:natom
    end
    com = centerofmass(ref[1, :], isweight=isweight, index=index)
    ref_x = ref[1, :].x .- com.x
    ref_y = ref[1, :].y .- com.y
    ref_z = ref[1, :].z .- com.z
    com = centerofmass(ta, isweight=isweight, index=index)
    x = ta.x .- com.x
    y = ta.y .- com.y
    z = ta.z .- com.z
    TrjArray(x, y, z, ta)
end

###### distance, angle, dihedral #################
"""
calcdistance

calculate distances between two atoms or groups of atoms
"""
function distance(ta1::TrjArray, ta2::TrjArray)::Vector{Float64}
    # TODO: support for PBC
    @assert ta1.nframe == ta2.nframe
    nframe = ta1.nframe
    com1 = centerofmass(ta1)
    com2 = centerofmass(ta2)
    dist = Vector{Float64}(undef, nframe)
    @threads for iframe in 1:nframe
        d = 0.0
        d += (com1.x[iframe] - com2.x[iframe])^2
        d += (com1.y[iframe] - com2.y[iframe])^2
        d += (com1.z[iframe] - com2.z[iframe])^2
        d = sqrt(d)
        dist[iframe] = d
    end
    dist
end

###### superimpose #################

