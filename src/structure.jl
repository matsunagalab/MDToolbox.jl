###### center of mass #################
"""
centerofmass

calculate center of mass of given trajectory
"""
function centerofmass(ta::TrjArray; isweight::Bool=true, index=nothing)
    nframe = ta.nframe
    natom = ta.natom
    if !isweight || isempty(ta.mass)
        weight = ones(Float64, natom)
    else
        weight = ta.mass
    end
    if index == nothing
        index = 1:natom
    end
    com_x = Vector{Float64}(undef, nframe)
    com_y = Vector{Float64}(undef, nframe)
    com_z = Vector{Float64}(undef, nframe)
    weight_sum = sum(weight)
    for iframe = 1:nframe
        com_x[iframe] = Statistics.mean(ta.x[iframe, index] .* weight) / weight_sum
        com_y[iframe] = Statistics.mean(ta.y[iframe, index] .* weight) / weight_sum
        com_z[iframe] = Statistics.mean(ta.z[iframe, index] .* weight) / weight_sum
    end
    com = TrjArray(x=com_x, y=com_y, z=com_z)
end

###### distance, angle, dihedral #################
"""
calcdistance

calculate distances between two atoms or groups of atoms
"""
function distance(ta1::TrjArray, ta2::TrjArray)
    # TODO: support for PBC
    @assert ta1.nframe == ta2.nframe
    nframe = ta1.nframe
    com1 = centerofmass(ta1)
    com2 = centerofmass(ta2)
    dist = Vector{Float64}(undef, nframe)
    for iframe = 1:nframe
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

