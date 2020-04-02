"""
afmize

simulate afm image from a given structure

example::
pdb = readpbd("1ake.pdb");
decenter!(pdb);
afm_image = afmize(pdb, (150.0, 150.0), (16, 16));
heatmap(afm_image) # using Plots
"""
function afmize(ta::TrjArray, box, npixel)::Matrix{Float64}
    rpixel = box ./ npixel
    nx = floor.(Int, (ta.x[1, :] .+ 0.5*box[1]) ./ rpixel[1]);
    ny = floor.(Int, (ta.y[1, :] .+ 0.5*box[2]) ./ rpixel[2]);
    #@show nx
    #@show ny
    afm_image = zeros(Float64, npixel)
    z_min = minimum(ta.z[1, :])
    z = ta.z[1, :] .- z_min
    for iatom in 1:ta.natom
        if (1 <= nx[iatom] <= npixel[1]) && (1 <= ny[iatom] <= npixel[2])
            if z[iatom] > afm_image[nx[iatom], ny[iatom]]
                afm_image[nx[iatom], ny[iatom]] = z[iatom]
            end
        end
    end
    return afm_image
end
