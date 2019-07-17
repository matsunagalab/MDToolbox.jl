include("Density.jl")

function calc_design_matrix(umbrella_centers, data, sigma)
    design_matrix = zeros(Float64, length(umbrella_centers) * length(data[1]), length(umbrella_centers))
    for i = 1:length(umbrella_centers)
        for j = 1:length(data[i])
            for k = 1:length(umbrella_centers)
                diff = get_gaussian(data[i][j], umbrella_centers[k], sigma) - get_gaussian(umbrella_centers[i], umbrella_centers[k], sigma)
                design_matrix[j + length(data[i]) * (i - 1), k] = diff
            end
        end
    end
    
    return design_matrix
end