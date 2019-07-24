include("Density.jl")
using Base.Threads

function calc_design_matrix(umbrella_centers, sim_all_datas, sigma_matrix)
    design_matrix = zeros(Float64, size(umbrella_centers, 1) * size(sim_all_datas[1], 1), size(umbrella_centers, 1))
    
    # i, jのループで全シミュレーションデータ点にアクセス
    Threads.@threads for i = 1:size(umbrella_centers, 1)
        for j = 1:size(sim_all_datas[i], 1)
            for k = 1:size(umbrella_centers, 1)
                diff = calc_gaussian(sim_all_datas[i][j, :], umbrella_centers[k, :], sigma_matrix)
                diff -= calc_gaussian(umbrella_centers[i, :], umbrella_centers[k, :], sigma_matrix)
                
                design_matrix[j + size(sim_all_datas[i], 1) * (i - 1), k] = diff
            end
        end
    end
    
    return design_matrix
end