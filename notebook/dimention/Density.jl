# 各点の密度を求める
function calc_density(input_point, sim_umbrella_datas, sigma_matrix)
    density = 0.0
    for i = 1:size(sim_umbrella_datas, 1)
        density += calc_gaussian(input_point, sim_umbrella_datas[i, :], sigma_matrix)
    end
    
    density /= size(sim_umbrella_datas, 1)
    return density
end

# 多変量ガウス
function calc_gaussian(input, mean, sigma_matrix)
    return pdf(MvNormal(mean, sigma_matrix), input) * (((2*pi)^(length(mean)/2)) * sqrt(det(sigma_matrix)))
end
