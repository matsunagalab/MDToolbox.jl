include("Density.jl");
include("BiasPotential.jl");

function calc_delta_potential_of_mean_force(umbrella_center, sim_umbrella_datas, sigma_matrix, kbt, input_point)
    delta_pmf = -kbt * log(calc_density(input_point,     sim_umbrella_datas, sigma_matrix))
    delta_pmf += kbt * log(calc_density(umbrella_center, sim_umbrella_datas, sigma_matrix))
    delta_pmf -= calc_bias_potential(umbrella_center, input_point, kbt)
    return delta_pmf
end

function calc_potential_of_mean_force(umbrella_centers, sample_datas, weights, sigma_matrix)
    pmf = []
    for i = 1:size(sample_datas, 1)
        sum_gaussian = 0.0
        for j = 1:size(umbrella_centers, 1)
            sum_gaussian += weights[j] * calc_gaussian(sample_datas[i], umbrella_centers[j, :], sigma_matrix)
        end
        push!(pmf, sum_gaussian)
    end
    
    return pmf
end