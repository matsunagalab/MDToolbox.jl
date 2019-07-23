include("Density.jl");
include("BiasPotential.jl");

function calc_delta_potential_of_mean_force(umbrella_center, data, sigma, kbt, input)
    delta_pmf = -kbt * log(calc_density(input, data, sigma))
    delta_pmf += kbt * log(calc_density(umbrella_center, data, sigma))
    delta_pmf -= calc_bias_potential(umbrella_center, input) - calc_bias_potential(umbrella_center, umbrella_center)
    return delta_pmf
end

function calc_potential_of_mean_force(umbrella_centers, data, weight, sigma)
    pmf = []
    for i = 1:length(data)
        tmp = 0.0
        for j = 1:length(umbrella_centers)
              tmp += weight[j] * get_gaussian(data[i], umbrella_centers[j], sigma)
        end
        push!(pmf, tmp)
    end
    
    return pmf
end