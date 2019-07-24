function calc_bias_potential(umbrella_center, sim_point_data, kbt)
    spring_constant = 50.0 * (pi/180.0).^2
    
    bias_potential = norm(sim_point_data - umbrella_center)
    #for i = 1:length(umbrella_center)
    #    bias_potential += (minimum_image(umbrella_center[i], sim_point_data[i])).^2
    #end
    
    bias_potential = bias_potential .- round.(bias_potential./360.0).*360.0
    
    return (spring_constant ./ kbt) * bias_potential.^2
end

function minimum_image(center, x)
  dx = x .- center
  dx = dx .- round.(dx./360.0).*360.0;
  dx
end