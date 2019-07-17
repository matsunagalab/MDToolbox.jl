function calc_bias_potential(umbrella_center, data)
    spring_constant = 200.0 * (pi/180.0)^2
    return spring_constant * (data - umbrella_center)^2
end