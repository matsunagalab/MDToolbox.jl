function sp_delta_pmf(umbrella_centers, sim_all_datas, kbt, spring_constant)
    nframe = size(sim_all_datas[1], 1)
    numbrella = size(umbrella_centers, 1)
    delta_pmf      = CuArrays.zeros(Float64, nframe, numbrella)
    bias_potential = CuArrays.zeros(Float64, nframe, numbrella)

    for i = 1:numbrella
        dd_mean = mean(sim_all_datas[i], dims=1)
        dd_cov = cov(sim_all_datas[i])
        dd_inv = dd_cov \ CuArray{Float64}(I, size(dd_cov))
        bias_potential = spring_constant * sum((sim_all_datas[i] .- umbrella_centers[i:i, :]).^2, dims=2)
        d = sim_all_datas[i] .- dd_mean
        tmp_pmf = kbt .* (-0.5 .* sum((d * dd_inv) .* d, dims=2))
        delta_pmf[:, i:i] .= .- tmp_pmf .+ tmp_pmf[1, 1] .- bias_potential .+ bias_potential[1, 1]
    end

    return delta_pmf
end

#######################################
function calc_design_matrix(umbrella_centers, sim_all_datas, stride, sigma_g)
    nframe = size(sim_all_datas[1], 1)
    subframes = range(1, stop=nframe, step=stride)
    design_matrix = CuArrays.zeros(Float64, length(subframes) * size(umbrella_centers, 1), size(umbrella_centers, 1))

    for i = 1:size(umbrella_centers, 1)
        for k = 1:size(umbrella_centers, 1)
            diff = exp.(- 0.5 .* sum((sim_all_datas[i] .- umbrella_centers[k:k, :]).^2, dims=2) ./ sigma_g.^2)
            diff .-= exp(- 0.5 .* sum((sim_all_datas[i][1, :] .- umbrella_centers[k, :]).^2) ./ sigma_g.^2)
            design_matrix[(length(subframes)*(i-1)+1):(i*length(subframes)), k] = diff
        end
    end

    return design_matrix
end

#######################################
function sp_design_matrix_xyz(umbrella_centers, sim_all_datas, sigma_g)
    nframe = size(sim_all_datas[1], 1)
    numbrella = size(umbrella_centers, 1)
    natom = Int(size(umbrella_centers, 2) / 3)
    design_matrix = CuArrays.zeros(Float64, nframe * numbrella, numbrella * natom)

    for i = 1:size(umbrella_centers, 1)
        for k = 1:size(umbrella_centers, 1)
            diff_x1 = exp.(- 0.5 .* ((sim_all_datas[i][:, 1:3:end]   .- umbrella_centers[k:k, 1:3:end]).^2) ./ sigma_g.^2)
            diff_y1 = exp.(- 0.5 .* ((sim_all_datas[i][:, 2:3:end]   .- umbrella_centers[k:k, 2:3:end]).^2) ./ sigma_g.^2)
            diff_z1 = exp.(- 0.5 .* ((sim_all_datas[i][:, 3:3:end]   .- umbrella_centers[k:k, 3:3:end]).^2) ./ sigma_g.^2)
            diff_x2 = exp.(- 0.5 .* ((sim_all_datas[i][1:1, 1:3:end] .- umbrella_centers[k:k, 1:3:end]).^2) ./ sigma_g.^2)
            diff_y2 = exp.(- 0.5 .* ((sim_all_datas[i][1:1, 2:3:end] .- umbrella_centers[k:k, 2:3:end]).^2) ./ sigma_g.^2)
            diff_z2 = exp.(- 0.5 .* ((sim_all_datas[i][1:1, 3:3:end] .- umbrella_centers[k:k, 3:3:end]).^2) ./ sigma_g.^2)
            design_matrix[(nframe*(i-1)+1):(nframe*i), (natom*(k-1)+1):(natom*k)] .= (diff_x1 .* diff_y1 .* diff_z1) .- (diff_x2 .* diff_y2 .* diff_z2)
        end
    end
    return design_matrix
end

#######################################
function funcS(cov, lambda)
    return sign(cov) * max((abs(cov) - lambda), 0)
end

#######################################
"""
Alternating Direction Method of Multipliers (ADMM) for solving lasso
"""
function sp_admm(y, X, lambda, rho=1.0, condition=1e-5, iter_max=10000)
    ncolumn = size(X, 2);
    beta = CuArrays.zeros(Float64, ncolumn)
    gamma = CuArrays.zeros(Float64, ncolumn)
    my = CuArrays.zeros(Float64, ncolumn)

    old_beta = copy(beta)
    old_gamma = copy(gamma)
    old_my = copy(my)

    U, S, V = svd(X' * X + ((CuArray{Float64}(I, (ncolumn, ncolumn)) .* (ncolumn .* rho))))
    inverse_M = V * inv(Diagonal(S)) * U'
    const_num = inverse_M

    max_diff = 100
    cnt::Int = 1
    X_y = X' * y
    while max_diff > condition
        beta .= const_num * (X_y .+ ((gamma .- (my .* (1.0 ./ rho))) .* (ncolumn .* rho)))
        gamma .= funcS.(beta .+ (my .* (1.0 ./ rho)), lambda ./ rho)
        my .= my .+ ((beta .- gamma) .* rho)

        max_diff = maximum(abs.(gamma .- old_gamma))
        old_beta .= beta
        old_gamma .= gamma
        old_my .= my
        cnt += 1
        if cnt > iter_max
            break
        end
    end

    println("[ Cycle Count = ", cnt, " ]")
    println("[ Complete Condition ]")
    println("  Max Differ = ", max_diff)
    println("\n")

    return gamma
end

#######################################
function sp_standardize!(M)
    mean_M = mean(M, dims=1)
    std_M = std(M, dims=1)
    M .= (M .- mean_M) ./ std_M
    return mean_M, std_M
end

#######################################
function sp_standardize(M)
    mean_M = mean(M, dims=1)
    std_M = std(M, dims=1)
    M_standardized = (M .- mean_M) ./ std_M
    return M_standardized, mean_M, std_M
end

#######################################
function sp_cumulate_pmf(M)
    pmf = []
    ndim = size(sample_datas, 2)
    ncenters = size(umbrella_centers, 1)
    nframe = size(sample_datas, 1)
    for i = 1:nframe #size(sample_datas, 1)
        sum_gaussian = 0.0
        for j = 1:ncenters #size(umbrella_centers, 1)
            sum_gaussian += weights[j] * (calc_gaussian(sample_datas[i, :], umbrella_centers[j, :], sigma2) - mean_M[j]) / std_M[j]
        end
        push!(pmf, sum_gaussian)
    end
    return pmf
end
