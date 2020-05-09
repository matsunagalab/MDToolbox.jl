function sp_delta_pmf(umbrella_center, data_k, kbt, spring_constant)
    nframe = size(data_k[1], 1)
    K = length(data_k)
    delta_pmf      = similar(umbrella_center, nframe, K)
    bias_potential = similar(umbrella_center, nframe, K)
    delta_pmf .= zero(umbrella_center[1])
    bias_potential .= zero(umbrella_center[1])

    for k = 1:K
        dd_mean = mean(data_k[k], dims=1)
        dd_cov = cov(data_k[k])
        #dd_inv = dd_cov \ CuArray{Float64}(I, size(dd_cov))
        dd_inv = dd_cov \ copyto!(similar(dd_cov, size(dd_cov,1), size(dd_cov,1)), I)
        bias_potential = spring_constant * sum((data_k[k] .- umbrella_center[k:k, :]).^2, dims=2)
        d = data_k[k] .- dd_mean
        tmp_pmf = kbt .* (-0.5 .* sum((d * dd_inv) .* d, dims=2))
        delta_pmf[:, k:k] .= .- tmp_pmf .+ tmp_pmf[1, 1] .- bias_potential .+ bias_potential[1, 1]
    end

    return delta_pmf
end

#######################################
function sp_design_matrix(umbrella_center, sim_all_datas, sigma_g; stride=1)
    nframe = size(sim_all_datas[1], 1)
    subframes = range(1, stop=nframe, step=stride)
    #design_matrix = CuArrays.zeros(Float64, length(subframes) * size(umbrella_center, 1), size(umbrella_center, 1))
    design_matrix = similar(umbrella_center, length(subframes) * size(umbrella_center, 1), size(umbrella_center, 1))
    design_matrix .= zero(umbrella_center[1])

    for i = 1:size(umbrella_center, 1)
        for k = 1:size(umbrella_center, 1)
            diff = exp.(- 0.5 .* sum((sim_all_datas[i] .- umbrella_center[k:k, :]).^2, dims=2) ./ sigma_g.^2)
            diff .-= exp(- 0.5 .* sum((sim_all_datas[i][1, :] .- umbrella_center[k, :]).^2) ./ sigma_g.^2)
            design_matrix[(length(subframes)*(i-1)+1):(i*length(subframes)), k] = diff
        end
    end

    return design_matrix
end

#######################################
function sp_design_matrix_atom(umbrella_center, sim_all_datas, sigma_g)
    nframe = size(sim_all_datas[1], 1)
    numbrella = size(umbrella_center, 1)
    natom = Int(size(umbrella_center, 2) / 3)
    #design_matrix = CuArrays.zeros(Float64, nframe * numbrella, numbrella * natom)
    design_matrix = similar(umbrella_center, nframe * numbrella, numbrella * natom)
    design_matrix .= zero(umbrella_center[1])

    for i = 1:size(umbrella_center, 1)
        for k = 1:size(umbrella_center, 1)
            diff_x1 = exp.(- 0.5 .* ((sim_all_datas[i][:, 1:3:end]   .- umbrella_center[k:k, 1:3:end]).^2) ./ sigma_g.^2)
            diff_y1 = exp.(- 0.5 .* ((sim_all_datas[i][:, 2:3:end]   .- umbrella_center[k:k, 2:3:end]).^2) ./ sigma_g.^2)
            diff_z1 = exp.(- 0.5 .* ((sim_all_datas[i][:, 3:3:end]   .- umbrella_center[k:k, 3:3:end]).^2) ./ sigma_g.^2)
            diff_x2 = exp.(- 0.5 .* ((sim_all_datas[i][1:1, 1:3:end] .- umbrella_center[k:k, 1:3:end]).^2) ./ sigma_g.^2)
            diff_y2 = exp.(- 0.5 .* ((sim_all_datas[i][1:1, 2:3:end] .- umbrella_center[k:k, 2:3:end]).^2) ./ sigma_g.^2)
            diff_z2 = exp.(- 0.5 .* ((sim_all_datas[i][1:1, 3:3:end] .- umbrella_center[k:k, 3:3:end]).^2) ./ sigma_g.^2)
            design_matrix[(nframe*(i-1)+1):(nframe*i), (natom*(k-1)+1):(natom*k)] .= (diff_x1 .* diff_y1 .* diff_z1) .- (diff_x2 .* diff_y2 .* diff_z2)
        end
    end
    return design_matrix
end

#######################################
function funcS(x, lambda)
    return sign(x) * max(abs(x) - lambda, 0.0)
end

#######################################
"""
Least squares
"""
function sp_lsquares(y, X)
    # y = X * weight
    U, S, V = svd(X) # svd
    inverse_M = V * inv(Diagonal(S)) * U' # pseudo-inverse
    beta = inverse_M * y

    residual = sum((y - X*beta).^2)
    return beta, residual
end

"""
Alternating Direction Method of Multipliers (ADMM) for solving lasso
"""
function sp_admm(y, X, lambda=0.1; rho=1.0, condition=1e-5, iter_max=10000)
    ncolumn = size(X, 2);
    #beta = CuArrays.zeros(Float64, ncolumn)
    #gamma = CuArrays.zeros(Float64, ncolumn)
    #my = CuArrays.zeros(Float64, ncolumn)
    beta = similar(y, ncolumn)
    gamma = similar(y, ncolumn)
    my = similar(y, ncolumn)
    beta .= zero(y[1])
    gamma .= zero(y[1])
    my .= zero(y[1])

    old_beta = copy(beta)
    old_gamma = copy(gamma)
    old_my = copy(my)

    B = X' * X
    #B += (I * (ncolumn * rho))
    B += copyto!(similar(y, size(X, 2), size(X, 2)), I) .* (ncolumn .* rho)
    #U, S, V = svd(X' * X + ((CuArray{Float64}(I, (ncolumn, ncolumn)) .* (ncolumn .* rho))))
    U, S, V = svd(B)
    #U, S, V = svd(B, CuArrays.CUSOLVER.JacobiAlgorithm)
    #U, S, V = svd(B, CuArrays.CUSOLVER.QRAlgorithm)
    pinverse_M = V * inv(Diagonal(S)) * U'
    const_num = pinverse_M
    #const_num = B \ CuArray{Float64}(I, size(B))

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

    residual = sum((y - X*beta).^2)

    return gamma, residual
end

#"""
#Coordinate descent
#"""
#function sp_descent(y, X, lambda=0.1; condition=1e-5, iter_max=10000)
#    nframe = size(X, 1)
#    nfeature = size(X, 2)
#    rhs1 = X' * y
#    B = X' * X
#
#    beta = similar(y, nfeature)
#    beta_old = similar(y, nfeature)
#    beta_old .= one(y[1])
#    rhs2 = similar(y, nfeature)
#    rhs3 = similar(y, nfeature)
#    rhs = similar(y, nfeature)
#    max_diff = 100.0
#    cnt = 1
#    while max_diff > condition
#        rhs2 .= B * beta_old
#        rhs3 .= nframe .* beta_old
#        rhs .= rhs1 .- rhs2 .+ rhs3
#        #@show rhs./nframe
#        beta .= funcS.(rhs./nframe, lambda)
#        #@show beta
#        max_diff = maximum(abs.(beta .- beta_old))
#        beta_old .= beta
#        cnt += 1
#        if cnt > iter_max
#            break
#        end
#    end
#
#    println("[ Cycle Count = ", cnt, " ]")
#    println("[ Complete Condition ]")
#    println("  Max Differ = ", max_diff)
#    println("\n")
#
#    return beta
#end

"""
Coordinate descent
"""
function sp_descent(y, X, lambda=0.1; condition=1e-5, iter_max=10000)
    nframe = size(X, 1)
    nfeature = size(X, 2)
    Xty = X' * y
    XtX = X' * X

    beta = similar(y, nfeature)
    beta_old = similar(y, nfeature)
    beta .= one(y[1])
    beta_old .= beta
    #rhs2 = similar(y, nfeature)
    #rhs3 = similar(y, nfeature)
    rhs = similar(y, nfeature)
    max_diff = 100.0
    cnt = 1
    while max_diff > condition
        cnt += 1
        #@show cnt, max_diff
        for j = 1:nfeature
            rhs1 = Xty[j]
            rhs2 =  sum(XtX[j, :] .* beta)
            rhs3 = nframe * beta[j]
            rhs = rhs1 - rhs2 + rhs3
            beta[j] = funcS(rhs/nframe, lambda)
        end
        max_diff = maximum(abs.(beta .- beta_old))
        beta_old .= beta
        if cnt > iter_max
            break
        end
    end

    println("[ Cycle Count = ", cnt, " ]")
    println("[ Complete Condition ]")
    println("  Max Differ = ", max_diff)
    println("\n")

    residual = sum((y - X*beta).^2)

    return beta, residual
end


#######################################
function sp_standardize!(M)
    mean_M = mean(M, dims=1)
    std_M = sqrt.(sum((M .- mean_M).^2, dims=1)) ./ sqrt(size(M, 1))
    M .= (M .- mean_M) ./ std_M
    return mean_M, std_M
end

#######################################
function sp_standardize(M)
    mean_M = mean(M, dims=1)
    #std_M = std(M, dims=1)
    std_M = sqrt.(sum((M .- mean_M).^2, dims=1)) ./ sqrt(size(M, 1))
    M_standardized = (M .- mean_M) ./ std_M
    return M_standardized, mean_M, std_M
end

#######################################
function sp_cumulate_pmf(x, weight, umbrella_center, sigma_rdf, mean_M, std_M)
    ndim = size(umbrella_center, 2)
    K = size(umbrella_center, 1)
    nframe = size(x, 1)

    pmf = zeros(Float64, nframe);
    for iframe = 1:nframe
        sum_rdf = 0.0
        for k = 1:K
            tmp = (exp(sum(- 0.5 .* (x[iframe, :] .- umbrella_center[k, :]).^2 ./ sigma_rdf.^2))  - mean_M[k]) / std_M[k]
            sum_rdf += tmp * weight[k]
        end
        pmf[iframe] = sum_rdf
    end

    return pmf .- minimum(pmf)
end

#######################################
function sp_cumulate_pmf_atom(data, weight, umbrella_center, sigma_rdf, mean_M, std_M)
    ndim = size(umbrella_center, 2)
    natom = Int(ndim / 3)
    natom3 = natom * 3
    K = size(umbrella_center, 1)
    nframe = size(data, 1)

    pmf = zeros(Float64, nframe);
    for iframe = 1:nframe
        sum_rdf = 0.0
        for k = 1:K
            index = ((k-1)*natom+1):(k*natom)
            index_x = 1:3:natom3
            index_y = 2:3:natom3
            index_z = 3:3:natom3
            tmp   = (data[iframe:iframe, index_x] .- umbrella_center[k:k, index_x]).^2
            tmp .+= (data[iframe:iframe, index_y] .- umbrella_center[k:k, index_y]).^2
            tmp .+= (data[iframe:iframe, index_z] .- umbrella_center[k:k, index_z]).^2
            tmp ./= sigma_rdf.^2
            tmp .*= -0.5
            tmp .= (exp.(tmp) .- mean_M[1:1, index]) ./ std_M[1:1, index]
            sum_rdf += sum(tmp[:] .* weight[index])
            #for iatom = 1:natom
            #    index = (k-1)*natom + iatom
            #    index3 = ((iatom-1)*3 + 1):(iatom*3)
            #    tmp = (exp(sum(- 0.5 .* (data[iframe, index3] .- umbrella_center[k, index3]).^2 ./ sigma_rdf.^2))  - mean_M[index]) / std_M[index]
            #    sum_rdf += tmp * weight[index]
            #end
        end
        pmf[iframe] = sum_rdf
    end

    return pmf .- minimum(pmf)
end

#######################################
function sp_lasso(M, delta_pmf, sigma_rdf)
end
