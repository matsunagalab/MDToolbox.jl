function logsumexp_dim1(x_km)
  max_m = maximum(x_km, dims=1)
  exp_km = exp.(x_km .- max_m)
  s_m = log.(sum(exp_km, dims=1)) .+ max_m
end

function logsumexp_dim2(x_km)
  max_k = maximum(x_km, dims=2)
  exp_km = exp.(x_km .- max_k)
  s_k = log.(sum(exp_km, dims=2)) .+ max_k
end

"""
wham

estimates (reduced) free energies of umbrella-windows and potential of mean force in data-bins by using the WHAM equations

References
[1] S. Kumar, D. Bouzida, R. H. Swendsen, P. A. Kollman, and
    J. M. Rosenberg, J. Comput. Chem. 13, 1011 (1992).
"""
function wham(h_km, bias_km)
    TOLERANCE = 10^(-8)

    # K: number of umbrella-windows
    K = size(h_km, 1)

    # M: number of bins
    M = size(h_km, 2)

    # conversion from Int64 to Float64
    h_km = map(Float64, h_km)

    # calculate counts (N_k)
    # N_k: number of samples counted in the histogram from umbrella-window (k)
    N_k = sum(h_km, dims=2)
    log_N_k = zeros(Float64, K, 1)
    log_N_k[N_k .> 0.0] .= log.(N_k[N_k .> 0.0])

    # calculate effective histogram (heff_m)
    # heff_m: effective number of independent samples in data-bin (m)
    g_km = ones(Float64, K, M)
    heff_m = sum(h_km./g_km, dims=1)
    log_heff_m = zeros(Float64, 1, M)
    log_heff_m[:, :] .= -10.0 #double(log(eps(Float32)))
    log_heff_m[heff_m .> 0.0] .= log.(heff_m[heff_m .> 0.0])

    # calculate effective counts (Neff_km)
    # Neff_km: effective number of independent samples for umbrella-window (k) in data-bin (m)
    Neff_km = N_k ./ g_km
    log_Neff_km = zeros(Float64, K, M)
    log_Neff_km[:, :] .= -10.0 #double(log(eps(Float32)))
    log_Neff_km[Neff_km .> 0] .= log.(Neff_km[Neff_km .> 0])

    # solve the WHAM equations by self-consistent iteration
    f_k = zeros(Float64, K, 1)
    check_convergence = Inf

    log_prob_m = zeros(Float64, 1, M)
    count_iteration = 0
    while check_convergence > TOLERANCE

      # 1st equation
      log_denominator_km = f_k .- bias_km
      log_denominator_km = log_denominator_km .+ log_Neff_km
      log_denominator_m  = logsumexp_dim1(log_denominator_km)

      log_numerator_m = log_heff_m

      log_prob_m = log_numerator_m .- log_denominator_m

      # 2nd equation
      f_k_new = - logsumexp_dim2(log_prob_m .- bias_km)
      f_k_new = f_k_new .- f_k_new[1]

      # check convergence
      check_convergence = maximum(abs.(f_k_new - f_k))./std(f_k_new)
      f_k = f_k_new;

      count_iteration = count_iteration + 1
      if mod(count_iteration, 100) == 0
        @printf("%dth iteration  delta = %e  tolerance = %e\n", count_iteration, check_convergence, TOLERANCE)
        @printf("free energies = ")
        for k = 1:K
          @printf("%f ", f_k[k])
        end
        @printf("\n")
        @printf("\n")
      end
    end

    pmf_m = - log_prob_m
    f_k, pmf_m
end
