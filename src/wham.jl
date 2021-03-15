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

function wham_equation!(F, x, K, M, bias_km, log_Neff_km, log_heff_m)
    log_prob_m = reshape(x[1:M], 1, M)
    f_k = reshape(x[(M+1):(M+K)], K, 1)

    # 1st equation
    log_denominator_km = f_k .- bias_km
    log_denominator_km = log_denominator_km .+ log_Neff_km
    log_denominator_m  = logsumexp_dim1(log_denominator_km)
    log_numerator_m = log_heff_m
    F[1:M] = - log_prob_m .+ log_numerator_m .- log_denominator_m

    # 2nd equation
    F[(M+1):(M+K)] = f_k .+ logsumexp_dim2(log_prob_m .- bias_km)
    F[(M+1):(M+K)] = F[(M+1):(M+K)] .- F[M+1]
end

"""
    wham(h_km, bias_km; ftol=1e-10, maxiterations=10^5) -> F

Estimates (reduced) free energies of umbrella-windows and potential of mean force in data-bins by using the WHAM equations.
K is # of umbrellas and M is # of bins. `h_km` is a K x M Array which is a histogram (data counts) of k-th umbrella data in m-th bin.
`bias_km` is also a K x M Array which is a bias-factor of k-th umbrella-window evaluated at m-th bin-center. 

Returns a NamedTuple object `F` whose members are `F.f_k` reduced relative free energies of umbrella-windows, 
and pmf_m reduced potential of mean force in data-bins under unbiased condition. 

# References
```
S. Kumar, D. Bouzida, R. H. Swendsen, P. A. Kollman, and J. M. Rosenberg, 
J. Comput. Chem. 13, 1011 (1992).
```
"""
function wham(h_km, bias_km; ftol=1e-10, iterations=10^5)
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

    #F_tmp = zeros(Float64, M+K)
    #initial_x = zeros(Float64, M+K)
    #wham_equation!(F_tmp, initial_x, K, M, bias_km, log_Neff_km, log_heff_m)

    f!(F, x) = wham_equation!(F, x, K, M, bias_km, log_Neff_km, log_heff_m)
    initial_x = zeros(Float64, M+K)
    sol = nlsolve(f!, initial_x, autodiff=:forward, ftol=ftol, iterations=iterations)

    log_prob_m = reshape(sol.zero[1:M], 1, M)
    f_k = reshape(sol.zero[(M+1):(M+K)], K, 1)
    pmf_m = - log_prob_m
    return (f_k=f_k .- f_k[1, 1], pmf=pmf_m .- pmf_m[1, 1])
end

"""
    wham_iteration(h_km, bias_km)

Old version of wham(). This function iteratively solve the wham equations.

# References
```
S. Kumar, D. Bouzida, R. H. Swendsen, P. A. Kollman, and J. M. Rosenberg, 
J. Comput. Chem. 13, 1011 (1992).
```
"""
function wham_iteration(h_km, bias_km)
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
    f_k .- f_k[1, 1], pmf_m .- pmf_m[1, 1]
end
