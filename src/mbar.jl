#function logsumexp(x)
#    max_x = maximum(x)
#    exp_x = exp.(x .- max_x)
#    log(sum(exp_x)) .+ max_x
#end

function mbar_equation!(F, f_k, K, N_max, N_k, u_kln, idx, log_wi_jn, log_term, FlogN)
    for k = 1:K
        FlogN[k] = log(Float64(N_k[k])) + f_k[k]
    end
    for i = 1:K
        for j = 1:K
            for n = 1:N_k[j]
                ui_xjn = u_kln[j, i, n];
                for k = 1:K
                    uk_xjn = u_kln[j, k, n];
                    log_term[k] = FlogN[k] - (uk_xjn - ui_xjn);
                end
                log_wi_jn[j, n] = -logsumexp(log_term);
            end
        end
        F[i] = f_k[i] + logsumexp(log_wi_jn[idx])
    end
    F[:] = F[:] .- F[1]
end

"""
    mbar(u_kl; tol=1e-8, iterations_self=10) -> F

Estimates the free energy differences between different thermodynamnic states by using the Multistate Bennet Acceptance Ratio Method (MBAR).
Let K be # of thermodynamic states. `u_kl` is a K x K Array whose elements are reduced bias-factors or 
potential energies of state k evaluated at state l. Solvs the MBAR equations by applying self-consistent iterations `iterations_self` times, 
followed by the Newton-Raphson iterations. The iteration stops when the difference from the previous falls below the threshold `tol`. 

Returns dimensionless free energies `f_k` whose elements corresponds to thermodynamics state k. 

# Example
```julia-repl
see Juypter notebooks in MDToolbox.jl/examples/
```

# References
```
M. R. Shirts and J. D. Chodera, J. Chem. Phys. 129, 124105 (2008).
```
"""
function mbar(u_kl; tol=1e-8, iterations_self=10)
    # K: number of umbrella windows
    K, L = size(u_kl)

    # N_k: number of data in k-th umbrella window
    N_k = zeros(Int64, K)
    for k = 1:K
        N_k[k] = length(u_kl[k, 1])
    end
    N_max = maximum(N_k)
    
    # conversion from array of array (u_kl) to array (u_kln)
    u_kln = zeros(Float64, K, K, N_max)
    for k = 1:K
        for l = 1:K
            u_kln[k, l, 1:N_k[k]] .= u_kl[k, l]
        end
    end

    dummy_kn = zeros(Int64, K, N_max)
    for k = 1:K
      dummy_kn[k, 1:N_k[k]] .= 1
    end
    idx = dummy_kn .> 0;

    # solve the MBAR equation
    #log_wi_jn = zeros(Float64, K, N_max)
    #log_term = zeros(Float64, K)
    #FlogN = zeros(Float64, K)
    ## f!(F, x) = mbar_equation!(F, x, K, N_max, N_k, u_kln, idx, log_wi_jn, log_term, FlogN)
    ## x_init = zeros(Float64, K)
    ## sol = nlsolve(f!, x_init, autodiff=:forward, ftol=tol, iterations=iterations)
    ## #sol = nlsolve(f!, x_init, method = :anderson, ftol=tol, iterations=iterations)
    ## f_k = sol.zero
    ## f_k .- f_k[1]

    #F = zeros(Float64, K)
    #x_init = zeros(Float64, K)
    #for n = 1:1000
    #    mbar_equation!(F, x_init, K, N_max, N_k, u_kln, idx, log_wi_jn, log_term, FlogN)
    #    x_init = x_init - F
    #    x_init = x_init .- x_init[1]
    #end
    #f_k = x_init

    # solve the MBAR equation by self-consistent iterations
    f_k = zeros(Float64, K)
    f_k_new = zeros(Float64, K)

    log_wi_jn = zeros(Float64, (K, N_max))
    for j = 1:K
      log_wi_jn[j, 1:N_k[j]] .= 1
    end
    index = log_wi_jn .> 0.5

    for count_iteration = 1:iterations_self
        for i = 1:K
            log_wi_jn = mbar_log_wi_jn(N_k, f_k, u_kln, u_kln[:, i, :], K, N_max)
            f_k_new[i] = - logsumexp_1d(log_wi_jn[index])
        end
          
        f_k_new .= f_k_new .- f_k_new[1]
        check_convergence = maximum(abs.(f_k_new .- f_k)) ./ std(f_k_new)
        f_k .= f_k_new

        @printf "iteration = %d  delta = %e  tol = %e\n" count_iteration check_convergence tol
        @printf "free energies = "
        for k = 1:K
            @printf " %f" f_k[k]
        end
        @printf "\n"
        @printf "\n"
    end

    first_gamma = 0.1;
    gamma = 1.0;
    
    check_convergence = Inf64
    W_nk = zeros(Float64, sum(N_k), K)
    
    count_iteration = iterations_self
    while check_convergence > tol
        f_k_new .= f_k
        for i = 1:K
          log_wi_jn = mbar_log_wi_jn(N_k, f_k, u_kln, u_kln[:, i, :], K, N_max)
          W_nk[:, i] .= exp.(log_wi_jn[index] .+ f_k[i])
        end
 
        g = zeros(Float64, K-1, 1) # gradient
        H = zeros(Float64, K-1, K-1) # hessian
        for i = 1:(K-1)
            g[i] = N_k[i+1] - N_k[i+1] * sum(W_nk[:, i+1])
            H[i, i] = - sum(N_k[i+1] .* W_nk[:, i+1] .* (1.0 .- N_k[i+1] .* W_nk[:, i+1]))
            for j = 1:(i-1)
                H[i, j] = sum((N_k[i+1] .* W_nk[:, i+1]) .* (N_k[j+1] .* W_nk[:, j+1]))
                H[j, i] = H[i, j]
            end
        end
    
        Hinvg = pinv(H, atol=1.0e-10) * g
        for k = 1:(K-1)
            if count_iteration == iterations_self
                f_k_new[k+1] = f_k_new[k+1] - first_gamma * Hinvg[k]
            else
                f_k_new[k+1] = f_k_new[k+1] - gamma * Hinvg[k]
            end
        end

        check_convergence = maximum(abs.(f_k_new .- f_k)) ./ std(f_k_new)
        f_k .= f_k_new
    
        count_iteration = count_iteration + 1
        @printf "iteration =%d  delta = %e  tol = %e\n" count_iteration check_convergence tol
        @printf "free energies = "
        for k = 1:K
            @printf " %f" f_k[k]
        end
        @printf "\n"
        @printf "\n"
    end
    
    return f_k
end

"""
    mbar_f(u_kl, f_k, u_k=nothing)

Calculate the free energy value at the target thermodynamics state (ensemble) 
by using the Multistate Bennet Acceptance Ratio Method (MBAR). 
Let K be # of thermodynamic states. `u_kl` is a K x K Array whose elements are reduced bias-factors or 
potential energies of state k evaluated at state l. `f_k` is dimension dimensionless free energies `f_k` 
obtained by prior MBAR calculation. If `u_k` (reduced potential energies of the target state) is given, 
weights for that state are computed. By default `u_k = nothing`, in thise case zero potential 
energies are used for computing weights (useful for umbrella sampling data). 

Returns an estimated free energy value `f` computed at the target thermodynamic state. 

# Example
```julia-repl
see Juypter notebooks in MDToolbox.jl/examples/
```

# References
```
M. R. Shirts and J. D. Chodera, J. Chem. Phys. 129, 124105 (2008).
```
"""
function mbar_f(u_kl, f_k, u_k=nothing)
    # K: number of umbrella windows
    K, L = size(u_kl)

    # N_k: number of data in k-th umbrella window
    N_k = zeros(Int64, K)
    for k = 1:K
        N_k[k] = length(u_kl[k, 1])
    end
    N_max = maximum(N_k)
    
    # conversion from array of array (u_kl) to array (u_kln)
    u_kln = zeros(Float64, K, K, N_max)
    for k = 1:K
        for l = 1:K
            u_kln[k, l, 1:N_k[k]] .= u_kl[k, l]
        end
    end

    # conversion from cell (u_k) to array (u_kn)
    u_kn = zeros(Float64, K, N_max)
    for k = 1:K
        if u_k === nothing
            u_kn[1, 1:N_k[k]] .= zero(Float64)
        else
            u_kn[k, 1:N_k[k]] .= u_k[k]
        end
    end

    log_w_kn = zeros(Float64, K, N_max)
    for k = 1:K
      log_w_kn[k, 1:N_k[k]] .= 1.0
    end
    idx = log_w_kn .> 0.5;

    log_w_kn = mbar_log_wi_jn(N_k, f_k, u_kln, u_kn, K, N_max)
    f_hat = - logsumexp_1d(log_wi_jn[idx])

    return f_hat
end

function ChainRulesCore.rrule(::typeof(mbar_f), u_kl, f_k, u_k)
    f_hat = mbar_f(u_kl, f_k, u_k)
    w_k = mbar_weight(u_kl, f_k, u_k)
    function mbar_f_pullback(df)
        du_k = deepcopy(w_k)
        exp_minus_f_hat = exp(-f_hat)
        for k = 1:length(w_k)
            for n = 1:length(w_k[k])
                du_k[k][n] = exp_minus_f_hat * w_k[k][n] * df
            end
        end
        return NoTangent(), ZeroTangent(), NoTangent(), du_k
    end
    return f_hat, mbar_f_pullback
end

"""
    mbar_weight(u_kl, f_k, u_k=nothing)

Calculate weights for computing expectations at the target thermodynamics state (ensemble) 
by using the Multistate Bennet Acceptance Ratio Method (MBAR). 
Let K be # of thermodynamic states. `u_kl` is a K x K Array whose elements are reduced bias-factors or 
potential energies of state k evaluated at state l. `f_k` is dimension dimensionless free energies `f_k` 
obtained by prior MBAR calculation. If `u_k` (reduced potential energies of the target state) is given, 
weights for that state are computed. By default `u_k = nothing`, in thise case zero potential 
energies are used for computing weights (useful for umbrella sampling data). 

Returns weights `w_k` computed at the target thermodynamic state. 

# Example
```julia-repl
see Juypter notebooks in MDToolbox.jl/examples/
```

# References
```
M. R. Shirts and J. D. Chodera, J. Chem. Phys. 129, 124105 (2008).
```
"""
function mbar_weight(u_kl, f_k, u_k=nothing)
    # K: number of umbrella windows
    K, L = size(u_kl)

    # N_k: number of data in k-th umbrella window
    N_k = zeros(Int64, K)
    for k = 1:K
        N_k[k] = length(u_kl[k, 1])
    end
    N_max = maximum(N_k)
    
    # conversion from array of array (u_kl) to array (u_kln)
    u_kln = zeros(Float64, K, K, N_max)
    for k = 1:K
        for l = 1:K
            u_kln[k, l, 1:N_k[k]] .= u_kl[k, l]
        end
    end

    # conversion from cell (u_k) to array (u_kn)
    u_kn = zeros(Float64, K, N_max)
    for k = 1:K
        if u_k === nothing
            u_kn[1, 1:N_k[k]] .= zero(Float64)
        else
            u_kn[k, 1:N_k[k]] .= u_k[k]
        end
    end

    log_w_kn = zeros(Float64, K, N_max)
    for k = 1:K
      log_w_kn[k, 1:N_k[k]] .= 1.0
    end
    idx = log_w_kn .> 0.5;

    log_w_kn = mbar_log_wi_jn(N_k, f_k, u_kln, u_kn, K, N_max)
    log_w_n  = log_w_kn[idx]

    s = logsumexp_1d(log_w_n)
    w_k = Vector{Vector{Float64}}(undef, K)
    for k = 1:K
      w_k[k] = exp.((log_w_kn[k, 1:N_k[k]] .- s))
    end

    return w_k
end

function ChainRulesCore.rrule(::typeof(mbar_weight), u_kl, f_k, u_k)
    w_k = mbar_weight(u_kl, f_k, u_k)
    function mbar_weight_pullback(dw_k)
        du_k = deepcopy(w_k)
        for k = 1:length(w_k)
            for n = 1:length(w_k[k])
                #du_k[k][n] = dw_k[k][n] * (- w_k[k][n])
                du_k[k][n] = dw_k[k][n] * (- w_k[k][n] + w_k[k][n]^2)
            end
        end
        return NoTangent(), ZeroTangent(), NoTangent(), du_k
    end
    return w_k, mbar_weight_pullback
end

# MATLAB-style coding
function mbar_log_wi_jn(N_k, f_k, u_kln, u_kn, K, N_max)
    log_wi_jn = zeros(Float64, (K, N_max))
    for k = 1:K
        x = repeat(log.(N_k), 1, N_k[k]) .+ repeat(f_k, 1, N_k[k]) .- (u_kln[k, :, 1:N_k[k]] .- repeat(u_kn[k:k, 1:N_k[k]], K, 1))
        log_wi_jn[k:k, 1:N_k[k]] .= - logsumexp_over_row(x)
    end
    return log_wi_jn
end

# FORTRAN-style coding
function mbar_log_wi_jn2(N_k, f_k, u_kln, u_kn, K, N_max)
    log_wi_jn = zeros(Float64, (K, N_max))
    FlogN = zeros(Float64, K)
    for k = 1:K
        FlogN[k] = log(N_k[k]) + f_k[k]
    end

    log_term = zeros(Float64, K)
    #Threads.@threads for k = 1:K
    for k = 1:K
        for n = 1:N_k[k]
            max_log_term = -Inf64

            u_k = u_kn[k, n]
            for l = 1:K
                u_l = u_kln[k, l, n]
                log_term[l] = FlogN[l] - (u_l - u_k)
                if log_term[l] > max_log_term
                    max_log_term = log_term[l]
                end
            end

            term_sum = 0.0
            for l = 1:K
                term_sum += exp(log_term[l] - max_log_term)
            end
            log_sum = log(term_sum) + max_log_term
            log_wi_jn[k, n] = - log_sum
        end
    end
    return log_wi_jn
end

function logsumexp_1d(x)
    max_x = maximum(x)
    exp_x = exp.(x .- max_x)
    s = log(sum(exp_x)) + max_x
    return s
end

function logsumexp_over_row(x)
    max_x = maximum(x, dims=1)
    exp_x = exp.(x .- max_x)
    s = log.(sum(exp_x, dims=1)) .+ max_x
    return s
end

