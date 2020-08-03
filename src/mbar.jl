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
mbar

estimates the free energy differences of umbrella-windowed systems by using the Multistate Bennet Acceptance Ratio Method (MBAR).

References
[1] M. R. Shirts and J. D. Chodera, J. Chem. Phys. 129, 124105 (2008).
"""
function mbar(u_kl; ftol=1e-8, iterations=10^2)
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
            u_kln[k, l, 1:N_k[k]] = u_kl[k, l]
        end
    end

    dummy_kn = zeros(Int64, K, N_max)
    for k = 1:K
      dummy_kn[k, 1:N_k[k]] .= 1
    end
    idx = dummy_kn .> 0;

    # solve the MBAR equation
    log_wi_jn = zeros(Float64, K, N_max)
    log_term = zeros(Float64, K)
    FlogN = zeros(Float64, K)
    # f!(F, x) = mbar_equation!(F, x, K, N_max, N_k, u_kln, idx, log_wi_jn, log_term, FlogN)
    # x_init = zeros(Float64, K)
    # sol = nlsolve(f!, x_init, autodiff=:forward, ftol=ftol, iterations=iterations)
    # #sol = nlsolve(f!, x_init, method = :anderson, ftol=ftol, iterations=iterations)
    # f_k = sol.zero
    # f_k .- f_k[1]

    F = zeros(Float64, K)
    x_init = zeros(Float64, K)
    for n = 1:1000
        mbar_equation!(F, x_init, K, N_max, N_k, u_kln, idx, log_wi_jn, log_term, FlogN)
        x_init = x_init - F
        x_init = x_init .- x_init[1]
    end
    f_k = x_init

end
