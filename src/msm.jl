function msmsample(p)
  p_cum = cumsum(p) ./ sum(p)
  r = rand()
  for i = 1:length(p_cum)
    if r <= p_cum[i]
      return i
    end
  end
end

"""
msmgenerate(nframe, T, pi_i) -> states
generate a discrete-state trajectory from given transition matrix, and equilibrium probabilities of states.
Examples
≡≡≡≡≡≡≡≡≡≡
julia> T, pi_i = msmtransitionmatrix(C)
julia> states = msmgenerate(1000, T, pi_i)
"""
function msmgenerate(nframe::Int, T, pi_i)
  states = zeros(typeof(nframe), nframe)

  states[1] = msmsample(pi_i)
  for iframe = 2:nframe
    states[iframe] = msmsample(T[states[iframe-1], :])
  end
  return states
end

"""
msmgenerate(nframe, T, pi_i, emission) -> states, observations
generate a discrete-state trajectory and observations from given transition matrix, equilibrium probabilities of states, and emissions.
Examples
≡≡≡≡≡≡≡≡≡≡
julia> T, pi_i = msmtransitionmatrix(C)
julia> states, observations = msmgenerate(1000, T, pi_i, emission)
"""
function msmgenerate(nframe::Int, T, pi_i, emission)
  states = zeros(typeof(nframe), nframe)
  observations = zeros(typeof(nframe), nframe)

  states[1] = msmsample(pi_i)
  observations[1] = msmsample(emission[states[1], :])
  for iframe = 2:nframe
    states[iframe] = msmsample(T[states[iframe-1], :])
    observations[iframe] = msmsample(emission[states[iframe], :])
  end
  return states, observations
end

function msmforward(data_list, T, pi_i, emission)
    ndata = length(data_list)
    nstate = length(T[1, :])
    logL = zeros(Float64, ndata)
    alpha_list = []
    factor_list = []
    for idata = 1:ndata
        data = data_list[idata]
        nframe = length(data)
        alpha  = zeros(Float64, (nframe, nstate))
        factor = zeros(Float64, nframe)
        alpha[1, :] = pi_i.*emission[:, data[1]]
        factor[1] = sum(alpha[1, :])
        alpha[1, :] = alpha[1, :]./factor[1]
        for iframe = 2:nframe
            alpha[iframe, :] = sum(alpha[iframe-1, :] .* T, dims=1)' .* emission[:, data[iframe]]
            factor[iframe] = sum(alpha[iframe, :])
            alpha[iframe, :] = alpha[iframe, :]./factor[iframe]
        end
        logL[idata] = sum(log.(factor))
        push!(alpha_list, alpha)
        push!(factor_list, factor)
    end
    logL, alpha_list, factor_list
end

function msmbackward(data_list, factor_list, T, pi_i, emission)
    ndata = length(data_list)
    nstate = length(T[1, :])
    logL = zeros(Float64, ndata)
    beta_list = []
    for idata = 1:ndata
        data   = data_list[idata]
        factor = factor_list[idata]
        nframe = length(data)
        beta   = zeros(Float64, (nframe, nstate))
        beta[nframe, :] .= 1.0
        for iframe = (nframe-1):-1:1
            beta[iframe, :] = sum((T .* (emission[:, data[iframe+1]] .* beta[iframe+1, :])'), dims=2) ./ factor[iframe+1]
        end
        logL[idata] = sum(log.(factor))
        push!(beta_list, beta)
    end
    logL, beta_list
end

function msmbaumwelch(data_list, T0, pi_i0, emission0)
    ## setup
    TOLERANCE = 10.0^(-4)
    check_convergence = Inf64
    count_iteration = 0
    logL_old = 1.0
    #if not isinstance(data_list, list):
    #    data_list = [data_list]
    ndata = length(data_list)
    nobs = length(emission0[1, :])
    nstate = length(T0[1, :])
    T = similar(T0)
    emission = similar(emission0)
    pi_i = similar(pi_i0)
    while check_convergence > TOLERANCE
        ## E-step
        logL, alpha_list, factor_list = msmforward(data_list, T0, pi_i0, emission0)
        #print("1"); println(logL)
        logL2, beta_list = msmbackward(data_list, factor_list, T0, pi_i0, emission0)
        #print("2"); println(logL2)
        log_alpha_list = []
        for a in alpha_list
            push!(log_alpha_list, log.(a))
        end
        log_beta_list = []
        for b in beta_list
            push!(log_beta_list, log.(b))
        end
        log_T0 = log.(T0)
        log_emission0 = log.(emission0)
        ## M-step
        # pi
        # pi = np.zeros(nstate, dtype=np.float64)
        # log_gamma_list = []
        # for idata in range(ndata):
        #     log_gamma_list.append(log_alpha_list[idata] + log_beta_list[idata])
        #     pi = pi + np.exp(log_gamma_list[idata][0, :])
        # pi = pi/np.sum(pi)
        pi_i = pi_i0
        # emission
        # emission = np.zeros((nstate, nobs), dtype=np.float64)
        # for idata in range(ndata):
        #     data = data_list[idata]
        #     for istate in range(nstate):
        #         for iobs in range(nobs):
        #             id = (data == iobs)
        #             if np.any(id):
        #                 emission[istate, iobs] = emission[istate, iobs] + np.sum(np.exp(log_gamma_list[idata][id, istate]))
        # emission[np.isnan(emission)] = 0.0
        # emission = emission / np.sum(emission, axis=1)[:, None]
        emission = emission0
        # T
        T = zeros(Float64, (nstate, nstate))
        for idata = 1:ndata
          data = data_list[idata]
          nframe = length(data)
          for iframe = 2:nframe
            #log_xi = bsxfun(@plus, log_alpha{idata}(iframe-1, :)', log_beta{idata}(iframe, :));
            log_xi = log_alpha_list[idata][iframe-1, :] .+ log_beta_list[idata][iframe, :]'
            #T = T .+ exp(bsxfun(@plus, log_xi, log_emission0(:, data(iframe))') + log_T0)./factor{idata}(iframe);
            T = T .+ exp.((log_xi .+ log_emission0[:, data[iframe]]') .+ log_T0) ./ factor_list[idata][iframe]
          end
        end
        #T[np.isnan(T)] = 0.0
        T = T ./ sum(T, dims=2)
        ## Check convergence
        count_iteration += 1
        logL = sum(logL)
        check_convergence = abs(logL_old - logL)
        if mod(count_iteration, 100) == 0
            Printf.@printf("%d iteration LogLikelihood = %e  delta = %e  tolerance = %e\n" , count_iteration, logL, check_convergence, TOLERANCE)
        end
        logL_old = logL
        pi_i0 = pi_i
        emission0 = emission
        T0 = T
    end
    T, pi_i, emission
end

"""
msmviterbi(T, pi_i, emission, observation) -> states
estimate most probable hidden state sequence from observation
Examples
≡≡≡≡≡≡≡≡≡≡
julia> states, observations = msmgenerate(1000, T, pi_i, emission)
julia> states_estimated = msmviterbi(T, pi_i, emission, observation)
"""
function msmviterbi(observation, T, pi_i, emission)
    nframe = size(observation, 1)
    nstate = size(T, 1)
    P = zeros(eltype(T), nstate, nframe)
    I = zeros(eltype(T), nstate, nframe)
    state_estimated = zeros(eltype(observation), nframe)

    # initialization
    P[:, 1] .= log.(pi_i) .+ log.(emission[:, observation[1]])
    I[:, 1] .= zeros(eltype(T), nstate)

    # argmax forward
    Z = zeros(eltype(T), nstate, nstate)
    for t = 2:nframe
        Z .= P[:, t-1] .+ log.(T)
        I[:, t] .= getindex.(argmax(Z, dims=1), 1)[:]
        P[:, t] .= maximum(Z, dims=1)[:] .+ log.(emission[:, observation[t]])
    end

    # termination
    P_star = maximum(P[:, nframe])
    state_estimated[nframe] = argmax(P[:, nframe])
    #@show P

    # decoding
    for t = (nframe-1):-1:1
        state_estimated[t] = I[state_estimated[t+1], t+1]
    end

    return state_estimated
end

function msmviterbi_original(observation, T, pi_i, emission)
    nframe = size(observation, 1)
    nstate = size(T, 1)
    P = zeros(eltype(T), nstate, nframe)
    I = zeros(eltype(T), nstate, nframe)
    state_estimated = zeros(eltype(observation), nframe)

    # initialization
    P[:, 1] .= pi_i .* emission[:, observation[1]]
    I[:, 1] .= zeros(eltype(T), nstate)

    # argmax forward
    Z = zeros(eltype(T), nstate, nstate)
    for t = 2:nframe
        Z .= P[:, t-1] .* T
        I[:, t] .= getindex.(argmax(Z, dims=1), 1)[:]
        P[:, t] .= maximum(Z, dims=1)[:] .* emission[:, observation[t]]
    end

    # termination
    P_star = maximum(P[:, nframe])
    state_estimated[nframe] = argmax(P[:, nframe])
    #@show P

    # decoding
    for t = (nframe-1):-1:1
        state_estimated[t] = I[state_estimated[t+1], t+1]
    end

    return state_estimated
end

####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

mutable struct Edge
    to::Int64
    value::Float64
end

function msmgenerate_edge(nframe, state_num, E, pi_i, emission)
    states = zeros(typeof(nframe), nframe)
    observations = zeros(typeof(nframe), nframe)

    states[1] = msmsample(pi_i)
    observations[1] = msmsample(emission[states[1], :])
    for iframe = 2:nframe
        tmp = []
        for edge in E[states[iframe - 1]]
            push!(tmp, edge.value)
        end
        states[iframe] = E[states[iframe - 1]][msmsample(tmp)].to
        observations[iframe] = msmsample(emission[states[iframe], :])
    end
    return states, observations
end

function msmforward_edge(data_list, state_num, E, pi_i, emission)
    alpha_list = []
    alpha_sum_list = []
    for data in data_list
        alpha_sum = []
        data_size = size(data, 1)
        alpha = zeros(Float64, (data_size, state_num))
        alpha[1, :] = pi_i .* emission[:, data[1]]
        push!(alpha_sum, sum(alpha[1, :]))
        alpha[1, :] ./= alpha_sum[1]
        
        for frame in 2:data_size
            for state in 1:state_num
                for edge in E[state]
                    to = edge.to
                    trans_prob = edge.value
                    alpha[frame, to] += alpha[frame - 1, state] * trans_prob * emission[to, data[frame]]
                end
            end
            push!(alpha_sum, sum(alpha[frame, :]))
            alpha[frame, :] ./= alpha_sum[frame]
        end
        push!(alpha_list, alpha)
        push!(alpha_sum_list, alpha_sum)
    end
    
    return alpha_list, alpha_sum_list
end

function msmbackward_edge(data_list, state_num, E, pi_i, emission, alpha_sum_list)
    beta_list = []
    data_num = length(data_list)
    for (data, data_id) in zip(data_list, 1:data_num)
        data_size = size(data, 1)
        beta = zeros(Float64, (data_size, state_num))
        # 全ての状態が最終状態となることを仮定
        beta[data_size, :] .= 1
        
        for frame in (data_size - 1):-1:1
            for state in 1:state_num
                for edge in E[state]
                    to = edge.to
                    trans_prob = edge.value
                    beta[frame, state] += beta[frame + 1, to] * trans_prob * emission[to, data[frame + 1]]
                end
            end
            beta[frame, :] ./= alpha_sum_list[data_id][frame + 1]
        end
        
        push!(beta_list, beta)
    end
    
    return beta_list
end

function msmbaumwelch_edge(data_list, state_num, E, pi_i, emission)
    ## setup
    tolerance = 10.0^(-4)
    check_convergence = Inf64
    count_iteration = 0
    data_num = length(data_list)
    obs_num = length(emission[1, :])
    while check_convergence > tolerance
        alpha_list, alpha_sum_list = msmforward_edge(data_list, state_num, E, pi_i, emission)
        beta_list = msmbackward_edge(data_list, state_num, E, pi_i, emission, alpha_sum_list)
        
        E_nxt = deepcopy(E)
        pi_nxt = deepcopy(pi_i)
        emission_nxt = deepcopy(emission)
        
        diff_sum = 0
        for state in 1:state_num
            state_sum = 0
            for edge_id in 1:size(E_nxt[state], 1)
                to = E_nxt[state][edge_id].to
                E_nxt[state][edge_id].value = 0
                for (data, data_id) in zip(data_list, 1:data_num)
                    data_size = size(data, 1)
                    for frame in 1:(data_size - 1)
                        E_nxt[state][edge_id].value += (alpha_list[data_id][frame, state]
                                                    * E[state][edge_id].value
                                                    * emission[to, data[frame + 1]]
                                                    * beta_list[data_id][frame + 1, to]
                                                    / alpha_sum_list[data_id][frame + 1])
                    end
                end
                
                state_sum += E_nxt[state][edge_id].value
            end
            
            for edge_id in 1:size(E_nxt[state], 1)
                # println("$(state) $(edge_id) $(E_nxt[state][edge_id].value) $(state_sum)")
                E_nxt[state][edge_id].value /= state_sum
                # println("$(state) $(edge_id) $(E_nxt[state][edge_id].value) $(state_sum)")
                diff_sum += abs(E[state][edge_id].value - E_nxt[state][edge_id].value)
            end
            
        end
        
        E = E_nxt
        pi_i = pi_nxt
        emission = emission_nxt
        check_convergence = diff_sum
        count_iteration += 1
        
        if count_iteration % 20 == 0
            println("count: ", count_iteration, " diff_sum: ", diff_sum)
        end
    end
    
    return E, pi_i, emission
end

function msmviterbi_edge(observation, state_num, E, pi_i, emission)
    # initialization
    frame_num = size(observation, 1)
    prob = zeros(Float64, frame_num, state_num)
    fill!(prob, -1e8)
    prob[1, :] .= log.(pi_i) .+ log.(emission[:, observation[1]])
    from = zeros(Int64, frame_num, state_num)

    for frame in 1:(frame_num - 1)
        for state in 1:state_num
            for edge in E[state]
                tmp = prob[frame, state] + log(edge.value) + log(emission[edge.to, observation[frame + 1]])
                if prob[frame + 1, edge.to] < tmp
                    prob[frame + 1, edge.to] = tmp
                    from[frame + 1, edge.to] = state
                end
            end
        end
    end

    # termination
    state_estimated = zeros(Int64, frame_num)
    state_estimated[frame_num] = argmax(prob[frame_num, :])

    # decoding
    for frame = (frame_num - 1):-1:1
        state_estimated[frame] = from[frame + 1, state_estimated[frame + 1]]
    end

    return state_estimated
end