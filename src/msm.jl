"""
    msmplot(T; pi_i=nothing, x=nothing, y=nothing, filename=nothing, 
            edgewidth_scale=10.0, arrow_scale=0.1, nodesize=0.5, fontsize=10, names=[], dpi=100)

Visualize the graphical structure of the given Markov state model parameters. 
`T` is a transition probability matrix whose elements Tij represents the probablity of transition from state i to j. 
`T` should satisfy the detailed balance condition. `pi_i` is a vector whose elements are equilibrium probabilities of states. 
`x` and `y` are X and Y coordinates of states, respectively. 

# Examples
```julia-repl
julia> nstate = 5
julia> T, pi_i = msmtransitionmatrix(rand(nstate, nstate))
julia> x = rand(nstate); y = rand(nstate)
julia> msmplot(T, pi_i=pi_i, x=x, y=y)
```
"""
function msmplot(T; pi_i=nothing, x=nothing, y=nothing, filename=nothing, 
                 edgewidth_scale=3.0, arrow_scale=0.0001, nodesize=0.5, fontsize=10, names=[], dpi=100)
  n = size(T, 1)
  g = MetaDiGraph(n, 0.0)
  U = typeof(T[1, 1])
  
  for i = 1:n
    for j = 1:n
      if !iszero(T[i, j])
        add_edge!(g, i, j)
      end
     end
  end
  
  for i = 1:n
    for j = 1:n
      if !iszero(T[i, j])
        set_prop!(g, Edge(i, j), :weight, T[i, j])
      end
    end
  end

  if isnothing(pi_i)
    pi_i = ones(1, n)./n
    for i = 1:(n*100)
      pi_i .= pi_i * T
    end
    pi_i = pi_i./sum(pi_i)
    pi_i = pi_i[:]
  end

  if !isnothing(x)
    std_x = std(x)
    std_y = std(y)
    max_xy = maximum([maximum(x),maximum(y)])
    xx = (x[:].-mean(x))./std_x .- std_x
    yy = (y[:].-mean(y))./std_y .- std_y
  else
    xx = nothing
    yy = nothing
  end

  gr()
#  pyplot()
  p = graphplot(g, edge_width=(s,d,w) -> get_prop(g, Edge(s,d), :weight)*edgewidth_scale, 
            #arrow=arrow(:closed, :head, 0.1, 0.1),
            arrow=arrow(:closed, :head, arrow_scale, arrow_scale),
            #markersize = markersize,
            node_weights = pi_i[:],
            #node_weights = 1:n,  
            #markercolor = range(colorant"yellow", stop=colorant"red", length=n),
            markercolor = :white,
            names = names,
            nodesize = nodesize, 
            fontsize = fontsize,
            fontcolor = :white, 
            nodeshape = :circle,
            linecolor = :darkgrey,
            #shorten=0.95,
            x = xx,
            y = yy,
            curves = true, 
            curvature_scalar = 0.8,
            dpi = dpi
            )

  if !isnothing(filename)
    savefig(p, filename)
  end

  return p
end

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
    msmgenerate(nframe::Int, T, pi_i)
                states = zeros(typeof(nframe), nframe)

Randomly samples a state trajectory from the the given transition matrix, and equilibrium probabilities of states. 

# Examples
```julia-repl
julia> T, pi_i = msmtransitionmatrix(C)
julia> states = msmgenerate(1000, T, pi_i)
```
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

Randomly samples a state trajectory and observations from the given transition matrix, equilibrium probabilities of states, and emissions.

# Examples
```julia-repl
julia> T, pi_i = msmtransitionmatrix(C)
julia> states, observations = msmgenerate(1000, T, pi_i, emission)
```
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

"""
    msmcountmatrix(indexOfCluster; tau=1) -> C::Matrix

Transition count matrix from a sinlgle binned trajectory or a set of binned trajectories. 
`indexOfCluster` is a vector or a set of vectors containing binned trajectorie(s). 
Returns count matrix for transitions from state i to state j in a time step of tau. 

# Examples
```julia-repl
julia> ta = mdload("ak.dcd")
julia> X = compute_distancemap(ta["atomname CA"])
julia> F = clusterkcenters(X)
julia> c = msmcountmatrix(F.indexOfCluster, tau=10)
```
"""
function msmcountmatrix(indexOfCluster; tau=1)
  nstate = maximum(indexOfCluster)
  nframe = size(indexOfCluster, 1)

  if typeof(indexOfCluster) <: AbstractArray
    indexOfCluster2 = [indexOfCluster]
  else
    indexOfCluster2 = indexOfCluster
  end

  ntraj = length(indexOfCluster2)
  U = typeof(indexOfCluster2[1][1])
  c = spzeros(U, nstate, nstate)
  for itraj = 1:ntraj
    nframe = length(indexOfCluster2[itraj])

    index_from = 1:(nframe-tau)
    index_to   = (1+tau):nframe
    indexOfCluster_from = indexOfCluster2[itraj][index_from]
    indexOfCluster_to   = indexOfCluster2[itraj][index_to]

    s = ones(U, length(indexOfCluster_from))
    c_1 = sparse(indexOfCluster_from, indexOfCluster_to, s, nstate, nstate)
    c .= c .+ c_1
  end

  return Array(c)
end

"""
    msmtransitionmatrix(C; TOLERANCE=10^(-4), verbose=true) -> T::Matrix, p::Vector

Estimate the transition probability matrix from count matrix `C`. 
Detailed balance is implicitly imposed in the estimation. 

Returns the transition probability matrix and the equilibrium probabilities of states. 

# Examples
```julia-repl
julia> ta = mdload("ak.dcd")
julia> X = compute_distancemap(ta["atomname CA"])
julia> F = clusterkcenters(X)
julia> C = msmcountmatrix(F.indexOfCluster, tau=10)
julia> T, p = msmtransitionmatrix(C)
```
# References
```
This routines uses the reversible maximum likelihood estimator described in 
K. A. Beauchamp, G. R. Bowman, T. J. Lane, L. Maibaum, I. S. Haque, and V. S. Pande, 
MSMBuilder2: Modeling Conformational Dynamics on the Picosecond to Millisecond Scale, 
J. Chem. Theory Comput. 7, 3412 (2011).
```
"""
function msmtransitionmatrix(C; TOLERANCE=10^(-4), verbose=true)
  c = Matrix{Float64}(C)
  nstate = size(c, 1)

  c_sym = c + c'
  x = c_sym

  c_rs = sum(c, dims=2)
  x_rs = sum(x, dims=2)

  logL_old = 2.0*TOLERANCE
  logL = 0.0
  count_iteration = 0

  x = zeros(Float64, nstate, nstate)
  x_new = zeros(Float64, nstate, nstate)

  while abs(logL_old - logL) > TOLERANCE
    count_iteration = count_iteration + 1
    logL_old = logL
  
    # fixed-point method
    for i = 1:nstate
      for j = 1:nstate
        denom = (c_rs[i]/x_rs[i]) + (c_rs[j]/x_rs[j])
        x_new[i, j] = (c[i, j] + c[j, i]) / denom
      end
    end
    
    # update
    x_rs .= sum(x_new, dims=2)
    x .= x_new
    logL = 0.0
    for i = 1:nstate
      for j = 1:nstate
        if !iszero(x[i, j]) & !iszero(x_rs[i])
          logL = logL + c[i, j] * log(x[i, j] / x_rs[i]);
        end
      end
    end
    
    if verbose & (mod(count_iteration, 10) == 0)
      Printf.@printf("%d iteration  LogLikelihood = %8.5e  delta = %8.5e  tolerance = %8.5e\n", count_iteration, logL, abs(logL_old-logL), TOLERANCE)
    end
  end
  
  pi_i = x_rs./sum(x_rs)
  pi_i = pi_i[:]
  t = x ./ x_rs

  return t, pi_i
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

function msmforward_missing(data_list, T, pi_i, emission)
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
      if data[1] === missing
        alpha[1, :] = pi_i
      else
        alpha[1, :] = pi_i.*emission[:, data[1]]
      end
      factor[1] = sum(alpha[1, :])
      alpha[1, :] = alpha[1, :]./factor[1]
      for iframe = 2:nframe
          if data[iframe] === missing
            alpha[iframe, :] = sum(alpha[iframe-1, :] .* T, dims=1)'
          else
            alpha[iframe, :] = sum(alpha[iframe-1, :] .* T, dims=1)' .* emission[:, data[iframe]]
          end
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

function msmbackward_missing(data_list, factor_list, T, pi_i, emission)
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
            if data[iframe+1] === missing
                beta[iframe, :] = sum((T .* beta[iframe+1, :]'), dims=2) ./ factor[iframe+1]
            else
                beta[iframe, :] = sum((T .* (emission[:, data[iframe+1]] .* beta[iframe+1, :])'), dims=2) ./ factor[iframe+1]
            end
        end
        logL[idata] = sum(log.(factor))
        push!(beta_list, beta)
    end
    logL, beta_list
end

"""
    msmbaumwelch(observations, T_init, p_init, emission_init; TOLERANCE=10.0^(-4), MAXITERATION=Inf64) -> T::Matrix, p::Vector, emission::Matrix

Baum-Welch algorithm estimates the most probable transition probabilities from the given observation data. 
In this function, detailed balance is implicitly imposed in the estimation, 
so the equilibrium probabilities can be determined from the estimated transition probabilities. 
Also, unlike the original Baum-Welch algorithm, the emission is NOT estimated in this function, 
because the emission probabilites are usually known a priori in cases of molecular experiments. 
`observations` is a set of observed vectors. `T_init`, `p_init` are initial transition probabilities and equilibrium probabilities, respectively. 
`emission_init` is a matrix whose rows correspond to states, and columns correspond to observations. 

Returns the transition probability matrix, the equilibrium probabilities of states, 
and the emission probabilities (though the emission does not change). 

# Examples
```julia-repl
julia> ta = mdload("ak.dcd")
julia> X = compute_distancemap(ta["atomname CA"])
julia> F = clusterkcenters(X)
julia> C = msmcountmatrix(F.indexOfCluster, tau=10)
julia> T, p = msmtransitionmatrix(C)
```
# References
```
The algorithm of this routines is based on the descriptions in PRML book by C. Bishop. 
```
"""
function msmbaumwelch(data_list, T0, pi_i0, emission0; TOLERANCE=10.0^(-4), MAXITERATION=Inf64)
    ## setup
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
    while (check_convergence > TOLERANCE) & (count_iteration <= MAXITERATION)
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
        # pi_i = pi_i0
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

        ## reversible T
        T, pi_i = msmtransitionmatrix(T, TOLERANCE=10^(-8), verbose=false)

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

function msmbaumwelch_missing(data_list, T0, pi_i0, emission0; TOLERANCE=10.0^(-4), MAXITERATION=Inf64)
  ## setup
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
  while (check_convergence > TOLERANCE) & (count_iteration <= MAXITERATION)
      ## E-step
      logL, alpha_list, factor_list = msmforward_missing(data_list, T0, pi_i0, emission0)
      #print("1"); println(logL)
      logL2, beta_list = msmbackward_missing(data_list, factor_list, T0, pi_i0, emission0)
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
      # pi_i = pi_i0
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
          log_xi = log_alpha_list[idata][iframe-1, :] .+ log_beta_list[idata][iframe, :]'
          if data[iframe] === missing
            T = T .+ exp.(log_xi .+ log_T0) ./ factor_list[idata][iframe]
          else
            T = T .+ exp.((log_xi .+ log_emission0[:, data[iframe]]') .+ log_T0) ./ factor_list[idata][iframe]
          end
        end
      end
      #T[np.isnan(T)] = 0.0
      T = T ./ sum(T, dims=2)

      ## reversible T
      T, pi_i = msmtransitionmatrix(T, TOLERANCE=10^(-8), verbose=false)

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
    msmviterbi(observation, T, p, emission) -> states::Vector

Viterbi algorithm estimates the most probable hidden state sequence from the observation data.
`observations` is a set of observed vectors. `T`, `p` are the transition probabilities and equilibrium probabilities, respectively. 
`emission` is a matrix whose rows correspond to states, and columns correspond to observations. 

Returns the transition probability matrix, the equilibrium probabilities of states, 
and the emission probabilities (though the emission does not change). 

# Examples
```julia-repl
julia> nframe = 1000
julia> states, observations = msmgenerate(nframe, T, pi_i, emission)
julia> states_estimated = msmviterbi(T, pi_i, emission, observation)
```
# References
```
The algorithm of this routines is based on the descriptions in PRML book by C. Bishop. 
```
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

function msmimpliedtime(indexOfCluster, tau)
  n = length(tau)
  k = min(10, n)
  nstate = Int(maximum(indexOfCluster))
  #implied_timescale = zeros(Float64, n, k)
  implied_timescale = Matrix{Union{Float64, Missing}}(undef, n, nstate-1)
  for i in 1:n
    t = tau[i]
    #@printf "calculating tau = %d\n" t
    C = msmcountmatrix(indexOfCluster, tau=t)
    # TODO fixme; tarjan algorithm to remove non-ergodic graphs
    T, p = msmtransitionmatrix(C, TOLERANCE=10^(-8), verbose=false)
    F = eigen(T, sortby = x -> -real(x))
    #@show F.values
    for istate = 2:nstate
      r = real(F.values[istate])
      if r >  0.0 
        implied_timescale[i, istate-1] = - t / log(r)
      else
        implied_timescale[i, istate-1] = missing
      end
    end
  end
  return implied_timescale
end


