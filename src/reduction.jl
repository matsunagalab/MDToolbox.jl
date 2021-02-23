"""
    clusterkcenters(t, kcluster; nReplicates=10)

Clustering with K-center algorithm

#  Example
```julia> scatter(t[:, 1], t[:, 2], c = r.indexOfCluster)
julia> t = rand(1000, 2)
julia> r = clusterkcenters(t, 3)
julia> scatter(t[:, 1], t[:, 2], c = r.indexOfCluster)
```
# References
```This function uses the method described in
[1] S. Dasgupta and P. M. Long, J. Comput. Syst. Sci. 70, 555 (2005).
[2] J. Sun, Y. Yao, X. Huang, V. Pande, G. Carlsson, and L. J. Guibas, Learning 24, 2 (2009).
```
"""
function clusterkcenters(t, kcluster; nReplicates=10)
    nframe = size(t, 1)
    #ta2, com = decenter(ta)
    dim = size(t, 2)

    indexOfCenter = zeros(Int64, kcluster)
    indexOfCluster = zeros(Int64, nframe)
    distanceFromCenter = zeros(Float64, nframe, 1)

    indexOfCluster_out = similar(indexOfCluster)
    distanceFromCenter_out = similar(distanceFromCenter)
    center_out = similar(t, kcluster, dim)
    distanceMax_out = Inf64

    for ireplica = 1:nReplicates
        # create the center of the 1st cluster
        indexOfCenter[1] = rand(1:nframe)
        # at first, all points belong to the 1st cluster
        indexOfCluster .= ones(Int64, nframe)
        # distance of the points from the 1st center
        ref = t[indexOfCluster[1]:indexOfCluster[1], :]
        distanceFromCenter .= sum((ref .- t).^2, dims=2)

        for i = 2:kcluster
            indexOfCenter[i] = argmax(distanceFromCenter[:])
            ref = t[indexOfCenter[i]:indexOfCenter[i], :]
            dist = sum((ref .- t).^2, dims=2)
            index = dist[:] .< distanceFromCenter[:]
            if any(index)
                # updated if the dist to a new cluster is smaller than the previous one
                distanceFromCenter[index] .= dist[index]
                indexOfCluster[index] .= i
            end
        end

        distanceMax = maximum(distanceFromCenter)
        if (ireplica == 1) | (distanceMax < distanceMax_out)
          distanceMax_out = distanceMax
          indexOfCluster_out = indexOfCluster
          center_out = t[indexOfCenter, :]
          distanceFromCenter_out .= distanceFromCenter
        end
        Printf.@printf("%d iteration  distance_max = %f  kcluster = %d\n", ireplica, sqrt(distanceMax_out), kcluster)
    end
    distanceFromCenter_out .= sqrt.(distanceFromCenter_out)

    return (indexOfCluster=indexOfCluster_out, center=center_out, distanceFromCenter=distanceFromCenter)
end

function compute_cov(X::AbstractMatrix; lagtime=0)
    nframe = size(X, 1)
    X_centerized = X .- mean(X, dims=1)
    cov = X_centerized[1:(end-lagtime), :]' * X_centerized[(1+lagtime):end, :] ./ (nframe - 1)
end

function pca(X)
    # eigendecomposition of covariance matrix
    covar = compute_cov(X)
    F = eigen(covar, sortby = x -> -x)
    #lambda = F.values[end:-1:1]
    #W = F.vectors[:, end:-1:1]
    variance = F.values
    mode = F.vectors

    # projection
    X_centerized = X .- mean(X, dims=1)
    projection = X_centerized * mode

    return (projection=projection, mode=mode, variance=variance)
end

function tica(X::AbstractMatrix, lagtime=1)
    # standard and time-lagged covariance matrices
    covar0 = compute_cov(X, lagtime=0)
    covar  = compute_cov(X, lagtime=lagtime)

    # symmetrize the time-lagged covariance matrix
    covar .= 0.5 .* (covar .+ covar')

    # calc pseudo-inverse of covar0
    #covar0_inv = pinv(covar0)

    # remove degeneracy for solving the generalized eigenvalue problem
    F = eigen(covar0, sortby = x -> -x)
    pvariance = F.values
    pmode = F.vectors
    index = pvariance .> 10^(-6)
    pmode = pmode[:, index]
    covar0 = pmode'*covar0*pmode
    covar = pmode'*covar*pmode

    # solve the generalized eigenvalue problem
    F = eigen(covar, covar0, sortby = x -> -x)
    variance = F.values
    mode = F.vectors

    # projection
    mode = pmode * mode
    X_centerized = X .- mean(X, dims=1)
    projection = X_centerized * mode

    # normalize mode vectors
    fac = sqrt.(sum(mode.^2, dims=1))
    mode .= mode ./ fac

    return (projection=projection, mode=mode, variance=variance)
end
