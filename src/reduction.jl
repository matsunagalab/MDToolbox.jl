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

function pca(X)
    nframe = size(X, 1)
    X_centerized = X .- mean(X, dims=1)
    e = eigen(X_centerized' * X_centerized ./ nframe)
    lambda = e.values[end:-1:1]
    W = e.vectors[:, end:-1:1]
    P = X_centerized * W
    return (p=P, mode=W, lambda=lambda)
end

function tica()
end

