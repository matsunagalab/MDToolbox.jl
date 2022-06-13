"""
    clusterkcenters(X::AbstractMatrix, kcluster::Int; nReplicates::Int=10) -> F

Perform clustering with K-center algorithm. Input data `X` should belong to AbstractMatrix type 
and its columns corresponds to variables, rows are frames. Also, the number of clusters `kcluster` should be specified. 

Returns a NamedTuple object `F` which contains cluster index for each sample in `F.indexOfCluster`,
the coordinates of cluster centers in `F.center`, the indices of cluster centers in `F.indexOfCenter`, 
and distances of samples from the nearest centers in `F.distanceFromCenter`.

#  Example
```julia-repl
julia> using MDToolbox, Plots
julia> X = rand(1000, 2)
julia> F = clusterkcenters(X, 3)
julia> scatter(X[:, 1], X[:, 2], c=F.indexOfCluster)
```
# References
```
This function uses the method described in
[1] S. Dasgupta and P. M. Long, J. Comput. Syst. Sci. 70, 555 (2005).
[2] J. Sun, Y. Yao, X. Huang, V. Pande, G. Carlsson, and L. J. Guibas, Learning 24, 2 (2009).
```
"""
function clusterkcenters(t::AbstractMatrix, kcluster::Int; nReplicates::Int=10)
    nframe = size(t, 1)
    #ta2, com = decenter(ta)
    dim = size(t, 2)

    indexOfCenter = zeros(Int64, kcluster)
    indexOfCluster = zeros(Int64, nframe)
    distanceFromCenter = zeros(Float64, nframe, 1)

    indexOfCenter_out = similar(indexOfCenter)
    indexOfCluster_out = similar(indexOfCluster)
    distanceFromCenter_out = similar(distanceFromCenter)
    center_out = similar(t, kcluster, dim)
    distanceMax_out = Inf64

    for ireplica = 1:nReplicates
        # create the center of the 1st cluster
        indexOfCenter[1] = rand(1:nframe)
        # first, all points belong to the 1st cluster
        indexOfCluster .= ones(Int64, nframe)
        # distances of the points from the 1st center
        ref = t[indexOfCenter[1]:indexOfCenter[1], :]
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
          indexOfCenter_out = indexOfCenter
          indexOfCluster_out = indexOfCluster
          center_out = t[indexOfCenter, :]
          distanceFromCenter_out .= distanceFromCenter
        end
        Printf.@printf("%d iteration  distance_max = %f  kcluster = %d\n", ireplica, sqrt(distanceMax_out), kcluster)
    end
    distanceFromCenter_out .= sqrt.(distanceFromCenter_out)

    return (indexOfCluster=indexOfCluster_out, center=center_out, distanceFromCenter=distanceFromCenter, indexOfCenter=indexOfCenter_out)
end

"""
    clusterkcenters(ta::AbstractMatrix, kcluster::Int; nReplicates::Int=10) -> F

Perform clustering with K-center algorithm for a TrjArray variable `ta`.
"""
function clusterkcenters(ta::TrjArray, kcluster::Int; nReplicates::Int=10)
    X = ta.xyz
    F = clusterkcenters(X, kcluster, nReplicates=nReplicates)
    center_ta = ta[F.indexOfCenter, :]
    return (indexOfCluster=F.indexOfCluster, center=center_ta, distanceFromCenter=F.distanceFromCenter, indexOfCenter=F.indexOfCenter)
end

"""
    clusterkmeans(X::AbstractMatrix, kcluster::Int; nReplicates::Int=10) -> F

Perform clustering with K-means algorithm. Input data `X` should belong to AbstractMatrix type 
and its columns corresponds to variables, rows are frames. Also, the number of clusters `kcluster` should be specified. 

Returns a NamedTuple object `F` which contains cluster index for each sample in `F.indexOfCluster`,
the coordinates of centroid in `F.center`, and distances of samples from the nearest centroid in `F.distanceFromCenter`.

#  Example
```julia-repl
julia> using MDToolbox, Plots
julia> X = rand(1000, 2)
julia> F = clusterkmeans(X, 3)
julia> scatter(X[:, 1], X[:, 2], c=F.indexOfCluster)
```
"""
function clusterkmeans(x::AbstractMatrix, kcluster::Int; nReplicates::Int=10)
    N = size(x, 1)
    P = size(x, 2)
    K = kcluster

    # return variables
    class_out = rand(1:K, N)
    centers_out = zeros(Float64, K, P)
    d_min_out = zeros(Float64, N)
    score_out = Inf64

    # initialization
    for iter = 1:nReplicates
        class = rand(1:K, N)
        centers = zeros(Float64, K, P)
        d_min = zeros(Float64, N)
        for k = 1:K
            centers[k:k, :] .= mean(x[class .== k, :], dims=1)
        end

        class_old = deepcopy(class)
        while true
            d_min .= Inf64
            class_old .= class

            # update clustering
            for k = 1:K
              d = sqrt.(sum((x .- centers[k:k, :]).^2, dims=2))
              idx = d .< d_min
              class[idx[:]] .= k
              d_min[idx[:]] .= d[idx]
            end

            # update centroid
            for k = 1:K
                centers[k:k, :] .= mean(x[class .== k, :], dims=1)
            end

            # check convergence
            hammingDistance = sum(class_old .!= class) / N
            if hammingDistance < 10^(-6)
                break
            end
        end
        score = sum(d_min.^2)

        if score < score_out
            class_out .= class
            centers_out .= centers_out
            d_min_out .= d_min
            score_out = score
        end
    end

    return (indexOfCluster=class_out, center=centers_out, distanceFromCenter=d_min_out)
end

"""
    clustercutoff(X::AbstractMatrix, rcut=10.0) -> F

Perform clustering using the given cutoff radius.

#  Example
```julia-repl
julia> using MDToolbox, Plots
julia> X = rand(1000, 2)
julia> F = clustercutoff(X, 0.5)
julia> scatter(X[:, 1], X[:, 2], c=F.indexOfCluster)
```
# References
```
This function uses the method described in
[1] X. Daura, K. Gademann, B. Jaun, D. Seebach, W. F. van Gunsteren, and A. E. Mark, Angewandte 38, 236–240 (1999).
```
"""
function clustercutoff(t::AbstractMatrix, rcut=10.0)
    nframe = size(t, 1)
    rcut2 = rcut^2

    indexOfCluster = zeros(Int64, nframe)
    id = indexOfCluster .== 0

    center = similar(t)
    k = 1

    while any(id)
        x = view(t, id, :)
        index = view(indexOfCluster, id)

        n_max = 0
        iframe_max = 0
        for iframe = 1:size(x, 1)
            dist = sum((x[iframe:iframe, :] .- x).^2, dims=2)
            n = sum(dist .< rcut2)
            if n_max < n
                n_max = n
                iframe_max = iframe
            end
        end

        center[k:k, :] .= x[iframe_max:iframe_max, :]
        dist = sum((center[k:k, :] .- x).^2, dims=2)
        index[dist[:] .< rcut2] .= k
        k += 1

        id = indexOfCluster .== 0
    end

    return (indexOfCluster=indexOfCluster, center=center[1:(k-1), :])
end


"""
    compute_cov(ta::AbstractMatrix, lagtime::Int=0) -> cov

Compute a variance-covariance or time-lagged covariance matrix from input data `X`
Input data `X` should belong to AbstractMatrix type and its columns corresponds to variables, rows are frames.
Optional input is the lagtime=`lagtime` for the calculation of the covariance matrix (default is `lagtime=0`). 

#  Example
```julia-repl
julia> X = rand(1000, 100)
julia> cov = compute_cov(X)
```
"""
function compute_cov(X::AbstractMatrix; lagtime::Int=0)
    nframe = size(X, 1)
    X_centerized = X .- mean(X, dims=1)
    cov = X_centerized[1:(end-lagtime), :]' * X_centerized[(1+lagtime):end, :] ./ (nframe - 1)
end

"""
    rsvd(X::AbstractMatrix; k::Number=10) -> F

Perform the randomized SVD for input data `X`.
Input data `X` should belong to AbstractMatrix type and its columns corresponds to variables, rows are frames. 
Users can specify the dimension `k` of subspace onto which the data is randomly projected (by default `k=10`). 

Returns a NamedTuple object `F` whose members are same as usual SVD, `F.V`, `F.S`, `F.U`. 

#  Example
```julia-repl
julia> X = randn(rand(1000, 10))
julia> F = rsvd(X)
```

# References
```
Halko, Nathan, Per-Gunnar Martinsson, and Joel A. Tropp. 
"Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions." 
SIAM review 53.2 (2011): 217-288.
```
"""
function rsvd(A::AbstractMatrix, k::Number=10)
    T = eltype(A)
    m, n = size(A)
    l = k + 2
    @assert l < n
    G = randn(T, n, l)
    #G .= G ./ sqrt.(sum(G.^2, dims=1))
    H = A * G
    @assert m > l
    Q = Matrix(qr!(H).Q)
    @assert m == size(Q, 1)
    @assert l == size(Q, 2)
    T = A' * Q
    V, σ, W = svd(T)
    U = Q * W
    λ = σ .* σ ./ m
    return (V=V[:,1:k], S=λ[1:k], U=U[:,1:k])
end

"""
    pca(X::AbstractMatrix; k=dimension) -> F

Perform principal component analysis (PCA). PCA captures degrees of freedom which have the largest variances.
Input data `X` should belong to AbstractMatrix type and its columns corresponds to variables, rows are frames. 

Returns a NamedTuple object `F` which contains the principal components in `F.projection`, 
the prncipal modes in the columns of the matrix `F.mode`, and the variances of principal components in `F.variance`. 

If `k=dimension` is specified, the randomized SVD is used for reducing memory. 
This algorithm first project the data into a randomly selected k+2-dimensional space, 
then PCA is performed in the projected data. See the references for details. 
Note that if the dimension of `X` is larger than 5000, the randomized SVD is forcibly used with `k=1000`. 

#  Example
```julia-repl
julia> using MDToolbox, Plots, Statistics
julia> X = cumsum(rand(1000, 10))
julia> F = pca(X)
julia> plot(F.projection[:, 1], F.projection[:, 2])
```

# References
```
Halko, N., Martinsson, P.-G., Shkolnisky, Y. & Tygert, M. 
An Algorithm for the Principal Component Analysis of Large Data Sets. 
SIAM J. Sci. Comput. 33, 2580–2594 (2011).
```
"""
function pca(X::AbstractMatrix; k=nothing)
    is_randomized_pca = false
    if !isnothing(k)
        is_randomized_pca = true
        @printf "Warning: a randomized SVD approximation is used for PCA with reduced dimension k = %d\n" k
    end
    dim = size(X, 2)
    if !is_randomized_pca & (dim > 5000)
        k = 1000
        is_randomized_pca = true
        @printf "Warning: the dimension of input data is too large (dim > 5000)\n"
        @printf "Warning: randomized SVD approximation is used for PCA with reduced dimension k = %d\n" k
    end
    if is_randomized_pca
        nframe = size(X, 1)
        X_centerized = X .- mean(X, dims=1)
        F = rsvd(X_centerized, k)
        variance = F.S.^2 ./ (nframe - 1)
        mode = F.V
        projection = X_centerized * mode
    else
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
    end

    return (projection=projection, mode=mode, variance=variance)
end

"""
    tica(X::AbstractMatrix, lagtime::Int=1) -> F

This routine performs time-structure based Independent Component Analysis (tICA).
tICA captures degrees of freedom which are most important in the sense that their motions are very slow.
`X` belongs to AbstractMatrix type and its columns corresponds to variables, rows are frames. 
User should specify the lagtime for the calculation of the time-lagged covariance matrix. 
Returns a NamedTuple object `F` which contains independent components in `F.projection`, 
the independent modes in the columns of the matrix `F.mode`, and the eigenvalues in `F.variance`. 

If the dimension of `X` is larger than 5000, the randomized SVD approximation is forcibly used with `k=1000` to reduce memory. 

#  Example
```julia-repl
julia> using MDToolbox, Plots, Statistics
julia> X = cumsum(rand(1000, 10))
julia> F = tica(X, 10)
julia> plot(F.projection[:, 1], F.projection[:, 2])
```

# References
```
Naritomi, Y. & Fuchigami, S. 
Slow dynamics in protein fluctuations revealed by time-structure based independent component analysis: 
The case of domain motions. 
The Journal of Chemical Physics 134, 065101 (2011).
```
"""
function tica(X::AbstractMatrix, lagtime::Int=1)
    # standard and time-lagged covariance matrices
    covar0 = compute_cov(X, lagtime=0)
    covar  = compute_cov(X, lagtime=lagtime)

    # symmetrize the time-lagged covariance matrix
    covar .= 0.5 .* (covar .+ covar')

    # calc pseudo-inverse of covar0
    #covar0_inv = pinv(covar0)

    # remove degeneracy for solving the generalized eigenvalue problem
    #F = eigen(covar0, sortby = x -> -x)
    F = pca(X)
    pvariance = F.variance
    pmode = F.mode
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
