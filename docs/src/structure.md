# Structural computations

```@docs
centerofmass(ta::TrjArray{T, U}; isweight::Bool=true, index::AbstractVector=Vector{U}(undef, 0)) where {T, U}
```

```@docs
decenter(ta::TrjArray{T, U}; isweight::Bool=true, index::AbstractVector=Vector{Int64}(undef, 0)) where {T, U}
```

```@docs
decenter!(ta::TrjArray{T, U}; isweight::Bool=true, index::AbstractVector=Vector{Int64}(undef, 0)) where {T, U}
```

```@docs
orient!(ta::TrjArray{T, U}) where {T, U}
```

```@docs
rotate(ta::TrjArray{T, U}, quater::AbstractVector{T}) where {T, U}
```

```@docs
rotate!(ta::TrjArray{T, U}, quater::AbstractVector{T}) where {T, U}
```

```@docs
rotate(ta_single::TrjArray{T, U}, quater::AbstractMatrix{T}) where {T, U}
```

```@docs
rotate_with_matrix(ta::TrjArray{T, U}, R::AbstractMatrix{T}) where {T, U}
```

```@docs
superimpose(ref::TrjArray{T, U}, ta::TrjArray{T, U}; isweight::Bool=true, index::AbstractVector=Vector{U}(undef, 0), isdecenter::Bool=false) where {T, U}
```

```@docs
compute_rmsd(ref::TrjArray{T, U}, ta::TrjArray{T, U}; isweight::Bool=true, index::AbstractVector=Vector{U}(undef, 0)) where {T, U}
```

```@docs
meanstructure(ta::TrjArray{T, U}; isweight::Bool=true, index::Vector{U}=Vector{U}(undef, 0)) where {T, U}
```

```@docs
compute_rmsf(ta::TrjArray{T, U}; isweight::Bool=true) where {T, U}
```

```@docs
compute_distance(ta::TrjArray{T, U}, index=[1 2]::Matrix{U}) where {T, U}
```

```@docs
compute_distance(ta1::TrjArray{T, U}, ta2::TrjArray{T, U}, index=[1 1]::Matrix{U}) where {T, U}
```

```@docs
compute_distancemap(ta::TrjArray{T, U}; kneighbor=3) where {T, U}
```

```@docs
compute_contactmap(ta::TrjArray{T, U}; rcut=8.5, kneighbor=3) where {T, U}
```

```@docs
compute_angle(ta1::TrjArray{T, U}, ta2::TrjArray{T, U}, ta3::TrjArray{T, U}) where {T, U}
```

```@docs
compute_dihedral(ta1::TrjArray{T, U}, ta2::TrjArray{T, U}, ta3::TrjArray{T, U}, ta4::TrjArray{T, U}) where {T, U}
```

```@docs
compute_dihedral(ta::TrjArray{T, U}, array_index) where {T, U}
```

```@docs
compute_qscore(native::TrjArray{T, U}, ta::TrjArray{T, U}) where {T, U}
```

```@docs
compute_drms(native1::TrjArray{T, U}, native2::TrjArray{T, U}, ta1::TrjArray{T, U}, ta2::TrjArray{T, U}) where {T, U}
```

```@docs
compute_pairlist(ta::TrjArray{T, U}, rcut::T; iframe=1::Int) where {T, U}
```

```@docs
compute_pairlist_bruteforce(ta::TrjArray{T, U}, rcut::T; iframe=1::Int) where {T, U}
```