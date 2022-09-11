#using .Bio3DView

function split_into_residues(ta::TrjArray)
    natom = ta.natom

    properties = Dict()
    if isempty(ta.chainid)
        properties["chainid"] = [nothing]
    else
        properties["chainid"] = unique(ta.chainid)
    end
    if isempty(ta.chainname)
        properties["chainname"] = [nothing]
    else
        properties["chainname"] = unique(ta.chainname)
    end
    if isempty(ta.resid)
        properties["resid"] = [nothing]
    else
        properties["resid"] = unique(ta.resid)
    end
    if isempty(ta.resname)
        properties["resname"] = [nothing]
    else
        properties["resname"] = unique(ta.resname)
    end

    res_array = []
    id = BitArray(undef, natom)
    id_c1 = BitArray(undef, natom)
    id_c2 = BitArray(undef, natom)
    id_r1 = BitArray(undef, natom)
    id_r2 = BitArray(undef, natom)
    for c1 in properties["chainid"]
        id .= id_c1 .= id_c2 .= id_r1 .= id_r2 .= 1
        if !isnothing(c1)
            id_c1 .= (id .& Bool(ta.chainid .== c1))
        end
        for c2 in properties["chainname"]
            if !isnothing(c2)
                id_c2 .= (id_c1 .& (ta.chainname .== c2))
            end
            for r1 in properties["resid"]
                if !isnothing(r1)
                    id_r1 .= (id_c2 .& (ta.resid .== r1))
                end
                for r2 in properties["resname"]
                    if !isnothing(r2)
                        id_r2 .= (id_r1 .& (ta.resname .== r2))
                    end
                    if any(id_r2)
                        push!(res_array, ta[:, id_r2])
                    end
                end
            end
        end
    end

    return res_array
end

function split_into_chains(ta::TrjArray)
    natom = ta.natom

    properties = Dict()
    if isempty(ta.chainid)
        properties["chainid"] = [nothing]
    else
        properties["chainid"] = unique(ta.chainid)
    end
    if isempty(ta.chainname)
        properties["chainname"] = [nothing]
    else
        properties["chainname"] = unique(ta.chainname)
    end

    chain_array = []
    id = BitArray(undef, natom)
    id_c1 = BitArray(undef, natom)
    id_c2 = BitArray(undef, natom)
    for c1 in properties["chainid"]
        id .= id_c1 .= id_c2 .= 1
        if !isnothing(c1)
            id_c1 .= (id .& Bool(ta.chainid .== c1))
        end
        for c2 in properties["chainname"]
            if !isnothing(c2)
                id_c2 .= (id_c1 .& (ta.chainname .== c2))
            end
            push!(chain_array, ta[:, id_c2])
        end
    end

    return chain_array
end

#"""
#viewstruct
#
#based on Bio3Dview.jl
#"""
#function viewstruc(ta::TrjArray; kwargs...)
#    io = IOBuffer()
#    writepdb(io, ta)
#    return Bio3DView.view("data-type='pdb'", String(take!(io)); style=Bio3DView.defaultstyle("pdb"), kwargs...)
#end

#cpu(m) = fmap(x -> adapt(Array, x), m)
#gpu(x) = use_cuda[] ? fmap(CUDA.cu, x) : x

function CuArray(ta::TrjArray{T, U}) where {T, U}
    TrjArray{T, U}(CuArray(ta.xyz), CuArray(ta.boxsize),
      ta.chainname, CuArray(ta.chainid),
      ta.resname, CuArray(ta.resid),
      ta.atomname, ta.atomtype, CuArray(ta.atomtype_id), CuArray(ta.atomid),
      CuArray(ta.mass), CuArray(ta.radius), CuArray(ta.charge), CuArray(ta.sasa),
      CuArray(ta.nnb), 
      CuArray(ta.list_bond), CuArray(ta.list_angle), CuArray(ta.list_dihedral),
      CuArray(ta.list_improper), CuArray(ta.list_cmap), 
      CuArray(ta.list_donor), CuArray(ta.list_acceptor))
end

function gpu(ta::TrjArray{T, U}) where {T, U}
    TrjArray{Float32, Int32}(CUDA.cu(ta.xyz), CUDA.cu(ta.boxsize),
      ta.chainname, CUDA.cu(ta.chainid),
      ta.resname, CUDA.cu(ta.resid),
      ta.atomname, ta.atomtype, CUDA.cu(ta.atomtype_id), CUDA.cu(ta.atomid),
      CUDA.cu(ta.mass), CUDA.cu(ta.radius), CUDA.cu(ta.charge), CUDA.cu(ta.sasa),
      CUDA.cu(ta.nnb), 
      CUDA.cu(ta.list_bond), CUDA.cu(ta.list_angle), CUDA.cu(ta.list_dihedral),
      CUDA.cu(ta.list_improper), CUDA.cu(ta.list_cmap), 
      CUDA.cu(ta.list_donor), CUDA.cu(ta.list_acceptor))
end

function cu(ta::TrjArray{T, U}) where {T, U}
    return gpu(ta)
end

function Array(ta::TrjArray{T, U}) where {T, U}
    TrjArray{Float64, Int64}(Array{Float64}(ta.xyz), Array{Float64}(ta.boxsize),
      ta.chainname, Array{Float64}(ta.chainid),
      ta.resname, Array{Float64}(ta.resid),
      ta.atomname, ta.atomtype, Array{Float64}(ta.atomtype_id), Array{Float64}(ta.atomid),
      Array{Float64}(ta.mass), Array{Float64}(ta.radius), Array{Float64}(ta.charge), Array{Float64}(ta.sasa),
      Array{Float64}(ta.nnb), 
      Array{Int64}(ta.list_bond), Array{Int64}(ta.list_angle), Array{Int64}(ta.list_dihedral),
      Array{Int64}(ta.list_improper), Array{Int64}(ta.list_cmap), 
      Array{Int64}(ta.list_donor), Array{Int64}(ta.list_acceptor))
end

"""
    logsumexp(X)

Compute `log(sum(exp, X))` avoiding possible overflows/undeflows.
`X` should be an iterator of real numbers.
"""
function logsumexp(X)
    isempty(X) && return log(sum(X))
    reduce(logaddexp, X)
end

function logsumexp(X::AbstractArray{T}; dims=:) where {T<:Real}
    # Do not use log(zero(T)) directly to avoid issues with ForwardDiff (#82)
    u = reduce(max, X, dims=dims, init=oftype(log(zero(T)), -Inf))
    u isa AbstractArray || isfinite(u) || return float(u)
    let u=u # avoid https://github.com/JuliaLang/julia/issues/15276
        # TODO: remove the branch when JuliaLang/julia#31020 is merged.
        if u isa AbstractArray
            u .+ log.(sum(exp.(X .- u); dims=dims))
        else
            u + log(sum(x -> exp(x-u), X))
        end
    end
end
