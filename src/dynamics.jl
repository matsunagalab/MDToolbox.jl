"""
dummy

dummy function
"""
function dummy(x::Matrix{Float64}, y::Matrix{Float64}, z::Matrix{Float64}, ps_per_frame)::Matrix{Float64}
    nframe = size(x, 1)
    natom = size(x, 2)

    nframe_of_msd = floor(Int64, nframe/10)
    msd = zeros(Float64, nframe_of_msd)
    time = (1:nframe_of_msd).*ps_per_frame

    for iframe = 1:nframe_of_msd
        @show iframe
        nblock = floor(Int64, nframe/iframe)
        d = zeros(Float64, 1, natom)
        d = d .+ sum((x[iframe:iframe:end, :] .- x[1:iframe:(end-iframe+1), :]).^2, dims=1)
        d = d .+ sum((y[iframe:iframe:end, :] .- y[1:iframe:(end-iframe+1), :]).^2, dims=1)
        d = d .+ sum((z[iframe:iframe:end, :] .- z[1:iframe:(end-iframe+1), :]).^2, dims=1)
        d = d ./ nblock
        msd[iframe] = mean(d)
    end
end


############################################################################
