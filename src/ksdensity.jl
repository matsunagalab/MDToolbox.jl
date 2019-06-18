function convolution1d!(tid, iframe, x, nx, rx, gridsize_x, grid_x, bandwidth, f_x, w)
    dx = x[iframe] - grid_x[1]
    ind_min = floor(Int, (dx - rx)/gridsize_x) + 1
    ind_min = ind_min >= 1 ? ind_min : 1
    ind_max = floor(Int, (dx + rx)/gridsize_x) + 2
    ind_max = ind_max <= nx ? ind_max : nx
    for ind in ind_min:ind_max
        dx = (grid_x[ind] - x[iframe])/bandwidth
        f_x[ind, tid] += w*exp(-0.5*dx*dx)/(sqrt(2.0*pi)*bandwidth)
    end
end

function convolution1d_box!(tid, iframe, boxsize, x, nx, rx, gridsize_x, grid_x, bandwidth, f_x, w)
    for ind in 1:nx
        dx = x[iframe] - grid_x[ind]
        dx = dx - round(dx/boxsize)*boxsize;
        if abs(dx) < rx
          dx = dx/bandwidth
          f_x[ind, tid] += w*exp(-0.5*dx*dx)/(sqrt(2.0*pi)*bandwidth)
        end
    end
end

function convolution2d!(tid, iframe,
    x, nx, rx, gridsize_x, grid_x,
    y, ny, ry, gridsize_y, grid_y,
    bandwidth, f_x, f_y, f_private, w)

    f_x[:, tid] .= 0.0
    f_y[:, tid] .= 0.0

    dx = x[iframe] - grid_x[1]
    ix_min = floor(Int, (dx - rx)/gridsize_x) + 1
    ix_min = ix_min >= 1 ? ix_min : 1
    ix_max = floor(Int, (dx + rx)/gridsize_x) + 2
    ix_max = ix_max <= nx ? ix_max : nx
    for ix in ix_min:ix_max
        dx = (grid_x[ix] - x[iframe])/bandwidth[1]
        f_x[ix, tid] += exp(-0.5*dx*dx)/(sqrt(2.0*pi)*bandwidth[1])
    end

    dy = y[iframe] - grid_y[1]
    iy_min = floor(Int, (dy - ry)/gridsize_y) + 1
    iy_min = iy_min >= 1 ? iy_min : 1
    iy_max = floor(Int, (dy + ry)/gridsize_y) + 2
    iy_max = iy_max <= ny ? iy_max : ny
    for iy in iy_min:iy_max
        dy = (grid_y[iy] - y[iframe])/bandwidth[2]
        f_y[iy, tid] += exp(-0.5*dy*dy)/(sqrt(2.0*pi)*bandwidth[2])
    end

    for ix in ix_min:ix_max
        for iy in iy_min:iy_max
            f_private[ix, iy, tid] += w*f_x[ix, tid]*f_y[iy, tid]
        end
    end
end

function convolution2d_box!(tid, iframe, boxsize,
    x, nx, rx, gridsize_x, grid_x,
    y, ny, ry, gridsize_y, grid_y,
    bandwidth, f_x, f_y, f_private, w,
    ix_array, iy_array)

    f_x[:, tid] .= 0.0
    f_y[:, tid] .= 0.0

    ix_count = 0
    for ix in 1:nx
        dx = x[iframe] - grid_x[ix]
        dx = dx - round(dx/boxsize[1])*boxsize[1];
        if abs(dx) < rx
          dx = dx/bandwidth[1]
          f_x[ix, tid] += exp(-0.5*dx*dx)/(sqrt(2.0*pi)*bandwidth[1])
          ix_count += 1
          ix_array[ix_count, tid] = ix
        end
    end

    iy_count = 0
    for iy in 1:ny
        dy = y[iframe] - grid_y[iy]
        dy = dy - round(dy/boxsize[2])*boxsize[2];
        if abs(dy) < ry
          dy = dy/bandwidth[2]
          f_y[iy, tid] += exp(-0.5*dy*dy)/(sqrt(2.0*pi)*bandwidth[2])
          iy_count += 1
          iy_array[iy_count, tid] = iy
        end
    end

    for i in 1:ix_count
        for j in 1:iy_count
            ix = ix_array[i, tid]
            iy = iy_array[j, tid]
            f_private[ix, iy, tid] += w*f_x[ix, tid]*f_y[iy, tid]
        end
    end
end

function convolution3d!(tid, iframe,
    x, nx, rx, gridsize_x, grid_x,
    y, ny, ry, gridsize_y, grid_y,
    z, nz, rz, gridsize_z, grid_z,
    bandwidth, f_x, f_y, f_z, f_private, w)

    f_x[:, tid] .= 0.0
    f_y[:, tid] .= 0.0
    f_z[:, tid] .= 0.0

    dx = x[iframe] - grid_x[1]
    ix_min = floor(Int, (dx - rx)/gridsize_x) + 1
    ix_min = ix_min >= 1 ? ix_min : 1
    ix_max = floor(Int, (dx + rx)/gridsize_x) + 2
    ix_max = ix_max <= nx ? ix_max : nx
    for ix in ix_min:ix_max
        dx = (grid_x[ix] - x[iframe])/bandwidth[1]
        f_x[ix, tid] += exp(-0.5*dx*dx)/(sqrt(2.0*pi)*bandwidth[1])
    end

    dy = y[iframe] - grid_y[1]
    iy_min = floor(Int, (dy - ry)/gridsize_y) + 1
    iy_min = iy_min >= 1 ? iy_min : 1
    iy_max = floor(Int, (dy + ry)/gridsize_y) + 2
    iy_max = iy_max <= ny ? iy_max : ny
    for iy in iy_min:iy_max
        dy = (grid_y[iy] - y[iframe])/bandwidth[2]
        f_y[iy, tid] += exp(-0.5*dy*dy)/(sqrt(2.0*pi)*bandwidth[2])
    end

    dz = z[iframe] - grid_z[1]
    iz_min = floor(Int, (dz - rz)/gridsize_z) + 1
    iz_min = iz_min >= 1 ? iz_min : 1
    iz_max = floor(Int, (dz + rz)/gridsize_z) + 2
    iz_max = iz_max <= nz ? iz_max : nz
    for iz in iz_min:iz_max
        dz = (grid_z[iz] - z[iframe])/bandwidth[3]
        f_z[iz, tid] += exp(-0.5*dz*dz)/(sqrt(2.0*pi)*bandwidth[3])
    end

    for ix in ix_min:ix_max
        for iy in iy_min:iy_max
            for iz in iz_min:iz_max
                f_private[ix, iy, iz, tid] += w*f_x[ix, tid]*f_y[iy, tid]*f_z[iz, tid]
            end
        end
    end
end

function convolution3d_box!(tid, iframe, boxsize,
    x, nx, rx, gridsize_x, grid_x,
    y, ny, ry, gridsize_y, grid_y,
    z, nz, rz, gridsize_z, grid_z,
    bandwidth, f_x, f_y, f_z, f_private, w,
    ix_array, iy_array, iz_array)

    f_x[:, tid] .= 0.0
    f_y[:, tid] .= 0.0
    f_z[:, tid] .= 0.0

    ix_count = 0
    for ix in 1:nx
        dx = x[iframe] - grid_x[ix]
        dx = dx - round(dx/boxsize[1])*boxsize[1];
        if abs(dx) < rx
          dx = dx/bandwidth[1]
          f_x[ix, tid] += exp(-0.5*dx*dx)/(sqrt(2.0*pi)*bandwidth[1])
          ix_count += 1
          ix_array[ix_count, tid] = ix
        end
    end

    iy_count = 0
    for iy in 1:ny
        dy = y[iframe] - grid_y[iy]
        dy = dy - round(dy/boxsize[2])*boxsize[2];
        if abs(dy) < ry
          dy = dy/bandwidth[2]
          f_y[iy, tid] += exp(-0.5*dy*dy)/(sqrt(2.0*pi)*bandwidth[2])
          iy_count += 1
          iy_array[iy_count, tid] = iy
        end
    end

    iz_count = 0
    for iz in 1:nz
        dz = z[iframe] - grid_z[iz]
        dz = dz - round(dz/boxsize[3])*boxsize[3];
        if abs(dz) < rz
          dz = dz/bandwidth[3]
          f_z[iz, tid] += exp(-0.5*dz*dz)/(sqrt(2.0*pi)*bandwidth[3])
          iz_count += 1
          iz_array[iz_count, tid] = iz
        end
    end

    for i in 1:ix_count
        for j in 1:iy_count
            for k in 1:iz_count
                ix = ix_array[i, tid]
                iy = iy_array[j, tid]
                iz = iz_array[k, tid]
                f_private[ix, iy, iz, tid] += w*f_x[ix, tid]*f_y[iy, tid]*f_z[iz, tid]
            end
        end
    end
end

"""
ksdensity

kernel density estimator with gaussian kernel
"""
function ksdensity(x::Vector{Float64};
    grid_x,
    weight::Vector{Float64}=Vector{Float64}(undef, 0),
    boxsize::Union{Float64,Missing}=missing,
    bandwidth::Union{Float64,Missing}=missing)
    nframe = length(x)
    if isempty(grid_x)
        grid_x2 = range(minimum(x), stop=maximum(x), length=100)
    else typeof(grid_x) != Vector{Float64}
        grid_x2 = map(Float64, grid_x)
    end
    nx = length(grid_x2)
    gridsize_x = abs(grid_x2[2]-grid_x2[1])
    if isempty(weight)
        weight2 = ones(Float64, nframe)/nframe
    else
        weight2 = weight
    end
    if typeof(bandwidth) == Missing
        bandwidth = gridsize_x
    end
    rx = bandwidth*5.0
    nth = Threads.nthreads()
    f_x = zeros(nx, nth)
    d = div(nframe, nth)
    if typeof(boxsize) == Missing
        Threads.@threads for id in 1:nth
            tid = Threads.threadid()
            frames = (tid-1)*d+1:tid*d
            for iframe in frames
                convolution1d!(tid, iframe, x, nx, rx, gridsize_x, grid_x2, bandwidth, f_x, weight2[iframe])
            end
        end
    else
        Threads.@threads for id in 1:nth
            tid = Threads.threadid()
            frames = (tid-1)*d+1:tid*d
            for iframe in frames
                convolution1d_box!(tid, iframe, boxsize, x, nx, rx, gridsize_x, grid_x2, bandwidth, f_x, weight2[iframe])
            end
        end
    end
    f = zeros(nx)
    for ind in 1:nx
        f[ind] = sum(f_x[ind, :])
    end
    f, grid_x2
end

function ksdensity_serial(x::Vector{Float64};
    grid_x,
    weight::Vector{Float64}=Vector{Float64}(undef, 0),
    boxsize::Union{Float64,Missing}=missing,
    bandwidth::Union{Float64,Missing}=missing)
    nframe = length(x)
    if isempty(grid_x)
        grid_x2 = range(minimum(x), stop=maximum(x), length=100)
    else typeof(grid_x) != Vector{Float64}
        grid_x2 = map(Float64, grid_x)
    end
    nx = length(grid_x2)
    gridsize_x = abs(grid_x2[2]-grid_x2[1])
    if isempty(weight)
        weight2 = ones(Float64, nframe)/nframe
    else
        weight2 = weight
    end
    if typeof(bandwidth) == Missing
        bandwidth = gridsize_x
    end
    rx = bandwidth*5.0
    nth = 1
    f_x = zeros(nx, nth)
    if typeof(boxsize) == Missing
        for iframe in 1:nframe
            convolution1d!(1, iframe, x, nx, rx, gridsize_x, grid_x2, bandwidth, f_x, weight2[iframe])
        end
    else
        for iframe in 1:nframe
            convolution1d_box!(1, iframe, boxsize, x, nx, rx, gridsize_x, grid_x2, bandwidth, f_x, weight2[iframe])
        end
    end
    f = zeros(nx)
    for ind in 1:nx
        f[ind] = sum(f_x[ind, :])
    end
    f, grid_x2
end

function ksdensity(x::Vector{Float64}, y::Vector{Float64};
    grid_x, grid_y,
    weight::Vector{Float64}=Vector{Float64}(undef, 0),
    boxsize::Vector{Float64}=Vector{Float64}(undef, 0),
    bandwidth::Vector{Float64}=Vector{Float64}(undef, 0))
    nframe = length(x)
    if isempty(grid_x)
        grid_x2 = range(minimum(x), stop=maximum(x), length=100)
    else typeof(grid_x) != Vector{Float64}
        grid_x2 = map(Float64, grid_x)
    end
    if isempty(grid_y)
        grid_y2 = range(minimum(y), stop=maximum(y), length=100)
    else typeof(grid_y) != Vector{Float64}
        grid_y2 = map(Float64, grid_y)
    end
    nx = length(grid_x2)
    ny = length(grid_y2)
    gridsize_x = abs(grid_x2[2]-grid_x2[1])
    gridsize_y = abs(grid_y2[2]-grid_y2[1])
    if isempty(weight)
        weight2 = ones(Float64, nframe)/nframe
    else
        weight2 = weight
    end
    if isempty(bandwidth)
        bandwidth = [gridsize_x, gridsize_y]
    end
    rx = bandwidth[1]*5.0
    ry = bandwidth[2]*5.0
    nth = Threads.nthreads()
    f_x = zeros(nx, nth)
    f_y = zeros(ny, nth)
    f_private = zeros(nx, ny, nth)
    d = div(nframe, nth)
    if isempty(boxsize)
        Threads.@threads for id in 1:nth
            tid = Threads.threadid()
            frames = (tid-1)*d+1:tid*d
            for iframe in frames
                convolution2d!(tid, iframe,
                    x, nx, rx, gridsize_x, grid_x2,
                    y, ny, ry, gridsize_y, grid_y2,
                    bandwidth, f_x, f_y, f_private, weight2[iframe])
            end
        end
    else
        ix_array = zeros(Int, nx, nth)
        iy_array = zeros(Int, ny, nth)
        Threads.@threads for id in 1:nth
            tid = Threads.threadid()
            frames = (tid-1)*d+1:tid*d
            for iframe in frames
                convolution2d_box!(tid, iframe, boxsize,
                    x, nx, rx, gridsize_x, grid_x2,
                    y, ny, ry, gridsize_y, grid_y2,
                    bandwidth, f_x, f_y, f_private, weight2[iframe],
                    ix_array, iy_array)
            end
        end

    end
    f = zeros(nx, ny)
    for ix in 1:nx
        for iy in 1:ny
            f[ix, iy] = sum(f_private[ix, iy, :])
        end
    end
    transpose(f), grid_x2, grid_y2
end

function ksdensity(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64};
    grid_x, grid_y, grid_z,
    weight::Vector{Float64}=Vector{Float64}(undef, 0),
    boxsize::Vector{Float64}=Vector{Float64}(undef, 0),
    bandwidth::Vector{Float64}=Vector{Float64}(undef, 0))
    nframe = length(x)
    if isempty(grid_x)
        grid_x2 = range(minimum(x), stop=maximum(x), length=100)
    else typeof(grid_x) != Vector{Float64}
        grid_x2 = map(Float64, grid_x)
    end
    if isempty(grid_y)
        grid_y2 = range(minimum(y), stop=maximum(y), length=100)
    else typeof(grid_y) != Vector{Float64}
        grid_y2 = map(Float64, grid_y)
    end
    if isempty(grid_z)
        grid_z2 = range(minimum(z), stop=maximum(z), length=100)
    else typeof(grid_z) != Vector{Float64}
        grid_z2 = map(Float64, grid_z)
    end
    nx = length(grid_x2)
    ny = length(grid_y2)
    nz = length(grid_z2)
    gridsize_x = abs(grid_x2[2]-grid_x2[1])
    gridsize_y = abs(grid_y2[2]-grid_y2[1])
    gridsize_z = abs(grid_z2[2]-grid_z2[1])
    if isempty(weight)
        weight2 = ones(Float64, nframe)/nframe
    else
        weight2 = weight
    end
    if isempty(bandwidth)
        bandwidth = [gridsize_x, gridsize_y, gridsize_z]
    end
    rx = bandwidth[1]*5.0
    ry = bandwidth[2]*5.0
    rz = bandwidth[3]*5.0
    nth = Threads.nthreads()
    f_x = zeros(nx, nth)
    f_y = zeros(ny, nth)
    f_z = zeros(nz, nth)
    f_private = zeros(nx, ny, nz, nth)
    d = div(nframe, nth)
    if isempty(boxsize)
        Threads.@threads for id in 1:nth
            tid = Threads.threadid()
            frames = (tid-1)*d+1:tid*d
            for iframe in frames
                convolution3d!(tid, iframe,
                    x, nx, rx, gridsize_x, grid_x2,
                    y, ny, ry, gridsize_y, grid_y2,
                    z, nz, rz, gridsize_z, grid_z2,
                    bandwidth, f_x, f_y, f_z, f_private, weight2[iframe])
            end
        end
    else
        ix_array = zeros(Int, nx, nth)
        iy_array = zeros(Int, ny, nth)
        iz_array = zeros(Int, nz, nth)
        Threads.@threads for id in 1:nth
            tid = Threads.threadid()
            frames = (tid-1)*d+1:tid*d
            for iframe in frames
                convolution3d_box!(tid, iframe, boxsize,
                    x, nx, rx, gridsize_x, grid_x2,
                    y, ny, ry, gridsize_y, grid_y2,
                    z, nz, rz, gridsize_z, grid_z2,
                    bandwidth, f_x, f_y, f_z, f_private, weight2[iframe],
                    ix_array, iy_array, iz_array)
            end
        end
    end
    f = zeros(nx, ny, nz)
    for ix in 1:nx
        for iy in 1:ny
            for iz in 1:nz
                f[ix, iy, iz] = sum(f_private[ix, iy, iz, :])
            end
        end
    end
    f, grid_x2, grid_y2, grid_z2
end

"""
getpmf

Potential of mean force estimator by a kernel density estimator
"""
function getpmf(x::Vector{Float64};
    grid_x,
    weight::Vector{Float64}=Vector{Float64}(undef, 0),
    boxsize::Union{Float64,Missing}=missing,
    bandwidth::Union{Float64,Missing}=missing)

    pdf, grid_x = ksdensity(x, grid_x=grid_x, weight=weight, boxsize=boxsize, bandwidth=bandwidth)
    #pdf[pdf .< eps()] .= missing
    pmf = -log.(pdf)
    pmf = pmf .- minimum(pmf)
    pmf, grid_x
end

function getpmf(x::Vector{Float64}, y::Vector{Float64};
    grid_x, grid_y,
    weight::Vector{Float64}=Vector{Float64}(undef, 0),
    boxsize::Vector{Float64}=Vector{Float64}(undef, 0),
    bandwidth::Vector{Float64}=Vector{Float64}(undef, 0))

    pdf, grid_x, grid_y = ksdensity(x, y, grid_x=grid_x, grid_y=grid_y, weight=weight, boxsize=boxsize, bandwidth=bandwidth)
    #pdf[pdf .< eps()] .= NaN
    pmf = -log.(pdf)
    pmf = pmf .- minimum(pmf)
    pmf, grid_x, grid_y
end

function getpmf(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64};
    grid_x, grid_y, grid_z,
    weight::Vector{Float64}=Vector{Float64}(undef, 0),
    boxsize::Vector{Float64}=Vector{Float64}(undef, 0),
    bandwidth::Vector{Float64}=Vector{Float64}(undef, 0))

    pdf, grid_x, grid_y, grid_z = ksdensity(x, y, z, grid_x=grid_x, grid_y=grid_y, grid_z=grid_z, weight=weight, boxsize=boxsize, bandwidth=bandwidth)
    #pdf[pdf .< eps()] .= NaN
    pmf = -log.(pdf)
    pmf = pmf .- minimum(pmf)
    pmf, grid_x, grid_y, grid_z
end
