###### center of mass #################
"""
centerofmass

calculate center of mass of given trajectory
"""
function centerofmass(ta::TrjArray; isweight::Bool=true, index::Vector{Int64}=Vector{Int64}(undef, 0))::TrjArray
    nframe = ta.nframe
    natom = ta.natom
    if isweight
        @assert length(ta.mass) == natom
        weight = ta.mass
    else
        weight = ones(Float64, natom)
    end
    if isempty(index)
        index2 = collect(1:natom)
    else
        index2 = index
    end
    natom_sub = length(index2)
    if natom_sub == 1
        return TrjArray(x=ta.x, y=ta.y, z=ta.z)
    else
        weight2 = reshape(weight[index2], 1, natom_sub)
        wsum_inv = 1.0 / sum(weight2)
        x = sum(weight2 .* view(ta.x, :, index2), dims=2) .* wsum_inv
        y = sum(weight2 .* view(ta.y, :, index2), dims=2) .* wsum_inv
        z = sum(weight2 .* view(ta.z, :, index2), dims=2) .* wsum_inv
        return TrjArray(x=x, y=y, z=z)
    end
end


###### superimpose #################
function innerproduct!(ref_x::Vector{Float64}, ref_y::Vector{Float64}, ref_z::Vector{Float64}, 
                       x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, 
                       weight::Vector{Float64}, index::Vector{Int64}, A::Vector{Float64})
    A[:] .= 0.0
    G1 = G2 = 0.0
    for i in index
        x1 = weight[i] * ref_x[i]
        y1 = weight[i] * ref_y[i]
        z1 = weight[i] * ref_z[i]

        G1 += x1 * ref_x[i] + y1 * ref_y[i] + z1 * ref_z[i]

        x2 = x[i]
        y2 = y[i]
        z2 = z[i]

        G2 += weight[i] * (x2 * x2 + y2 * y2 + z2 * z2)

        A[1] +=  (x1 * x2)
        A[2] +=  (x1 * y2)
        A[3] +=  (x1 * z2)

        A[4] +=  (y1 * x2)
        A[5] +=  (y1 * y2)
        A[6] +=  (y1 * z2)

        A[7] +=  (z1 * x2)
        A[8] +=  (z1 * y2)
        A[9] +=  (z1 * z2)
    end
    return (G1 + G2) * 0.5
end

function fastCalcRMSDAndRotation!(A::Vector{Float64}, E0::Float64, wsum_inv::Float64, rot::Vector{Float64}, C::Vector{Float64})
    oldg = 0.0
    #C = zeros(Float64, 3)
    #rot = zeros(Float64, 9)
    evecprec = 1e-6
    evalprec = 1e-11

    Sxx = A[1]; Sxy = A[2]; Sxz = A[3]
    Syx = A[4]; Syy = A[5]; Syz = A[6]
    Szx = A[7]; Szy = A[8]; Szz = A[9]

    Sxx2 = Sxx * Sxx
    Syy2 = Syy * Syy
    Szz2 = Szz * Szz

    Sxy2 = Sxy * Sxy
    Syz2 = Syz * Syz
    Sxz2 = Sxz * Sxz

    Syx2 = Syx * Syx
    Szy2 = Szy * Szy
    Szx2 = Szx * Szx

    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz)
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2

    C[3] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2)
    C[2] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz)

    SxzpSzx = Sxz + Szx
    SyzpSzy = Syz + Szy
    SxypSyx = Sxy + Syx
    SyzmSzy = Syz - Szy
    SxzmSzx = Sxz - Szx
    SxymSyx = Sxy - Syx
    SxxpSyy = Sxx + Syy
    SxxmSyy = Sxx - Syy
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2

    C[1] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2 + 
         (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2) + 
         (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz)) + 
         (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz)) + 
         (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz)) + 
         (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz))

    # Newton-Raphson
    mxEigenV = E0
    icount = 0
    for i in 1:50
        icount += 1
        oldg = mxEigenV
        x2 = mxEigenV*mxEigenV
        b = (x2 + C[3])*mxEigenV
        a = b + C[2]
        delta = ((a*mxEigenV + C[1])/(2.0*x2*mxEigenV + b + a));
        mxEigenV -= delta;
        # printf("\n diff[%3d]: %16g %16g %16g", i, mxEigenV - oldg, evalprec*mxEigenV, mxEigenV)
        if abs(mxEigenV - oldg) < abs(evalprec*mxEigenV)
            break
        end
    end

    if icount == 50
        print("More than 50 iterations needed!")
    end

    # the abs() is to guard against extremely small, but *negative* numbers due to floating point error */
    rmsd = sqrt(abs(2.0 * (E0 - mxEigenV) * wsum_inv))
    # printf("\n\n %16g %16g %16g \n", rmsd, E0, 2.0 * (E0 - mxEigenV)/len)

    a11 = SxxpSyy + Szz-mxEigenV; a12 = SyzmSzy; a13 = - SxzmSzx; a14 = SxymSyx
    a21 = SyzmSzy; a22 = SxxmSyy - Szz-mxEigenV; a23 = SxypSyx; a24= SxzpSzx
    a31 = a13; a32 = a23; a33 = Syy-Sxx-Szz - mxEigenV; a34 = SyzpSzy
    a41 = a14; a42 = a24; a43 = a34; a44 = Szz - SxxpSyy - mxEigenV
    a3344_4334 = a33 * a44 - a43 * a34; a3244_4234 = a32 * a44-a42*a34
    a3243_4233 = a32 * a43 - a42 * a33; a3143_4133 = a31 * a43-a41*a33
    a3144_4134 = a31 * a44 - a41 * a34; a3142_4132 = a31 * a42-a41*a32
    q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233
    q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133
    q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132
    q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132

    qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4

    # The following code tries to calculate another column in the adjoint matrix when the norm of the
    # current column is too small.
    # Usually this block will never be activated.  To be absolutely safe this should be
    # uncommented, but it is most likely unnecessary.
    while qsqr < evecprec
        q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233
        q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133
        q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132
        q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132
        qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4

        if qsqr < evecprec
            a1324_1423 = a13 * a24 - a14 * a23; a1224_1422 = a12 * a24 - a14 * a22
            a1223_1322 = a12 * a23 - a13 * a22; a1124_1421 = a11 * a24 - a14 * a21
            a1123_1321 = a11 * a23 - a13 * a21; a1122_1221 = a11 * a22 - a12 * a21

            q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322
            q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321
            q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221
            q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221
            qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4

            if qsqr < evecprec
                q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322
                q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321
                q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221
                q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221
                qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4;

                if qsqr < evecprec
                    # if qsqr is still too small, return the identity matrix.
                    rot[1] = rot[5] = rot[9] = 1.0
                    rot[2] = rot[3] = rot[4] = rot[6] = rot[7] = rot[8] = 0.0
                    return rmsd
                end
            end
        end
    end

    normq = sqrt(qsqr)
    q1 /= normq
    q2 /= normq
    q3 /= normq
    q4 /= normq

    a2 = q1 * q1
    x2 = q2 * q2
    y2 = q3 * q3
    z2 = q4 * q4

    xy = q2 * q3
    az = q1 * q4
    zx = q4 * q2
    ay = q1 * q3
    yz = q3 * q4
    ax = q1 * q2

    rot[1] = a2 + x2 - y2 - z2
    rot[2] = 2 * (xy + az)
    rot[3] = 2 * (zx - ay)
    rot[4] = 2 * (xy - az)
    rot[5] = a2 - x2 + y2 - z2
    rot[6] = 2 * (yz + ax)
    rot[7] = 2 * (zx + ay)
    rot[8] = 2 * (yz - ax)
    rot[9] = a2 - x2 - y2 + z2
    return rmsd
end

"""
superimpose

superimpose ta::TrjArray to ref:TrjArray

this code is licensed under the BSD license (Copyright 2009-2012 Pu Liu and Douglas L. Theobald), see LICENSE.md 
"""
function superimpose(ref::TrjArray, ta::TrjArray; isweight::Bool=true, index::Vector{Int64}=Vector{Int64}(undef, 0), iscomremoved::Bool=false)::Tuple{Array{Float64,1},TrjArray}
    nframe = ta.nframe
    natom = ta.natom
    if isempty(index)
        index2 = collect(1:natom)
    else
        index2 = index
    end
    if isweight
        @assert length(ta.mass) == natom
        weight = ta.mass
    else
        weight = ones(Float64, natom)
    end

    if iscomremoved
        ta_x = ta.x
        ta_y = ta.y
        ta_z = ta.z
        ref_x = ref.x[1, :]
        ref_y = ref.y[1, :]
        ref_z = ref.z[1, :]
        weight2 = reshape(weight[index2], length(index2))
        wsum_inv = 1.0 / sum(weight2)
    else
        com = centerofmass(ta, isweight=isweight, index=index2)
        ta_x = ta.x .- com.x
        ta_y = ta.y .- com.y
        ta_z = ta.z .- com.z
        com = centerofmass(ref[1, :], isweight=isweight, index=index2)
        ref_x = ref.x[1, :] .- com.x[1]
        ref_y = ref.y[1, :] .- com.y[1]
        ref_z = ref.z[1, :] .- com.z[1]
        weight2 = reshape(weight[index2], length(index2))
        wsum_inv = 1.0 / sum(weight2)
    end

    x = Matrix{Float64}(undef, nframe, natom)
    y = Matrix{Float64}(undef, nframe, natom)
    z = Matrix{Float64}(undef, nframe, natom)
    rmsd = Vector{Float64}(undef, nframe)
    A = Vector{Float64}(undef, 9)
    rot = Vector{Float64}(undef, 9)
    C = Vector{Float64}(undef, 3)
    for iframe in 1:nframe
        frame_x = ta_x[iframe, :]
        frame_y = ta_y[iframe, :]
        frame_z = ta_z[iframe, :]
        E0 = innerproduct!(ref_x, ref_y, ref_z, frame_x, frame_y, frame_z, weight2, index2, A)
        r = fastCalcRMSDAndRotation!(A, E0, wsum_inv, rot, C)
        rmsd[iframe] = r
        x[iframe, :] .= rot[1] .* frame_x .+ rot[2] .* frame_y .+ rot[3] .* frame_z
        y[iframe, :] .= rot[4] .* frame_x .+ rot[5] .* frame_y .+ rot[6] .* frame_z
        z[iframe, :] .= rot[7] .* frame_x .+ rot[8] .* frame_y .+ rot[9] .* frame_z
    end
    if !iscomremoved
        x = x .+ com.x[1]
        y = y .+ com.y[1]
        z = z .+ com.z[1]
    end
    ta_fit = TrjArray(x, y, z, ta)
    rmsd, ta_fit
end

###### rmsd #################
"""
rmsd

root mean square deviation of ta::TrjArray from ref:TrjArray
"""
function calcrmsd(ref::TrjArray, ta::TrjArray; isweight::Bool=true, index::Vector{Int64}=Vector{Int64}(undef, 0))::Vector{Float64}
    nframe = ta.nframe
    natom = ta.natom
    if isweight
        @assert length(ta.mass) == natom
        weight = ta.mass
    else
        weight = ones(Float64, natom)
    end
    if isempty(index)
        index2 = cumsum(ones(Int64, natom))
    else
        index2 = index
    end
    wsum_inv = 1.0 / sum(weight[index2])
    weight2 = reshape(weight[index2], 1, length(index2))
    ref_x = ref.x[1:1, index2]
    ref_y = ref.y[1:1, index2]
    ref_z = ref.z[1:1, index2]
    ta_x = ta.x[:, index2]
    ta_y = ta.y[:, index2]
    ta_z = ta.z[:, index2]
    d =  sum(weight2 .* ((ta_x .- ref_x).^2 .+ (ta_y .- ref_y).^2 .+ (ta_z .- ref_z).^2), dims=2)
    d = d .* wsum_inv
    d = sqrt.(d)
    rmsd = reshape(d, nframe)
end


###### distance, angle, dihedral #################
"""
calcdistance

calculate distances between two atoms or groups of atoms
"""
function calcbond(ta1::TrjArray, ta2::TrjArray)::Vector{Float64}
    # TODO: support for PBC
    @assert ta1.nframe == ta2.nframe
    nframe = ta1.nframe
    com1 = centerofmass(ta1)
    com2 = centerofmass(ta2)
    dist = Vector{Float64}(undef, nframe)
    @threads for iframe in 1:nframe
        d = 0.0
        d += (com1.x[iframe] - com2.x[iframe])^2
        d += (com1.y[iframe] - com2.y[iframe])^2
        d += (com1.z[iframe] - com2.z[iframe])^2
        d = sqrt(d)
        dist[iframe] = d
    end
    dist
end


