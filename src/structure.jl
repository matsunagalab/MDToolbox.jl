"""
centerofmass

calculate center of mass of given trajectory
"""
function centerofmass(ta::TrjArray{T, U};
    isweight::Bool=true, index::Vector{Int64}=Vector{Int64}(undef, 0))::TrjArray{T, U} where {T, U}
    nframe = ta.nframe
    natom = ta.natom
    if isempty(index)
        index2 = Colon()
        natom_sub = natom
    else
        index2 = index
        natom_sub = length(index2)
    end
    if natom_sub == 1
        return ta[:, index2]
    elseif isweight && length(ta.mass) == natom
        weight = reshape(ta.mass, 1, natom_sub)
        wsum_inv = one(T) / sum(weight)
        x = sum(weight .* view(ta.x, :, index2), dims=2) .* wsum_inv
        y = sum(weight .* view(ta.y, :, index2), dims=2) .* wsum_inv
        z = sum(weight .* view(ta.z, :, index2), dims=2) .* wsum_inv
        return TrjArray{T, U}(x=x, y=y, z=z)
    else
        wsum_inv = one(T) / Float64(natom_sub)
        x = sum(view(ta.x, :, index2), dims=2) .* wsum_inv
        y = sum(view(ta.y, :, index2), dims=2) .* wsum_inv
        z = sum(view(ta.z, :, index2), dims=2) .* wsum_inv
        return TrjArray{T, U}(x=x, y=y, z=z)
    end
end


############################################################################
"""
decenter

remove center of mass
"""
function decenter(ta::TrjArray{T, U};
    isweight::Bool=true, index::Vector{Int64}=Vector{Int64}(undef, 0))::Tuple{TrjArray{T, U}, TrjArray{T, U}} where {T, U}
    com = centerofmass(ta, isweight=isweight, index=index)
    TrjArray(ta.x .- com.x, ta.y .- com.y, ta.z .- com.z, ta), com
end

"""
decenter!

remove center of mass
"""
function decenter!(ta::TrjArray{T, U};
    isweight::Bool=true, index::Vector{Int64}=Vector{Int64}(undef, 0)) where {T, U}
    com = centerofmass(ta, isweight=isweight, index=index)
    ta.x .= ta.x .- com.x
    ta.y .= ta.y .- com.y
    ta.z .= ta.z .- com.z
    nothing
end


###### superimpose #################
function innerproduct(iframe::U, ref::TrjArray{T, U}, ta::TrjArray{T, U},
                       index::Vector{U}, isweight::Bool) where {T, U}
    A = zeros(T, 9)
    # A = Vector{Float64}(undef, 9)
    # A[:] .= 0.0
    G1 = G2 = zero(T)
    if isweight
        for i in index
            x1 = ta.mass[i] * ref.x[i]
            y1 = ta.mass[i] * ref.y[i]
            z1 = ta.mass[i] * ref.z[i]
            G1 += x1 * ref.x[i] + y1 * ref.y[i] + z1 * ref.z[i]
            x2 = ta.x[iframe, i]
            y2 = ta.y[iframe, i]
            z2 = ta.z[iframe, i]
            G2 += ta.mass[i] * (x2 * x2 + y2 * y2 + z2 * z2)
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
    else
        for i in index
            x1 = ref.x[i]
            y1 = ref.y[i]
            z1 = ref.z[i]
            G1 += x1 * ref.x[i] + y1 * ref.y[i] + z1 * ref.z[i]
            x2 = ta.x[iframe, i]
            y2 = ta.y[iframe, i]
            z2 = ta.z[iframe, i]
            G2 += (x2 * x2 + y2 * y2 + z2 * z2)
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
    end
    return A, (G1 + G2) * T(0.5)
end


############################################################################
"""
this code is licensed under the BSD license (Copyright (c) 2009-2016 Pu Liu and Douglas L. Theobald), see LICENSE.md
"""
function fastCalcRMSDAndRotation(A::Vector{T}, E0::T, wsum_inv::T) where {T}
    oldg = zero(T)
    C = zeros(T, 3)
    rot = zeros(T, 9)
    evecprec = T(1e-6)
    evalprec = T(1e-11)

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

    SyzSzymSyySzz2 = 2*(Syz*Szy - Syy*Szz)
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2

    C[3] = -2 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2)
    C[2] = 8 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz)

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
        delta = ((a*mxEigenV + C[1])/(T(2)*x2*mxEigenV + b + a));
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
    r1 = sqrt(abs(2 * (E0 - mxEigenV) * wsum_inv))
    # printf("\n\n %16g %16g %16g \n", r1, E0, 2.0 * (E0 - mxEigenV)/len)

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
                    rot[1] = rot[5] = rot[9] = one(T)
                    rot[2] = rot[3] = rot[4] = rot[6] = rot[7] = rot[8] = zero(T)
                    return r1, rot
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
    return r1, rot
end

function applyrotation!(iframe, x, y, z, ta2, rot)
    for iatom in 1:ta2.natom
        @inbounds x[iframe, iatom] = rot[1] * ta2.x[iframe, iatom] + rot[2] * ta2.y[iframe, iatom] + rot[3] * ta2.z[iframe, iatom]
        @inbounds y[iframe, iatom] = rot[4] * ta2.x[iframe, iatom] + rot[5] * ta2.y[iframe, iatom] + rot[6] * ta2.z[iframe, iatom]
        @inbounds z[iframe, iatom] = rot[7] * ta2.x[iframe, iatom] + rot[8] * ta2.y[iframe, iatom] + rot[9] * ta2.z[iframe, iatom]
    end
end

############################################################################
function rotate(ta::TrjArray{T, U}, quater::Vector{T})::TrjArray{T, U} where {T, U}
    x = zeros(T, ta.nframe, ta.natom)
    y = zeros(T, ta.nframe, ta.natom)
    z = zeros(T, ta.nframe, ta.natom)
    rot = zeros(T, 9)
    rot[1] = 1 - 2 * quater[2] * quater[2] - 2 * quater[3] * quater[3];
    rot[2] = 2 * (quater[1] * quater[2] + quater[3] * quater[4]);
    rot[3] = 2 * (quater[1] * quater[3] - quater[2] * quater[4]);
    rot[4] = 2 * (quater[1] * quater[2] - quater[3] * quater[4]);
    rot[5] = 1 - 2 * quater[1] * quater[1] - 2 * quater[3] * quater[3];
    rot[6] = 2 * (quater[2] * quater[3] + quater[1] * quater[4]);
    rot[7] = 2 * (quater[1] * quater[3] + quater[2] * quater[4]);
    rot[8] = 2 * (quater[2] * quater[3] - quater[1] * quater[4]);
    rot[9] = 1 - 2 * quater[1] * quater[1] - 2 * quater[2] * quater[2];
    #@show det(reshape(rot, (3, 3)))
    for iframe in 1:ta.nframe
        for iatom in 1:ta.natom
            @inbounds x[iframe, iatom] = rot[1] * ta.x[iframe, iatom] + rot[2] * ta.y[iframe, iatom] + rot[3] * ta.z[iframe, iatom]
            @inbounds y[iframe, iatom] = rot[4] * ta.x[iframe, iatom] + rot[5] * ta.y[iframe, iatom] + rot[6] * ta.z[iframe, iatom]
            @inbounds z[iframe, iatom] = rot[7] * ta.x[iframe, iatom] + rot[8] * ta.y[iframe, iatom] + rot[9] * ta.z[iframe, iatom]
        end
    end
    return TrjArray(x, y, z, ta)
end

function rotate!(ta::TrjArray{T, U}, quater::Vector{T}) where {T, U}
    rot = zeros(T, 9)
    rot[1] = 1 - 2 * quater[2] * quater[2] - 2 * quater[3] * quater[3];
    rot[2] = 2 * (quater[1] * quater[2] + quater[3] * quater[4]);
    rot[3] = 2 * (quater[1] * quater[3] - quater[2] * quater[4]);
    rot[4] = 2 * (quater[1] * quater[2] - quater[3] * quater[4]);
    rot[5] = 1 - 2 * quater[1] * quater[1] - 2 * quater[3] * quater[3];
    rot[6] = 2 * (quater[2] * quater[3] + quater[1] * quater[4]);
    rot[7] = 2 * (quater[1] * quater[3] + quater[2] * quater[4]);
    rot[8] = 2 * (quater[2] * quater[3] - quater[1] * quater[4]);
    rot[9] = 1 - 2 * quater[1] * quater[1] - 2 * quater[2] * quater[2];
    #@show det(reshape(rot, (3, 3)))
    for iframe in 1:ta.nframe
        for iatom in 1:ta.natom
            @inbounds x = rot[1] * ta.x[iframe, iatom] + rot[2] * ta.y[iframe, iatom] + rot[3] * ta.z[iframe, iatom]
            @inbounds y = rot[4] * ta.x[iframe, iatom] + rot[5] * ta.y[iframe, iatom] + rot[6] * ta.z[iframe, iatom]
            @inbounds z = rot[7] * ta.x[iframe, iatom] + rot[8] * ta.y[iframe, iatom] + rot[9] * ta.z[iframe, iatom]
            @inbounds ta.x[iframe, iatom] = x
            @inbounds ta.y[iframe, iatom] = y
            @inbounds ta.z[iframe, iatom] = z
        end
    end
end

############################################################################
"""
superimpose

superimpose ta to ref

this code is licensed under the BSD license (Copyright (c) 2009-2016 Pu Liu and Douglas L. Theobald), see LICENSE.md
"""
function superimpose(ref::TrjArray{T, U}, ta::TrjArray{T, U};
    isweight::Bool=true, index::Vector{U}=Vector{U}(undef, 0), isdecenter::Bool=false)::TrjArray{T, U} where {T, U}
    nframe = ta.nframe
    natom = ta.natom

    if isempty(index)
        index2 = collect(1:natom)
    else
        index2 = index
    end

    if isweight && length(ref.mass) == natom && length(ta.mass) == natom
        isweight2 = true
        weight = ta.mass
    else
        isweight2 = false
    end

    if isweight2
        weight2 = reshape(weight[index2], length(index2))
        wsum_inv = one(T) / sum(weight2)
    else
        wsum_inv = one(T) / T(length(index2))
    end

    if isdecenter
        ta2 = copy(ta)
        ref2 = ref[1, :]
    else
        ta2, = decenter(ta, isweight=isweight2, index=index)
        ref2, com = decenter(ref[1, :], isweight=isweight2, index=index)
    end

    x = Matrix{T}(undef, nframe, natom)
    y = Matrix{T}(undef, nframe, natom)
    z = Matrix{T}(undef, nframe, natom)
    Threads.@threads for iframe in 1:nframe
        A, E0 = innerproduct(iframe, ref2, ta2, index2, isweight2)
        rmsd, rot = fastCalcRMSDAndRotation(A, E0, wsum_inv)
        applyrotation!(iframe, x, y, z, ta2, rot)
    end

    if !isdecenter
        x = x .+ com.x
        y = y .+ com.y
        z = z .+ com.z
    end
    ta_fit = TrjArray(x, y, z, ta)
    #rmsd, ta_fit
end


############################################################################
function superimpose_serial(ref::TrjArray{T, U}, ta::TrjArray{T, U};
    isweight::Bool=true, index::Vector{U}=Vector{U}(undef, 0), isdecenter::Bool=false)::TrjArray{T, U} where {T, U}
    nframe = ta.nframe
    natom = ta.natom

    if isempty(index)
        index2 = collect(1:natom)
    else
        index2 = index
    end

    if isweight && length(ref.mass) == natom && length(ta.mass) == natom
        isweight2 = true
        weight = ta.mass
    else
        isweight2 = false
    end

    if isweight2
        weight2 = reshape(weight[index2], length(index2))
        wsum_inv = one(T) / sum(weight2)
    else
        wsum_inv = one(T) / T(length(index2))
    end

    if isdecenter
        ta2 = copy(ta)
        ref2 = ref[1, :]
    else
        ta2, = decenter(ta, isweight=isweight2, index=index)
        ref2, com = decenter(ref[1, :], isweight=isweight2, index=index)
    end

    x = Matrix{T}(undef, nframe, natom)
    y = Matrix{T}(undef, nframe, natom)
    z = Matrix{T}(undef, nframe, natom)
    for iframe in 1:nframe
        A, E0 = innerproduct(iframe, ref2, ta2, index2, isweight2)
        rmsd, rot = fastCalcRMSDAndRotation(A, E0, wsum_inv)
        applyrotation!(iframe, x, y, z, ta2, rot)
    end

    if !isdecenter
        x = x .+ com.x
        y = y .+ com.y
        z = z .+ com.z
    end
    ta_fit = TrjArray(x, y, z, ta)
    #rmsd, ta_fit
end


############################################################################
"""
calcrmsd

rmsd (root mean square deviation)
"""
function getrmsd(ref::TrjArray{T, U}, ta::TrjArray{T, U};
    isweight::Bool=true, index::Vector{U}=Vector{U}(undef, 0))::Vector{T} where {T, U}
    nframe = ta.nframe
    natom = ta.natom
    if isweight && length(ta.mass) == natom
        weight = ta.mass
    else
        weight = ones(T, natom)
    end
    if isempty(index)
        index2 = cumsum(ones(U, natom))
    else
        index2 = index
    end
    wsum_inv = one(T) / sum(weight[index2])
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
    reshape(d, nframe)
end


############################################################################
"""
meanstructure

compute average structure by iterative superimposes
"""
function meanstructure(ta::TrjArray{T, U};
    isweight::Bool=true, index::Vector{U}=Vector{U}(undef, 0))::Tuple{TrjArray{T, U},TrjArray{T, U}} where {T, U}
    nframe = ta.nframe
    natom = ta.natom
    ta2 = ta
    ref = ta2[1, :]

    r = [one(T)]
    tolerance = tolerance = T(10^(-6))
    while r[1] > tolerance
        ref_old = ref;
        ta2 = superimpose(ref, ta2, isweight=isweight, index=index)
        ref = TrjArray{T, U}(x=mean(ta2.x, dims=1), y=mean(ta2.y, dims=1), z=mean(ta2.z, dims=1)) # TODO: mean(ta) should be available in the futre
        r .= getrmsd(ref_old, ref, isweight=isweight, index=index)
        println("rmsd from the previous mean structure: ", r[1])
    end

    ta2 = superimpose(ref, ta2, isweight=isweight, index=index)
    ref, ta2
end


############################################################################
"""
calcrmsf

rmsf (root mean square fluctuation)
"""
function getrmsf(ta::TrjArray{T, U};
    isweight::Bool=true)::Vector{T} where {T, U}
    nframe = ta.nframe
    natom = ta.natom
    if isweight && length(ta.mass) == natom
        weight = ta.mass
    else
        weight = ones(T, natom)
    end
    wsum_inv = 1.0 / sum(weight[index2])
    weight2 = reshape(weight[index2], 1, length(index2))
    mean_x = (ta.x .* weight2)

    ref_x = ref.x[1:1, index2]
    ref_y = ref.y[1:1, index2]
    ref_z = ref.z[1:1, index2]
    ta_x = ta.x[:, index2]
    ta_y = ta.y[:, index2]
    ta_z = ta.z[:, index2]
    d =  sum(weight2 .* ((ta_x .- ref_x).^2 .+ (ta_y .- ref_y).^2 .+ (ta_z .- ref_z).^2), dims=2)
    d = d .* wsum_inv
    d = sqrt.(d)
    reshape(d, nframe)
end


############################################################################
"""
calcbond

distance between two atoms or groups of atoms
"""
function getdistance(ta1::TrjArray{T, U}, ta2::TrjArray{T, U})::Vector{T} where {T, U}
    # TODO: support for PBC
    # TODO: hypot
    nframe = ta1.nframe
    com1 = centerofmass(ta1, isweight=true)
    com2 = centerofmass(ta2, isweight=true)
    dist = sqrt.((com1.x .- com2.x).^2 .+ (com1.y .- com2.y).^2 .+ (com1.z .- com2.z).^2)
    reshape(dist, nframe)
end

############################################################################
"""
calcangle

angle of three atoms or groups of atoms
"""
function getangle(ta1::TrjArray{T, U}, ta2::TrjArray{T, U}, ta3::TrjArray{T, U})::Vector{T} where {T, U}
    nframe = ta1.nframe
    com1 = centerofmass(ta1, isweight=true)
    com2 = centerofmass(ta2, isweight=true)
    com3 = centerofmass(ta3, isweight=true)
    a = zeros(T, nframe)
    for iframe in 1:nframe
        d1 = [com1.x[iframe] - com2.x[iframe]; com1.y[iframe] - com2.y[iframe]; com1.z[iframe] - com2.z[iframe]]
        d2 = [com3.x[iframe] - com2.x[iframe]; com3.y[iframe] - com2.y[iframe]; com3.z[iframe] - com2.z[iframe]]
        a[iframe] = acos(dot(d1, d2)/(norm(d1)*norm(d2)))
    end
    a = (a ./ pi) .* T(180)
end


############################################################################
"""
calcdihedral

dihedral of four atoms or groups of atoms
"""
function getdihedral(ta1::TrjArray{T, U}, ta2::TrjArray{T, U}, ta3::TrjArray{T, U}, ta4::TrjArray{T, U})::Vector{T} where {T, U}
    nframe = ta1.nframe
    com1 = centerofmass(ta1, isweight=true)
    com2 = centerofmass(ta2, isweight=true)
    com3 = centerofmass(ta3, isweight=true)
    com4 = centerofmass(ta4, isweight=true)
    a = zeros(T, nframe)
    # Threads.@threads for iframe in 1:nframe
    for iframe in 1:nframe
        d1 = [com1.x[iframe] - com2.x[iframe]; com1.y[iframe] - com2.y[iframe]; com1.z[iframe] - com2.z[iframe]]
        d2 = [com3.x[iframe] - com2.x[iframe]; com3.y[iframe] - com2.y[iframe]; com3.z[iframe] - com2.z[iframe]]
        d3 = [com3.x[iframe] - com4.x[iframe]; com3.y[iframe] - com4.y[iframe]; com3.z[iframe] - com4.z[iframe]]
        m1 = cross(d1, d2)
        m2 = cross(d2, d3)
        a[iframe] = acos(dot(m1, m2)/(norm(m1)*norm(m2)))
        rotdirection = dot(d2,cross(m1,m2))
        if rotdirection < zero(T)
            a[iframe] = -a[iframe]
        end
    end
    a = (a ./ pi) .* T(180)
end
