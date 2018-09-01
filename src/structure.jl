###### center of mass #################
"""
centerofmass

calculate center of mass of given trajectory
"""
function centerofmass(ta::TrjArray; isweight::Bool=true, index::Vector{Int64}=Vector{Int64}(undef, 0))::TrjArray
    nframe = ta.nframe
    natom = ta.natom
    if isempty(index)
        index = Colon()
        natom_sub = natom
    else
        index = index
        natom_sub = length(index)
    end
    if natom_sub == 1
        return ta[:, index]
    elseif isweight && length(ta.mass) == natom
        weight = reshape(ta.mass, 1, natom_sub)
        wsum_inv = 1.0 / sum(weight)
        x = sum(weight .* view(ta.x, :, index), dims=2) .* wsum_inv
        y = sum(weight .* view(ta.y, :, index), dims=2) .* wsum_inv
        z = sum(weight .* view(ta.z, :, index), dims=2) .* wsum_inv
        return TrjArray(x=x, y=y, z=z)
    else
        wsum_inv = 1.0 / Float64(natom_sub)
        x = sum(view(ta.x, :, index), dims=2) .* wsum_inv
        y = sum(view(ta.y, :, index), dims=2) .* wsum_inv
        z = sum(view(ta.z, :, index), dims=2) .* wsum_inv
        return TrjArray(x=x, y=y, z=z)
    end
end


"""
decenter

remove center of mass
"""
function decenter(ta::TrjArray; isweight::Bool=true, index::Vector{Int64}=Vector{Int64}(undef, 0))::Tuple{TrjArray, TrjArray}
    nframe = ta.nframe
    natom = ta.natom
    com = centerofmass(ta, isweight=isweight, index=index)
    TrjArray(ta.x .- com.x, ta.y .- com.y, ta.z .- com.z, ta), com
end


###### superimpose #################
function innerproduct!(ref_x::Vector{Float64}, ref_y::Vector{Float64}, ref_z::Vector{Float64},
                       x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64},
                       index::Vector{Int64}, A::Vector{Float64}, isweight::Bool, ta::TrjArray)::Float64
    A[:] .= 0.0
    G1 = G2 = 0.0
    if isweight
        for i in index
            x1 = ta.mass[i] * ref_x[i]
            y1 = ta.mass[i] * ref_y[i]
            z1 = ta.mass[i] * ref_z[i]
            G1 += x1 * ref_x[i] + y1 * ref_y[i] + z1 * ref_z[i]
            x2 = x[i]
            y2 = y[i]
            z2 = z[i]
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
            x1 = ref_x[i]
            y1 = ref_y[i]
            z1 = ref_z[i]
            G1 += x1 * ref_x[i] + y1 * ref_y[i] + z1 * ref_z[i]
            x2 = x[i]
            y2 = y[i]
            z2 = z[i]
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
    return (G1 + G2) * 0.5
end

function innerproduct!(iframe::Int64, ref::TrjArray, ta::TrjArray,
                       index::Vector{Int64}, A::Vector{Float64}, isweight::Bool)::Float64
    A[:] .= 0.0
    G1 = G2 = 0.0
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
    return (G1 + G2) * 0.5
end

"""
this code is licensed under the BSD license (Copyright (c) 2009-2016 Pu Liu and Douglas L. Theobald), see LICENSE.md
"""
function fastCalcRMSDAndRotation!(A::Vector{Float64}, E0::Float64, wsum_inv::Float64, rot::Vector{Float64}, C::Vector{Float64})::Float64
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

this code is licensed under the BSD license (Copyright (c) 2009-2016 Pu Liu and Douglas L. Theobald), see LICENSE.md
"""
function superimpose(ref::TrjArray, ta::TrjArray; isweight::Bool=true, index::Vector{Int64}=Vector{Int64}(undef, 0), isdecenter::Bool=false)::TrjArray
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
        wsum_inv = 1.0 / sum(weight2)
    else
        wsum_inv = 1.0 / Float64(length(index2))
    end

    if isdecenter
        ta2 = copy(ta)
        ref2 = ref[1, :]
    else
        ta2, = decenter(ta, isweight=isweight2, index=index)
        ref2, com = decenter(ref[1, :], isweight=isweight2, index=index)
    end

    x = Matrix{Float64}(undef, nframe, natom)
    y = Matrix{Float64}(undef, nframe, natom)
    z = Matrix{Float64}(undef, nframe, natom)
    rmsd = Vector{Float64}(undef, nframe)
    A = Vector{Float64}(undef, 9)
    rot = Vector{Float64}(undef, 9)
    C = Vector{Float64}(undef, 3)
    for iframe in 1:nframe
        E0 = innerproduct!(iframe, ref2, ta2, index2, A, isweight2)
        r = fastCalcRMSDAndRotation!(A, E0, wsum_inv, rot, C)
        rmsd[iframe] = r
        for iatom in 1:natom
            @inbounds x[iframe, iatom] = rot[1] * ta2.x[iframe, iatom] + rot[2] * ta2.y[iframe, iatom] + rot[3] * ta2.z[iframe, iatom]
            @inbounds y[iframe, iatom] = rot[4] * ta2.x[iframe, iatom] + rot[5] * ta2.y[iframe, iatom] + rot[6] * ta2.z[iframe, iatom]
            @inbounds z[iframe, iatom] = rot[7] * ta2.x[iframe, iatom] + rot[8] * ta2.y[iframe, iatom] + rot[9] * ta2.z[iframe, iatom]
        end
    end

    if !isdecenter
        x = x .+ com.x
        y = y .+ com.y
        z = z .+ com.z
    end
    ta_fit = TrjArray(x, y, z, ta)
    #rmsd, ta_fit
end

###### rmsd #################
"""
rmsd

root mean square deviation of ta::TrjArray from ref:TrjArray
"""
function calcrmsd(ref::TrjArray, ta::TrjArray; isweight::Bool=true, index::Vector{Int64}=Vector{Int64}(undef, 0))::Vector{Float64}
    nframe = ta.nframe
    natom = ta.natom
    if isweight && length(ta.mass) == natom
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

###### mean structure #################
"""
meanstructure

compute average structure by iterative superimposes
"""
function meanstructure(ta::TrjArray; isweight::Bool=true, index::Vector{Int64}=Vector{Int64}(undef, 0))::Tuple{TrjArray,TrjArray}
    nframe = ta.nframe
    natom = ta.natom

    ref = ta[1, :]
    rmsd = [1.0]
    tolerance = tolerance = 10^(-6)
    while rmsd[1] > tolerance
        ref_old = ref;
        ta = superimpose(ref, ta, isweight=isweight, index=index)
        ref = TrjArray(x=mean(ta.x, dims=1), y=mean(ta.y, dims=1), z=mean(ta.z, dims=1)) # TODO: mean(ta) should be available in the futre
        rmsd = calcrmsd(ref_old, ref, isweight=isweight, index=index)
        println("rmsd from the previous mean structure: ", rmsd[1])
    end

    mean_crd = ref
    ta = superimpose(mean_crd, ta, isweight=isweight, index=index)
    mean_crd, ta
end

###### distance, angle, dihedral #################
"""
calcdistance

calculate distances between two atoms or groups of atoms
"""
function calcbond(ta::TrjArray, index1::Vector{Int64}=Vector{Int64}(undef, 0), index2::Vector{Int64}=Vector{Int64}(undef, 0))::Vector{Float64}
    # TODO: support for PBC
    nframe = ta.nframe
    if isempty(index1)
        index1 = [1]
    end
    if isempty(index2)
        index2 = [2]
    end
    com1 = centerofmass(ta, isweight=true, index=index1)
    com2 = centerofmass(ta, isweight=true, index=index2)
    dist = sqrt.((com1.x .- com2.x).^2 .+ (com1.y .- com2.y).^2 .+ (com1.z .- com2.z).^2)
    reshape(dist, nframe)
end
