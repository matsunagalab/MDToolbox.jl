"""
    centerofmass(ta::TrjArray; isweight=true, index=[]) -> com::TrjArray

Calculates the center of mass coordinates (COM) of the given trajectory, a `TrjArray` object `ta`. 
If `isweight` is `true` (default), coordinates are weighted by `ta.mass` (as long as `ta.mass` is not empty).
Uses can also specify atom IDs for the COM calculation by giving `index` which should be an Array of integers. 

Returns the center of mass coordinates as a TrjArray object, which have virtual single "atom" whose coordinates are COMs. 

# Example
```julia-repl
julia> ta = mdload("ak.pdb")
julia> com = centersofmass(ta)
```
"""
function centerofmass(ta::TrjArray{T, U};
    isweight::Bool=true, index::AbstractVector=Vector{Int64}(undef, 0))::TrjArray{T, U} where {T, U}
    nframe = ta.nframe
    natom = ta.natom
    if isempty(index)
        index = 1:natom
        index3 = to3(1:natom)
        natom_sub = natom
    else
        index3 = to3(index)
        if eltype(index3) == Bool
            natom_sub = sum(index)
        else
            natom_sub = length(index)
        end
    end

    xyz = similar(ta.xyz, nframe, 3)
    if natom_sub == 1
        return ta[:, index]
    elseif isweight && length(ta.mass) == natom
        weight = reshape(ta.mass, 1, natom_sub)
        wsum_inv = one(T) / sum(weight)
        xyz[:, 1:1] .= sum(weight .* view(ta.xyz, :, index3[1:3:end]), dims=2) .* wsum_inv
        xyz[:, 2:2] .= sum(weight .* view(ta.xyz, :, index3[2:3:end]), dims=2) .* wsum_inv
        xyz[:, 3:3] .= sum(weight .* view(ta.xyz, :, index3[3:3:end]), dims=2) .* wsum_inv
    else
        wsum_inv = one(T) / T(natom_sub)
        xyz[:, 1:1] .= sum(view(ta.xyz, :, index3[1:3:end]), dims=2) .* wsum_inv
        xyz[:, 2:2] .= sum(view(ta.xyz, :, index3[2:3:end]), dims=2) .* wsum_inv
        xyz[:, 3:3] .= sum(view(ta.xyz, :, index3[3:3:end]), dims=2) .* wsum_inv
    end

    TrjArray{T, U}(xyz=xyz)
end

"""
    decenter(ta::TrjArray; isweight=true, index=[]) -> ta_decentered::TrjArray, com::TrjArray

Calculates the centers of mass coordinates (COM) of the given trajectory, a TrjArray object ta.
And translates the coordinates so that the COM of the output is identical to the origin (x=y=z=0). 
If `isweight` is `true` (default), coordinates are weighted by `ta.mass` (as long as `ta.mass` is not empty).
Uses can also specify atom IDs for the COM calculation by giving `index` which should be an Array of integers. 

Returns a TrjArray object whose COM is the origin (x=y=z=0), and the COM of the given TrjArray `ta`. 

# Example
```julia-repl
julia> ta = mdload("ak.pdb")
julia> ta_decentered = decenter(ta)
```
"""
function decenter(ta::TrjArray{T, U}; isweight::Bool=true, index::Vector{Int64}=Vector{Int64}(undef, 0))::Tuple{TrjArray{T, U}, TrjArray{T, U}} where {T, U}
    com = centerofmass(ta, isweight=isweight, index=index)
    xyz = similar(ta.xyz)
    xyz[:, 1:3:end] .= ta.xyz[:, 1:3:end] .- com.xyz[:, 1:1]
    xyz[:, 2:3:end] .= ta.xyz[:, 2:3:end] .- com.xyz[:, 2:2]
    xyz[:, 3:3:end] .= ta.xyz[:, 3:3:end] .- com.xyz[:, 3:3]
    TrjArray(ta, xyz=xyz), com
end

"""
    decenter!(ta::TrjArray; isweight=true, index=[])

Translates the coordinates of the given trajectory, a TrjArray `ta`, 
so that the COM of it become identical to the origin (x=y=z=0). 
If `isweight` is `true` (default), coordinates are weighted by `ta.mass` (as long as `ta.mass` is not empty).
Uses can also specify atom IDs for the COM calculation by giving `index` which should be an Array of integers. 

# Example
```julia-repl
julia> ta = mdload("ak.pdb")
julia> decenter!(ta)
```
"""
function decenter!(ta::TrjArray{T, U}; isweight::Bool=true, index::Vector{Int64}=Vector{Int64}(undef, 0)) where {T, U}
    com = centerofmass(ta, isweight=isweight, index=index)
    ta.xyz[:, 1:3:end] .= ta.xyz[:, 1:3:end] .- com.x[:, 1]
    ta.xyz[:, 2:3:end] .= ta.xyz[:, 2:3:end] .- com.y[:, 2]
    ta.xyz[:, 3:3:end] .= ta.xyz[:, 3:3:end] .- com.z[:, 3]
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
            ref_x = ref.xyz[1, 3*(i-1)+1]
            ref_y = ref.xyz[1, 3*(i-1)+2]
            ref_z = ref.xyz[1, 3*(i-1)+3]
            x1 = ta.mass[i] * ref_x
            y1 = ta.mass[i] * ref_y
            z1 = ta.mass[i] * ref_z
            G1 += x1 * ref_x + y1 * ref_x + z1 * ref_z
            x2 = ta.xyz[iframe, 3*(i-1)+1]
            y2 = ta.xyz[iframe, 3*(i-1)+2]
            z2 = ta.xyz[iframe, 3*(i-1)+3]
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
            ref_x = ref.xyz[1, 3*(i-1)+1]
            ref_y = ref.xyz[1, 3*(i-1)+2]
            ref_z = ref.xyz[1, 3*(i-1)+3]
            x1 = ref_x
            y1 = ref_y
            z1 = ref_z
            G1 += x1 * ref_x + y1 * ref_y + z1 * ref_z
            x2 = ta.xyz[iframe, 3*(i-1)+1]
            y2 = ta.xyz[iframe, 3*(i-1)+2]
            z2 = ta.xyz[iframe, 3*(i-1)+3]
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
    for i in 1:ta2.natom
        @inbounds x[iframe, i] = rot[1] * ta2.xyz[iframe, 3*(i-1)+1] + rot[2] * ta2.xyz[iframe, 3*(i-1)+2] + rot[3] * ta2.xyz[iframe, 3*(i-1)+3]
        @inbounds y[iframe, i] = rot[4] * ta2.xyz[iframe, 3*(i-1)+1] + rot[5] * ta2.xyz[iframe, 3*(i-1)+2] + rot[6] * ta2.xyz[iframe, 3*(i-1)+3]
        @inbounds z[iframe, i] = rot[7] * ta2.xyz[iframe, 3*(i-1)+1] + rot[8] * ta2.xyz[iframe, 3*(i-1)+2] + rot[9] * ta2.xyz[iframe, 3*(i-1)+3]
    end
end

"""
    rotate(ta::TrjArray, quaternion::Vector) -> ta_rotated::TrjArray

Performs rotation of coordinates of the given TrjArray object `ta`. 
Rotation matrix is constructed fromt the given `quaternion` vector. 

Returns the trajecotry with rotated coordinates.

# Example
```julia-repl
julia> ta = mdload("ak.dcd")
julia> quaternion = [0.80951112,  0.10657412,  0.35146912,  0.45804312]
julia> ta_rotated = rotate(ta, quaternion)
```
"""
function rotate(ta::TrjArray{T, U}, quater::AbstractVector{T})::TrjArray{T, U} where {T, U}
    xyz = similar(ta.xyz)
    rot = similar(quater, 9)
    rot[1] = 1 - 2 * quater[2] * quater[2] - 2 * quater[3] * quater[3]
    rot[2] = 2 * (quater[1] * quater[2] + quater[3] * quater[4])
    rot[3] = 2 * (quater[1] * quater[3] - quater[2] * quater[4])
    rot[4] = 2 * (quater[1] * quater[2] - quater[3] * quater[4])
    rot[5] = 1 - 2 * quater[1] * quater[1] - 2 * quater[3] * quater[3]
    rot[6] = 2 * (quater[2] * quater[3] + quater[1] * quater[4])
    rot[7] = 2 * (quater[1] * quater[3] + quater[2] * quater[4])
    rot[8] = 2 * (quater[2] * quater[3] - quater[1] * quater[4])
    rot[9] = 1 - 2 * quater[1] * quater[1] - 2 * quater[2] * quater[2]
    x = rot[1] .* ta.xyz[:, 1:3:end] .+ rot[2] .* ta.xyz[:, 2:3:end] .+ rot[3] .* ta.xyz[:, 3:3:end]
    y = rot[4] .* ta.xyz[:, 1:3:end] .+ rot[5] .* ta.xyz[:, 2:3:end] .+ rot[6] .* ta.xyz[:, 3:3:end]
    z = rot[7] .* ta.xyz[:, 1:3:end] .+ rot[8] .* ta.xyz[:, 2:3:end] .+ rot[9] .* ta.xyz[:, 3:3:end]
    xyz[:, 1:3:end] .= x
    xyz[:, 2:3:end] .= y
    xyz[:, 3:3:end] .= z
    return TrjArray(ta, xyz=xyz)
end

"""
    rotate!(ta::TrjArray, quaternion::Vector)

Performs rotation of coordinates of the given TrjArray object `ta`. 
Rotation matrix is constructed fromt the given `quaternion` vector. 

The coordinates of `ta` is replaced with the rotated ones. 

# Example
```julia-repl
julia> ta = mdload("ak.dcd")
julia> quaternion = [0.80951112,  0.10657412,  0.35146912,  0.45804312]
julia> rotate!(ta, quaternion)
```
"""
function rotate!(ta::TrjArray{T, U}, quater::AbstractVector{T}) where {T, U}
    rot = similar(quater, 9)
    rot[1] = 1 - 2 * quater[2] * quater[2] - 2 * quater[3] * quater[3]
    rot[2] = 2 * (quater[1] * quater[2] + quater[3] * quater[4])
    rot[3] = 2 * (quater[1] * quater[3] - quater[2] * quater[4])
    rot[4] = 2 * (quater[1] * quater[2] - quater[3] * quater[4])
    rot[5] = 1 - 2 * quater[1] * quater[1] - 2 * quater[3] * quater[3]
    rot[6] = 2 * (quater[2] * quater[3] + quater[1] * quater[4])
    rot[7] = 2 * (quater[1] * quater[3] + quater[2] * quater[4])
    rot[8] = 2 * (quater[2] * quater[3] - quater[1] * quater[4])
    rot[9] = 1 - 2 * quater[1] * quater[1] - 2 * quater[2] * quater[2]
    x = rot[1] .* ta.xyz[:, 1:3:end] .+ rot[2] .* ta.xyz[:, 2:3:end] .+ rot[3] .* ta.xyz[:, 3:3:end]
    y = rot[4] .* ta.xyz[:, 1:3:end] .+ rot[5] .* ta.xyz[:, 2:3:end] .+ rot[6] .* ta.xyz[:, 3:3:end]
    z = rot[7] .* ta.xyz[:, 1:3:end] .+ rot[8] .* ta.xyz[:, 2:3:end] .+ rot[9] .* ta.xyz[:, 3:3:end]
    ta.xyz[:, 1:3:end] .= x
    ta.xyz[:, 2:3:end] .= y
    ta.xyz[:, 3:3:end] .= z
    return nothing
end

"""
    rotate(ta::TrjArray, quaternions::Matrix) -> ta_rotated::TrjArray

Performs multiple rotations of the coordinates of the 1st frame of the given TrjArray object `ta`. 
Rotation matrices are construcuted from the 2-dimensional Array `quaternions`.
One set of quaternions for rotaton should be given in a row vector in `quaternions`. 

Returns the trajecotry with rotated coordinates from the 1st frame of `ta`. 

# Example
```julia-repl
julia> ta = mdload("ak.pdb")
julia> quaternions = [0.80951112  0.10657412  0.35146912  0.45804312; 
>                     0.10657412  0.80951112  -0.35146912  0.45804312]
julia> ta_rotated = rotate(ta, quaternions)
```
"""
function rotate(ta_single::TrjArray{T, U}, quater::AbstractMatrix{T})::TrjArray{T, U} where {T, U}
    rot = similar(quater, size(quater, 1), 9)
    rot[:, 1:1] .= 1 .- 2 .* quater[:, 2:2] .* quater[:, 2:2] .- 2.0 .* quater[:, 3:3] .* quater[:, 3:3]
    rot[:, 2:2] .= 2 .* (quater[:, 1:1] .* quater[:, 2:2] .+ quater[:, 3:3] .* quater[:, 4:4])
    rot[:, 3:3] .= 2 .* (quater[:, 1:1] .* quater[:, 3:3] .- quater[:, 2:2] .* quater[:, 4:4])
    rot[:, 4:4] .= 2 .* (quater[:, 1:1] .* quater[:, 2:2] .- quater[:, 3:3] .* quater[:, 4:4])
    rot[:, 5:5] .= 1 .- 2 .* quater[:, 1:1] .* quater[:, 1:1] .- 2 .* quater[:, 3:3] .* quater[:, 3:3]
    rot[:, 6:6] .= 2 .* (quater[:, 2:2] .* quater[:, 3:3] .+ quater[:, 1:1] .* quater[:, 4:4])
    rot[:, 7:7] .= 2 .* (quater[:, 1:1] .* quater[:, 3:3] .+ quater[:, 2:2] .* quater[:, 4:4])
    rot[:, 8:8] .= 2 .* (quater[:, 2:2] .* quater[:, 3:3] .- quater[:, 1:1] .* quater[:, 4:4])
    rot[:, 9:9] .= 1 .- 2 .* quater[:, 1:1] .* quater[:, 1:1] .- 2 .* quater[:, 2:2] .* quater[:, 2:2]
    x = rot[:, 1:1] .* ta_single.xyz[1:1, 1:3:end] .+ rot[:, 2:2] .* ta_single.xyz[1:1, 2:3:end] .+ rot[:, 3:3] .* ta_single.xyz[1:1, 3:3:end]
    y = rot[:, 4:4] .* ta_single.xyz[1:1, 1:3:end] .+ rot[:, 5:5] .* ta_single.xyz[1:1, 2:3:end] .+ rot[:, 6:6] .* ta_single.xyz[1:1, 3:3:end]
    z = rot[:, 7:7] .* ta_single.xyz[1:1, 1:3:end] .+ rot[:, 8:8] .* ta_single.xyz[1:1, 2:3:end] .+ rot[:, 9:9] .* ta_single.xyz[1:1, 3:3:end]
    xyz[:, 1:3:end] .= x
    xyz[:, 2:3:end] .= y
    xyz[:, 3:3:end] .= z
    return TrjArray(ta_single, xyz=xyz)
end

"""
    rotate_with_matrix(ta::TrjArray, R::Matrix) -> ta_rotated::TrjArray

Performs rotation of coordinates of the given TrjArray object `ta` with the given rotation matrix. 

The coordinates of `ta` is replaced with the rotated ones. 

# Example
```julia-repl
julia> ta = mdload("ak.dcd")
julia> R = [0.36 0.48 -0.8;
            -0.8 0.60 0.0;
            -.48 0.64 0.60]
julia> ra_rotated = rotate(ta, R)
```
"""
function rotate_with_matrix(ta::TrjArray{T, U}, R::AbstractMatrix{T})::TrjArray{T, U} where {T, U}
    x = R[1, 1].*ta.xyz[:, 1:3:end] .+ R[1, 2].*ta.xyz[:, 2:3:end] .+ R[1, 3].*ta.xyz[:, 3:3:end]
    y = R[2, 1].*ta.xyz[:, 1:3:end] .+ R[2, 2].*ta.xyz[:, 2:3:end] .+ R[2, 3].*ta.xyz[:, 3:3:end]
    z = R[3, 1].*ta.xyz[:, 1:3:end] .+ R[3, 2].*ta.xyz[:, 2:3:end] .+ R[3, 3].*ta.xyz[:, 3:3:end]
    ta.xyz[:, 1:3:end] .= x
    ta.xyz[:, 2:3:end] .= y
    ta.xyz[:, 3:3:end] .= z
    return TrjArray(ta, xyz=xyz)
end

"""
    superimpose(ref::TrjArray, ta::TrjArray; isweight=true, index=[], isdecenter=false) -> ta_superimposed::TrjArray

Performs the least-squares fitting of the given trajectory (TrjArray object `ta`) 
to the reference structure (1st frame of TrjArray object `ref`). 
If `isweight` is `true` (default), coordinates are weighted by `ta.mass` (as long as `ta.mass` is not empty).
Uses can specify atom IDs for the COM calculation by giving `index` which should be an Array of integers. 
When `isdecenter` is `true`, the function assumes the COMs of both structures are located at the origin (default is `isdecenter=false`). 

Returns the superimposed trajecotry. 

The code of this function is licensed under the BSD license (Copyright (c) 2009-2016 Pu Liu and Douglas L. Theobald), see LICENSE.md

# Example
```julia-repl
julia> ref = mdload("ak.pdb")
julia> ta = mdload("ak.dcd", top=pdb)
julia> ta_superimposed = superimpose(ref, ta)
```
# References
```
The algorithm of this function is based on 
P. Liu, D. K. Agrafiotis, and D. L. Theobald, J. Comput. Chem. 31, 1561-1563 (2010).
```
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
        x = x .+ com.xyz[:, 1]
        y = y .+ com.xyz[:, 2]
        z = z .+ com.xyz[:, 3]
    end

    xyz = Matrix{T}(undef, nframe, natom*3)
    xyz[:, 1:3:end] .= x
    xyz[:, 2:3:end] .= y
    xyz[:, 3:3:end] .= z

    return TrjArray(ta, xyz=xyz)
end

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

    xyz = Matrix{T}(undef, nframe, natom*3)
    xyz[:, 1:3:end] .= x
    xyz[:, 2:3:end] .= y
    xyz[:, 3:3:end] .= z

    return TrjArray(ta, xyz=xyz)
end

"""
    compute_rmsd(ref::TrjArray, ta::TrjArray; isweight=true, index=[]) -> rmsd

Calculates the root mean square deviations (RMSD) of the given trajectory (TrjArray object `ta`)
from the reference structure (1st frame of TrjArray object `ref`). 
If `isweight` is `true` (default), coordinates are weighted by `ta.mass` (as long as `ta.mass` is not empty).
Uses can specify atom IDs for the RMSD calculation by giving `index` which should be an Array of integers. 

Returns RMSD values. 

# Example
```julia-repl
julia> ref = mdload("ak.pdb")
julia> ta = mdload("ak.nc", top=pdb)
julia> ta_superimposed = superimpose(ref, ta)
julia> rmsd = compute_rmsd(ref, ta_superimposed)
```
"""
function compute_rmsd(ref::TrjArray{T, U}, ta::TrjArray{T, U};
    isweight::Bool=true, index::Vector{U}=Vector{U}(undef, 0))::Vector{T} where {T, U}
    nframe = ta.nframe
    natom = ta.natom
    if isweight && length(ta.mass) == natom
        weight = ta.mass
    else
        weight = ones(T, natom)
    end
    if isempty(index)
        index2 = 1:ta.natom
        index3 = to3(index2)
    else
        index2 = index
        index3 = to3(index2)
    end
    wsum_inv = one(T) / sum(weight[index2])
    weight2 = reshape(weight[index2], 1, length(index2))
    ref_x = view(ref.xyz, 1:1, index3[1:3:end])
    ref_y = view(ref.xyz, 1:1, index3[2:3:end])
    ref_z = view(ref.xyz, 1:1, index3[3:3:end])
    ta_x = view(ta.xyz, :, index3[1:3:end])
    ta_y = view(ta.xyz, :, index3[2:3:end])
    ta_z = view(ta.xyz, :, index3[3:3:end])
    d =  sum(weight2 .* ((ta_x .- ref_x).^2 .+ (ta_y .- ref_y).^2 .+ (ta_z .- ref_z).^2), dims=2)
    d = d .* wsum_inv
    d = sqrt.(d)
    reshape(d, nframe)
end

"""
    meanstructure(ta::TrjArray; isweight=true, index=[]) -> ta_mean::TrjArray, ta_superimposed::TrjArray

Calculates the mean structure of the given trajectory (TrjArray object `ta`)
by iteratively fitting the trajecotry to tentative mean structures until the mean structure converges. 
If `isweight` is `true` (default), coordinates are weighted by `ta.mass` (as long as `ta.mass` is not empty).
Uses can specify atom IDs for the fitting (superimpose function) by giving `index` which should be an Array of integers. 

Returns the mean structure as TrjArray object, and the supserimposed trajectory to the mean structure. 

# Example
```julia-repl
julia> ta = mdload("ak.nc", top=pdb)
julia> ta_mean, ta_superimposed = meanstructure(ta)
```
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
        r .= compute_rmsd(ref_old, ref, isweight=isweight, index=index)
        println("rmsd from the previous mean structure: ", r[1])
    end

    ta2 = superimpose(ref, ta2, isweight=isweight, index=index)
    ref, ta2
end

"""
    compute_rmsf(ta::TrjArray{T, U}; isweight::Bool=true) -> rmsf

rmsf (root mean square fluctuation)

Calculates the root mean square fluctuations (RMSFs) of the atoms in the given trajectory (TrjArray object `ta`). 
from the average coordinates of the trajectory. 
If `isweight` is `true` (default), coordinates are weighted by `ta.mass` (as long as `ta.mass` is not empty).

Returns RMSF values. 

# Example
```julia-repl
julia> ta = mdload("ak.nc", top=pdb)
julia> ta_mean, ta_superimposed = meanstructure(ta)
julia> rmsf = compute_rmsf(ta)
```
"""
function compute_rmsf(ta::TrjArray{T, U}; isweight::Bool=true)::Vector{T} where {T, U}
    nframe = ta.nframe
    natom = ta.natom
    if isweight && length(ta.mass) == natom
        weight = ta.mass
    else
        weight = ones(T, natom)
    end
    wsum_inv = 1.0 / sum(weight)
    weight_row = reshape(weight, 1, length(weight))

    ref_x = sum(ta.xyz[1:1, 1:3:end] .* weight_row * wsum_inv, dims=1)
    ref_y = sum(ta.xyz[1:1, 2:3:end] .* weight_row * wsum_inv, dims=1)
    ref_z = sum(ta.xyz[1:1, 3:3:end] .* weight_row * wsum_inv, dims=1)
    ta_x = ta.xyz[:, 1:3:end]
    ta_y = ta.xyz[:, 2:3:end]
    ta_z = ta.xyz[:, 3:3:end]
    d =  sum(weight_row .* ((ta_x .- ref_x).^2 .+ (ta_y .- ref_y).^2 .+ (ta_z .- ref_z).^2), dims=2)
    d = d .* wsum_inv
    d = sqrt.(d)
    rmsf = mean(d, dims=1)
    reshape(rmsf, natom)
end

"""
    compute_distance(ta::TrjArray, index::Matrix)

Calculates distances between two atom pairs specified by the Matrix object `index`. 
Each row vector in `index` contains two atom IDs for calculating distance. 
Default of `index` is `[1 2]`. 

Returns distances specified pairs in `index`. 

# Example
```julia-repl
julia> ta = mdload("ak.dcd")
julia> d = compute_distance(ta, [1 2; 1; 3])
```
"""
function compute_distance(ta::TrjArray{T, U}, index=[1 2]::Matrix{U})::Matrix{T} where {T, U}
    xyz1 = view(ta.xyz, :, to3(index[:, 1]))
    xyz2 = view(ta.xyz, :, to3(index[:, 2]))
    dx = xyz1[:, 1:3:end] .- xyz2[:, 1:3:end]
    dy = xyz1[:, 2:3:end] .- xyz2[:, 2:3:end]
    dz = xyz1[:, 3:3:end] .- xyz2[:, 3:3:end]
    if !isempty(ta.boxsize)
        wrap_dx!(dx, ta.boxsize[:, 1])
        wrap_dx!(dy, ta.boxsize[:, 2])
        wrap_dx!(dz, ta.boxsize[:, 3])
    end
    d = sqrt.(dx.^2 .+ dy.^2 .+ dz.^2)
    return d
end

"""
    compute_distance(ta1::TrjArray, ta2::TrjArray, index::Matrix)

Calculates distances between two atom pairs in `ta1` and `ta2` specified by the Matrix object `index`. 
Each row vector in `index` contains two atom id in `ta1` and atom id in `ta2` for calculating distance. 
Default of `index` is `[1 1]`. 

Returns distances specified pairs in `index`. 

# Example
```julia-repl
julia> ta = mdload("ak.pdb")
julia> ta = mdload("ak.dcd", top=ta)
julia> d = compute_distance(ta[:, "atomid 1"], ta[:, "atomid 9"])
```
"""
function compute_distance(ta1::TrjArray{T, U}, ta2::TrjArray{T, U}, index=[1 1]::Matrix{U})::Matrix{T} where {T, U}
    xyz1 = view(ta1.xyz, :, to3(index[:, 1]))
    xyz2 = view(ta2.xyz, :, to3(index[:, 2]))
    dx = xyz1[:, 1:3:end] .- xyz2[:, 1:3:end]
    dy = xyz1[:, 2:3:end] .- xyz2[:, 2:3:end]
    dz = xyz1[:, 3:3:end] .- xyz2[:, 3:3:end]
    if !isempty(ta1.boxsize)
        wrap_dx!(dx, ta1.boxsize[:, 1])
        wrap_dx!(dy, ta1.boxsize[:, 2])
        wrap_dx!(dz, ta1.boxsize[:, 3])
    end
    d = sqrt.(dx.^2 .+ dy.^2 .+ dz.^2)
    return d
end

"""
    compute_distancemap(ta::TrjArray; kneighbor=3)

Calculates distance-map vectors from the given TrjArray `ta`. 
If `atomid=j` and `atomid=i` are closer to each other than the specified `kneighbor`, i.e., `abs(i-j) < kneighbor`, 
the calculation of the distance is ignored, considering that they are not nonbonded pairs. 
The default value of `kneighbor` is `3`. 

Returns distance-map vectors of all frames. The distance-map of each frame is contained in a row vector of the output.

# Example
```julia-repl
julia> ta = mdload("ak.pdb")
julia> ta = mdload("ak.dcd", top=ta)
julia> d = compute_distancemap(ta["atomname CA"])
```
"""
function compute_distancemap(ta::TrjArray{T, U}; kneighbor=3)::Matrix{T} where {T, U}
    natom = ta.natom
    nframe = ta.nframe
    count = 0
    for i = 1:natom
        for j = (i+kneighbor):natom
            count += 1
        end
    end
    index_all = zeros(U, count, 2)
    count = 0
    for i = 1:natom
        for j = (i+kneighbor):natom
            count += 1
            index_all[count, 1] = i
            index_all[count, 2] = j
        end
    end
    d = compute_distance(ta, index_all)
    return d
end

"""
    compute_contactmap(ta::TrjArray; rcut=8.5, kneighbor=3)

Calculates contact-map vectors from the given TrjArray `ta`. 
Cut-off value for contact can be given in `rcut` (default is `rcut=8.5`). 
If `atomid=j` and `atomid=i` are closer to each other than the specified `kneighbor`, i.e., `abs(i-j) < kneighbor`, 
the calculation of the distance is ignored, considering that they are not nonbonded pairs. 
The default value of `kneighbor` is `3`. 

Returns contact-map vectors of all frames. The contact-map of each frame is contained in a row vector of the output. 

# Example
```julia-repl
julia> ta = mdload("ak.pdb")
julia> ta = mdload("ak.dcd", top=ta)
julia> d = compute_distancemap(ta["atomname CA"])
```
"""
function compute_contactmap(ta::TrjArray{T, U}; rcut=8.5, kneighbor=3)::Matrix{T} where {T, U}
    d = compute_distancemap(ta, kneighbor=kneighbor)
    c = T.(d .< rcut)
    return c
end

"""
    compute_angle(ta1::TrjArray, ta2::TrjArray, ta3::TrjArray)

Calculates angles for the atom coordinates or the COMs of groups from the triplet of `ta1`, `ta2` and `ta3`. 

Returns angles.

# Example
```julia-repl
julia> ta = mdload("ak.pdb")
julia> ta = mdload("ak.dcd", top=ta)
julia> a = compute_distance(ta[:, "atomid 1"], ta[:, "atomid 9"], ta[:, "atomid 11"])
```
"""
function compute_angle(ta1::TrjArray{T, U}, ta2::TrjArray{T, U}, ta3::TrjArray{T, U})::Vector{T} where {T, U}
    nframe = ta1.nframe
    com1 = centerofmass(ta1, isweight=true)
    com2 = centerofmass(ta2, isweight=true)
    com3 = centerofmass(ta3, isweight=true)
    a = zeros(T, nframe)
    for iframe in 1:nframe
        d1 = [com1.xyz[iframe, 1] - com2.xyz[iframe, 1]; com1.xyz[iframe, 2] - com2.xyz[iframe, 2]; com1.xyz[iframe, 3] - com2.xyz[iframe, 3]]
        d2 = [com3.xyz[iframe, 1] - com2.xyz[iframe, 1]; com3.xyz[iframe, 2] - com2.xyz[iframe, 2]; com3.xyz[iframe, 3] - com2.xyz[iframe, 3]]
        a[iframe] = acos(dot(d1, d2)/(norm(d1)*norm(d2)))
    end
    a = (a ./ pi) .* T(180)
end

"""
    compute_dihedral(ta1::TrjArray, ta2::TrjArray, ta3::TrjArray, ta4::TrjArray)

Calculates dihedral angles for the atom coordinates or the COMs of groups from the quadruplet of `ta1`, `ta2`, `ta3` and `ta4`. 

Returns dihedral angles.

# Example
```julia-repl
julia> ta = mdload("ak.pdb")
julia> ta = mdload("ak.dcd", top=ta)
julia> d = compute_dihedral(ta[:, "atomid 1"], ta[:, "atomid 9"], ta[:, "atomid 11", ta[:, "atom 13"]])
```
"""
function compute_dihedral(ta1::TrjArray{T, U}, ta2::TrjArray{T, U}, ta3::TrjArray{T, U}, ta4::TrjArray{T, U})::Vector{T} where {T, U}
    nframe = ta1.nframe
    com1 = centerofmass(ta1, isweight=true)
    com2 = centerofmass(ta2, isweight=true)
    com3 = centerofmass(ta3, isweight=true)
    com4 = centerofmass(ta4, isweight=true)
    a = zeros(T, nframe)
    # Threads.@threads for iframe in 1:nframe
    for iframe in 1:nframe
        d1 = [com1.xyz[iframe, 1] - com2.xyz[iframe, 1]; com1.xyz[iframe, 2] - com2.xyz[iframe, 2]; com1.xyz[iframe, 3] - com2.xyz[iframe, 3]]
        d2 = [com3.xyz[iframe, 1] - com2.xyz[iframe, 1]; com3.xyz[iframe, 2] - com2.xyz[iframe, 2]; com3.xyz[iframe, 3] - com2.xyz[iframe, 3]]
        d3 = [com3.xyz[iframe, 1] - com4.xyz[iframe, 1]; com3.xyz[iframe, 2] - com4.xyz[iframe, 2]; com3.xyz[iframe, 3] - com4.xyz[iframe, 3]]
        m1 = cross(d1, d2)
        m2 = cross(d2, d3)
        a[iframe] = acos(dot(m1, m2)/(norm(m1)*norm(m2)))
        rotdirection = dot(d2,cross(m1,m2))
        if rotdirection < zero(T)
            a[iframe] = -a[iframe]
        end
    end
    a .= (a ./ pi) .* T(180)
end

"""
    compute_dihedral(ta1::TrjArray, ta2::TrjArray, ta3::TrjArray, ta4::TrjArray)

Calculates dihedral angles for the atom coordinates or the COMs of groups from the quadruplet of `ta1`, `ta2`, `ta3` and `ta4`. 

Returns dihedral angles.

# Example
```julia-repl
julia> ta = mdload("ak.pdb")
julia> ta = mdload("ak.dcd", top=ta)
julia> d = compute_dihedral(ta[:, "atomid 1"], ta[:, "atomid 9"], ta[:, "atomid 11", ta[:, "atom 13"]])
```
"""
function compute_dihedral(ta::TrjArray{T, U}, array_index) where {T, U}
    nframe = ta.nframe
    ntorsion = length(array_index)
    natom3 = ta.natom*3
    a = zeros(T, nframe, ntorsion)
    x = view(ta.xyz, :, 1:3:natom3)
    y = view(ta.xyz, :, 2:3:natom3)
    z = view(ta.xyz, :, 3:3:natom3)
    for itorsion = 1:ntorsion
        id = array_index[itorsion]
        for iframe in 1:nframe
            d1 = [x[iframe, id[1]] - x[iframe, id[2]]; y[iframe, id[1]] - y[iframe, id[2]]; z[iframe, id[1]] - z[iframe, id[2]]]
            d2 = [x[iframe, id[3]] - x[iframe, id[2]]; y[iframe, id[3]] - y[iframe, id[2]]; z[iframe, id[3]] - z[iframe, id[2]]]
            d3 = [x[iframe, id[3]] - x[iframe, id[4]]; y[iframe, id[3]] - y[iframe, id[4]]; z[iframe, id[3]] - z[iframe, id[4]]]
            m1 = cross(d1, d2)
            m2 = cross(d2, d3)
            a[iframe, itorsion] = acos(dot(m1, m2)/(norm(m1)*norm(m2)))
            rotdirection = dot(d2,cross(m1,m2))
            if rotdirection < zero(T)
                a[iframe, itorsion] = -a[iframe, itorsion]
            end
        end
    end
    a .= (a ./ pi) .* T(180)
end

function compute_phi(ta::TrjArray{T, U}) where {T, U}
    return 0
end

"""
    compute_qscore(native::TrjArray, ta::TrjArray) -> qscore::Array

Calculates all-atom based Q-score from given heavy atom coordinates. 

Returns Q-scores. 

# Example
```julia-repl
julia> ta = mdload("ak.pdb")
julia> qscore = compute_qscore(ta["not hydrogen"])
```

# References
```
Definition of Q-score from heavy atoms is given in 
R. B. Best, G. Hummer, and W. A. Eaton, PNAS 110, 17874 (2013).
```
"""
function compute_qscore(native::TrjArray{T, U}, ta::TrjArray{T, U})::Vector{T} where {T, U}
    natom = native.natom
    index_all = zeros(U, U(natom*(natom-1)/2), 2)
    count = 1
    for i = 1:natom
        for j = (i+1):natom
            index_all[count, 1] = i
            index_all[count, 2] = j
            count += 1
        end
    end
    id = abs.(native.resid[index_all[:, 1]] .- native.resid[index_all[:, 2]]) .> 3
    index_all = index_all[id, :]
    d = compute_distance(native[1, :], index_all)
    id = (d .< 4.5)[:]
    index_contact = index_all[id, :]
    d_contact = d[:, id]

    d = compute_distance(ta, index_contact)
    BETA = 5.0
    LAMBDA = 1.8
    # qscore = mean(1.0./(1.0 + exp(BETA * bsxfun(@minus, dist, LAMBDA*dist0))), 2);
    qscore = mean(1.0 ./ (1.0 .+ exp.(BETA .* (d .- (LAMBDA .* d_contact)))), dims=2)[:]
    return qscore
end

"""
    compute_drms(native1::TrjArray, native2::TrjArray, ta1::TrjArray, ta2::TrjArray) -> drms

Calculates distance-based root mean square displacements (DRMS) from the given trajectories. 
First, native contacts between molecule1 and molecule2 are identified from the two native structures `native1` and `native2` where cutoff distance of 6.0 Angstrom is used. 
Then, the DRMS are calculated from the trajectories of molecule1 `ta1` and molecule2 `ta2`. 

Returns distance-based root mean square displacements. 

# Example
```julia-repl
julia> native1 = mdload("protein.pdb")["not hydrogen"]
julia> native2 = mdload("ligand.pdb")["not hydrogen"]
julia> ta1 = mdload("protein.dcd")["not hydrogen"]
julia> ta2 = mdload("ligand.dcd")["not hydrogen"]
julia> drms = compute_drms(native1, native2, ta1, ta2)
```

# References
```
Definition of DRMS is given in 
J. Domański, G. Hedger, R. B. Best, P. J. Stansfeld, and M. S. P. Sansom, 
Convergence and Sampling in Determining Free Energy Landscapes for Membrane Protein Association, 
J. Phys. Chem. B 121, 3364 (2017).
```
"""
function compute_drms(native1::TrjArray{T, U}, native2::TrjArray{T, U}, ta1::TrjArray{T, U}, ta2::TrjArray{T, U})::Vector{T} where {T, U}
    natom1 = native1.natom
    natom2 = native2.natom
    index_all = zeros(U, natom1*natom2, 2)
    count = 1
    for i = 1:natom1
        for j = 1:natom2
            index_all[count, 1] = i
            index_all[count, 2] = j
            count += 1
        end
    end
    d = compute_distance(native1[1, :], native2[1, :], index_all)
    id = (1.0 .< d)[:] .& (d .< 6.0)[:]
    index_contact = index_all[id, :]
    d_contact = d[:, id]

    d = compute_distance(ta1, ta2, index_contact)
    ncontact = T(size(index_contact, 1))
    drms = sum((d_contact .- d).^2, dims=2)[:]
    drms .= drms ./ ncontact
    return drms
end

function wrap_x!(x, boxsize)
    x .= x .- floor.(x./boxsize).*boxsize
end

function wrap_dx!(dx, boxsize)
    dx .= dx .- round.(dx./boxsize).*boxsize
end

minimum_image(n::Int, N::Int) = n - N * sign(n) * floor(typeof(n), (Float64(abs(n)) + 0.5*Float64(N)) / Float64(N))

function pairwise_distance(x::AbstractVector, y::AbstractVector, z::AbstractVector, index::AbstractVector, rcut2::Real)
    n = length(index)
    npair = Int(n*(n-1)/2)
    pair = zeros(eltype(index), npair, 2)
    dist = zeros(eltype(x), npair)
    count = 0
    @inbounds for i = 1:n
        ii = index[i]
        for j = (i + 1):n
            jj = index[j]
            r2 = (x[ii] - x[jj])^2 + (y[ii] - y[jj])^2 + (z[ii] - z[jj])^2
            count += 1
            pair[count, 1] = ii
            pair[count, 2] = jj
            dist[count] = r2
        end
    end
    logical_index = dist .< rcut2
    pair[logical_index, :], sqrt.(dist[logical_index])
end

function pairwise_distance(x::AbstractVector, y::AbstractVector, z::AbstractVector, index::AbstractVector, rcut2::Real, boxsize::AbstractVector)
    n = length(index)
    npair = Int(n*(n-1)/2)
    pair = zeros(eltype(index), npair, 2)
    dist = zeros(eltype(x), npair)
    count = 0
    for i = 1:n
        ii = index[i]
        for j = (i + 1):n
            jj = index[j]
            dx = x[ii] - x[jj]
            dy = y[ii] - y[jj]
            dz = z[ii] - z[jj]
            dx = dx - round(dx/boxsize[1])*boxsize[1]
            dy = dy - round(dy/boxsize[2])*boxsize[2]
            dz = dz - round(dz/boxsize[3])*boxsize[3]
            r2 = dx^2 + dy^2 + dz^2
            count += 1
            pair[count, 1] = ii
            pair[count, 2] = jj
            dist[count] = r2
        end
    end
    logical_index = dist .< rcut2
    pair[logical_index, :], sqrt.(dist[logical_index])
end

function pairwise_distance(x::AbstractVector, y::AbstractVector, z::AbstractVector, index1::AbstractVector, index2::AbstractVector, rcut2::Real)
    n1 = length(index1)
    n2 = length(index2)
    npair = n1*n2
    pair = zeros(eltype(index1), npair, 2)
    dist = zeros(eltype(x), npair)
    count = 0
    @inbounds for i = 1:n1
        ii = index1[i]
        for j = 1:n2
            jj = index2[j]
            r2 = (x[ii] - x[jj])^2 + (y[ii] - y[jj])^2 + (z[ii] - z[jj])^2
            count += 1
            pair[count, 1] = ii
            pair[count, 2] = jj
            dist[count] = r2
        end
    end
    logical_index = dist .< rcut2
    pair[logical_index, :], sqrt.(dist[logical_index])
end

function pairwise_distance(x::AbstractVector, y::AbstractVector, z::AbstractVector, index1::AbstractVector, index2::AbstractVector, rcut2::Real, boxsize::AbstractVector)
    n1 = length(index1)
    n2 = length(index2)
    npair = n1*n2
    pair = zeros(eltype(index1), npair, 2)
    dist = zeros(eltype(x), npair)
    count = 0
    @inbounds for i = 1:n1
        ii = index1[i]
        for j = 1:n2
            jj = index2[j]
            dx = x[ii] - x[jj]
            dy = y[ii] - y[jj]
            dz = z[ii] - z[jj]
            dx = dx - round(dx/boxsize[1])*boxsize[1]
            dy = dy - round(dy/boxsize[2])*boxsize[2]
            dz = dz - round(dz/boxsize[3])*boxsize[3]
            r2 = dx^2 + dy^2 + dz^2
            count += 1
            pair[count, 1] = ii
            pair[count, 2] = jj
            dist[count] = r2
        end
    end
    logical_index = dist .< rcut2
    pair[logical_index, :], sqrt.(dist[logical_index])
end

"""
    compute_pairlist(ta::TrjArray, rcut; iframe=1::Int) -> F

Make a pairlist for the given strcuture `ta1` by searching pairs within a cutoff distance `rcut`. 
By default, the 1st frame of `ta1` is used. Users can specify the frame by `iframe`. 

Returns a NamedTuple object `F` which contains the pair lists in `F.pair`, 
and the distances of corresponding pairs in `F.dist`. 

# Example
```julia-repl
julia> ta = mdload("ak.pdb")
julia> F = compute_pairlist(ta, 8.0)
```

# References
```
The algorithm of this function is based on 
T.N. Heinz, and P.H. Hünenberger, 
A fast pairlist-construction algorithm for molecular simulations under periodic boundary conditions. 
J Comput Chem 25, 1474–1486 (2004). 
```
"""
function compute_pairlist(ta::TrjArray{T, U}, rcut::T; iframe=1::Int) where {T, U}
    # TODO: PBC support
    natom = ta.natom
    rcut2 = rcut^2
    is_pbc = !isempty(ta.boxsize)
    x = ta.xyz[iframe, 1:3:end]
    y = ta.xyz[iframe, 2:3:end]
    z = ta.xyz[iframe, 3:3:end]
    boxsize = zeros(T, 3)

    if is_pbc
        boxsize .= ta.boxsize[iframe, :]
        wrap_x!(x, boxsize[1])
        wrap_x!(y, boxsize[2])
        wrap_x!(z, boxsize[3])
    else
        x .= x .- minimum(x)
        y .= y .- minimum(y)
        z .= z .- minimum(z)
        boxsize[1] = maximum(x) + one(T)
        boxsize[2] = maximum(y) + one(T)
        boxsize[3] = maximum(z) + one(T)
    end

    if any(boxsize .< (2.0*rcut))
        println("working with brute-force search...")
        return compute_pairlist_bruteforce(ta, rcut, iframe=iframe)
    end

    rcell = map(x -> max(x, rcut), T.([8.0, 8.0, 8.0]))
    ncell = floor.(U, boxsize./T.(rcell))
    rcell .= boxsize./T.(ncell)

    nx = floor.(U, x./rcell[1])
    ny = floor.(U, y./rcell[2])
    nz = floor.(U, z./rcell[3])

    m = (ncell[1]*ncell[2]).*nz .+ ncell[1].*ny .+ nx .+ one(U)
    M = prod(ncell)
    
    # calculate cell mask
    dm = collect(U, 1:(M-1))
    dmx = mod.(dm, ncell[1])
    dmy = floor.(U, mod.(dm, ncell[1]*ncell[2]) ./ ncell[1])
    dmz = floor.(U, dm ./ (ncell[1]*ncell[2]))
    dnx = abs.(minimum_image.(dmx, ncell[1]))
    dny = similar(dmy)
    logical_index = iszero.(dmx) .| (isequal.(dmy, ncell[2] .- one(U)) .& isequal.(dmz, ncell[3] .- one(U)))
    dny[logical_index] .= abs.(minimum_image.(dmy[logical_index], ncell[2]))
    dny[.~logical_index] .= min.(abs.(minimum_image.(dmy[.~logical_index], ncell[2])), abs.(minimum_image.(dmy[.~logical_index] .+ 1, ncell[2])))
    dnz = similar(dmz)
    logical_index = isequal.(dmz, ncell[3] .- one(U)) .| (iszero.(dmx) .& iszero.(dmy))
    dnz[logical_index] .= abs.(minimum_image.(dmz[logical_index], ncell[3]))
    dnz[.~logical_index] .= min.(abs.(minimum_image.(dmz[.~logical_index], ncell[3])), abs.(minimum_image.(dmz[.~logical_index] .+ 1, ncell[3])))
    mask = T.((max.(dnx, one(U)) .- one(U)).^2) .* rcell[1].^2 .+ 
           T.((max.(dny, one(U)) .- one(U)).^2) .* rcell[2].^2 .+ 
           T.((max.(dnz, one(U)) .- one(U)).^2) .* rcell[3].^2 .< rcut2
    mask_index = findall(mask)

    # calculate cell mask pointer
    K = round(U, (natom/M) * 20) ######### TODO, 20 is really OK?
    cell_pointer = zeros(U, M)
    cell_atom = zeros(U, K, M)
    for iatom = 1:natom
        cell_pointer[m[iatom]] += 1
        cell_atom[cell_pointer[m[iatom]], m[iatom]] = iatom
    end
    if any(cell_pointer .> K)
        error("# of atoms in a cell greater than expected")
    end
    
    # calculate distances of atoms in masked cells
    pair = Matrix{U}(undef, natom*1000, 2)
    dist = Vector{T}(undef, natom*1000)
    count = 0
    @inbounds for m1 = 1:M
        if iszero(cell_pointer[m1])
            continue
        end
        m1_index = cell_atom[1:cell_pointer[m1], m1]
        if isempty(m1_index)
            continue
        elseif length(m1_index) > 1
            if is_pbc
                pair_local, dist_local = pairwise_distance(x, y, z, m1_index, rcut2, boxsize)
            else
                pair_local, dist_local = pairwise_distance(x, y, z, m1_index, rcut2)
            end
            num = length(dist_local)
            if (num > 0)
                pair[(count+1):(count+num), :] .= pair_local
                dist[(count+1):(count+num)] .= dist_local
                count += num
            end
        end

        for m2 = (m1 .+ mask_index)
            if m2 > M
                break
            end
            if iszero(cell_pointer[m2])
                continue
            end
            m2_index = cell_atom[1:cell_pointer[m2], m2]
            if is_pbc
                pair_local, dist_local = pairwise_distance(x, y, z, m1_index, m2_index, rcut2, boxsize)
            else
                pair_local, dist_local = pairwise_distance(x, y, z, m1_index, m2_index, rcut2)
            end
            num = length(dist_local)
            if (num > 0)
                pair[(count+1):(count+num), :] .= pair_local
                dist[(count+1):(count+num)] .= dist_local
                count += num    
            end
        end
    end
    return (pair=pair[1:count, :], dist=dist[1:count])
end

"""
    compute_pairlist_bruteforce(ta::TrjArray, rcut; iframe=1::Int) -> F

Bruteforce search version of `compute_pairlist` function. Mainly used for debugging. 

Returns a NamedTuple object `F` which contains the pair lists in `F.pair`, 
and the distances of corresponding pairs in `F.dist`. 

# Example
```julia-repl
julia> ta = mdload("ak.pdb")
julia> F = compute_pairlist(ta, 8.0)
```
"""
function compute_pairlist_bruteforce(ta::TrjArray{T, U}, rcut::T; iframe=1::Int) where {T, U}
    rcut2 = rcut^2
    natom3 = ta.natom*3
    x = view(ta.xyz, iframe, 1:3:natom3)
    y = view(ta.xyz, iframe, 2:3:natom3)
    z = view(ta.xyz, iframe, 3:3:natom3)
    index = collect(U, 1:ta.natom)
    is_pbc = !isempty(ta.boxsize)
    boxsize = zeros(T, 3)
    if is_pbc
        boxsize .= ta.boxsize[iframe, :]
        pair, dist = pairwise_distance(x, y, z, index, rcut2, boxsize)
    else
        pair, dist = pairwise_distance(x, y, z, index, rcut2)
    end
    return (pair=pair, dist=dist)
end

