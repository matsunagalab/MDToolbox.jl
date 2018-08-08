function InnerProduct!(A::Array{Float64, 1}, coords1::Array{Float64, 1}, double **coords2, const size_t len, const double *weight)
    double          x1, x2, y1, y2, z1, z2;
    size_t          i;
    const double   *fx1 = coords1[0], *fy1 = coords1[1], *fz1 = coords1[2];
    const double   *fx2 = coords2[0], *fy2 = coords2[1], *fz2 = coords2[2];
    double          G1 = 0.0, G2 = 0.0;

    A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = 0.0;

    if (weight != NULL)
    {
        for (i = 0; i < len; ++i)
        {
            x1 = weight[i] * fx1[i];
            y1 = weight[i] * fy1[i];
            z1 = weight[i] * fz1[i];

            G1 += x1 * fx1[i] + y1 * fy1[i] + z1 * fz1[i];

            x2 = fx2[i];
            y2 = fy2[i];
            z2 = fz2[i];

            G2 += weight[i] * (x2 * x2 + y2 * y2 + z2 * z2);

            A[0] +=  (x1 * x2);
            A[1] +=  (x1 * y2);
            A[2] +=  (x1 * z2);

            A[3] +=  (y1 * x2);
            A[4] +=  (y1 * y2);
            A[5] +=  (y1 * z2);

            A[6] +=  (z1 * x2);
            A[7] +=  (z1 * y2);
            A[8] +=  (z1 * z2);   
        }
    }
    else
    {
        for (i = 0; i < len; ++i)
        {
            x1 = fx1[i];
            y1 = fy1[i];
            z1 = fz1[i];

            G1 += x1 * x1 + y1 * y1 + z1 * z1;

            x2 = fx2[i];
            y2 = fy2[i];
            z2 = fz2[i];

            G2 += (x2 * x2 + y2 * y2 + z2 * z2);

            A[0] +=  (x1 * x2);
            A[1] +=  (x1 * y2);
            A[2] +=  (x1 * z2);

            A[3] +=  (y1 * x2);
            A[4] +=  (y1 * y2);
            A[5] +=  (y1 * z2);

            A[6] +=  (z1 * x2);
            A[7] +=  (z1 * y2);
            A[8] +=  (z1 * z2);  
        }
    }

    return (G1 + G2) * 0.5
end


function superimpose()

    ## initialization
    @show "called"
    nframe = size(trj, 1)
    npair = size(pair, 1)

    ## calculation
    bond = zeros(Float64, (nframe, npair))
    for ipair = 1:npair
        index1 = to3(pair[ipair, 1])
        index2 = to3(pair[ipair, 2])
        for iframe = 1:nframe
            d = sum((trj[iframe, index1] - trj[iframe, index2]).^2)
            d = sqrt(d)
            bond[iframe, ipair] = d
        end
    end

    bond
end
