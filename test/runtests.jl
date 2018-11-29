using MDToolbox
using Test

@testset "trjarray.jl" begin
    @testset "null object" begin
        @test TrjArray().natom == 0
        @test TrjArray().nframe == 0
    end

    @testset "constructor" begin
        nframe = 5
        natom = 5
        ta = TrjArray(x=rand(nframe, natom), y=rand(nframe, natom), z=rand(nframe, natom))

        @test ta.natom == natom
        @test ta.nframe == nframe
    end

    @testset "slicing" begin
        nframe = 5
        natom = 5
        ta = TrjArray(x=rand(nframe, natom), y=rand(nframe, natom), z=rand(nframe, natom))

        @test ta[1, :].x ≈ ta.x[1:1, :]
        @test ta[1:3, :].x ≈ ta.x[1:3, :]
        @test ta[[1, 3, 5], :].x ≈ ta.x[[1, 3, 5], :]
        @test ta[[true, true, true, false, false], :].x ≈ ta.x[[true, true, true, false, false], :]

        @test ta[:, 1].x ≈ ta.x[:, 1:1]
        @test ta[:, 1:3].x ≈ ta.x[:, 1:3]
        @test ta[:, [1, 3, 5]].x ≈ ta.x[:, [1, 3, 5]]
        @test ta[:, [true, true, true, false, false]].x ≈ ta.x[:, [true, true, true, false, false]]
    end

    @testset "slicing with boxsize" begin
        nframe = 5
        natom = 5
        ta = TrjArray(x=rand(nframe, natom), y=rand(nframe, natom), z=rand(nframe, natom), boxsize=rand(nframe, 3))

        @test ta[1, :].boxsize ≈ ta.boxsize[1:1, :]
        @test ta[1:3, :].boxsize ≈ ta.boxsize[1:3, :]
        @test ta[[1, 3, 5], :].boxsize ≈ ta.boxsize[[1, 3, 5], :]
        @test ta[[true, true, true, false, false], :].boxsize ≈ ta.boxsize[[true, true, true, false, false], :]
    end

    @testset "atom selection" begin
        nframe = 5
        natom = 5
        ta = TrjArray(x=rand(nframe, natom), y=rand(nframe, natom), z=rand(nframe, natom), atomid=1:natom, chainname=fill("A", natom))

        @test ta["atomid 1:3"].x ≈ ta.x[:, 1:3]
        @test ta["chainname A"].x ≈ ta.x[:, :]
        @test ta["atomid 3:5 and chainname A"].x ≈ ta.x[:, 3:5]
        @test ta["atomid 3:5 or chainname A"].x ≈ ta.x[:, :]
        @test ta["(atomid 3:5 and chainname A) and atomid 3"].x ≈ ta.x[:, 3]
        @test ta["(atomid 3:5 and chainname A) and (atomid 3 and chainname A)"].x ≈ ta.x[:, 3]
        @test ta["not atomid 1:3"].x ≈ ta.x[:, 4:5]
        @test ta["not (atomid 1:2 or atomid 4:5)"].x ≈ ta.x[:, 3]
    end

    @testset "merging" begin
        nframe = 5
        natom = 5
        topology = TrjArray(atomid=1:natom)
        ta = TrjArray(x=rand(nframe, natom), y=rand(nframe, natom), z=rand(nframe, natom))

        ta_new = [topology; ta]
        ta_new["atomid 1:3"].x ≈ ta.x[:, 1:3]
    end
end

@testset "structure.jl" begin
    @testset "centerofmass" begin
        c = [-1.0 0.0 1.0;
             -1.0 0.0 1.0]
        ta = TrjArray(x=c, y=c, z=c)
        com = centerofmass(ta)
        @test com.x ≈ com.y ≈ com.z ≈ [0.0, 0.0]

        c = [-2.0 0.0 1.0;
             -2.0 0.0 1.0]
        ta = TrjArray(x=c, y=c, z=c, mass=[0.5, 1.0, 1.0])
        com = centerofmass(ta)
        @test com.x ≈ com.y ≈ com.z ≈ [0.0, 0.0]

        c = [rand(Float64)]
        c = reshape(c, 1, 1)
        ta = TrjArray(x=c, y=c, z=c)
        com = centerofmass(ta)
        @test com.x ≈ com.y ≈ com.z ≈ c
    end

    @testset "decenter" begin
        c = [1.0 1.0 1.0;
             2.0 2.0 2.0]
        ta = TrjArray(x=c, y=c, z=c)
        ta2, com = decenter(ta)
        @test ta2.x ≈ ta2.y ≈ ta2.z ≈ [0.0 0.0 0.0; 0.0 0.0 0.0]

        c = [rand(Float64)]
        c = reshape(c, 1, 1)
        ta = TrjArray(x=c, y=c, z=c)
        ta2, com = decenter(ta)
        @test ta2.x ≈ ta2.y ≈ ta2.z ≈ [0.0]
    end

    @testset "superimpose and get_rmsd" begin
        A = [-2.803  -15.373   24.556;
             0.893  -16.062   25.147;
             1.368  -12.371   25.885;
             -1.651  -12.153   28.177;
             -0.440  -15.218   30.068;
             2.551  -13.273   31.372;
             0.105  -11.330   33.567]
        x = A[:, 1]
        y = A[:, 2]
        z = A[:, 3]
        ref = TrjArray(x=x', y=y', z=z', mass=1:7)
        A = [-14.739  -18.673   15.040;
             -12.473  -15.810   16.074;
             -14.802  -13.307   14.408;
             -17.782  -14.852   16.171;
             -16.124  -14.617   19.584;
             -15.029  -11.037   18.902;
             -18.577  -10.001   17.996]
        x = A[:, 1]
        y = A[:, 2]
        z = A[:, 3]
        ta = TrjArray(x=x', y=y', z=z', mass=1:7)
        # r, ta_fit = superimpose(ref, ta, isweight=true)
        ta_fit = superimpose(ref, ta, isweight=true)
        #@test r[1] ≈ 0.7450161471
        r2 = get_rmsd(ref, ta_fit)
        #@test r2[1] ≈ r[1]
        @test r2[1] ≈ 0.7450161471
    end

    @testset "get_distance" begin
        x = [1.0 0.0]
        y = [0.0 0.0]
        z = [0.0 0.0]
        ta = TrjArray(x=x, y=y, z=z)
        d = get_distance(ta[:, 1], ta[:, 2])
        @test d ≈ [1.0]
    end

    @testset "get_angle" begin
        x = [1.0 0.0 0.0]
        y = [0.0 0.0 1.0]
        z = [0.0 0.0 0.0]
        ta = TrjArray(x=x, y=y, z=z)
        a = get_angle(ta[:, 1], ta[:, 2], ta[:, 3])
        @test a ≈ [90.0]
    end

    @testset "get_dihedral" begin
        x = [-1.0 -1.0 1.0 1.0]
        y = [-1.0  0.0 0.0 1.0]
        z = [ 0.0  0.0 0.0 0.0]
        ta = TrjArray(x=x, y=y, z=z)
        a = get_dihedral(ta[:, 1], ta[:, 2], ta[:, 3], ta[:, 4])
        @test a ≈ [180.0]
    end
end

@testset "ksdensity.jl" begin
    @testset "ksdensity 1d" begin
        grid_x=collect(-3:0.1:3)
        f, _ = ksdensity(randn(10^6), grid_x=grid_x);
        sigma = 1.0
        f_true = (1.0/(sqrt(2.0*pi)*sigma)).*exp.(-grid_x.*grid_x./(2.0*sigma^2));
        @test maximum(abs.(f - f_true)) < 0.01
    end

    @testset "ksdensity 2d" begin
        grid_x=collect(-2:0.1:2)
        grid_y=collect(-2:0.1:2)
        f, _, _ = ksdensity(randn(10^6), randn(10^6), grid_x=grid_x, grid_y=grid_y);
        sigma = 1.0
        f_true = zeros(Float64, length(grid_x), length(grid_y))
        for i in 1:length(grid_x)
            for j in 1:length(grid_y)
                f_true[i, j] = (1.0/(sqrt(2.0*pi)*sigma))*exp(-grid_x[i]*grid_x[i]/(2.0*sigma^2))*
                            (1.0/(sqrt(2.0*pi)*sigma))*exp(-grid_y[j]*grid_y[j]/(2.0*sigma^2));
            end
        end
        @test maximum(abs.(f - f_true)) < 0.01
    end

    @testset "ksdensity 3d" begin
        grid_x=collect(-1:0.1:1)
        grid_y=collect(-1:0.1:1)
        grid_z=collect(-1:0.1:1)
        f, _, _, _ = ksdensity(randn(10^6), randn(10^6), randn(10^6), grid_x=grid_x, grid_y=grid_y, grid_z=grid_z);
        sigma = 1.0
        f_true = zeros(Float64, length(grid_x), length(grid_y), length(grid_z))
        for i in 1:length(grid_x)
            for j in 1:length(grid_y)
                for k in 1:length(grid_z)
                    f_true[i, j, k] = (1.0/(sqrt(2.0*pi)*sigma))*exp(-grid_x[i]*grid_x[i]/(2.0*sigma^2))*
                            (1.0/(sqrt(2.0*pi)*sigma))*exp(-grid_y[j]*grid_y[j]/(2.0*sigma^2))*
                            (1.0/(sqrt(2.0*pi)*sigma))*exp(-grid_z[j]*grid_z[j]/(2.0*sigma^2));
                end
            end
        end
        @test maximum(abs.(f - f_true)) < 0.05
    end
end
