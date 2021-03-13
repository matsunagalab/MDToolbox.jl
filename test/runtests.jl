using MDToolbox
using Test

@testset "trjarray.jl" begin
    @testset "null object" begin
        @test TrjArray{Float64, Int64}().natom == 0
        @test TrjArray{Float64, Int64}().nframe == 0
    end

    @testset "constructor" begin
        nframe = 5
        natom = 5
        ta = TrjArray{Float64, Int64}(xyz=rand(nframe, natom*3))

        @test ta.natom == natom
        @test ta.nframe == nframe
    end

    @testset "slicing" begin
        nframe = 5
        natom = 5
        ta = TrjArray{Float64, Int64}(xyz=rand(nframe, natom*3))

        @test ta[1, :].xyz ≈ ta.xyz[1:1, :]
        @test ta[1:3, :].xyz ≈ ta.xyz[1:3, :]
        @test ta[[1, 3, 5], :].xyz ≈ ta.xyz[[1, 3, 5], :]
        @test ta[[true, true, true, false, false], :].xyz ≈ ta.xyz[[true, true, true, false, false], :]

        @test ta[:, 1].xyz ≈ ta.xyz[:, MDToolbox.to3(1:1)]
        @test ta[:, 1:3].xyz ≈ ta.xyz[:, MDToolbox.to3(1:3)]
        @test ta[:, [1, 3, 5]].xyz ≈ ta.xyz[:, MDToolbox.to3([1, 3, 5])]
        @test ta[:, [true, true, true, false, false]].xyz ≈ ta.xyz[:, MDToolbox.to3([true, true, true, false, false])]
    end

    @testset "slicing with boxsize" begin
        nframe = 5
        natom = 5
        ta = TrjArray{Float64, Int64}(xyz=rand(nframe, natom*3), boxsize=rand(nframe, 3))

        @test ta[1, :].boxsize ≈ ta.boxsize[1:1, :]
        @test ta[1:3, :].boxsize ≈ ta.boxsize[1:3, :]
        @test ta[[1, 3, 5], :].boxsize ≈ ta.boxsize[[1, 3, 5], :]
        @test ta[[true, true, true, false, false], :].boxsize ≈ ta.boxsize[[true, true, true, false, false], :]
    end

    @testset "atom selection" begin
        nframe = 5
        natom = 5
        ta = TrjArray{Float64, Int64}(xyz=rand(nframe, natom*3), atomid=1:natom, chainname=fill("A", natom))

        @test ta["atomid 1:3"].xyz ≈ ta.xyz[:, MDToolbox.to3(1:3)]
        @test ta["chainname A"].xyz ≈ ta.xyz[:, :]
        @test ta["atomid 3:5 and chainname A"].xyz ≈ ta.xyz[:, MDToolbox.to3(3:5)]
        @test ta["atomid 3:5 or chainname A"].xyz ≈ ta.xyz[:, :]
        @test ta["(atomid 3:5 and chainname A) and atomid 3"].xyz ≈ ta.xyz[:, MDToolbox.to3(3)]
        @test ta["(atomid 3:5 and chainname A) and (atomid 3 and chainname A)"].xyz ≈ ta.xyz[:, MDToolbox.to3(3)]
        @test ta["not atomid 1:3"].xyz ≈ ta.xyz[:, MDToolbox.to3(4:5)]
        @test ta["not (atomid 1:2 or atomid 4:5)"].xyz ≈ ta.xyz[:, MDToolbox.to3(3)]
    end

    @testset "merging" begin
        nframe = 5
        natom = 5
        topology = TrjArray{Float64, Int64}(atomid=1:natom)
        ta = TrjArray{Float64, Int64}(xyz=rand(nframe, natom*3))

        ta_new = [topology; ta]
        @test ta_new["atomid 1:3"].xyz ≈ ta.xyz[:, MDToolbox.to3(1:3)]
    end
end

@testset "structure.jl" begin
    @testset "centerofmass" begin
        c = [-1.0 0.0 1.0;
             -1.0 0.0 1.0]
        xyz = zeros(Float64, 2, 9)
        xyz[:, 1:3:end] .= c
        xyz[:, 2:3:end] .= c
        xyz[:, 3:3:end] .= c
        ta = TrjArray{Float64, Int64}(xyz=xyz)
        com = centerofmass(ta)
        @test com.xyz ≈ [0.0 0.0 0.0; 0.0 0.0 0.0]

        c = [-2.0 0.0 1.0;
             -2.0 0.0 1.0]
        xyz = zeros(Float64, 2, 9)
        ta = TrjArray{Float64, Int64}(xyz=xyz, mass=[0.5, 1.0, 1.0])
        com = centerofmass(ta)
        @test com.xyz ≈ [0.0 0.0 0.0; 0.0 0.0 0.0]

        c = rand(Float64)
        c = [c c c]
        ta = TrjArray{Float64, Int64}(xyz=c)
        com = centerofmass(ta)
        @test com.xyz ≈ c
    end

    @testset "decenter" begin
        c = [1.0 1.0 1.0;
             2.0 2.0 2.0]
        xyz = zeros(Float64, 2, 9)
        xyz[:, 1:3:end] .= c
        xyz[:, 2:3:end] .= c
        xyz[:, 3:3:end] .= c
        ta = TrjArray{Float64, Int64}(xyz=xyz)
        ta2, com = decenter(ta)
        @test ta2.xyz ≈ zeros(Float64, 2, 9)

        c = rand(Float64)
        c = [c c c]
        ta = TrjArray{Float64, Int64}(xyz=c)
        ta2, com = decenter(ta)
        @test ta2.xyz ≈ zeros(Float64, 1, 3)
    end

    @testset "superimpose and getrmsd" begin
        A = [-2.803  -15.373   24.556;
             0.893  -16.062   25.147;
             1.368  -12.371   25.885;
             -1.651  -12.153   28.177;
             -0.440  -15.218   30.068;
             2.551  -13.273   31.372;
             0.105  -11.330   33.567]
        ref = TrjArray{Float64, Int64}(xyz=A'[:]', mass=collect(1:7))
        A = [-14.739  -18.673   15.040;
             -12.473  -15.810   16.074;
             -14.802  -13.307   14.408;
             -17.782  -14.852   16.171;
             -16.124  -14.617   19.584;
             -15.029  -11.037   18.902;
             -18.577  -10.001   17.996]
        ta = TrjArray{Float64, Int64}(xyz=A'[:]', mass=collect(1:7))
        # r, ta_fit = superimpose(ref, ta, isweight=true)
        ta_fit = superimpose(ref, ta, isweight=true)
        #@test r[1] ≈ 0.7450161471
        r2 = compute_rmsd(ref, ta_fit)
        #@test r2[1] ≈ r[1]
        @test r2[1] ≈ 0.7450161471
    end

    @testset "getdistance" begin
        x = [1.0 0.0]
        y = [0.0 0.0]
        z = [0.0 0.0]
        xyz = rand(Float64, 1, 6)
        xyz[1:1, 1:3:end] .= x
        xyz[1:1, 2:3:end] .= y
        xyz[1:1, 3:3:end] .= z
        ta = TrjArray{Float64, Int64}(xyz=xyz)
        d = compute_distance(ta[:, 1], ta[:, 2])
        @test d ≈ [1.0]
    end

    @testset "getangle" begin
        x = [1.0 0.0 0.0]
        y = [0.0 0.0 1.0]
        z = [0.0 0.0 0.0]
        xyz = rand(Float64, 1, 9)
        xyz[1:1, 1:3:end] .= x
        xyz[1:1, 2:3:end] .= y
        xyz[1:1, 3:3:end] .= z
        ta = TrjArray{Float64, Int64}(xyz=xyz)
        a = compute_angle(ta[:, 1], ta[:, 2], ta[:, 3])
        @test a ≈ [90.0]
    end

    @testset "getdihedral" begin
        x = [-1.0 -1.0 1.0 1.0]
        y = [-1.0  0.0 0.0 1.0]
        z = [ 0.0  0.0 0.0 0.0]
        xyz = rand(Float64, 1, 12)
        xyz[1:1, 1:3:end] .= x
        xyz[1:1, 2:3:end] .= y
        xyz[1:1, 3:3:end] .= z
        ta = TrjArray{Float64, Int64}(xyz=xyz)
        a = compute_dihedral(ta[:, 1], ta[:, 2], ta[:, 3], ta[:, 4])
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
