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
    @test 1 == 1
end

