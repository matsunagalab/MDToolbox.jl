using MDToolbox
using Test

@testset "TrjArray operations" begin
    @testset "null object" begin
        @test TrjArray().natom == 0
        @test TrjArray().nframe == 0
    end

    @testset "slicing" begin
        natom = 5
        nframe = 5
        ta = TrjArray(x=rand(nframe, natom), y=rand(nframe, natom), z=rand(nframe, natom))

        @test ta.natom == natom
        @test ta.nframe == nframe

        @test ta[1, :].x ≈ ta.x[1:1, :]
        @test ta[1:3, :].x ≈ ta.x[1:3, :]
        @test ta[[1, 3, 5], :].x ≈ ta.x[[1, 3, 5], :]
        @test ta[[true, true, true, false, false], :].x ≈ ta.x[[true, true, true, false, false], :]

        @test ta[:, 1].x ≈ ta.x[:, 1:1]
        @test ta[:, 1:3].x ≈ ta.x[:, 1:3]
        @test ta[:, [1, 3, 5]].x ≈ ta.x[:, [1, 3, 5]]
        @test ta[:, [true, true, true, false, false]].x ≈ ta.x[:, [true, true, true, false, false]]
    end

    @testset "atom selection" begin
        natom = 5
        nframe = 5
        ta = TrjArray(x=rand(nframe, natom), y=rand(nframe, natom), z=rand(nframe, natom), atomid=1:natom, chainname=fill("A", natom))
        @test ta["atomid 1:3"].x ≈ ta.x[:, 1:3]
        @test ta["chainname A"].x ≈ ta.x[:, :]
        @test ta["atomid 3:5 and chainname A"].x ≈ ta.x[:, 3:5]
        @test ta["atomid 3:5 or chainname A"].x ≈ ta.x[:, :]
        @test ta["(atomid 3:5 and chainname A) and atomid 3"].x ≈ ta.x[:, 3]
        @test ta["(atomid 3:5 and chainname A) and (atomid 3 and chainname A)"].x ≈ ta.x[:, 3]
    end

    @testset "merging" begin
        natom = 5
        nframe = 5
        topology = TrjArray(atomid=1:natom)
        ta = TrjArray(x=rand(nframe, natom), y=rand(nframe, natom), z=rand(nframe, natom))

        ta_new = [topology; ta]
        ta_new["atomid 1:3"].x ≈ ta.x[:, 1:3]
    end
end

