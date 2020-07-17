@testset "Bacteria" begin
    bact = Bacterium(30, (2, 3), 3, 10.0, nothing)

    @testset "AbstractBacteria" begin
        @test bact isa AbstractAgent
        @test bact isa AbstractBacterium

        @test bacteria(bact)
        @test species(bact) == 3
        @test energy(bact) ≈ 10
        energy!(bact, 9.0)
        @test energy(bact) ≈ 9.0
        @test !haslatent(bact)
        prophage!(bact, 4)
        @test haslatent(bact)
        @test prophage(bact) == 4

        energyupdate(bact, ConstantEnergyUpdate(1.0, 15.0))
        @test energy(bact) ≈ 10

        energyupdate(bact, RandomEnergyUpdate(50.0, 15.0, 0.1))
        @test energy(bact) ≈ 15
    end

    @testset "BacteriaRules" begin

        bactrule = BacteriaRules(0.1, [3, 4, 5], ConstantEnergyUpdate(1.0, 15.0))

        @test pmove(bact, bactrule) ≈ 0.1
        @test Ediv(bact, bactrule) ≈ 5

    end
end
