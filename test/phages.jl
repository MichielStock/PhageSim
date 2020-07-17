
@testset "Phages" begin

    phage = Phage(1, (3, 4), 2)

    @testset "AbstractPhages" begin
        @test phage isa AbstractAgent
        @test phage isa AbstractPhage
        @test phage isa Phage

        @test phages(phage)
        @test species(phage) == 2

    end

    @testset "PhageRules" begin
        phagerules = PhageRules(0.3)
        @test phagerules.pmove ≈ 1.0

        @test_throws AssertionError PhageRules(0.3, 1.2)
        @test_throws AssertionError PhageRules(-0.3, 0.1)
    end

    @testset "InteractionRules" begin
        bact = Bacterium(1, (3, 4), 1)
        phage = Phage(1, (3, 4), 2)

        @testset "homogeneous" begin
            interactrules = InteractionRules(1.0, 2.1, false, 0.01, 2)

            @test interactrules isa InteractionRules
            @test !lysogenic(phage, bact, interactrules)
            @test infects(phage, bact, interactrules)

        end

        @testset "heterogeneous" begin
            interactrules = InteractionRules([true false; true true], 1.1, [false, true], 0.01, 2)

            @test interactrules isa InteractionRules
            @test lysogenic(phage, bact, interactrules)
            @test !infects(phage, bact, interactrules)

        end
    end

end
