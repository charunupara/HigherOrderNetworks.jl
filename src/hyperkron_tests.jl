include("hyperkron.jl")
using Test

@testset "hyperkron_tests" begin
    @testset "num_multiset_permutations" begin
        @test num_multiset_permutations(Integer[]) == 1
        @test num_multiset_permutations([0]) == 1
        @test num_multiset_permutations([0,0,0]) == 1
        @test num_multiset_permutations([0,0,1]) == 3
        @test num_multiset_permutations([1,2,3,3]) == 12
        @test num_multiset_permutations([1,2,3,4]) == 24
        @test num_multiset_permutations([1,1,5,6,7,7,8]) == 1260
    end

    @testset "_next_distinct_character" begin
        @test _next_distinct_character(Integer[], 0, 0) == 1
        @test _next_distinct_character([0], 1, 1) == 2
        @test _next_distinct_character([0], 4, 1) == 2
        @test _next_distinct_character([0,0,0], 1, 3) == 4
        @test _next_distinct_character([0,0,0], 6, 3) == 4
        @test _next_distinct_character([0,0,1], 1, 3) == 3
        @test _next_distinct_character([0,0,1], 2, 3) == 3
        @test _next_distinct_character([0,0,1], 3, 3) == 4
        @test _next_distinct_character([1,1,1,3,5], 2, 5) == 4
        @test _next_distinct_character([4,5,7,7,7], 3, 5) == 6
        @test _next_distinct_character([4,5,7,7,7,8], 1, 6) == 2
    end

    @testset "unrank5!" begin
        @test unrank5!(Integer[], 1, 1, 0) == Integer[]
        @test unrank5!([0], 1, 1, 1) == [0]
        @test unrank5!([0], 2, 1, 1) == [0]
        @test unrank5!([0], 0, 1, 1) == [0]
        @test unrank5!([0,0,1], 1, 3, 3) == [0,0,1]
        @test unrank5!([0,0,1], 2, 3, 3) == [0,1,0]
        @test unrank5!([0,1,2,2], 7, 12, 4) == [2,0,1,2]
    end

    @testset "_morton_decode3" begin
        #=@test _morton_decode3(Integer[], 2, 8) == (1,1,1)
        @test _morton_decode3([0], 2, (2,2,2)) == (1,1,1)
        @test _morton_decode3([1], 2, (2,2,2)) == (1,1,2)
        @test _morton_decode3([7], 2, (2,2,2)) == (2,2,2)
        @test _morton_decode3([0,0], 2, (4,4,4)) == (1,1,1)
        @test _morton_decode3([3,4], 2, (4,4,4)) == (4,1,2)=#
        @test _morton_decode3([7,1,2], 2, (8,8,8)) == (4,5,7)
    end

    @testset "hyperkron_graph" begin
        G, _ = hyperkron_graph(kron_params(0.4, 0.8, 0.2, 0.5), 3)
        @test all([G.colptr[i] != G.rowval[i] for i = 1:size(G,1)]) # no self-loops
    end
end
