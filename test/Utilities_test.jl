using Test
include("Utilities.jl")

@testset "Utilities" begin
    G1 = MatrixNetwork(load_matrix_network("cores_example")) # Test other small networks
    G2 = make_undirected(MatrixNetwork([(1,2), (2,4), (2,5), (3,4), (3,5)], 5))
    @testset "Degree Distribution" begin
        for i = 2:10 # Testing common subgraphs: path, cycle, clique
            @test deg_distr(make_undirected(path(i))) == [1; [2 for j = 1:i-1]; 1]
            @test deg_distr(make_undirected(cycle(i+1))) == [2 for j = 1:i+1]
            @test deg_distr(clique(i)) == [i-1 for j = 1:i]
        end

        @test deg_distr(G1) == [1,4,2,6,5,3,3,4,3,5,4,1,1,3,3,0,2,2,2,3,1]
        @test deg_distr(G2) == [1,3,2,2,2]
    end

    @testset "Nodes by Degree" begin # Getting sorted list of nodes by their degree
        @test nodes_by_deg(G1) == [16, 1, 12, 13, 21, 3, 17, 18, 19, 6, 7, 9, 14, 15, 20, 2, 8, 11, 5, 10, 4]
        @test nodes_by_deg(G1; descending=true) == [4, 5, 10, 2, 8, 11, 6, 7, 9, 14, 15, 20, 3, 17, 18, 19, 1, 12, 13, 21, 16]

        @test nodes_by_deg(G2) == [1, 3, 4, 5, 2]
        @test nodes_by_deg(G2; descending=true) == [2, 3, 4, 5, 1]
    end

    @testset "Outweight and Inweight Distributions" begin # Weight coming out of and going in to each node
        @test
    end
end
