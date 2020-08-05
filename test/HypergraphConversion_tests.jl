include("hypergraphs_conversions.jl")
using Test

@testset "hyper_to_matrix_tests" begin
    # Trivial graphs
    @test to_adjacency_matrix(Hypergraphs(Array{Int64,1}[], 1, 0))[1,:] == [0]
    @test to_adjacency_matrix(Hypergraphs([[1]], 1, 1))[1,:] == [0]
    @test to_adjacency_matrix(Hypergraphs([[1,1]], 1, 1))[1,:] == [2]

    # Only self-loops
    @test to_adjacency_matrix(Hypergraphs([[1,1,1], [2,2], [3]], 3, 3)) == [6 0 0; 0 2 0; 0 0 0]

    # Multi-edges and self-loops
    @test to_adjacency_matrix(Hypergraphs([[1,4], [2,3], [1,2,3]], 4, 3)) == [0 1 1 1; 1 0 2 0; 1 2 0 0; 1 0 0 0]
    @test to_adjacency_matrix(Hypergraphs([[1,2], [1,3], [3,3]], 3, 3)) == [0 1 1; 1 0 0; 1 0 2]

    # Conversion to simple
    @test to_adjacency_matrix(Hypergraphs([[1,4], [2,3], [1,2,3]], 4, 3); simple=true) == [0 1 1 1; 1 0 1 0; 1 1 0 0; 1 0 0 0]  
    @test to_adjacency_matrix(Hypergraphs([[1,2,3], [2,3,4], [1,4], [5,6]], 6, 4); simple=true) == [0 1 1 1 0 0; 1 0 1 1 0 0; 1 1 0 1 0 0; 1 1 1 0 0 0; 0 0 0 0 0 1; 0 0 0 0 1 0]
end

@testset "matrix_to_hyper_tests" begin # STUB
end

