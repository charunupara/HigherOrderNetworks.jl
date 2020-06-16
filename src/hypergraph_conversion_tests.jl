include("hypergraphs_conversions.jl")
using Test

@testset "hyper_to_matrix_tests" begin
    @test to_adjacency_matrix(VertexHypergraph(Array{Int64,1}[], 1, 0))[1,:] == [0]
    @test to_adjacency_matrix(StubHypergraph([[1.5]], 1, 1))[1,:] == [1] # loopy multigraphs as well?

    @test to_adjacency_matrix(VertexHypergraph([[1,2,3], [2,3,4], [1,4], [5,6]], 6, 4)) == [0 1 1 1 0 0; 1 0 1 1 0 0; 1 1 0 1 0 0; 1 1 1 0 0 0; 0 0 0 0 0 1; 0 0 0 0 1 0] 
    @test to_adjacency_matrix(StubHypergraph([[1.5,4.5], [2.5,3.5], [4/3,7/3,10/3]], 4, 3)) == [0 1 1 1; 1 0 1 0; 1 1 0 0; 1 0 0 0] 

    stub = StubHypergraph([[1.5,10/3,4.5], [2.5,3.25,6.25], [4/3,5.5,19/3], [3.5,13/3,6.5]], 6, 4)
    @test to_adjacency_matrix(stub) == to_adjacency_matrix(VertexHypergraph(stub))
end

@testset "matrix_to_hyper_tests" begin # STUB
end

