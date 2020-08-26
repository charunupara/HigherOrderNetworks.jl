using Test
include("RandomHypergraphs.jl")
include("load_data.jl")

@testset "Stub Matching" begin
    # Degree and edge-dimension sequences are correct
    DK(H::Hypergraph) = H.D, H.K
    sequence_match(D::Vector{Int64}, K::Vector{Int64}) = return DK(stub_matching(D,K)) == (D, sort(K, rev=true))

    # Realizable hypergraphs
    @test sequence_match([1,1,1,1], [2,2])
    @test sequence_match([1,2,3,4,5,6],[1,2,3,4,5,6])
    @test sequence_match([6,3,3,3], [4,3,3,1,1,1,1,1])
    @test sequence_match([1,2,3,1,1,1],[3,2,3,1])

    # Unrealizable hypergraphs
    @test_throws ArgumentError sequence_match([0],[5])
    @test_throws ArgumentError sequence_match([5,5,5,5,5],[5])
    @test_throws ArgumentError sequence_match([1,1,1,1],[4,4,4])

end
@testset "MCMC" begin
    @testset "Vertex" begin
        @test_throws AssertionError MCMC(Hypergraph(Array{Array{Int64,1},1}(),4,0))
        @test_throws AssertionError MCMC(Hypergraph([[1,2]],2,1))
        @test typeof(MCMC(Hypergraph([[1,2],[3,4],[1,4]], 4, 3))) <: Hypergraph
        @test typeof(MCMC(stub_matching([1,2,3,1,1,1],[3,2,3,1]); samples=1000)) <: Hypergraph
        @test typeof(MCMC(load_hypergraph("Enron"); samples=10000)) <: Hypergraph
    end

    @testset "Stub" begin
        @test_throws AssertionError typeof(MCMC(Hypergraph(Array{Array{Int64,1},1}(),4,0); type="stub")) <: Hypergraph
        @test typeof(MCMC(Hypergraph([[1,2],[3,4],[1,4]], 4, 3); type="stub")) <: Hypergraph
        @test typeof(MCMC(stub_matching([1,2,3,1,1,1],[3,2,3,1]); type="stub")) <: Hypergraph
        @test typeof(MCMC(load_hypergraph("Enron"); samples=10000, type="stub")) <: Hypergraph
    end
end
