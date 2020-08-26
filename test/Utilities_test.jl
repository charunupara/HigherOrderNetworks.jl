using Test
include("Utilities.jl")

@testset "Utilities" begin
    G1 = MatrixNetwork(load_matrix_network("cores_example")) # Simple networks
    G2 = MatrixNetwork([1,2,2,3,3], [2,4,5,4,5])
    uG2 = make_undirected(G2)

    M1 = MatrixNetwork([1,1,3,3,3,3,3,4,4,6,6,6], [3,1,1,2,2,4,6,6,6,1,2,7])
    M2 = MatrixNetwork([1,1,2,2,2,3,3,4,4,5,5,6,6,6,7,7,7],
                       [2,6,1,1,4,2,5,1,3,3,4,1,6,7,4,4,5])

    @testset "Degree Distribution" begin
        for i = 2:10 # Testing common subgraphs: path, cycle, clique
            @test deg_distr(make_undirected(path(i))) == [1; [2 for j = 1:i-1]; 1]
            @test deg_distr(make_undirected(cycle(i+1))) == [2 for j = 1:i+1]
            @test deg_distr(clique(i)) == [i-1 for j = 1:i]
        end

        @test deg_distr(G1) == [1,4,2,6,5,3,3,4,3,5,4,1,1,3,3,0,2,2,2,3,1]
        @test deg_distr(uG2) == [1,3,2,2,2]
    end

    @testset "Nodes by Degree" begin # Getting sorted list of nodes by their degree
        @test nodes_by_deg(G1) == [16, 1, 12, 13, 21, 3, 17, 18, 19, 6, 7, 9, 14, 15, 20, 2, 8, 11, 5, 10, 4]
        @test nodes_by_deg(G1; descending=true) == [4, 5, 10, 2, 8, 11, 6, 7, 9, 14, 15, 20, 3, 17, 18, 19, 1, 12, 13, 21, 16]

        @test nodes_by_deg(uG2) == [1, 3, 4, 5, 2]
        @test nodes_by_deg(uG2; descending=true) == [2, 3, 4, 5, 1]
    end

    @testset "Get Edges" begin
        function edge_compare(A::MatrixNetwork)
            ei, ej = directed_edges(A)
            return get_edges(A) == [[ei[x],ej[x]] for x = 1:length(ei)]
        end
        @test edge_compare(G1)
        @test edge_compare(G2)
        @test edge_compare(uG2)
        @test edge_compare(M1)
        @test edge_compare(M2)
    end

    @testset "Simplify" begin
        simplify_test_simple(A::MatrixNetwork) = get_edges(simplify(A)) == get_edges(A) # Checks that simplifying a simple graph keeps it simple
        @test simplify_test_simple(G1)
        @test simplify_test_simple(G2)

        function simplify_test(A::MatrixNetwork)
            M = Array(sparse(simplify(A)))
            return Diagonal(M) == zeros(Int64, (A.n, A.n)) && minimum(M) == 0 && maximum(M) == 1
        end

        @test simplify_test(M1)
        @test simplify_test(M2)
    end

    @testset "Directed -> Undirected" begin
        function undirected_test(A::MatrixNetwork)
            U = sparse(make_undirected(A))
            return U == U' == sparse(A) + sparse(A)'
        end

        @test undirected_test(G1)
        @test undirected_test(G2)
        @test undirected_test(uG2)

        @test undirected_test(M1)
        @test undirected_test(M2)
    end

    @testset "Clustering Coefficients" begin
        @test clustering(G1) == (0.3555555555555555, 0.44)
        @test_throws ErrorException clustering(G2)

        @test clustering(make_undirected(M1)) == (0.5898904756952629, 0.585695488070208)
        @test clustering(make_undirected(M2)) == (0.4950303252150675, 0.48426879640604437)
    end

end
