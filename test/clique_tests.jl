@testset "k_cliques tests" begin

    @testset "triangle tests" begin
        @test length(kcliques(MatrixNetwork(sparse([0 1; 1 0])),3)[2]) == 0
        @test length(kcliques(MatrixNetwork(sparse(fill(0,10,10))),3)[2]) == 0
        @test length(kcliques(MatrixNetwork(sparse(fill(0,10,10))),3)[2]) == 0
        @test length(kcliques(MatrixNetwork(sparse(fill(0,10,10))),3)[2]) == 0 
        @test length(kcliques(MatrixNetwork(sparse([0 1 1; 1 0 1; 1 1 0])), 3)[2]) == 1  
        @test length(kcliques(load_data("facebook_data"), 3)[2]) == 1612010
        @test length(kcliques(MatrixNetwork(load_matrix_network("clique-10")), 3)[2]) == 120
    end

end

