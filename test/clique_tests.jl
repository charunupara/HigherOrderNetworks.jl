@testset "triangles tests" begin

    @testset "emptiness tests" begin
        @test length(triangles(MatrixNetwork(sparse([0 1; 1 0])))[1]) == 0
        @test length(triangles(MatrixNetwork(sparse(fill(0,10,10))))[1]) == 0
        @test length(triangles(MatrixNetwork(sparse(fill(0,10,10))))[1]) == 0
        @test length(triangles(MatrixNetwork(sparse(fill(0,10,10))))[1]) == 0 
    end

    @testset "small cases" begin
        @test length(triangles(MatrixNetwork(sparse([0 1 1; 1 0 1; 1 1 0])))[1]) == 1     
    end

    @testset "data sets tests" begin
       @test length(triangles(load_data("facebook_data"))[1]) == 1612010
       @test length(triangles(MatrixNetwork(load_matrix_network("clique-10")))[1]) == 120
    end

end

@testset "k_cliques tests" begin

    @testset "emptiness tests" begin
        @test length(k_cliques(MatrixNetwork(sparse([0 1; 1 0])), 3)) == 0
        @test length(k_cliques(MatrixNetwork(sparse(fill(0,10,10))), 2)) == 0
        @test length(k_cliques(MatrixNetwork(sparse(fill(0,10,10))), 3)) == 0
        @test length(k_cliques(MatrixNetwork(sparse(fill(0,10,10))), 10)) == 0 
    end
    
    
end

