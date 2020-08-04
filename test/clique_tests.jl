# @testset "triangles tests" begin
#     @testset "emptiness tests" begin
#        @test  
#     end
# end

@testset "emptiness tests" begin
    @test length(k_cliques(MatrixNetwork(sparse(fill(0,10,10))), 2)) == 0
    @test length(k_cliques(MatrixNetwork(sparse(fill(0,10,10))), 3)) == 0
    @test length(k_cliques(MatrixNetwork(sparse(fill(0,10,10))), 10)) == 0 
end

