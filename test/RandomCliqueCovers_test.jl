using Test
include("RandomCliqueCovers.jl")

@testset "rcc_test" begin
    A = [0.1, 0.5, 3., 15.]
    S = collect(0:0.1:0.9)
    C = [0.01, 0.1, 0.5, 1., 10.]
    N = [1, 5, 10, 100]
    for a in A
        for s in S
            for c in C
                for n in N
                    @test typeof(random_clique_cover(a,s,c; N=n)) <: MatrixNetwork
                end
            end
        end
    end
end
