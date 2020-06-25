using SpecialFunctions
using MatrixNetworks
using LinearAlgebra
using SparseArrays
using Distributions: Poisson
using Statistics
using PyCall
using PyPlot
include("Utilities.jl")
#include("PlotGraph.jl")

"""
`random_clique_cover`
=====================

Generate a random graph exhibiting global sparsity and local density by randomly generating unbounded cliques.
Based on https://arxiv.org/pdf/1810.06738.pdf.

Arguments
---------
    - `a::Float64 (> 0)`: Controls the expected number of vertices per clique
    - `s::Float64 (∈ [0,1))`: Controls the skewness of the degree distribution
    - `c::Float64 (> -s)`: Controls the expected number of cliques each vertex is a part of. Larger value = less overlap
    - `N (=rand(Poisson(τ))`: The total number of cliques

"""
function random_clique_cover(a::Float64, s::Float64, c::Float64; N=rand(Poisson(2*pi)))
    @assert a > 0
    @assert s >= 0 && s < 1
    @assert c > -s
    P(x) = rand(Poisson(x)) # Sample from Poisson distribution with λ = x

    function gamma_ratio_approx(n) # Permits large N, since it avoids using large numbers
        radical = ((ℯ^(1-s)) * ((n+c+s)^(s+0.5)) / ((n+c+1)^1.5)) ^ (1/(n+c))
        return (radical * (n+c+s)/(n+c+1)) ^ (n+c) # Derived from Stirling's approximation for gamma
    end

    agammr = a*gamma(c+1) / gamma(c+s)

    cliques = [[i for i = 1:P(a)]] # Stores cliques
    clique_count = [1 for i = 1:size(cliques[1],1)] # Stores how many cliques each node is part of

    for n = 1:N-1
        push!(cliques, []) # New empty clique
        for k = 1:size(clique_count,1) # Decide which of existing nodes to place in new clique
            if rand(Float64) <= (clique_count[k]-s) / (n+c)
                push!(cliques[n+1], k)
                clique_count[k] += 1
            end
        end

        for j = 1:P(agammr*gamma_ratio_approx(n)) # Add new nodes to graph and clique
            push!(cliques[n+1], size(clique_count,1)+1)
            push!(clique_count, 1)
        end
    end

    Z = zeros(Int64, N, size(clique_count,1)) # Convert clique list to edge clique cover. Z[i,j] = 1 iff clique i contains node j
    for i = 1:N
        for j in cliques[i]
            Z[i,j] = 1
        end
    end

    mat = Z' * Z # Convert edge clique cover to corresponding adjacency matrix
    for i = 1:size(mat,1)
        mat[i,i] = 0
    end
    mat = min.(mat, 1)
    return MatrixNetwork(sparse(mat))
end

pygui(true)
a = 10.
s = 0.9
c = 0.1
N = 100
@time rcc = random_clique_cover(a,s,c;N=N)
#graphplot(sparse(rcc), igraph_layout(sparse(rcc)))
#plot(1:rcc.n, sort(deg_distr(rcc), rev=true))
#boxplot([deg_distr(random_clique_cover(a,s,c;N=N)) for s in 0.:0.1:0.9])
