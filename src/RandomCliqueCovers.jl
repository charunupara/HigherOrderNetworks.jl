using SpecialFunctions
using MatrixNetworks
using LinearAlgebra
using SparseArrays
using Distributions: Poisson
include("PlotGraph.jl")

"""
`random_clique_cover`
=====================

Generate a random graph exhibiting global sparsity and local density by randomly generating unbounded cliques.
Based on https://arxiv.org/pdf/1810.06738.pdf.

Arguments
---------
    - `a::Float64`: Controls the expected number of vertices per clique
    - `s::Float64`: Controls how close the clique size distribution is to a power law
    - `c::Float64`: Controls the expected number of cliques each vertex is a part of. Larger value = less overlap
    - `N (=rand(Poisson(τ))`: The total number of cliques

"""
function random_clique_cover(a::Float64, s::Float64, c::Float64; N=rand(Poisson(2*pi)))
    @assert a > 0
    @assert s >= 0 && s < 1
    @assert c > -s
    p(x) = rand(Poisson(x)) # Sample from Poisson distribution

    cliques = [[i for i = 1:p(a)]] # Stores cliques
    clique_count = [1 for i = 1:size(cliques[1],1)] # Stores how many cliques each node is part of
    #println(cliques)
    #println(clique_count)

    for n = 1:N-1
        push!(cliques, []) # New empty clique
        for k = 1:size(clique_count,1) # Decide which of existing nodes to place in new clique
            if rand(Float64) <= (clique_count[k]-s) / (n+c)
                push!(cliques[n+1], k)
                clique_count[k] += 1
            end
        end

        for j = 1:p(a*((gamma(c+1) * gamma(n+c+s) / (gamma(c+s) * gamma(n+c+1))))) # Add new nodes to graph and clique
            push!(cliques[n+1], size(clique_count,1)+1)
            push!(clique_count, 1)
        end
    end
    #println(cliques)
    #println(clique_count)

    Z = zeros(Int64, N, size(clique_count,1)) # Convert clique list to edge clique cover. Z[i,j] = 1 iff clique i contains node j
    for i = 1:N
        for j in cliques[i]
            Z[i,j] = 1
        end
    end
    #println(Z)

    mat = Z' * Z # Convert edge clique cover to corresponding adjacency matrix
    for i = 1:size(mat,1)
        mat[i,i] = 0
    end
    mat = min.(mat, 1)
    sp = sparse(mat)
    #graphplot(sp, igraph_layout(sp)) uncomment this line and comment below line to see graph
    return MatrixNetwork(sp)
end

random_clique_cover(5.2, 0.5, 0.4, N=5)
