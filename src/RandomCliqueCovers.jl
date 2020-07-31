"""
Random graph model exhibiting low GCC and high ALCC. Creates a fixed number
of cliques each with a random number of nodes. Nodes can partake in multiple
cliques.

Based on "Random clique covers for graphs with local density and global sparsity"
by Sinead A. Williamson and Mauricio Tec. The original paper may be found
at https://arxiv.org/pdf/1704.03913.pdf.
"""

using MatrixNetworks
using SpecialFunctions # Used for gamma
using LinearAlgebra # Matrix transpose
using SparseArrays
using Distributions: Poisson

"""
`rcc`
=====================

Generate a random graph by creating random clique covers.

Arguments
---------
    - `a::Float64 (> 0)`: The expected number of vertices per clique
    - `s::Float64 (∈ [0,1))`: Controls how close the degree distribution is to a power law
    - `c::Float64 (> -s)`: Controls (but not equal to) the expected number of cliques each vertex is a part of.
                           Larger c = less overlap
    - `N::Int64 (=10)`: The total number of cliques

Preconditions
-------------
    - `a > 0`
    - `s ∈ [0,1)`
    - `c > -s`
    - `N >= 0`

Examples
--------
~~~~
random_clique_cover(10.,0.9,0.1;N=100) # 100-clique graph with ~10 nodes per clique with high degree skew and medium clique overlap
random_clique_cover(50.,0.1,5.) # 10-clique graph with ~50 nodes per clique with more uniform degree distribution and high clique overlap
~~~~

The random cliques are generated using the following process. Imagine a fixed
number of restaurant customers coming up to an (infinitely long) buffet one by one.
The first customer takes the first n dishes, where n is random. Each customer
thereafter
    1) samples each already-tried dish at random, with the probability
    of taking a dish decreasing with the number of times it has been tried
    2) selects the next k (also random) dishes that have not yet been tried.

This graph model works by substituting customers for cliques and dishes for nodes. Bon appétit!
"""
function rcc(a::Float64, s::Float64, c::Float64; N::Int64=10)
    @assert a > 0
    @assert s >= 0 && s < 1
    @assert c > -s
    @assert N >= 0

    P(x) = rand(Poisson(x)) # Sample from Poisson distribution with λ = x

    ### Approximate Γ(n+c+s−1)/Γ(n+c) to avoid overflow with large n ###
    function gamma_ratio_approx(n) # Permits large N, since it avoids using large numbers
        radical = ((ℯ^(1-s)) * ((n+c+s)^(s+0.5)) / ((n+c+1)^1.5)) ^ (1/(n+c))
        return (radical * (n+c+s)/(n+c+1)) ^ (n+c) # Derived from Stirling's approximation for gamma
    end

    agammr = a*gamma(c+1) / gamma(c+s) # Constant term for determining the number of new nodes to add to a clique

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

    Z = zeros(Int64, N, size(clique_count,1)) # Convert clique list to edge clique cover.
                                              # Z[i,j] = 1 iff clique i contains node j.
    for i = 1:N
        for j in cliques[i]
            Z[i,j] = 1
        end
    end

    mat = Z' * Z # Convert edge clique cover to corresponding adjacency matrix
    mat -= Diagonal(mat) # Simplify matrix
    mat = min.(mat, 1)

    return MatrixNetwork(sparse(mat))
end
