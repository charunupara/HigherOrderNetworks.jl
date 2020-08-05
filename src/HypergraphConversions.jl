"""
Convert from a hypergraph to a dyadic graph by turning hyperedges into cliques.
Convert between hypergraph and bipartite representations.
"""

using MatrixNetworks
using SparseArrays
using LinearAlgebra
include("Hypergraphs.jl")
include("RandomHypergraphs.jl")


"""
`to_adjacency_matrix`
=====================

Convert a hypergraph to an adjacency matrix

Arguments
----------
    - `H::Hypergraphs`: The hypergraph that will be converted into an adjacency matrix
    - `simple::Bool (=false)`: Whether to exclude multi-edges and self-loops from the matrix (they are included by default)

Examples
----------
~~~~
to_adjacency_matrix(Hypergraphs([[1,2,3,4], [3,4,1], [1,4,5]], ["One", "Two", "Three", "Four", "Five"], 5, 3); simple=true)
~~~~
"""
function to_adjacency_matrix(H::Hypergraphs; simple::Bool=false)
    A = zeros(Int64, H.n, H.n) # Create adjacency matrix of size n x n

    # Loop through each hyperedge and create edges in the new adjacency matrix
    for i = 1:H.m
        for j = 1:length(H.edges[i])-1
            for k = j+1:length(H.edges[i])
                index1 = h.edges[i][j]
                index2 = h.edges[i][k]
                A[index1,index2] += 1
                A[index2,index1] += 1
            end
        end
    end

    if simple
        A -= Diagonal(A) # Remove self-loops
        A = min.(A, 1) # Remove multi-edges
    end

    return A
end



"""
`hypergraph_to_matrixnetworks`
=====================

Convert a hypergraph to a MatrixNetwork

Arguments
----------
    - `H::Hypergraphs`: The hypergraph that will be converted into a MatrixNetwork
    - `simple::Bool (=false)`: Whether to exclude multi-edges and self-loops from the MatrixNetwork (they are included by default)

Examples
----------
~~~~
hypergraph_to_matrixnetwork(Hypergraphs([[1,2,3,4], [3,4,1], [1,4,5]], ["One", "Two", "Three", "Four", "Five"], 5, 3); simple=true)
~~~~
"""
function hypergraph_to_matrixnetwork(H::Hypergraphs; simple=false)
    return MatrixNetwork(sparse(to_adjacency_matrix(H; simple=simple)))
end
