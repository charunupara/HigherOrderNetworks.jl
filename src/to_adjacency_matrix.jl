include("Hypergraphs.jl")
include("RandomHypergraphs.jl")
using MatrixNetworks

# TODO: Convert to MatrixNetworks type

"""
`to_adjacency_matrix`
=====================

Convert a VertexHypergraph to an adjacency matrix

Arguments
----------
    - `A::VertexHypergraph : The VertexHypergraph that will be converted into an adjacency matrix`

Examples
----------
~~~~
to_adjacency_matrix(VertexHyperA([[1,2,3,4], [3,4,1], [1,4,5]], ["One", "Two", "Three", "Four", "Five"], 5, 3))
~~~~
"""
function to_adjacency_matrix(A::VertexHypergraph)
    adjacency_matrix = zeros(Int64, A.n, A.n) # create adjacency matrix of size n x n

    # loop through each hyperedge and create edges in the new adjacency matrix
    for i in 1:A.m
        for j in 1:size(A.edges[i],1)
            for k in 1:size(A.edges[i],1)
                if j == k
                    continue # pass because we don't want self-loops
                else
                    index1 = A.edges[i][j]
                    index2 = A.edges[i][k]  
                    adjacency_matrix[index1,index2] = 1
                end
            end
        end
    end
    return adjacency_matrix
end


