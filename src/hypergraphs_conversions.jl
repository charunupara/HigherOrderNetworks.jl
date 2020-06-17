include("Hypergraphs.jl")
include("RandomHypergraphs.jl")
using MatrixNetworks
using SparseArrays


"""
`to_adjacency_matrix`
=====================

Convert a VertexHypergraph to an adjacency matrix

Arguments
----------
    - `h::Hypergraphs`: The hypergraph that will be converted into an adjacency matrix
    - `simple::Bool (=false)`: Whether to exclude multi-edges and self-loops from the matrix (they are included by default)

Examples
----------
~~~~
to_adjacency_matrix(VertexHypergraph([[1,2,3,4], [3,4,1], [1,4,5]], ["One", "Two", "Three", "Four", "Five"], 5, 3); simple=true)
~~~~
"""
function to_adjacency_matrix(h::Hypergraphs; simple::Bool=false)
    adjacency_matrix = zeros(Int64, h.n, h.n) # Create adjacency matrix of size n x n

    # Loop through each hyperedge and create edges in the new adjacency matrix
    for i = 1:h.m
        for j = 1:size(h.edges[i],1)-1
            for k = j+1:size(h.edges[i],1)
                index1 = vertex_floor(h.edges[i][j])
                index2 = vertex_floor(h.edges[i][k])
                if simple
                    adjacency_matrix[index1,index2] = adjacency_matrix[index2,index1] = 1 # Setting to 1 excludes multi-edges
                else
                    adjacency_matrix[index1,index2] += 1 # Incrementing tracks edge multiplicity
                    adjacency_matrix[index2,index1] += 1
                end
            end
        end
    end

    if simple # Remove self-loops
        for i = 1:h.n
            adjacency_matrix[i,i] = 0
        end
    end

    return adjacency_matrix
end



"""
`hypergraph_to_matrixnetworks`
=====================

Convert a VertexHypergraph to a MatrixNetwork

Arguments
----------
    - `h::Hypergraphs`: The hypergraph that will be converted into a MatrixNetwork
    - `simple::Bool (=false)`: Whether to exclude multi-edges and self-loops from the MatrixNetwork (they are included by default)

Examples
----------
~~~~
hypergraph_to_matrixnetwork(VertexHypergraph([[1,2,3,4], [3,4,1], [1,4,5]], ["One", "Two", "Three", "Four", "Five"], 5, 3); simple=true)
~~~~
"""
function hypergraph_to_matrixnetwork(h::Hypergraphs; simple=false)
    return MatrixNetwork(sparse(to_adjacency_matrix(h; simple=simple)))
end
