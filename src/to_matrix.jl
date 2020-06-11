include("Hypergraphs.jl")
include("RandomHypergraphs.jl")
using MatrixNetworks

A = VertexHypergraph([[1,2,3,4], [3,4,1], [1,4,5]], ["One", "Two", "Three", "Four", "Five"], 5, 3)



function to_matrix(A::VertexHypergraph)
    B = zeros(Int64, graph.n, graph.n) # create adjacency matrix of size n x n

    # loop through each hyperedge and create edges in the new adjacency matrix
    for i in 1:graph.m
        for j in 1:size(graph.edges[i],1)
            for k in 1:size(graph.edges[i],1)
                if j == k
                    continue # pass because we don't want self-loops
                else
                    index1 = graph.edges[i][j]
                    index2 = graph.edges[i][k]  
                    B[index1,index2] = 1
                end
            end
        end
    end
    
    return B
end


