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
                index1 = H.edges[i][j]
                index2 = H.edges[i][k]
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


"""
`hyper_to_bipartite`
====================

Converts a hypergraph into its bipartite representation. For all vertices v and
edges e in the hypergraph, an edge (v, e) exists in the bipartite if and only if
v ∈ e in the hypergraph.
"""
function hyper_to_bipartite(H::Hypergraphs)
   di_edges::Vector{Vector{Int64}} = []

   for i = 1:h.m
      for j in h.edges[i]
         push!(di_edges, [i + h.n, j])
      end
   end

   return MatrixNetwork(di_edges, h.n + h.m)
end

"""
`is_bipartite`
==============

Checks whether a graph is bipartite by checking whether it is 2-colorable.
"""
function is_bipartite(A::MatrixNetwork)
   group1::Set{Int64} = Set()
   group2::Set{Int64} = Set()

   function color_neighbors(v::Int64, is_group1::Bool)
      for n in A.ci[A.rp[v]:A.rp[v+1]-1]
         if is_group1 # v is in group 2, so all neighbors should be in group 1
            if n in group2
               return false
            end
            push!(group1, n)
         else
            if n in group1
               return false
            end
            push!(group2, n)
         end
      end
      return true
   end

   is_group1 = true

   for i = 1:A.n
      is_group1 = !is_group1
      if !color_neighbors(i, is_group1)
         return false
      end
   end
   return (group1, group2)
end

"""
`bipartite_to_hyper`
====================

Returns the left and right projections of a bipartite graph into a hypergraph.
In the left projection, all vertices in group two that are adjacent to the same
vertex in group one are put together in a hyperedge. The right projection is the
same except groups one and two are swapped in their roles.
"""
function bipartite_to_hyper(A::MatrixNetwork)
   b = is_bipartite(A)
   if typeof(b) <: Bool
      error("ArgumentError: The network must be bipartite to be projected to a hypergraph.")
   end

   g1, g2 = collect.(b)
   left_edges = [map(x -> findfirst(y -> x == y, g2), A.ci[A.rp[i]:A.rp[i+1]-1]) for i in g1]
   right_edges = [map(x -> findfirst(y -> x == y, g1), A.ci[A.rp[j]:A.rp[j+1]-1]) for j in g2]

   return Hypergraphs(left_edges, A.vals[g2], length(g2), length(left_edges)),
          Hypergraphs(right_edges, A.vals[g1], length(g1), length(right_edges))
end
