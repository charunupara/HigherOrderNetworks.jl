"""
Hypergraph creation and manipulation.
"""

using MatrixNetworks
using StatsBase
include("Utilities.jl")

"""
`Hypergraphs{T}`
=====================

Represents a vertex-labeled Hypergraphs as a list of edges.

Fields
------
   - `vals::Vector{T}`: The values corresponding to each node (or hyperedge, depending on use case)
   - `n::Int64`: The number of nodes
   - `m::Int64`: The number of hyperedges
   - `D::Vector{Int64}`: The degree sequence, where `D[i] = degree of node i`
   - `K::Vector{Int64}`: The edge dimension sequence, where `K[j] = size of edge j`
   - `edges::Vector{Vector{Int64}}`: The hyperedges and their members

A hypergraph is a generalization of a graph in which edges can contain any number
of nodes. So, for example, if we wanted to denote a three-way relationship between
nodes 1, 2, and 3, the hypergraph would contain the edge `{1, 2, 3}`.

Hypergraphs are useful because they allow us to more accurately indicate relationships
among groups of things. Consider a coauthorship network. In a dyadic graph,
we might represent a five-way collaboration as a 5-clique. While this does capture
the fact that all pairs of authors have now appeared in a paper together, it doesn't seem
right to categorize the event as 15 pairwise interactions. With a hypergraph, we
can more succinctly and intuitively represent the event as a single hyperedge of five nodes.
"""
mutable struct Hypergraphs{T}
   vals::Vector{T}
   n::Int64
   m::Int64
   D::Vector{Int64}
   K::Vector{Int64}
   edges::Vector{Vector{Int64}}
end

"""
`Hypergraphs_kernel`
===================

Verifies that:
   - The number of edges == `m`
   - The number of vals == `n`
   - All nodes are between 0 and `n`

If all conditions are met, returns degree and edge dimension sequences
"""
function Hypergraphs_kernel(edges::Vector{Vector{Int64}}, vals::Vector{T},
                           n::Int64, m::Int64) where T
   @assert length(edges) == m # m is the number of edges
   @assert length(vals) == n # Each node has a val associated with it
   @assert all([0 < edges[i][j] < n + 1 for i = 1:m for j = 1:length(edges[i])]) # No node exceeds n

   D = zeros(Int64, n)
   K = zeros(Int64, m)
   for e = 1:m
      K[e] = size(edges[e],1)
      for v in edges[e]
         D[v] += 1
      end
   end
   edges = sort!.(edges, by=v -> D[v], rev=true) # Sort edges by descending node degree

   return D, K
end

"""
`Hypergraphs` constructors
=============================

Functions
---------
   - `Hypergraphs(edges, vals, n, m)`: Produces a hypergraph with the given edges, values, and size
   - `Hypergraphs(edges, n, m)`: Produces a hypergraph with the given edge set and size, with all values as 1.0
   - `Hypergraphs(edges)`: Produces a hypergraph with the given set of edges

Examples
--------
~~~~
Hypergraphs([[1,2], [3,4], [1,4]], ["One", "Two", "Three", "Four"], 4, 3)
Hypergraphs([[i,i+2] for i=1:3], 4, 3)
Hypergraphs([[1,2,3], [2,4,6], [1,5]])
~~~~
"""
function Hypergraphs(edges::Vector{Vector{Int64}}, vals::Vector{T},
                          n::Int64, m::Int64) where T
   D, K= Hypergraphs_kernel(edges, vals, n, m)
   return Hypergraphs(vals, n, m, D, K, edges)
end

function Hypergraphs(edges::Vector{Vector{Int64}}, n::Int64, m::Int64)
   return Hypergraphs(edges, ones(n), n, m)
end

function Hypergraphs(edges::Vector{Vector{Int64}})
   return Hypergraphs(edges, maximum([e[i] for e in edges for i = 1:length(e)]), length(edges))
end

Base.copy(H::Hypergraphs) = Hypergraphs(deepcopy(h.edges), deepcopy(h.vals), h.n, h.m)

"""
`add_node!`
===========

Add a node with a specific value to the hypergraph.

Arguments
---------
   - `H::Hypergraphs`: The hypergraph to add to
   - `val::T`: The value to associate with the new node
"""
function add_node!(H::Hypergraphs, val::T) where T
   H.n += 1
   push!(H.D, 0)
   push!(H.vals, val)
end

"""
`add_node!`
===========

Add a node to the hypergraph.
"""
function add_node!(H::Hypergraphs)
   add_node!(H, nothing)
end

"""
`add_edge!`
===========

Add a hyperedge to a hypergraph.

Arguments
---------
   - `V::Hypergraphs`: The hypergraph to add an edge to
   - `e::Vector{Int64}`: The edge to add
"""
function add_edge!(V::Hypergraphs, e::Vector{Int64})
   push!(V.edges, e)
   push!(V.K, length(e))

   for v in e
      if v > V.n
         add_node!(V)
      end
      V.D[v] += 1
   end
end

"""
`add_value_edge!`
=================

Add a hyperedge to a hypergraph by specifying the corresponding values to be linked.

Arguments
---------
   - `V::Hypergraphs`: The hypergraph to add an edge to
   - `e::Vector{T}`: The values to connect

Example
-------
~~~~
add_value_edge!(V, ["Kim", "Jim", "Tim"])
~~~~
"""
function add_edge_by_value!(V::Hypergraphs, e::Vector{T}) where T
   node_edge = []

   for v in e
      i = findfirst(x -> x==v, V.vals)
      if i == -1
         add_node!(V, v)
         i = V.n + 1
      end
      V.D[i] += 1
      push!(node_edge, i)
   end

   push!(V.edges, node_edge)
   push!(V.K, length(node_edge))
end

"""
`remove_node!`
==============

Delete the node at index `n` from a hypergraph.
"""
function remove_node!(H::Hypergraphs, n::Int64)
   for i = 1:H.m
      filter!(x -> x == n, H.edges[i]) # Remove from edge
      map!(x -> x > n ? x - 1 : x, H.edges[i]) # Adjust node identities
      if H.edges[i] == [] # Clear empty edges
         deleteat!(H.edges, i)
         deleteat!(H.K, i)
      end
      H.K[i] -= count(x -> x==n, H.edges[i]) # Edit edge dimension sequence
   end

   deleteat!(H.D, n) # Remove from degree sequence

   H.n -= 1
   H.m = length(H.edges)
end

"""
`remove_edge!`
==============

Delete the edge at index `e` from a hypergraph.
"""
function remove_edge!(H::Hypergraphs, e::Int64)
   for n in H.edges[e]
      H.D[n] -= 1 # Edit degree sequence
   end
   deleteat!(H.edges, e) # Remove edge
   deleteat!(H.K, e)
   H.m -= 1
end

"""
`edge_intersect`
================

Get the set of intersection between two hyperdges. Returns a set of integers.

Arguments
---------
   - `H::Hypergraphs`: The hypergraph where the edges live
   - `e1::Int64`: The index of the first edge
   - `e2::Int64`: The index of the second edge

Examples
--------
~~~~
edge_intersect(Hypergraphs([[1,3,5], [2,3,5], [3,5]], 1, 2)) -> {3, 5}
~~~~
"""
function edge_intersect(H::Hypergraphs, e1::Int64, e2::Int64)
   return intersect(H.edges[e1], H.edges[e2])
end

"""
`num_parallel`
==============

Get the number of hyperedges parallel (multiset equal) to a given edge.

Arguments
---------
   - `H::Hypergraphs`: The hypergraph to be analyzed
   - `e::Int64`: The index of the edge
"""
function num_parallel(H::Hypergraphs, e::Int64)
   return count(x -> x == H.edges[e], H.edges)
end

"""
`random_edges`
==============

Select `n` distinct random edges from a hypergraph.

Example
========
~~~~
random_edges(G, 3) -> [[e1], [e2], [e3]]
~~~~
"""
function random_edges(H::Hypergraphs, n::Int64)
   return H.edges[random_edge_indices(H, n)]
end

"""
`random_edge_indices`
==============

Select `n` distinct random edge indices from a hypergraph.

Example
========
~~~~
random_edge_indices(G, 3) -> [e1, e2, e3]
~~~~
"""
function random_edge_indices(H::Hypergraphs, n::Int64)
   return sample(1:H.m, n, replace=false)
end
