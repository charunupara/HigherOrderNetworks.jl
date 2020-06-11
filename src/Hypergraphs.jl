using MatrixNetworks
using StatsBase
include("Utilities.jl")

abstract type Hypergraphs end

"""
`VertexHypergraph{T}`
=====================
Represents a vertex-labeled hypergraph as a list of edges.

Fields
------
   - `vals::Vector{T}`: The values corresponding to each node
   - `n::Int64`: The number of nodes
   - `m::Int64`: The number of hyperedges
   - `D::Vector{Int64}`: The degree sequence
   - `K::Vector{Int64}`: The edge dimension sequence
   - `edges::Vector{Vector{Int64}}`: The hyperedges and their members
"""
mutable struct VertexHypergraph{T} <: Hypergraphs
   vals::Vector{T}
   n::Int64
   m::Int64
   D::Vector{Int64}
   K::Vector{Int64}
   edges::Vector{Vector{Int64}}
end

"""
`StubHypergraph{T}`
===================

Represents a stub-labeled hypergraph as a list of edges.

Fields
------
   - `vals::Vector{T}`: The values corresponding to each node
   - `n::Int64`: The number of nodes
   - `m::Int64`: The number of hyperedges
   - `D::Vector{Int64}`: The degree sequence
   - `K::Vector{Int64}`: The edge dimension sequence
   - `edges::Vector{Vector{Float64}}`: The hyperedges and their members
Note: The nth stub of node i is represented as `i + \\tfrac{1}{n+1}`.
E.g., `2_2 \\rightarrow \\tfrac{7}{3}` and `4_1 \\rightarrow \\tfrac{9}{2}`
"""
mutable struct StubHypergraph{T} <: Hypergraphs
   vals::Vector{T}
   n::Int64
   m::Int64
   D::Vector{Int64}
   K::Vector{Int64}
   edges::Vector{Vector{Float64}}
end

"""
`Hypergraph_kernel`
===================

Verifies that:
   - The number of nodes == `m`
   - The number of vals == `n`
   - All nodes are between 0 and `n`
   - If StubHypergraph, that all stub numberings are valid and none are skipped

If all conditions are met, returns degree and edge dimension sequences
"""
function Hypergraph_kernel(edges::Vector{Vector{Te}}, vals::Vector{T},
                           n::Int64, m::Int64) where Te <: Real where T
   @assert size(edges, 1) == m # m is the number of edges
   @assert size(vals, 1) == n # Each node has a val associated with it
   @assert all([0 < edges[i][j] < n + 1 for i = 1:m for j = 1:size(edges[i], 1)]) # No node exceeds n

   count(v, e) = sum([i == e ? 1 : 0 for k = 1:m for i in vertex_floor.(v[k])])

   if Te == Float64
      possible_stubs = [i + 1/(j+1) for i = 1:n for j = 1:count(edges, i)]
      @assert all([i in possible_stubs for k = 1:m for i in edges[k]]) # Valid stub numberings
   end

   D::Vector{Int64} = [count(edges, v) for v = 1:n]
   K::Vector{Int64} = [size(edges[e], 1) for e = 1:m]
   edges = map(sort!, edges) # Sort edge multisets

   return D, K
end

"""
VertexHypergraph constructors
=============================

Functions
---------
   - `VertexHypergraph(edges, vals, n, m)`: Produces a vertex-labeled hypergraph with the given edge set, values, and size
   - `VertexHypergraph(edges, n, m)`: Produces a vertex-labeled hypergraph with the given edge set and size, with all values as 1.0
   - `VertexHypergraph(s)`: Converts a stub-labeled hypergraph into a vertex-labeled hypergraph

Examples
--------
~~~~
VertexHypergraph([[1,2], [3,4], [1,4]], ["One", "Two", "Three", "Four"], 4, 3)
VertexHypergraph([[i,i+2] for i=1:3], 4, 3)
VertexHypergraph(StubMatching([3,2,2,2], [4,3,2]))
~~~~
"""
function VertexHypergraph(edges::Vector{Vector{Int64}}, vals::Vector{T},
                          n::Int64, m::Int64) where T
   D, K = Hypergraph_kernel(edges, vals, n, m)
   return VertexHypergraph(vals, n, m, D, K, edges)
end

VertexHypergraph(edges::Vector{Vector{Int64}}, n::Int64, m::Int64) = VertexHypergraph(edges, ones(n), n, m)
VertexHypergraph(s::StubHypergraph) = VertexHypergraph([vertex_floor.(s.edges[i]) for i = 1:s.m], s.vals, s.n, s.m)

"""
StubHypergraph constructors
===========================

Functions
---------
   - `StubHypergraph(edges, vals, n, m)`: Produces a stub-labeled hypergraph with the given edge set, values, and size
   - `StubHypergraph(edges, n, m)`: Produces a stub-labeled hypergraph with the given edge set and size, with all values as 1.0

Examples
--------
~~~~
StubHypergraph([[1.33333333,2.5], [3.5,4.5], [1.5,4.33333333]], ["One", "Two", "Three", "Four"], 4, 3)
StubHypergraph([[1.5,2.33333333], [1.33333333,2.5]], 2, 2)
~~~~
"""
function StubHypergraph(edges::Vector{Vector{Float64}}, vals::Vector{T},
                        n::Int64, m::Int64) where T
   D, K = Hypergraph_kernel(edges, vals, n, m)
   return StubHypergraph(vals, n, m, D, K, edges)
end

StubHypergraph(edges::Vector{Vector{Float64}}, n::Int64, m::Int64) = StubHypergraph(edges, ones(n), n, m)

"""
`remove!`
=========

Functions
---------
   - `remove!(h, n)`: Delete the node at index n from a hypergraph
   - `remove!(h, e)`: Delete the edge at index e from a hypergraph
"""
function remove!(h::Hypergraphs, n::Int64)
   for i = 1:m
      filter!(x -> vertex_floor(x) == n, h.edges[i]) # Remove from edge
      map!(x -> x > n ? x - 1 : x, h.edges[i]) # Adjust node identities
      if h.edges[i] == [] # Clear empty edges
         deleteat!(h.edges, i)
         deleteat!(h.K, i)
      end
      h.K[i] -= count(x -> x == n, vertex_floor.(h.edges[i])) # Edit edge dimension sequence
   end

   deleteat!(h.D, n) # Remove from degree sequence

   h.n -= 1
   h.m = size(h.edges,1)
end

function remove!(h::Hypergraphs, e::Int64)
   for n in h.edges[e]
      h.D[vertex_floor(n)] -= 1 # Edit degree sequence
   end
   deleteat!(h.edges, e) # Remove edge
   h.m -= 1
end

"""
`edge_intersect`
================

Get the set of intersection between two hyperdges. Returns a set of integers.

Arguments
---------
   - `h::Hypergraphs`: The hypergraph where the edges live
   - `e1::Int64`: The index of the first edge
   - `e2::Int64`: The index of the second edge

Examples
--------
~~~~
edge_intersect(VertexHypergraph([[1,3,5], [2,3,5], [3,5]], 1, 2)) -> {3, 5}
edge_intersect(StubHypergraph([[1.5,3.33333333], [2.5,4.5,5.5], [3.5,5.33333333]], 1, 3)) -> {3}
~~~~
"""
edge_intersect(h::Hypergraphs, e1::Int64, e2::Int64) = Set(filter(x -> x in vertex_floor.(h.edges[e2]), vertex_floor.(h.edges[e1])))

"""
`num_parallel`
==============

Get the number of hyperedges parallel (multiset equal) to a given edge.

Arguments
---------
   - `h::Hypergraphs`: The hypergraph to be analyzed
   - `e::Int64`: The index of the edge
"""
num_parallel(h::Hypergraphs, e::Int64) = count(x -> vertex_floor.(x) == vertex_floor.(h[e]), h.edges)

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
function random_edges(h::Hypergraphs, n::Int64)
   @assert 0 < n <= h.m # The number to be sampled is valid

   edges_copy = copy(h.edges)
   for i = 1:h.m - n
      deleteat!(edges_copy, rand(1:size(edges_copy,1)))
   end

   return edges_copy
end

"""
`random_edges`
==============

Select `n` distinct random edge indices from a hypergraph.

Example
========
~~~~
random_edge_indices(G, 3) -> [e1, e2, e3]
~~~~
"""
function random_edge_indices(h::Hypergraphs, n::Int64)
   @assert 0 < n <= h.m

   indices = [i for i = 1:h.m]
   for i = 1:h.m - n
      deleteat!(indices, rand(1:size(indices,1)))
   end

   return indices
end

"""
`hyper_to_bipartite`
====================

Converts a hypergraph into its bipartite representation. For all vertices v and
edges e in the hypergraph, an edge (v, e) exists in the bipartite if and only if
v ∈ e in the hypergraph.
"""
function hyper_to_bipartite(h::Hypergraphs)
   di_edges::Vector{Vector{Int64}} = []

   for i = 1:h.m
      for j in h.edges[i]
         push!(di_edges, [i + h.n, vertex_floor(j)])
      end
   end

   return MatrixNetwork(di_edges, h.n + h.m)
end

"""
`bipartite_to_hyper`
====================

Converts a bipartite graph into a hypergraph.
"""
function bipartite_to_hyper(b::MatrixNetwork) end # STUB
