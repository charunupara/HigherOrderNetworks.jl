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
mutable struct VertexHypergraph{T}
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
mutable struct StubHypergraph{T}
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
                           n::Int64, m::Int64) where Te <: Number where T
   @assert size(edges, 1) == m # m is the number of edges
   @assert size(vals, 1) == n # Each node has a val associated with it
   @assert all([0 < edges[i][j] < n + 1 for i = 1:m for j = 1:size(edges[i], 1)]) # No node exceeds n

   if Te == Float64
      possible_stubs = [i + 1/(j+1) for i = 1:n for j = 1:count(edges, i)]
      @assert all([i in possible_stubs for k = 1:m for i in edges[k]]) # Valid stub numberings
   end

   count(v, e) = sum([i == e ? 1 : 0 for k = 1:m for i in v[k]])

   D::Vector{Int64} = [count(edges, v) for v = 1:n]
   K::Vector{Int64} = [size(edges[e], 1) for e = 1:m]

   return D, K
end

"""
VertexHypergraph constructors.

Functions
---------
   - VertexHypergraph(edges, vals, n, m): Produces a vertex-labeled hypergraph with the given edge set, values, and size
   - VertexHypergraph(edges, n, m): Produces a vertex-labeled hypergraph with the given edge set and size, with all values as 1.0
   - VertexHypergraph(s): Converts a stub-labeled hypergraph into a vertex-labeled hypergraph

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

VertexHypergraph(edges::Vector{Vector{Int64}}, n::Int64, m::Int64) = VertexHypergraph(edges, ones(m), n, m)
VertexHypergraph(s::StubHypergraph) = VertexHypergraph([Int64.(floor.(s.edges[i])) for i = 1:s.m], s.vals, s.n, s.m)

"""
StubHypergraph constructors.

Functions
---------
   - StubHypergraph(edges, vals, n, m): Produces a stub-labeled hypergraph with the given edge set, values, and size
   - StubHypergraph(edges, n, m): Produces a stub-labeled hypergraph with the given edge set and size, with all values as 1.0

Examples
--------
~~~~
StubHypergraph([[1.33333333,2.5], [3.5,4.5], [1.5,4.33333333]], ["One", "Two", "Three", "Four"], 4, 3)
StubHypergraph([[1.5,2.33333333], [1.33333333,2.5]], 2, 2)
~~~~
"""
function StubHypergraph(edges::Vector{Vector{T}}, vals::Vector{T},
                        n::Int64, m::Int64) where T <: Number
   D, K = Hypergraph_kernel(edges, vals, n, m)
   return StubHypergraph(vals, n, m, D, K, edges)
end

StubHypergraph(edges::Vector{Vector{Float64}}, n::Int64, m::Int64) = StubHypergraph(edges, ones(m), n, m)
