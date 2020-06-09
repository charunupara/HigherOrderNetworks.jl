"""
Represents a Hypergraph as an incidence matrix.

# Fields
- `ci::Vector{Int64}`: Column indices, correspond to nodes
- `ri::Vector{Int64}`: Row indices, correspond to hyperedges
- `vals::Vector{T}`: The nonzero values of the matrix
- `n::Int64`: The number of nodes
- `m::Int64`: The number of hyperedges
- `D::Vector{Int64}`: The degree sequence
- `K::Vector{Int64}`: The edge dimension sequence
"""

abstract type hg end

mutable struct VertexHypergraph <: hg
   vals::Vector{Float64}
   n::Int64
   m::Int64
   D::Vector{Int64}
   K::Vector{Int64}
   edges::Dict{Int64, Vector{Int64}}
end

"""
The nth stub for node `i` is stored as `i + 1/(n+1)`
"""
mutable struct StubHypergraph <: hg
   vals::Vector{Float64}
   n::Int64
   m::Int64
   D::Vector{Int64}
   K::Vector{Int64}
   edges::Dict{Int64, Vector{Float64}}
end

"""
Verifies that the number of nodes and edges match `n` and `m` and that
all nodes are less than `n`. If successful, returns the degree and edge dimension sequences.
For VertexHypergraphs
"""
function Hypergraph_kernel(edges::Dict{Int64, Vector{Int64}}, vals::Vector{Float64},
                           n::Int64, m::Int64)
   @assert size(collect(keys(edges)), 1) == m # m is the number of edges
   @assert size(vals, 1) == n # Each node has a val associated with it
   @assert all([0 < edges[e][i] < n + 1 for e = 1:m for i = 1:size(edges[e], 1)]) # No node exceeds n

   count(d, e) = sum([i == e ? 1 : 0 for k in keys(d) for i in d[k]])

   D::Vector{Int64} = [count(edges, v) for v = 1:n]
   K::Vector{Int64} = [size(edges[e], 1) for e in keys(edges)]

   return D, K
end

"""
Verifies that the number of nodes and edges match `n` and `m` and that
all nodes are less than `n`. If successful, returns the degree and edge dimension sequences.
For StubHypergraphs
"""
function Hypergraph_kernel(edges::Dict{Int64, Vector{Float64}}, vals::Vector{Float64},
                           n::Int64, m::Int64)
   @assert size(collect(keys(edges)), 1) == m # m is the number of edges
   @assert size(vals, 1) == n # Each node has a val associated with it
   @assert all([0 < edges[e][i] < n + 1 for e = 1:m for i = 1:size(edges[e], 1)]) # No node exceeds n

   count(d, e) = sum([i == e ? 1 : 0 for k in keys(d) for i in d[k]])

   D::Vector{Int64} = [count(edges, v) for v = 1:n]
   K::Vector{Int64} = [size(edges[e], 1) for e in keys(edges)]

   return D, K
end

function VertexHypergraph(edges::Dict{Int64, Vector{Int64}}, vals::Vector{Float64},
                          n::Int64, m::Int64)
   D, K = Hypergraph_kernel(edges, vals, n, m)
   return VertexHypergraph(vals, n, m, D, K, edges)
end

"""
Creates VertexHypergraph with all node values as ones
"""
VertexHypergraph(edges::Dict{Int64, Vector{Int64}}, n::Int64, m::Int64) = VertexHypergraph(edges, ones(m), n, m)
VertexHypergraph(s::StubHypergraph) = VertexHypergraph(Dict(i => Int64.(floor.(s.edges[i])) for i in keys(s.edges)), s.vals, s.n, s.m)

function StubHypergraph(edges::Dict{Int64, Vector{Float64}}, vals::Vector{Float64},
                        n::Int64, m::Int64)
   D, K = Hypergraph_kernel(edges, vals, n, m)
   return StubHypergraph(vals, n, m, D, K, edges)
end

"""
Creates StubHypergraph with all node values as ones
"""
StubHypergraph(edges::Dict{Int64, Vector{Float64}}, n::Int64, m::Int64) = StubHypergraph(edges, ones(m), n, m)

function StubMatching(D::Vector{Int64}, K::Vector{Int64})
   hypergraph = StubHypergraph(Dict{Int64, Vector{Float64}}(), 0, 0)
   stubs::Vector{Vector{Float64}} = [[i+(1/(n+1)) for n = 1:i] for i = 1:size(D, 1)]
   for k = 1:size(K,1)
      edge = zeros(K[k])
      for d = 1:K[k]
         node = rand(stubs)
         stub = rand(node)
         edge[d] = stub
         setdiff!(node, [stub])
         if isempty(node)
            setdiff!(stubs, [node])
         end
      end
      hypergraph.edges[k] = edge
   end
   hypergraph.n = size(D, 1)
   hypergraph.m = size(K, 1)
   hypergraph.D = D
   hypergraph.K = K
   hypergraph.vals = ones(hypergraph.n)
   return hypergraph
end

print(VertexHypergraph(StubMatching([2,3,2,2], [3,2,4])))



#print(VertexHypergraph(StubHypergraph(Dict(1 => [1.5,2.5,3.5,4.5], 2 => [1.33,2.33,4.33], 3 => [2.25,3.25]), ones(4), 4, 3)))
