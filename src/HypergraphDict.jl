"""
Represents a Hypergraph as an incidence matrix.

# Fields
- `vals::Vector{Float64}`: The nonzero values of the matrix
- `n::Int64`: The number of nodes
- `m::Int64`: The number of hyperedges
- `D::Vector{Int64}`: The degree sequence
- `K::Vector{Int64}`: The edge dimension sequence
- `edges::Vector{Vector{Int64 or Float64}}`: The hyperedges and their members
"""

# abstract type hg{T} <: Number end

mutable struct VertexHypergraphs{T}
   vals::Vector{T}
   n::T
   m::T
   D::Vector{T}
   K::Vector{T}
   edges::Vector{Vector{T}}
end

"""
The nth stub for node `i` is stored as `i + 1/(n+1)`
"""
mutable struct StubHypergraphs{T}
   vals::Vector{T}
   m::T
   n::T
   D::Vector{T}
   K::Vector{T}
   edges::Vector{Vector{T}}
end

"""
Verifies that the number of nodes and edges match `n` and `m` and that
all nodes are less than `n`. If successful, returns the degree and edge dimension sequences.
For VertexHypergraphs
"""
function Hypergraph_kernel(edges::Vector{Vector{T}}, vals::Vector{T},
                           n::T, m::T) where T <: Number
   @assert size(edges, 1) == m # m is the number of edges
   @assert size(vals, 1) == n # Each node has a val associated with it
   @assert all([0 < edges[i][j] < n + 1 for i = 1:m for j = 1:size(edges[i], 1)]) # No node exceeds n

   count(v, e) = sum([i == e ? 1 : 0 for k = 1:m for i in v[k]])

   D::Vector{T} = [count(edges, v) for v = 1:n]
   K::Vector{T} = [size(edges[e], 1) for e = 1:m]

   return D, K
end

"""
Verifies that the number of nodes and edges match `n` and `m` and that
all nodes are less than `n`. If successful, returns the degree and edge dimension sequences.
For StubHypergraphs
"""
function Hypergraph_kernel(edges::Vector{Vector{T}}, vals::Vector{T},
                           n::T, m::T) where T <: Number
   @assert size(edges, 1) == m # m is the number of edges
   @assert size(vals, 1) == n # Each node has a val associated with it
   @assert all([0 < edges[i][j] < n + 1 for i = 1:m for j = 1:size(edges[i], 1)]) # No node exceeds n

   count(v, e) = sum([i == e ? 1 : 0 for k = 1:m for i in v[k]])

   D::Vector{T} = [count(edges, v) for v = 1:n]
   K::Vector{T} = [size(edges[e], 1) for e = 1:m]

   return D, K
end

function VertexHypergraph(edges::Vector{Vector{T}}, vals::Vector{T},
                          n::T, m::T) where T <: Number
   D, K = Hypergraph_kernel(edges, vals, n, m)
   return VertexHypergraph(vals, n, m, D, K, edges)
end

"""
Creates VertexHypergraph with all node values as ones
"""
VertexHypergraph(edges::Vector{Vector{Int64}}, n::Int64, m::Int64) = VertexHypergraph(edges, ones(m), n, m)
VertexHypergraph(s::StubHypergraph) = VertexHypergraph([Int64.(floor.(s.edges[i])) for i = 1:s.m], s.vals, s.n, s.m)

function StubHypergraph(edges::Vector{Vector{T}}, vals::Vector{T},
                        n::Int64, m::Int64) where T <: Number
   D, K = Hypergraph_kernel(edges, vals, n, m)
   return StubHypergraph(vals, n, m, D, K, edges)
end

"""
Creates StubHypergraph with all node values as ones
"""
StubHypergraph(edges::Vector{Vector{Float64}}, n::Int64, m::Int64) = StubHypergraph(edges, ones(m), n, m)

function StubMatching(D::Vector{T}, K::Vector{T}; vals=ones(size(D,1))) where T <: Number
   hypergraph = StubHypergraph(Vector{Vector{T}}(), 0, 0)
   stubs::Vector{Vector{T}} = [[i+(1/(n+1)) for n = 1:D[i]] for i = 1:size(D,1)]
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
      push!(hypergraph.edges, edge)
   end
   hypergraph.n = size(D, 1)
   hypergraph.m = size(K, 1)
   hypergraph.D = D
   hypergraph.K = K
   hypergraph.vals = vals
   return hypergraph
end

print(VertexHypergraph(StubMatching([3,2,2,2], [4,3,2])))



#print(VertexHypergraph(StubHypergraph(Dict(1 => [1.5,2.5,3.5,4.5], 2 => [1.33,2.33,4.33], 3 => [2.25,3.25]), ones(4), 4, 3)))
