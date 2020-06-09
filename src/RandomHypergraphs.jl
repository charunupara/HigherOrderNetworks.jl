include("Hypergraphs.jl")

"""
`stub_matching`
===============

Produce a random (likely degenerate) hypergraph with given degree and edge
dimension sequences via stub matching.

Arguments
---------
   - `D::Vector{Int64}`: The desired degree sequence
   - `K::Vector{Int64}`: The desired edge dimension sequence
   - `vals::Any (=ones())`: The values located at each node

Examples
-------
~~~~
stub_matching([3,2,2,2], [4,3,2])
stub_matching([1], [2], ["foo"])
~~~~
"""
function stub_matching(D::Vector{Int64}, K::Vector{Int64}; vals=ones(size(D,1)))
   hypergraph = StubHypergraph(Vector{Vector{Float64}}(), 0, 0)
   stubs::Vector{Vector{Float64}} = [[i+(1/(n+1)) for n = 1:D[i]] for i = 1:size(D,1)]
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

"""
`pairwise_reshuffle!`
=====================

Performs a pairwise reshuffle on the edges numbered `e1` and `e2`.
"""
function pairwise_reshuffle!(h::Hypergraphs, e1::Int64, e2::Int64)
   @assert 0 < e1 <= h.m && 0 < e2 <= h.m # Edges to be reshuffled are valid

   if typeof(h) <: VertexHypergraph # Call proper function
      pairwise_reshuffle_v!(h, e1, e2)
   else
      pairwise_reshuffle_s!(h, e1, e2)
   end
end

function pairwise_reshuffle_v!(h::VertexHypergraph, e1::Int64, e2::Int64)
   inter = filter(x -> x in h.edges[e2], h.edges[e1]) # edge1 ⋂ edge2
   exclu = setdiff([h.edges[e1]; h.edges[e2]], inter) # (edge1 ⋃ edge2) \ inter

   h.edges[e1] = copy(inter) # Initialize e1 as intersection

   for i = 1:h.K[e1] - size(inter,1) # Randomly assign |e1 \ e2| random nodes to e1
      new_n = rand(exclu)
      push!(h.edges[e1], new_n)
      setdiff!(exclu, [new_n])
   end

   h.edges[e2] = [inter; exclu] # Assign the rest to e2
end

function pairwise_reshuffle_s!() end # STUB
