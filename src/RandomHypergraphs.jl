include("Hypergraphs.jl")
"""
The models in this file come from Chodrow (2019) "Configuration Models of Random Hypergraphs". https://arxiv.org/abs/1902.09302
Any occurrence of 'pp. X' indicates that a bit of code is inspired by material on page X of the paper.
"""

"""
`stub_matching`
===============

Produce a random (likely degenerate) hypergraph with given degree and edge
dimension sequences via stub matching. pp. 2, 6

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
   hypergraph = StubHypergraph(Vector{Vector{Float64}}(), 0, 0) # Blank hypergraph
   stubs::Vector{Vector{Float64}} = [[i+(1/(n+1)) for n = 1:D[i]] for i = 1:size(D,1)]
   for k = 1:size(K,1) # Iterate over edge dimension sequence
      edge = zeros(K[k]) # Empty edge of length K[k]
      for d = 1:K[k] # Choose K[k] random stubs for the edge
         node = rand(stubs) # Node to get the stub from
         stub = rand(node)
         edge[d] = stub
         setdiff!(node, [stub]) # Remove stub
         if isempty(node)
            setdiff!(stubs, [node])
         end
      end
      push!(hypergraph.edges, edge)
   end
   hypergraph.n = size(D,1) # Set hypergraph fields
   hypergraph.m = size(K,1)
   hypergraph.D = D
   hypergraph.K = K
   hypergraph.vals = vals
   return hypergraph
end

"""
`pairwise_reshuffle!`
=====================

Performs a pairwise reshuffle on the edges numbered `e1` and `e2`.

A pairwise shuffle is an operation that randomly swaps nodes between
two edges while leaving their intersection intact. Stubs within the
intersection may be exchanged. pp. 6, 7
"""
function pairwise_reshuffle!(h::Hypergraphs, e1::Int64, e2::Int64)
   @assert 0 < e1 <= h.m && 0 < e2 <= h.m # Edges to be reshuffled are valid

   if typeof(h) <: VertexHypergraph # Call proper function
      pairwise_reshuffle_v!(h, e1, e2)
   else
      pairwise_reshuffle_s!(h, e1, e2)
   end
end

"""
`pairwise_reshuffle_v!`
=======================

Perform a pairwise reshuffle on edges at indices `e1` and `e2` in a vertex-labeled hypergraph.
"""
function pairwise_reshuffle_v!(h::VertexHypergraph, e1::Int64, e2::Int64)
   inter = collect(edge_intersect(h, e1, e2)) # edge1 ⋂ edge2
   exclu = setdiff([h.edges[e1]; h.edges[e2]], inter) # (edge1 ⋃ edge2) \ inter

   h.edges[e1] = copy(inter) # Initialize e1 as intersection

   for i = 1:h.K[e1] - size(inter,1) # Randomly assign |e1 \ e2| nodes to e1
      new_n = rand(exclu)
      push!(h.edges[e1], new_n)
      setdiff!(exclu, [new_n])
   end

   h.edges[e2] = [inter; exclu] # Assign the rest to e2
end

"""
`pairwise_reshuffle_s!`
=======================

Perform a pairwise reshuffle on edges at indices `e1` and `e2` in a stub-labeled hypergraph.
"""
function pairwise_reshuffle_s!(h::StubHypergraph, e1::Int64, e2::Int64)
   inter = edge_intersect(h, e1, e2) # Get the set of intersection nodes
   exclu = filter(x -> !(Int(floor(x)) in inter), [h.edges[e1]; h.edges[e2]]) # Get stubs that are not in intersection
   filter!(x -> vertex_floor(x) in inter, h.edges[e1]) # Set e1 and e2 to only the stubs that are in the intersection
   filter!(x -> vertex_floor(x) in inter, h.edges[e2])
   sort!(h.edges[e1])
   sort!(h.edges[e2])

   for i = 1:length(inter) # Shuffle intersection stubs
      if rand(Float64) <= 0.5
         temp = h.edges[e1][i]
         h.edges[e1][i] = h.edges[e2][i]
         h.edges[e2][i] = temp
      end
   end

   for i = 1:h.K[e1] - length(inter)
      new_n = rand(exclu)
      push!(h.edges[e1], new_n)
      setdiff!(exclu, [new_n])
   end

   h.edges[e2] = [h.edges[e2]; exclu]
end

"""
`probability`
=============

Get the probability that a pairwise reshuffle between `e1` and `e2` occurs. pp. 7, 8
"""
function probability(h::Hypergraphs, e1::Int64, e2::Int64)
   intersect_size = length(edge_intersect(h, e1, e2))
   pstub_realization = 1/((2 ^ intersect_size) * C(h.K[e1] + h.K[e2] - 2*intersect_size, h.K[e1] - intersect_size)) # (2), pp. 7
   overall_stub = pstub_realization/C(h.m, 2) # (3), pp. 7
   if typeof(h) <: StubHypergraph
      return overall_stub
   else
      return (overall_stub * 2^intersect_size) / (num_parallel(h, e1)*num_parallel(h, e2)) # Statement of Theorem 2
   end
end

"""
`MCMC`
======

Randomly explore the space of hypergraphs (vertex- or stub-labeled) with the degree
and edge-dimension sequences of a given graph via random pairwise shuffling. pp. 7-9

Arguments
---------
   - `intiial::Hypergraphs`: The hypergraph from which exploration begins
   - `samples(=1000)`: The number of iterations to run 

Example
--------
~~~~
MCMC(G, 40, 50)
~~~~
"""
function MCMC(initial::Hypergraphs; samples=1000)
   @assert samples > 0 # Samples are positive

   current = copy(initial)
   for t = 1:samples
      if rand(Float64) <= 0.5#probability(current, sample_edges[1], sample_edges[2]) # Randomly decide whether to shuffle
         sample_edges = random_edge_indices(current, 2) # Random pair of edges
         pairwise_reshuffle!(current, sample_edges[1], sample_edges[2])
      end
   end
   return current
end

#=g = StubHypergraph([[1.5,2.5],[7/3,3.5],[4/3,4.5]], 4, 3)
print(g)
MCMC(g, 10, 10)
print(g)=#
