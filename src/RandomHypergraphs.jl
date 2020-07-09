include("Hypergraphs.jl")

"""
The models in this file come from Chodrow (2019) "Configuration Models of Random Hypergraphs". https://arxiv.org/abs/1902.09302
Any occurrence of 'pp. X' indicates that a bit of code is inspired by material on page X of the paper.
"""

"""
`stub_matching`
===============

Produce a random (likely degenerate) hypergraph with given degree and edge
dimension sequences via stub matching. (pp. 2, 6)

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
      push!(hypergraph.edges, sort(edge, by=x -> D[x], rev=true))
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
intersection may be exchanged. (pp. 6-7)
"""
function pairwise_reshuffle!(h::Hypergraphs, e1::Int64, e2::Int64)
   @assert 0 < e1 <= h.m && 0 < e2 <= h.m # Edges to be reshuffled are valid

   if typeof(h) <: VertexHypergraph # Call proper function
      pairwise_reshuffle_v!(h, e1, e2)
   else
      pairwise_reshuffle_s!(h, e1, e2)
   end

   sort!(h.edges[e1], by=x -> h.D[x], rev=true)
   sort!(h.edges[e2], by=x -> h.D[x], rev=true)
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
   exclu = filter(x -> !(vertex_floor(x) in inter), [h.edges[e1]; h.edges[e2]]) # Get stubs that are not in intersection
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
`MCMC`
======

Randomly explore the space of hypergraphs (vertex- or stub-labeled) with the degree
and edge-dimension sequences of a given graph via random pairwise shuffling. (pp. 7-9)

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

   if typeof(initial) <: StubHypergraph
      return MCMC_s(initial; samples=samples)
   else
      return MCMC_v(initial; samples=samples)
   end
end

function MCMC_s(initial::Hypergraphs; samples=1000)
   c = copy(initial)
   for t = 1:samples
      sample_edges = random_edge_indices(c, 2) # Random pair of edges
      pairwise_reshuffle_s!(c, sample_edges[1], sample_edges[2])
   end
   return c
end

function MCMC_v(initial::Hypergraphs; samples=1000)
   k = 0 # Number of iterations performed
   n_rejected = 0 # Number of iterations in which a swap was rejected

   c = copy(initial)
   parallels = Dict()
   for i = 1:c.m
      if !haskey(parallels, c.edges[i])
         parallels[c.edges[i]] = 0
      end
      parallels[c.edges[i]] += 1
   end

   while k - n_rejected < samples # Number of performed swaps < desired
      n_rand = 20000

      n_ = 1
      p_ = 1
      inds = zeros(Int64, n_rand) # Random edge indices
      for i = 0:2:n_rand-1
         pair = sample(1:initial.m, 2, replace=false)
         inds[i+1], inds[i+2] = pair[1], pair[2]
      end
      probs = rand(Float64, Int(n_rand / 2))

      while true
         if n_ >= n_rand / 2
            n_ = 1
            p_ = 1
            inds = zeros(Int64, n_rand)
            for i = 0:2:n_rand-1
               pair = sample(1:initial.m, 2, replace=false)
               inds[i+1], inds[i+2] = pair[1], pair[2]
            end
            probs = rand(Float64, n_rand)
         end 

         i, j = inds[n_], inds[n_+1]
         ei, ej = c.edges[i], c.edges[j]
         
         n_ += 2
         p_ += 1
         inter = 2.0^-length(intersect(ei, ej))

         if probs[p_] > inter / (parallels[ei] * parallels[ej]) # Randomly decide whether to shuffle
            n_rejected += 1
            k += 1
         else
            parallels[ei] -= 1
            parallels[ej] -= 1
            if parallels[ei] <= 0
               delete!(parallels, ei)
            end
            if parallels[ej] <= 0
               delete!(parallels, ej)
            end

            pairwise_reshuffle_v!(c, i, j)

            new_ei = c.edges[i]
            new_ej = c.edges[j]

            if !haskey(parallels, new_ei)
               parallels[new_ei] = 0
            end
            if !haskey(parallels, new_ej)
               parallels[new_ej] = 0
            end
            parallels[new_ei] += 1
            parallels[new_ej] += 1

            k += 1
            break
         end
      end
   end

   return c
end
