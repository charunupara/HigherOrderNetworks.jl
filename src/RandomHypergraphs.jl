"""
Hypergraph configuration model for fixed degree and edge-dimension sequences using
a Markov Chain Monte Carlo scheme.

Based on "Configuration Models of Random Hypergraphs" by Philip S. Chodrow.
The original paper may be found at https://arxiv.org/abs/1902.09302.
"""


"""
`stub_matching`
===============

Produce a random (likely degenerate) hypergraph with given degree and edge
dimension sequences via stub matching. (pp. 2,6)

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
function stub_matching(D::Vector{Int64}, K::Vector{Int64}; vals=ones(length(D)))
   hypergraph = Hypergraph(Vector{Vector{Float64}}(), 0, 0) # Blank hypergraph
   degs = copy(D)
   valid_nodes = 1:length(D)

   for k = 1:size(K,1) # Iterate over edge dimension sequence
      edge = zeros(K[k]) # Empty edge of length K[k]
      for d = 1:K[k] # Choose K[k] random stubs for the edge
         node = rand(1:length(degs))
         edge[d] = node
         degs[node] -= 1 # Remove stub
         if degs[node] == 0
            setdiff!(valid_nodes, [node])
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

Arguments
---------
   - `H::Hypergraph`: The hypergraph containing the edges to shuffle
   - `e1::Int64`: The index of the first edge
   - `e2::Int64`: The index of the second edge

Preconditions
-------------
   - `e1, e2 ∈ (0,m]` where `m` is the number of edges in the hypergraph
   - `e1 ≠ e2`

A pairwise reshuffle is an operation that randomly swaps nodes between
two edges while leaving their intersection intact. Stubs within the
intersection may be exchanged. (pp. 6-7)

Example
-------
~~~~
pairwise_reshuffle(H, 1, 2) # Shuffles edges 1 and 2
~~~~
"""
function pairwise_reshuffle!(H::Hypergraph, e1::Int64, e2::Int64)
   @assert 0 < e1 <= H.m && 0 < e2 <= H.m && e1 != e2 # Edges to be reshuffled are valid

   inter = collect(edge_intersect(h, e1, e2)) # edge1 ⋂ edge2
   exclu = setdiff([h.edges[e1]; h.edges[e2]], inter) # (edge1 ⋃ edge2) \ inter

   h.edges[e1] = copy(inter) # Initialize e1 as intersection
   for i = 1:h.K[e1] - size(inter,1) # Randomly assign |e1 \ e2| nodes to e1
      new_n = rand(exclu)
      push!(h.edges[e1], new_n)
      setdiff!(exclu, [new_n])
   end

   h.edges[e2] = [inter; exclu] # Assign the rest to e2

   sort!(h.edges[e1], by=x -> H.D[x], rev=true)
   sort!(h.edges[e2], by=x -> H.D[x], rev=true)
end

"""
`MCMC`
======

Randomly explore the space of hypergraphs (vertex- or stub-labeled) with the degree
and edge-dimension sequences of a given hypergraph via random pairwise shuffling. (pp. 7-9)

Arguments
---------
   - `intiial::Hypergraph`: The hypergraph from which exploration begins
   - `samples::Int64 (=1000)`: The number of pairwise reshuffles to perform
   - `type::String (="vertex")`: Which distribution to draw from. Choices are
                                 "vertex" or "stub".

Preconditions
-------------
   - `samples > 0`
"""
function MCMC(initial::Hypergraph; samples::Int64=1000, type::String="vertex")
   @assert samples > 0 # Positive number of samples
   if type == "stub"
      return MCMC_s(initial; samples=samples)
   elseif type == "vertex"
      return MCMC_v(initial; samples=samples)
   else
      error("Invalid type argument. Choices are 'vertex' and 'stub'.")
   end
end

"""
`MCMC_s`
========

See `MCMC`.

We generate a random stub-labeled hypergraph by performing a swap at every
iteration of the loop.
"""
function MCMC_s(initial::Hypergraph; samples=1000)
   c = copy(initial)
   for t = 1:samples
      sample_edges = random_edge_indices(c, 2) # Random pair of edges
      pairwise_reshuffle!(c, sample_edges[1], sample_edges[2])
   end
   return c
end

"""
`MCMC_v`
========

See `MCMC`.

To generate a random vertex-labeled hypergraph, we can't swap at every iteration;
otherwise, we effectively perform the stub-labeled process, with the effect of
biasing the end result towards vertex-labeled hypergraphs that have more equivalent
realizations as stub-labeled hypergraphs. Thus, we must sample from the below
distribution.
"""
function MCMC_v(initial::Hypergraph; samples=1000)
   k = 0 # Number of iterations performed
   n_rejected = 0 # Number of iterations in which a swap was rejected

   c = copy(initial)
   parallels = Dict() # Keeps track of the number of hyperedges parallel (multiset equal) to each distinct hyperedge
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
         pair = sample(1:initial.m, 2, replace=false) # Prohibit swaps between an edge and itself
         inds[i+1], inds[i+2] = pair[1], pair[2]
      end
      probs = rand(Float64, Int(n_rand / 2)) # Generate many random floats at once

      while true
         if n_ >= n_rand / 2 # Refresh random values
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

         n_ += 2 # Go to next pair of edges
         p_ += 1 # Next probability
         inter = 2.0^-length(intersect(ei, ej)) # Larger intersection = less likely to swap

         if probs[p_] > inter / (parallels[ei] * parallels[ej]) # Randomly decide whether to shuffle
            n_rejected += 1
            k += 1
         else
            if haskey(parallels, ei)
               parallels[ei] -= 1
               if parallels[ei] <= 0
                  delete!(parallels, ei)
               end
            end

            if haskey(parallels, ej)
               parallels[ej] -= 1
               if parallels[ej] <= 0
                  delete!(parallels, ej)
               end
            end

            pairwise_reshuffle!(c, i, j) # Perform shuffle

            new_ei = c.edges[i]
            new_ej = c.edges[j]

            if !haskey(parallels, new_ei)
               parallels[new_ei] = 0
            end
            if !haskey(parallels, new_ej)
               parallels[new_ej] = 0
            end
            parallels[new_ei] += 1 # Update parallels
            parallels[new_ej] += 1

            k += 1
            break
         end
      end
   end

   return c
end
