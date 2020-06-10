include("Hypergraphs.jl")
using Combinatorics

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

Performs a pairwise reshuffle on the edges numbered `e1` and `e2`. A pairwise
shuffle is an operation that randomly swaps nodes between two edges while
leaving their intersection untouched.
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

function pairwise_reshuffle_s!(h::StubHypergraph, e1::Int64, e2::Int64)
   inter = edge_intersect(h, e1, e2) # Get the set of intersection nodes
   exclu = filter(x -> !(Int(floor(x)) in inter), [h.edges[e1]; h.edges[e2]]) # Get stubs that are not in intersection
   filter!(x -> Int(floor(x)) in inter, h.edges[e1]) # Set e1 and e2 to only the stubs that are in the intersection
   filter!(x -> Int(floor(x)) in inter, h.edges[e2])
   sort!(h.edges[e1])
   sort!(h.edges[e2])

   for i = 1:length(inter) # Shuffle intersection stubs
      if rand(Float64) < 0.5
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

# n choose k
C(n,k) = factorial(n) / (factorial(k) * factorial(n-k))

function probability(h::Hypergraphs, e1::Int64, e2::Int64)
   intersect_size = length(edge_intersect(h, e1, e2))
   pstub_realization = 1/((2 ^ intersect_size) * C(h.K[e1] + h.K[e2] - 2*intersect_size, h.K[e1] - intersect_size))
   overall_stub = pstub_realization/C(h.m, 2)
   if typeof(h) <: StubHypergraph
      return overall_stub
   else
      return (overall_stub * 2^intersect_size) / (num_parallel(e1)*num_parallel(e2))
   end
end

function MCMC(initial::Hypergraphs, interval::Int64, sz::Int64) 
   current = initial
   for t = 1:interval*sz
      sample_edges = random_edge_indices(current, 2)
      if rand(Float64) <= probability(current, sample_edges[1], sample_edges[2])
         pairwise_reshuffle!(current, sample_edges[1], sample_edges[2])
      end
   end
   return current
end

#=g = StubHypergraph([[1.5,2.5],[7/3,3.5],[4/3,4.5]], 4, 3)
print(g)
MCMC(g, 10, 10)
print(g)=#