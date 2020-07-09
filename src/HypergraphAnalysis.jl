include("RandomHypergraphs.jl")
using Statistics
using Random
using PyCall
const phil = pyimport("hypergraph")

"""
The analyses in this file come from Chodrow (2019) "Configuration Models of Random Hypergraphs". https://arxiv.org/abs/1902.09302
Any occurrence of 'pp. X' indicates that a bit of code is inspired by material on page X of the paper.
"""

"""
`choose_pairs`
==============

Sample `n_pairs` pairs of nodes in the hypergraph that share an edge.
Returns an `n_pairs` x 2 matrix, where each row is a pair and each entry is the degree of the selected node.

Arguments
---------
    - `h::Hypergraphs`: The hypergraph from which to select the pairs
    - `n_pairs::Int64`: The number of pairs to sample
    - `choice_function::Function`: The method for choosing the nodes from each edge. Options are `uniform`, `top_2`, and `top_bottom`.
"""
function choose_pairs(h::Hypergraphs, n_pairs::Int64, choice_function::Function)
    candidates = filter(x -> size(h.edges[x],1) >= 2, 1:h.m) # Filters out singleton edges
    if n_pairs > size(candidates,1) # Adjust n_pairs if too large
        n_pairs = size(candidates,1)
    end
    
    edges = sample(candidates, n_pairs, replace=string(choice_function) == "uniform") # Randomly sample edges, and replace only if uniformly selecting nodes
    pairs = [h.D[v] for e in edges for v in choice_function(h, e)] # Choose from edges and convert to degrees

    return reshape(pairs, 2, size(pairs,1) ÷ 2)'
end

"""
`assortativity`
===============

Calculate the degree assortativity of a given hypergraph (pp. 13). Degree assortativity measures how much nodes tend to associate with other
nodes of similar degree.

Arguments
---------
    - `h::Hypergraphs`: The hypergraph to be analyzed
    - `choice_function::String (="uniform")`: The method used to select nodes from edges. Options are `"uniform"`, which uniformly selects two nodes from each edge,
      `"top_2"`, which selects the two highest-degree nodes from each edge, and `"top_bottom"`, which selects the highest and lowest degree nodes from each edge.
    - `method::String (="spearman")`: The type of correlation coefficient. Options are `"spearman"`, which measures monotonic non-decreasing correlation, and
      `"pearson"`, which measures linear correlation.
    - `samples::Int64 (=h.m)`: The number of edges to sample

Example
--------
~~~~
assortativity(G; choice_function="top_bottom", method="pearson", samples = 100)
~~~~
"""

function assortativity(h::Hypergraphs; choice_function::String = "uniform", method::String = "spearman", samples::Int64 = h.m)
    function uniform(h::Hypergraphs, e::Int64)
        return sample(h.edges[e], 2, replace=false)
    end
    
    function top_2(h::Hypergraphs, e::Int64)
        return shuffle(h.edges[e][1:2])
    end
    
    function top_bottom(h::Hypergraphs, e::Int64)
        return shuffle([h.edges[e][1], h.edges[e][h.K[e]]])
    end

    choice_functions = Dict("uniform" => uniform, "top_2" => top_2, "top_bottom" => top_bottom)

    pairs = choose_pairs(h, samples, choice_functions[choice_function])

    if method == "spearman"
        pairs = [sortperm(sortperm(pairs[:,1])) sortperm(sortperm(pairs[:,2]))] # Rank degrees within each column. sortperm(sortperm(x)) ranks the elements of x
    end

    return cor(pairs[:,1], pairs[:,2])
end

"""
`conditional_profile`
=====================

Computes the probability that two randomly selected edges of sizes `k` and `l` have an intersection of size `j`. (pp. 14-15)

Example
-------
~~~~
conditional_profile(G, 4, 6, 10)
~~~~
"""
function conditional_profile(h::Hypergraphs, j::Int64, k::Int64, l::Int64) 
    if (j > k || j > l) # Can't have intersection of size j if edge size is less than j
        return 0
    end

    k_i = findall(i -> h.K[i] == k, 1:h.m) # Find indices of edges of size k
    if size(k_i) == 0
        return 0
    end

    if k == l
        l_i = k_i
    else
        l_i = findall(i -> h.K[i] == l, 1:h.m) # Find indices of edges of size l
    end

    if size(l_i) == 0
        return 0
    end

    return mean([length(edge_intersect(h, x, y)) == j ? 1 : 0 for x in k_i for y in l_i if x != y])
end

"""
`conditional_average`
=====================

Compute the average intersection size of all edges of sizes `k` and `l`. (pp. 15, last paragraph)
"""
function conditional_average(h::Hypergraphs, k::Int64, l::Int64)
    k_i = findall(i -> h.K[i] == k, 1:h.m)
    if size(k_i) == 0
        return 0
    end

    if k == l
        l_i = k_i
    else
        l_i = findall(i -> h.K[i] == l, 1:h.m)
    end

    if size(l_i) == 0
        return 0
    end
    
    return mean([length(edge_intersect(h, x, y)) for x in k_i for y in l_i if x != y])
end

"""
`marginal_profile`
==================

Computes the probability that two randomly selected edges have an intersection of size `j`. (pp. 15, Definition 8)
"""
function marginal_profile(h::Hypergraphs, j::Int64)
    return mean([length(edge_intersect(h, x, y)) == j ? 1 : 0 for x in 1:h.m-1 for y in x+1:h.m])
end


swaps = 5000

g = VertexHypergraph("new-enron.txt")
m = MCMC(g, samples=swaps)
p = phil.hypergraph([[i - 1 for i in v] for v in g.edges])

samples = g.m
func = "uniform"
meth = "spearman"


#println(assortativity(g; choice_function=func, method=meth, samples=samples)) #Outputs a number, but the values are off from the paper
#println(p.assortativity(choice_function=func, method=meth, n_samples=samples))
println(assortativity(m; choice_function=func, method=meth, samples=samples))
p.MH(n_steps=swaps)
println(p.assortativity(choice_function=func, method=meth, n_samples=samples))