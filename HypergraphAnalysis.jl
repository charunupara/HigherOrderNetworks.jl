"""
Analyses for hypergraphs: degree assortativity and intersection profiles.

Based on "Configuration Models of Random Hypergraphs" by Philip S. Chodrow.
The original paper may be found at https://arxiv.org/abs/1902.09302.
"""

using Random # Used to randomly sample pairs of adjacent nodes

"""
`choose_pairs`
==============

Sample random pairs of nodes in the hypergraph that share an edge.

Arguments
---------
    - `H::Hypergraph`: The hypergraph from which to select the pairs
    - `n_pairs::Int64`: The number of pairs to sample
    - `choice_function::Function`: The method for choosing the nodes from each edge.
                                   Options are `uniform`, `top_2`, and `top_bottom`.

Returns
-------
    - An `n_pairs` x 2 matrix, where each row is a pair and each entry is the
      degree of one of the nodes.
"""
function choose_pairs(H::Hypergraph, n_pairs::Int64, choice_function::Function)
    candidates = filter(x -> H.K[x] >= 2, 1:H.m) # Filters out singleton edges
    if n_pairs > size(candidates,1) # Adjust n_pairs if too large
        n_pairs = size(candidates,1)
    end

    # Randomly sample edges, and replace only if uniformly selecting nodes
    edges = sample(candidates, n_pairs, replace=string(choice_function) == "uniform")
    pairs = [H.D[v] for e in edges for v in choice_function(H, e)] # Choose from edges and convert to degrees

    return reshape(pairs, 2, size(pairs,1) ÷ 2)'
end

"""
`assortativity`
===============

Calculate the degree assortativity of a hypergraph (pp. 13).
Degree assortativity measures how much nodes tend to associate with other
nodes of similar degree.

Arguments
---------
    - `H::Hypergraph`: The hypergraph to be analyzed
    - `choice_function::String (="uniform")`: The method used to select nodes from edges. Options are
                                                - `"uniform"`, which selects two nodes from each edge u.a.r.,
                                                - `"top_2"`, which selects the two highest-degree nodes from each edge, and
                                                - `"top_bottom"`, which selects the highest- and lowest-degree nodes from each edge.
    - `method::String (="spearman")`: The type of correlation coefficient. Options are
                                        - `"spearman"`, which measures monotonic non-decreasing correlation, and
                                        - `"pearson"`, which measures linear correlation.
    - `samples::Int64 (=H.m)`: The number of edges to sample

Returns
-------
    - The correlation coefficient `r` ∈ [-1,1]

Example
--------
~~~~
assortativity(H; choice_function="top_bottom", method="pearson", samples=100)
~~~~
"""

function assortativity(H::Hypergraph; choice_function::String="uniform", method::String="spearman", samples::Int64=H.m)
    function uniform(H::Hypergraph, e::Int64)
        return sample(H.edges[e], 2, replace=false)
    end

    function top_2(H::Hypergraph, e::Int64)
        return shuffle(H.edges[e][1:2])
    end

    function top_bottom(H::Hypergraph, e::Int64)
        return shuffle([H.edges[e][1], H.edges[e][H.K[e]]])
    end

    choice_functions = Dict("uniform" => uniform, "top_2" => top_2, "top_bottom" => top_bottom)

    pairs = choose_pairs(h, samples, choice_functions[choice_function])

    if method == "spearman"
        # Rank degrees within each column
        pairs = [sortperm(sortperm(pairs[:,1])) sortperm(sortperm(pairs[:,2]))]
    end

    return cor(pairs[:,1], pairs[:,2])
end

"""
`conditional_profile`
=====================

Computes the probability that two randomly selected edges of sizes `k` and `l`
have an intersection of size `j`. (pp. 14-15)

Arguments:
----------
    - `H::Hypergraph`: The hypergraph to be analyzed
    - `j::Int64`: The intersection size to look for
    - `k::Int64`: Edge size
    - `l::Int64`: Edge size

Example
-------
~~~~
conditional_profile(H, 4, 6, 10)
~~~~
"""
function conditional_profile(H::Hypergraph, j::Int64, k::Int64, l::Int64)
    if (j > k || j > l) # Can't have intersection of size j if edge size is less than j
        return 0
    end

    k_i = findall(i -> H.K[i] == k, 1:H.m) # Find indices of edges of size k
    if size(k_i) == 0 # No k, no way!
        return 0
    end

    if k == l
        l_i = k_i
    else
        l_i = findall(i -> H.K[i] == l, 1:H.m) # Find indices of edges of size l
    end

    if size(l_i) == 0
        return 0
    end

    return mean(
                [length(edge_intersect(h, x, y)) == j ? 1 : 0
                for x in k_i for y in l_i if x != y]
                )
end

"""
`conditional_average`
=====================

Compute the average intersection size of all edges of sizes `k` and `l`. (pp. 15, last paragraph)

Arguments:
----------
    - `H::Hypergraph`: The hypergraph to be analyzed
    - `k::Int64`: Edge size
    - `l::Int64`: Edge size
"""
function conditional_average(H::Hypergraph, k::Int64, l::Int64)
    k_i = findall(i -> H.K[i] == k, 1:H.m)
    if size(k_i) == 0
        return 0
    end

    if k == l
        l_i = k_i
    else
        l_i = findall(i -> H.K[i] == l, 1:H.m)
    end

    if size(l_i) == 0
        return 0
    end

    return mean(
                [length(edge_intersect(h, x, y))
                for x in k_i for y in l_i if x != y]
                )
end

"""
`marginal_profile`
==================

Computes the probability that two randomly selected edges have an intersection of size `j`. (pp. 15, Definition 8)

Arguments:
----------
    - `H::Hypergraph`: The hypergraph to be analyzed
    - `j::Int64`: The intersection size to look for
"""
function marginal_profile(H::Hypergraph, j::Int64)
    return mean(
                [length(edge_intersect(h, x, y)) == j ? 1 : 0
                for x in 1:H.m for y in x+1:H.m]
                )
end
