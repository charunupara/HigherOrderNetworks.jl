include("RandomHypergraphs.jl")
using Statistics

function triadic_closure(h::Hypergraphs) # Wait for Charun's code
    return 
end

function uniform(h::Hypergraphs, e::Int64)
    return sample(h.edges[e], 2, replace=false)
end

function top_2(h::Hypergraphs, e::Int64)
    return sort(h.edges[e], by=x -> h.D[x], rev=true)[1:2]
end

function top_bottom(h::Hypergraphs, e::Int64)
    sorted_by_deg = sort(h.edges[e], by=x -> h.D[x], rev=true)
    return [sorted_by_deg[1], sorted_by_deg[h.K[e]]]
end


function assortativity(h::Hypergraphs; choice_function = uniform)
    pairs = [choice_function(h, e) for e = 1:h.m if h.K[e] >= 2]
    
    sorted_by_deg = sorted_nodes(h)
    get_rank(n) = findfirst(isequal(n), sorted_by_deg)

    println(pairs[1:10])
    println([get_rank(x[1]) for x in pairs[1:10]])

    return cor([get_rank(x[1]) for x in pairs], [get_rank(x[2]) for x in pairs])
end

function conditional_profile(h::Hypergraphs, j::Int64, k::Int64, l::Int64) # Haven't tried
    k_i = findall(i -> h.K[i] == k, 1:h.m)
    if k == l
        l_i = k_i
    else
        l_i = findall(i -> h.K[i] == l, 1:h.m)
    end

    if size(k_i) == 0 || size(l_i) == 0
        return 0
    end

    return mean([length(edge_intersect(h, x, y)) == j ? 1 : 0 for x in k_i for y in l_i if x != y])
end

function marginal_profile(h::Hypergraphs, j::Int64) # Super slow
    return mean([length(edge_intersect(h, x, y)) == j ? 1 : 0 for x in 1:h.m-1 for y in x+1:h.m])
end

g = VertexHypergraph("new-enron.txt")
print(assortativity(g; choice_function=top_2)) #Outputs a number, but the values are off from the paper
#print(marginal_profile(g, 0))

"""
Here is what I've verified works correctly with the assortativity function:
    - uniform
    - pairs
    - sorted_nodes(h)
    - get_rank(n)
    - cor
    - get_rank mapping onto pairs

All of these parts work individually, so they should come together to output the correct results, but they're not.
So I must have some conceptual misunderstanding of the process. What do you guys think?
"""

