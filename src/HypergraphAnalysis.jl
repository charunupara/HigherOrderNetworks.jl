include("RandomHypergraphs.jl")
using Statistics
using Random
using PyCall
const phil = pyimport("hypergraph")

function triadic_closure(h::Hypergraphs) # Wait for Charun's code
    return 
end

function uniform(h::Hypergraphs, e::Int64)
    return sample(h.edges[e], 2, replace=false)
end

function top_2(h::Hypergraphs, e::Int64)
    return shuffle(sort(h.edges[e], by=x -> h.D[x], rev=true)[1:2])
end

function top_bottom(h::Hypergraphs, e::Int64)
    sorted_by_deg = sort(h.edges[e], by=x -> h.D[x], rev=true)
    return shuffle([sorted_by_deg[1], sorted_by_deg[h.K[e]]])
end


function assortativity(h::Hypergraphs; choice_function = uniform, method = "spearman", n_pairs = 10)
    pairs = [[h.D[v] for v in choice_function(h, e)] for e = 1:n_pairs if h.K[e] >= 2]
    node_degrees = collect(Set([n for p in pairs for n in p]))

    get_rank(n) = findfirst(isequal(n), sort(node_degrees, rev=true))

    if method == "spearman"
        pairs = [[get_rank(v) for v in p] for p in pairs]
    end
    #println([get_rank(x[1]) for x in pairs[1:10]])

    #print(pairs)

    #println([(get_rank(x[1]), x[1]) for x in pairs])
    #println([(get_rank(x[2]), x[2]) for x in pairs])

    return cor([x[1] for x in pairs], [x[2] for x in pairs])
end

function conditional_profile(h::Hypergraphs, j::Int64, k::Int64, l::Int64) # Haven't tried
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

    return mean([length(edge_intersect(h, x, y)) == j ? 1 : 0 for x in k_i for y in l_i if x != y])
end

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
    
    return mean([length(edge_intersect(h, x, y)) for x in k_i for y in l_i if x != y for j = 1:maximum(h.K)])
end

function marginal_profile(h::Hypergraphs, j::Int64) # Super slow
    return mean([length(edge_intersect(h, x, y)) == j ? 1 : 0 for x in 1:h.m-1 for y in x+1:h.m])
end

g = VertexHypergraph("new-enron.txt")
m = MCMC(g, samples=10000)
#p = phil.hypergraph([[i - 1 for i in v] for v in g.edges])

println(conditional_profile(m, 10, 10, 10))
#println(log10(conditional_average(g, 10, 10) / conditional_average(m, 10, 10)))

"""println(assortativity(g; choice_function=top_2, method="spearman")) #Outputs a number, but the values are off from the paper
println(assortativity(m; choice_function=top_2, method="spearman"))
println(p.assortativity(n_samples=g.m, method="spearman", choice_function="top_2"))
p.stub_edge_MH()
println(p.assortativity(n_samples=g.m, method="spearman", choice_function="top_2"))
#print(marginal_profile(g, 0))
"""

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

"""
The discrepancy arises from how we are creating the ranked arrays for the Spearman.
"""