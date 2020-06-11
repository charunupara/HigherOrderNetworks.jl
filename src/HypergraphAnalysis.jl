include("Hypergraphs.jl")
using Statistics

function triadic_closure(h::Hypergraphs) # Wait for Charun's code
    return 
end

function uniform(h::Hypergraphs, e::Int64)
    return sample(edge, 2, replace=false)
end

function top_2(h::Hypergraphs, e::Int64)
    return sort(map(x -> h.D[vertex_floor(x)], h.edges[e]))[1:2]
end

function top_bottom(h::Hypergraphs, e::Int64)
    sorted_by_deg = sort(map(x -> h.D[vertex_floor(x)], h.edges[e]))
    return [sorted_by_deg[1], sorted_by_deg[h.K[e]]]
end

function degree_assortativity(h::Hypergraphs; choice_function = uniform)
    pairs = [p for p in [choice_function(h, e) for e = 1:h.m]]
    return cor(map(first, pairs), map(x -> x[2], pairs))
end

function conditional_profile(h::Hypergraphs, j::Int64, k::Int64, l::Int64)
    k_i = findall(i -> h.K[i] == k, 1:h.m)
    if k == l
        l_i = k_i
    else
        l_i = findall(i -> h.K[i] == l, 1:h.m)
    end

    if size(k_i) == 0 || size(l_i) == 0
        return 0
    end

    return sum([length(edge_intersect(h, x, y)) == j ? 1 : 0 for x in k_i for y in l_i if x != y]) / (size(k_i,1) * size(l_i,1))
end

function marginal_profile(h::Hypergraphs, j::Int64)
    return sum([length(edge_intersect(h, x, y)) == j ? 1 : 0 for x in 1:h.m-1 for y in x+1:h.m]) / (factorial(h.m-1))
end

g = VertexHypergraph([[1,2,3]], 3, 1)
println(g)
println(triadic_closure(g))