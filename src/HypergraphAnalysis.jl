include("RandomHypergraphs.jl")
using Statistics
using Random
using PyCall
const phil = pyimport("hypergraph")

function choose_pairs(h::Hypergraphs, samples::Int64, choice_function::Function)
    candidates = filter(x -> size(h.edges[x],1) >= 2, 1:h.m)
    if samples > size(candidates,1)
        samples = size(candidates,1)
    end

    edges = sample(candidates, samples, replace=false)
    pairs = [h.D[v] for e in edges for v in choice_function(h, e)]

    return reshape(pairs, 2, size(pairs,1) ÷ 2)'
end

function assortativity(h::Hypergraphs; choice_function = "uniform", method = "spearman", samples = h.m)
    function uniform(h::Hypergraphs, e::Int64)
        return sample(h.edges[e], 2, replace=false)
    end
    
    function top_2(h::Hypergraphs, e::Int64)
        return shuffle(sort(h.edges[e], by=x -> h.D[x])[1:2])
    end
    
    function top_bottom(h::Hypergraphs, e::Int64)
        sorted_by_deg = sort(h.edges[e], by=x -> h.D[x])
        return shuffle([sorted_by_deg[1], sorted_by_deg[h.K[e]]])
    end

    choice_functions = Dict("uniform" => uniform, "top_2" => top_2, "top_bottom" => top_bottom)

    pairs = choose_pairs(h, samples, choice_functions[choice_function])

    if method == "spearman"
        pairs = [sortperm(sortperm(pairs[:,1])) sortperm(sortperm(pairs[:,2]))]
    end

    return cor(pairs[:,1], pairs[:,2])
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
    
    return mean([length(edge_intersect(h, x, y)) for x in k_i for y in l_i if x != y])
end

function marginal_profile(h::Hypergraphs, j::Int64) # Super slow
    return mean([length(edge_intersect(h, x, y)) == j ? 1 : 0 for x in 1:h.m-1 for y in x+1:h.m])
end

g = VertexHypergraph("new-enron.txt")
#m = MCMC(g, samples=1000)
p = phil.hypergraph([[i - 1 for i in v] for v in g.edges])
samples = g.m
func = "top_2"
meth = "spearman"
#println(conditional_profile(m, 10, 10, 10))
#println(log10(conditional_average(g, 10, 10) / conditional_average(m, 10, 10)))

println(assortativity(g; choice_function=func, method=meth, samples=samples)) #Outputs a number, but the values are off from the paper
#println(assortativity(m; choice_function=top_2, method="spearman", samples=samples))
println(p.assortativity(choice_function=func, method=meth, n_samples=samples))
#p.stub_edge_MH()
#println(p.assortativity(n_samples=samples, method="spearman", choice_function="top_2"))