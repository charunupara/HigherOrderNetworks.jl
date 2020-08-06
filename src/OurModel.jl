using MatrixNetworks
using LinearAlgebra
using Statistics
using DataFrames
using CSV
using Dates
include("hypergraphs_conversions.jl")
include("clique_m.jl")
include("Utilities.jl")
include("RandomCliqueCovers.jl")
include("TGPA_generator_v2.jl")
include("hyperkron.jl")
#include("PlotGraph.jl")

"""
`to_hypergraph`
===============

Converts a MatrixNetwork into a hypergraph by converting all cliques into hyperedges.
"""
function to_hypergraph(A::MatrixNetwork)
    cliques = find_triangle(AdjacencyListGraph(A))
    return VertexHypergraph(vcat(cliques...))
end

"""
`print_stats`
=============

Prints statistics about a MatrixNetwork. Currently supported analyses are:
    - Clustering coefficient (Average local, global, higher-order)
    - Degree distribution
    - Eigenvalue distribution
"""
function print_stats(A::MatrixNetwork)
    dd = deg_distr(A)
    cc, gcc = clust(A)

    function print_dist(D::Vector{T}) where T <: Real
        d = sort(D)
        println("- Size: $(size(d,1))")
        println("- Lowest: $(d[1])")
        println("- Highest: $(d[size(d,1)])")
        println("- Mean: $(mean(d))")
        println("- Variance: $(var(d))")
        println("- Median: $(median(d))")
        println("- IQR: $(iqr(d))")
    end

    println("CLUSTERING COEFFICIENTS\n-----------------------")
    println("- Average local: $cc")
    println("- Global: $gcc")
    println("- Higher-order: TODO")

    println("\nDEGREE DISTRIBUTION\n-------------------")
    print_dist(dd)

    #=println("\nEIGENVALUE DISTRIBUTION\n-----------------------")
    print_dist(eigvals(Array(sparse(m))))=#
end

"""
`clust`
=======

Computes the average local and global clustering coefficients for a MatrixNetwork.

Example
-------
~~~~
lcc, gcc = clust(A)
~~~~
"""
function clust(A::MatrixNetwork)
    dd = deg_distr(A)
    cc = clustercoeffs(A)
    wedge_count = [dd[i]*(dd[i]-1) for i = 1:size(dd,1)]
    gcc = sum(cc .* wedge_count) / sum(wedge_count)

    return mean(cc), gcc
end

function edge_MCMC(initial::VertexHypergraph; e_bound=1, samples=1000)
    k = 0 # Number of iterations performed
    n_rejected = 0 # Number of iterations in which a swap was rejected

    g = hypergraph_to_matrixnetwork(initial)
    c = copy(initial) # Hypergraph that keeps track of location in hypergraph-space
    node_edges = [[] for i = 1:c.n]

    for e = 1:c.m
        for v in c.edges[e]
            push!(node_edges[v], e)
        end
    end

    quit = false
    while k - n_rejected < samples # Number of performed swaps < desired
        n_rand = 20000
        inds = rand(1:initial.m, n_rand)  # Random edge indices
        n_ = 1 # Location in index list
        tried = deepcopy(node_edges)

        if quit
            println("Unable to find edge-count-preserving swap.")
            break
        end

        while true # Propose shuffles until one is accepted
            if n_ >= n_rand # Generate new random values if n_ gets large
                inds = rand(1:initial.m, n_rand)
                n_ = 1
            end

            ei = c.edges[inds[n_]] # Current edge
            if size(ei,1) < e_bound # Check ei is big enough
                n_rejected += 1
                k += 1
                n_ += 1
                continue
            end
            vi, vt = rand(ei), rand(setdiff(1:c.n, ei)) # Random vertex in edge, and random vertex outside

            e_vi = setdiff(node_edges[vi], [inds[n_]])
            e_vt = node_edges[vt]

            if length(e_vi) == 0
                inter_i = Set()
            else
                inter_i = Set(intersect(ei, reduce(vcat, c.edges[e_vi]))) # Intersection size between ei and other edges vi belongs to
            end

            if length(e_vt) == 0
                inter_t = Set()
            else
                inter_t = setdiff(Set(intersect(ei, reduce(vcat, c.edges[e_vt]))), [vi])
            end


            if length(inter_i)-1 != length(inter_t) && !(length(inter_i) == 0 && length(inter_t) == 0)
                n_rejected += 1
                k += 1
                #println("Intersections: $(length(inter_i)), $(length(inter_t)), Edge: $ei, Vertices: $vi, $vt") # Swap stats
                setdiff!(tried[vi], [inds[n_]])
                if all([isempty(t) for t in tried]) quit = true end
            else # Success!
                #println("Intersections: $(length(inter_i)), $(length(inter_t)), Edge: $ei, Vertices: $vi, $vt") # Swap stats
                setdiff!(ei, [vi])
                setdiff!(node_edges[vi], [inds[n_]])
                push!(ei, vt)
                push!(node_edges[vt], inds[n_])
                k += 1
                #=h = hypergraph_to_matrixnetwork(c; simple=true)
                if size(h.ci,1) - size(g.ci,1) != 0
                    #println(vcat(c.edges[e_vi], [ei]), c.edges[e_vt])
                    #println(inter_i, inter_t)
                    #println(size(h.ci,1) - size(g.ci,1))
                    #g = deepcopy(h)
                end=# # Check change in edges in projection for debugging
                break
            end
            n_ += 1 # Hop to next index
        end
    end
    return c
end

"""
`our_model`
===========

Generates a new sample of a MatrixNetwork by randomly shuffling nodes between cliques

Arguments
---------
    - `A::MatrixNetwork`: The network to be "resampled"
    - `samples::Int64 (=1000)`: The number of clique shuffles to perform
    - `analyses::Bool (=true)`: Whether to print analyses comparing the initial graph to the new graph
    - `min_clique::Int64 (=3)`: The minimum size a clique has to be to be selected for swapping

Examples
--------
~~~~
our_model(A; samples=10000) # Performs 10000 swaps, prints analyses, and only swaps 3-cliques or larger
our_model(A; analyses=false, min_clique=2) # Performs 1000 swaps, doesn't print, and swaps any clique
~~~~
"""
function our_model(A::MatrixNetwork; samples::Int64=1000, analyses::Bool=true, min_clique::Int64=3, maintain_edges=false)
    if analyses
        println("\nInitial graph:\n==================")
        print_stats(A)
    end
    if maintain_edges
        r = hypergraph_to_matrixnetwork(edge_MCMC(to_hypergraph(A); samples=samples, e_bound=min_clique); simple=true)
    else
        r = hypergraph_to_matrixnetwork(MCMC_v(to_hypergraph(A); samples=samples, e_bound=min_clique); simple=true)
    end
    if analyses
        println("\n\nRandomized graph:\n================")
        print_stats(r)
    end
    return r
end

experiments = false
if experiments
    CC_d = []
    GCC_d = []
    CC1 = []
    CC2 = []
    GCC1 = []
    GCC2 = []
    change = []
    @time for i = 1:10
        #rcc = random_clique_cover(15., 0.9, 1.; N=100)
        #rcc = pa_graph(1000,5,2)
        rcc = erdos_renyi_undirected(10, 0.5)
        #rcc = MatrixNetwork(load_matrix_network("cores_example"))
        g = our_model(rcc; analyses=false, samples=1000, min_clique=2, maintain_edges=false)
        #println(size(g.ci,1) - size(rcc.ci,1))
        cc1, gcc1 = clust(rcc)
        cc2, gcc2 = clust(g)
        push!(CC1, cc1)
        push!(CC2, cc2)
        push!(GCC1, gcc1)
        push!(GCC2, gcc2)
        push!(CC_d, cc2 - cc1)
        push!(GCC_d, gcc2 - gcc1)


        rcc_t = length(find_triangle(AdjacencyListGraph(rcc))[2])
        g_t = length(find_triangle(AdjacencyListGraph(g))[2])
        #push!(change, g_t - rcc_t)
        println("$rcc_t -> $g_t")#: $(g_t - rcc_t)")#, $rcc_e -> $g_e: $(g_e - rcc_e)")
    end
    println("LCC Before: $(mean(CC1)), LCC After $(mean(CC2)), LCC Change: $(mean(CC_d))")
    println("GCC Before: $(mean(GCC1)), GCC After $(mean(GCC2)), GCC Change: $(mean(GCC_d))")
    println("Triangle ratio: $(mean(change))")
end

params = Dict("PA" => (1000,5,2),
              "TGPA" => (1000,0.5,5,2),
              "RCC" => (15.,0.9,1.,100),
              "Hyperkron" => (0.4,0.8,0.2,0.5,5))
Models = Dict("PA" => pa_graph,
              "TGPA" => tgpa_graph_v2,
              "RCC" => (a,s,c,N) -> random_clique_cover(a,s,c;N=N),
              "Hyperkron" => (a,b,c,d,r) -> MatrixNetwork(hyperkron_graph(kron_params(a,b,c,d), r)[1]))
Stats = Dict("lcc/gcc" => x -> round.(clust(x), digits=4),
             "tris" => x -> length(find_triangle(AdjacencyListGraph(x))[1]),
             "edges" => x -> size(x.ci,1))

function run_experiments(models::Vector{String}, stats::Vector{String}, n::Int64; samples=1000, min_clique=2, maintain_edges=false, filename="C:\\Users\\Elizabeth Turner\\Documents\\Josh\\MAP\\HO Networks\\Code\\Analyses$(Dates.format(now(), "u-dd-YY-HH-MM-SS")).csv")
    df = DataFrame([[String],[Int64],(Base.return_types(Stats[stats[i]]) for i = 1:length(stats))...],[["Model", "Samples"]; stats])

    function add_graph!(name::String, A::MatrixNetwork)
        stat = (Stats[s](A) for s in stats)
        push!(df, (name, n, stat...))
    end

    for model in models
        t1 = now()
        if !haskey(Models, model)
            println("Invalid model '$model'.")
            continue
        end
        for i = 1:n # TODO: Multiple samples
            before = Models[model](params[model]...)
            after = our_model(before; samples=samples, analyses=false, min_clique=min_clique, maintain_edges=maintain_edges)
            add_graph!("$model before", before)
            add_graph!("$model after", after)
        end
        println("$model finished in $((now() - t1)).")
    end
    println("Analyses complete. Data written to $filename.")
    delete!(df, [1])
    CSV.write(filename, df)
end
println("schweef")
run_experiments(collect(keys(Models)), collect(keys(Stats)), 1; maintain_edges=true)
