module HigherOrderNetworks
using MatrixNetworks
using SparseArrays
using StatsBase

include("Utilities.jl")
include("Hypergraphs.jl")
export VertexHypergraph, StubHypergraph, add_node!, remove_node!, remove_edge!, hyper_to_bipartite

include("RandomHypergraphs.jl")
export stub_matching, pairwise_reshuffle!, pairwise_reshuffle_v!, pairwise_reshuffle_s!, MCMC, MCMC_s, MCMC_v

include("kcliques.jl")
export k_cliques, triangles

include("manage_data.jl")
export load_data

end # module
