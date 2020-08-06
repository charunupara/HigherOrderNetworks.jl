module HigherOrderNetworks
using MatrixNetworks
using SparseArrays
using StatsBase
using LinearAlgebra

include("Utilities.jl")
include("Hypergraphs.jl")
export VertexHypergraph, StubHypergraph, add_node!, remove_node!, remove_edge!, hyper_to_bipartite

include("Hyperkron.jl")
export num_multiset_permutations, _next_distinct_character, unrank5!, _morton_decode3, hyperkron_graph, kron_params

include("HypergraphConversions.jl")
export to_adjacency_matrix,hypergraph_to_matrixnetwork, hyper_to_bipartite, is_bipartite, bipartite_to_hyper, Hypergraphs

include("RandomCliqueCovers.jl")
export rcc

include("RandomHypergraphs.jl")
export stub_matching, pairwise_reshuffle!, pairwise_reshuffle_v!, pairwise_reshuffle_s!, MCMC, MCMC_s, MCMC_v

include("Cliques.jl")
export kcliques

include("ManageData.jl")
export load_data

end # module
