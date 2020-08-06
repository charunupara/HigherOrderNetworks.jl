module HigherOrderNetworks
using MatrixNetworks
using SparseArrays
using StatsBase
using LinearAlgebra

include("Utilities.jl")
include("Hypergraphs.jl")
export Hypergraphs, add_node!, remove_node!, add_edge!, remove_edge!

include("Hyperkron.jl")
export hyperkron_graph, kron_params

include("HypergraphConversions.jl")
export to_adjacency_matrix, hypergraph_to_matrixnetwork, hyper_to_bipartite, bipartite_to_hyper, is_bipartite

include("RandomCliqueCovers.jl")
export rcc

include("RandomHypergraphs.jl")
export stub_matching, pairwise_reshuffle!, MCMC, MCMC_s, MCMC_v

include("Cliques.jl")
export kcliques

include("ManageData.jl")
export load_data

end # module
