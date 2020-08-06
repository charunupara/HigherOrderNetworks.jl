"""
Higher-order clustering coefficient, which generalizes the traditional
clustering coefficients to closure of cliques larger than triangles.

Based on "Higher-order clustering in networks" by Hao Yin, Austin R. Benson, and
Jure Leskovec. The original paper may be found at https://arxiv.org/pdf/1704.03913.pdf.
Code from https://github.com/arbenson/HigherOrderClustering.jl.
"""

using LinearAlgebra
using SparseArrays
using StatsBase
using MatrixNetworks
include("cliques.jl")

"""
hoccf_data
----------

This is the data type returned by the function that computes higher-order clustering coefficients.

The field values are

order::Int64
    The order of the clustering coefficient (order = 2 is the classical
    clustering coefficient definition).

global_hoccf::Float64
    The global higher-order clustering coefficient.

avg_hoccf::Float64
    The average higher-order clustering coefficient (the mean is taken over
    nodes at the center of at least one wedge).

avg_hoccf2::Float64
     The average higher-order clustering coefficient, where the local clustering
    of a node not in any wedge is considered to be 0.

local_hoccfs::Vector{Float64}
    The vector of local higher-order clustering coefficients. If a node is not
    the center of at least one wedge, then the value is 0.

ho_wedge_counts::Vector{Int64}
    Higher-order wedge counts of nodes: ho_wedge_counts[v] is the number of higher-order
    wedges with node v at its center.

clique_counts::Vector{Int64}
    Clique counts of nodes: clique_counts[v] is the number of k-cliques containing
    node v where k = order + 1.
"""
struct hoccf_data
    order::Int64
    global_hoccf::Float64
    avg_hoccf::Float64
    avg_hoccf2::Float64
    local_hoccfs::Vector{Float64}
    ho_wedge_counts::Vector{Int64}
    clique_counts::Vector{Int64}
end

"""
`higher_order_ccfs`
======

Computes the global and average local l-th order clustering coefficients of a
network.

Arguments
---------
    - `A::MatrixNetwork`: The network to be analyzed
    - `l::Int64`: The order of the clustering coefficients

The authors define an l-wedge to be an l-clique with an extra node of degree 1 attached.
For example, this is a 3-wedge: https://i.imgur.com/Wc1pRzl.png
The clique node that is adjacent to the extra node is the "center" of the wedge.
The global l-th order clustering coefficient measures the ratio of (l+1)-cliques
to l-wedges: the probability that a randomly chosen l-wedge is "closed" into an
(l+1)-clique. The local l-th order clustering coefficient is defined similarly
as the ratio of the number of (l+1)-cliques a node is part of to the number of
l-wedges it is the center of.

The traditional clustering coefficients are a special case `hocc(A,2)` of the
above concept.
"""
function higher_order_ccfs(A::SparseMatrixCSC{T,Int64}, l::Int64) where T
    A = min.(A, 1)
    A -= Diagonal(A)
    n = size(A, 1)
    # Get clique counts
    clique_counts1 = kcliques(A, l)[1]
    clique_counts2 = kcliques(A, l + 1)[1]
    degs = vec(sum(A, dims=2))
    # normalize
    wedge_counts = (degs .- l .+ 1) .* clique_counts1
    nz_inds = findall(wedge_counts .> 0)
    local_hoccfs = zeros(Float64, n)
    local_hoccfs[nz_inds] = l .* clique_counts2[nz_inds] ./ wedge_counts[nz_inds]
    return hoccf_data(l, l .* sum(clique_counts2) / sum(wedge_counts),
                      mean(local_hoccfs[nz_inds]), mean(local_hoccfs),
                      local_hoccfs, wedge_counts, clique_counts2)
end
