"""
Miscellaneous functions.
"""

using SparseArrays
using LinearAlgebra
using MatrixNetworks

"""
`deg_distr`
===========

Gets the degree distribution of an undirected MatrixNetwork `A`.
"""
deg_distr(A::MatrixNetwork) = [A.rp[i+1]-A.rp[i] for i = 1:A.n]

"""
`nodes_by_deg`
==============

Returns the nodes of a MatrixNetwork sorted by degree.

Arguments
---------
    - `A::MatrixNetwork`: The network to get the sorted nodes of
    - `descending::Bool (=false)`: Whether the nodes are sorted by descending or
                                   or ascending degree
"""
function nodes_by_deg(A::MatrixNetwork; descending::Bool=false)
    distr = deg_distr(A)
    return sort(collect(1:A.n), by=x->distr[x], rev=descending)
end

"""
`outweight_distr`
================

Computes the weight going out of each node in a weighted MatrixNetwork.

Returns
-------
    - A vector of Reals, where `d[i]` is the weight going out of the ith node.
"""
function outweight_distr(A::MatrixNetwork) # For unweighted, equal to outdegree distribution
    deg_dist = zeros(Int64, A.n)
    if eltype(A.vals) <: Real
        for i = 1:size(A.rp,1)-1
            deg_dist[i] += sum(A.vals[A.rp[i]:A.rp[i+1]-1])
        end
    else
        for i = 1:size(A.rp,1)-1
            deg_dist[i] += (A.rp[i+1] - A.rp[i])
        end
    end
    return deg_dist
end

"""
`inweight_distr`
================

Computes the weight coming into each node in a weighted MatrixNetwork.

Returns
-------
    - A vector of Reals, where `d[i]` is the weight coming into the ith node.
"""
function inweight_distr(A::MatrixNetwork)
    deg_dist = zeros(Int64, A.n)
    if eltype(A.vals) <: Real
        for j = 1:size(A.ci,1)
            deg_dist[A.ci[j]] += A.vals[j]
        end
    else
        for j in A.ci
            deg_dist[j] += 1
        end
    end
    return deg_dist
end

"""
`total_edge_weights`
====================

Computes the total weight between all pairs of adjacent nodes in a MatrixNetwork.
Can be either directed or undirected. In the undirected case, weights[[i,j]] = weights[[j,i]].

Preconditions
-------------
    - The values of the MatrixNetwork are of type `Real`

Returns
-------
    - A dictionary of pairs to Reals, where weights[[i,j]] is the total directed weight
      from i to j.
"""
function total_edge_weights(A::MatrixNetwork)
    @assert eltype(A.vals) <: Real
    weights = Dict{Vector{Int64}, eltype(A.vals)}()
    i = 1
    for j = 1:length(A.ci)
        if j >= A.rp[i+1]
            i += 1
         end

        if !haskey(weights, [i,A.ci[j]])
            weights[[i,A.ci[j]]] = 0
        end
        weights[[i,A.ci[j]]] += A.vals[j]
    end
    return weights
end

"""
`get_edges`
===========

Returns the edges of a MatrixNetwork.
"""
function get_edges(A::MatrixNetwork)
    return [(i,j) for i = 1:length(A.rp)-1 for j in A.ci[A.rp[i]:A.rp[i+1]-1]]
end

"""
`simplify`
==========

Removes any self-loops and multi-edges from a MatrixNetwork.
"""
function simplify(A::MatrixNetwork)
    M = Array(sparse(A))
    M -= Diagonal(M)
    M = min.(M,1)
    return MatrixNetwork(sparse(M))
end

"""
`make_undirected`
=================

Symmetrizes a MatrixNetwork.
"""
function make_undirected(A::MatrixNetwork)
    arr = Array(sparse(A))
    arr = triu(arr) + triu(arr)'
    return MatrixNetwork(sparse(arr))
end

"""
`clustering`
============

Calculates the global and average local clustering coefficients for a MatrixNetwork `A`.

Output tuple is of the form `(lcc, gcc)`.
"""
function clustering(A::MatrixNetwork)
    dd = deg_distr(A)
    cc = clustercoeffs(A)
    cc[isnan.(cc)] .= 0.0

    wedge_count = [dd[i]*(dd[i]-1) for i = 1:length(dd)]
    gcc = sum(cc .* wedge_count) / sum(wedge_count)

    return mean(cc), gcc
end
