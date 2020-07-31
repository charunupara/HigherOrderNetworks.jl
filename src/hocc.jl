"""
Higher-order clustering coefficient, which generalizes the traditional
clustering coefficients to closure of cliques larger than triangles.

Based on "Higher-order clustering in networks" by Hao Yin, Austin R. Benson, and
Jure Leskovec. The original paper may be found at https://arxiv.org/pdf/1704.03913.pdf.
"""

using MatrixNetworks
include("kcliques.jl") # Used to find cliques
include("Utilities.jl") # Used to get degree distribution

"""
`hocc`
======

Computes the global and average local l-th order clustering coefficients of a
network.

Arguments
---------
    - `A::MatrixNetwork`: The network to be analyzed
    - `l::Int64`: The order of the clustering coefficients

Preconditions
-------------
    - `l > 1`

The authors define an l-wedge to be an l-clique with an extra node of degree 1.
The global l-th order clustering coefficient measures the ratio of (l+1)-cliques
to l-wedges: the probability that a randomly chosen l-wedge is "closed" into an
(l+1)-clique. The local l-th order clustering coefficient is defined similarly
as the ratio of the number of (l+1)-cliques a node is part of to the number of
l-wedges it is part of.

The traditional clustering coefficients are a special case `hocc(A,2)` of the
above concept.
"""
function hocc(A::MatrixNetwork, l::Int64)
    @assert l > 1

    dd = deg_distr(A)
    edges = Set(get_edges(A))

    Wl = k_cliques(A, l)
    Kl = Set()

    # Find all (l+1)-cliques by finding all nodes that are completely connected
    # to an l-clique
    for i = 1:A.n
        for w in Wl
            if all([[i,v] in edges for v in w])
                push!(Kl, sort([w; [i]]))
            end
        end
    end

    clique_membership = zeros(A.n)
    wedge_membership = zeros(A.n)

    for K in Kl
        for v in K
            clique_membership[v] += l # Each vertex in an (l+1)-clique closes
                                      # l l-wedges
        end
    end

    for W in Wl
        for v in W
            wedge_membership[v] += dd[v]-l+1 # Each vertex is the center (incident
                                             # to the extra vertex) of dv-l+1 l-wedges
        end
    end

    cc = [clique_membership[i] / wedge_membership[i] for i = 1:A.n] # HOLCC for each node
    cc[isnan.(cc)] .= 0.0

    return mean(cc), sum(clique_membership) / sum(wedge_membership)
end
