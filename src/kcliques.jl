"""
Algorithm for finding k-cliques (complete subgraphs of order k).

Based on "Arbority and Subgraph Listing Algorithms" by Norshige Chiba and Takao
Nishizeki. The original paper may be found at
https://pdfs.semanticscholar.org/0d19/245a27bc65a87a8014d5b8a66fb514c8ff0b.pdf?_ga=2.251120710.1988770043.1596204199-972144793.1596204199.
"""

using MatrixNetworks

"""
`AdjacencyListGraph`
====================

Represents a dyadic graph as an adjacency list.

Fields
------
    - `vals::Vector{Int64}`: The values corresponding to each node
    - `n::Int64`: The number of nodes
    - `D::Vector{Int64}`: The degree sequence
    - `list::Vector{Vector{Int64}}`: The adjacency list, where the ith vector
       contains all nodes adjacent to node i
"""
mutable struct AdjacencyListGraph
    vals::Vector{Int64}
    n::Int64
    D::Vector{Int64}
    list::Vector{Vector{Int64}}
end

"""
`AdjacencyListGraph`
====================

Constructor that creates an `AdjacencyListGraph` from a `MatrixNetwork`.

The `vals` field from the `MatrixNetwork` is not transferred to the new structure.
Each new `AdjacencyListGraph` has an array of zeros as its values.
"""
function AdjacencyListGraph(A::MatrixNetwork)
    list = [collect(A.ci[A.rp[i]:A.rp[i+1]-1]) for i = 1:A.n]
    return AdjacencyListGraph(
        zeros(Int64, A.n),
        A.n,
        length.(list), # Degree sequence[i] = length(list[i])
        list,
    )
end

"""
`delete_node!`
==============

Functionally remove a node from an `AdjacencyListGraph` without actually
discarding it.

Arguments
---------
    - `G::AdjacencyListGraph`: The graph from which to remove a node
    - `node::Int64`: The node to be removed

Preconditions
-------------
    - `node <= G.n`

Postconditions
--------------
    - No adjacency list contains `node`
    - `G.D[node] = -1`
    - `G.vals[node] = G.n + 1`
    - The degrees of all neighbors of `node` are decreased by 1
    - `G.list[node] = []`
    - `G'.n = G.n` (`n` stays the same)
"""
function delete_node!(G::AdjacencyListGraph, node::Int64)
    @assert node <= G.n

    for n_list in G.list[G.list[node]] # Undirected, so we only need to
        setdiff!(n_list, [node])       # remove from neighbor lists
    end

    G.D[node] = -1 # Mark as deleted
    G.vals[node] = G.n + 1
    for v in G.list[node]
        G.D[v] -= 1 # Account for degree loss in neighbors
    end
    G.list[node] = []
end

"""
`k_cliques`
===========

The Chiba-Nishizeki algorithm for finding all k-cliques in a graph. Found on
pp. 216 of the paper.

Arguments
---------
    - `A::MatrixNetwork`: The graph in which to search for k-cliques
    - `k::Int64`: The desired clique size

The algorithm takes advantage of the following observations:
    1) a node is in a k-clique iff the subgraph induced by its neighbors
       contains a (k-1)-clique.
    2) if v is in the subgraph "A" induced by the neighbors of u, and w is in
       the subgraph induced by the neighbors of v in "A", then [u,v,w] is a clique.
So, the algorithm recursively finds k-cliques by exploring subgraphs of subgraphs
of subgraphs ... of neighbors, taking note of which vertices remain after k
iterations.
"""
function k_cliques(A::MatrixNetwork, k::Int64)
    C = Vector{Int64}() # Global stack containing the vertices from 2) above
    cliques = Set()

    function K(U::AdjacencyListGraph, l::Int64)
        if l == 2
            for n = 1:length(U.list)
                for u in U.list[n]
                    c = [[n, u]; C] # All [n,u] pairs are adjacent to all
                                    # vertices in C, so they form a clique per
                                    # 2) above
                    if length(c) == length(Set(c)) # No duplicate vertices
                        push!(cliques, sort(c))
                    end
                end
            end
        else
            # Numbers indicate steps of the pseudocode on pp. 216

            if sum(U.D) == -U.n # No valid node, stop searching here
                return # This original step was not in the pseudocode
            end

            # 2: Sort nodes by degree in nonincreasing order
            sorted_nodes = sort(1:U.n, by=x -> length(U.list[x]), rev=true)

            for v in sorted_nodes
                # 3: Generate the subgraph of U induced by neighbors of v that have
                #    the label "l". Relabel all vertices in the subgraph "l-1"
                U_prime = deepcopy(U)
                for n in setdiff(1:U_prime.n, U_prime.list[v]) # Delete non-neighbors
                    delete_node!(U_prime, n)
                end

                for n in U_prime.list[v]
                    if U_prime.vals[n] != l # Delete neighbors with incorrect labels
                        delete_node!(U_prime, n)
                    else
                        U_prime.vals[n] = l - 1
                    end
                end
                # 4: Sort adjacency lists to put lower-labeled vertices first
                sort!.(U_prime.list, by=x->U_prime.vals[x])
                # 6 (Skip 5, unnecessary): v is a member of some clique, push it
                #                          to the stack
                push!(C, v)
                # 7: Recurse, searching for (l-1)-clique in neighbor subgraph
                K(U_prime, l-1)
                # 8: Remove v from top
                pop!(C)
                # 9: Restore labels
                U_prime.vals = l * ones(U_prime.n)
                # 10: Logically delete v from U_prime
                U_prime.vals[v] = l + 1
                # 11: See step 4
                sort!.(U_prime.list, by=x->U_prime.vals[x])
            end
        end
    end

    G = AdjacencyListGraph(A)
    G.vals = k * ones(G.n)

    K(G, k)
    return cliques
end
