using MatrixNetworks
include("Utilities.jl")
"""
`motif_modularity`
==================

Calculate the motif modularity of a graph given a motif and community partition.

Arguments
---------
    - `A::MatrixNetwork`: The graph to calculate motif modularity for
    - `C::Vector{Vector{Int64}}`: The community partition of `A`. The ith subarray of `C` contains the nodes of the ith community.
    - `M::MatrixNetwork (=MatrixNetwork([1],[2]))`: The given motif. A directed edge by default, will measure the standard modularity.
"""
function motif_modularity(A::MatrixNetwork, C::Vector{Vector{Int64}}; M::MatrixNetwork=MatrixNetwork([1],[2]))
    node_communities = zeros(Int64, A.n) # Easy lookup for which community a node is part of
    for i = 1:size(C,1)
        for n in C[i]
            node_communities[n] = i
        end
    end

    k_delt(n1, n2) = return node_communities[n1] == node_communities[n2] # Kronecker delta function, returns true if nodes are in same community

    a_e = get_edges(A) # Network edges
    m_e = get_edges(M) # Motif edges
    out_w = outweight_distr(A)
    in_w = inweight_distr(A)
    weights = total_edge_weights(A) # Total weights between each pair of nodes

    shared_strengths = Dict() # Relates sets of communities to the motif strength among them.

    function inc_inds!(inds::Vector{Int64}) # Simulate one iteration of an arbitrary number of nested for loops
        an = A.n
        sz = size(inds,1)
        if inds == an * ones(Int64, sz)
            inds .-= inds
            return
        end

        for i = sz:-1:1
            if inds[i] == an
                inds[i] = 1
                continue
            else
                inds[i] += 1
                return
            end
        end
        if any([inds[x] == inds[y] for x = 1:length(inds) for y = x+1:length(inds)])
            inc_inds!(inds)
        end
    end

    indices = ones(Int64, M.n)
    z = zeros(Int64, M.n)
    w,W,u,U = 0,0,0,0 # w: strength of motifs fully contained within communities, W: total strength of motifs
                      # u: w for a random, strength-preserving graph, U: u for the same

    while indices != z # While all indices haven't been checked
        w_inc = 1
        W_inc = 1
        weight = 1

        p = 1
        nf = 1

        motif_membership = Set() # Communities that the motif part of

        for e in m_e
            if !(indices[e] in a_e) # Motif edge not in graph
                w_inc = 0
                W_inc = 0
            end
            if !k_delt(indices[e]...) # Edge crosses community lines
                w_inc = 0
                nf = 0
                push!(motif_membership, node_communities[indices[e[1]]])
                push!(motif_membership, node_communities[indices[e[2]]])
            end
            if w_inc == 1
                weight *= weights[[indices[e[1]], indices[e[2]]]]
            end
            p *= out_w[indices[e[1]]] * in_w[indices[e[2]]]
        end
        w += w_inc * weight
        W += W_inc * weight
        u += p * nf
        U += p

        if !isempty(motif_membership) && length(reduce(vcat, C[collect(motif_membership)])) >= M.n
            motif_membership = sort(collect(motif_membership))
            if !haskey(shared_strengths, motif_membership)
                shared_strengths[motif_membership] = [0,0]
            end
            shared_strengths[motif_membership][1] += W_inc * weight
            shared_strengths[motif_membership][2] += p
        end

        inc_inds!(indices)
    end
    println(shared_strengths)
    #println("$C: $w/$W,$u/$U")
    return (w/W) - (u/U), shared_strengths, W, U
end

"""
`louvain_motif`
======================

Determine an approximation for the best community division using a modified Louvain algorithm (https://en.wikipedia.org/wiki/Louvain_modularity) on a given motif.

Arguments
---------
    - `A::MatrixNetwork`: The graph in which to find communities
    - `M::MatrixNetwork (=MatrixNetwork([1],[2]))`: The motif that forms the basis of a community
"""

function louvain_motif(A::MatrixNetwork; M::MatrixNetwork=MatrixNetwork([1],[2]))
    C = [[i] for i = 1:A.n] # Community partition
    C_last = deepcopy(C) # Copy of above, used to check if any merges occurred in previous round

    mod, strengths, W, U = motif_modularity(A,C; M=M) # Initial values for modularity shortcut

    ### Merge the communities at indices i and j in C. Updates both C and strengths. ###
    function merge_communities!(i,j) # TODO: For merges need to take into account 'imaginary' motifs for n_ij
        if i == j return end # Can't merge a community with itself!

        new_groups = collect(deepcopy(keys(strengths)))
        cpy = deepcopy(new_groups) # Remember old keys so we can get rid of ones that no longer exist at the end

        new_strengths = Dict()
        covered = Set() # Pairwise relationships that have already been used in updates

        for k = 1:size(new_groups,1)
            g = new_groups[k]
            g[g.==j] .= i # Replace j with i

            g = sort(collect(Set(g))) # If i already included, remove duplicate and sort for consistency

            if length(reduce(vcat, C[g])) < M.n  || g in keys(new_strengths) continue end # Don't do anything for degenerate groups or repeats
            #if length(g) < 2 || g in keys(new_strengths) continue end # Don't do anything for degenerate groups or repeats

            if i in g # Since g contains the new merged community, we have to update its dict values
                new_strengths[g] = [0,0]

                for c in g
                    if c == i continue end # No change needed for [i,i]
                    ic = sort([i,c])
                    jc = sort([j,c])

                    if !(ic in covered) && haskey(strengths, ic) # ic is new and i and c have a relationship TODO AMBIGUOUS!!
                        new_strengths[g] .+= strengths[ic] # Update values
                        push!(covered, ic)
                    end

                    if !(jc in covered) && haskey(strengths, jc)
                        new_strengths[g] .+= strengths[jc]
                        push!(covered, jc)
                    end
                end
            else
                new_strengths[g] = strengths[cpy[k]] # g doesn't contain the new community, so we leave its values alone
            end
        end

        for g in setdiff(cpy, new_groups) # Erase entries that no longer exist
            delete!(new_strengths,g)
        end

        old_keys = collect(deepcopy(keys(new_strengths)))
        new_keys = []
        vals = []
        for k = 1:length(old_keys) # Relabel community relationships to account for the loss of j
            key = old_keys[k]
            push!(vals, new_strengths[key])
            for r = 1:length(key)
                if key[r] > j
                    key[r] -= 1
                end
            end
            push!(new_keys, key)
        end

        # Update partition
        C[i] = [C[i]; C[j]]
        setdiff!(C, [C[j]])

        strengths = Dict(new_keys[i] => vals[i] for i = 1:length(new_keys)) # Ta da!
    end

    ### Perform one iteration of the modified Louvain algorithm. For each community, finds the merge that most increases modularity, then performs it. ###
    function louvain_iter!()
        for i = 1:size(C,1)
            if i > size(C,1) break end
            pre_mod = mod # Initial modularity, used to compute the the modularity resulting from a hypothetical merge
            best_dQ = typemin(Float64)
            best = i # Keep track of the community with which merging would maximally increase modularity
            for j = 1:size(C,1)
                if i == j continue end # Merging with self doesn't change modularity, so we can skip
                if j > size(C,1) break end

                new_C = deepcopy(C) # Copy of C used to simulate merging i and j
                ordered = sort([i,j])

                # Perform the merge
                new_C[i] = [C[i]; C[j]]
                setdiff!(new_C, [C[j]])

                if !(ordered in keys(strengths)) # i and j have no relationship, so merging them won't change modularity # TODO: This is where I assume edges
                    dQ = 0.0
                else
                    dQ = (strengths[ordered][1]/W) - (strengths[ordered][2]/U) # Change in modularity from performing merge
                end
                #real_mod, s = motif_modularity(A, new_C; M=M)[1:2]

                post_mod = pre_mod + dQ
                println("$(C[i]), $(C[j]), $dQ")
                if post_mod > mod
                    mod = post_mod
                    best = j
                    #println(strengths)
                end
            end
            merge_communities!(i, best) # Perform best merge
        end
    end

    louvain_iter!() # Get the ball rolling, make C different from C_last
    while C_last != C # Stop merging when no merge increases modularity
        C_last = deepcopy(C)
        louvain_iter!()
    end
    println(mod)
    return map(sort, C)
end

A = MatrixNetwork([(1,2), (1,3), (2,1), (2,3), (2,4), (3,1), (3,2), (4,2), (4,5), (4,8), (5,4), (5,6), (5,7), (6,5), (6,7), (7,5), (7,6), (8,4), (8,9), (8,10), (9,8), (9,10), (10,9), (10,8)], 10)
#A = MatrixNetwork([1,1,2,2,3,3,3,4],[2,3,1,3,1,2,4,3])
println("START")
println("")

@time println(louvain_motif(A))

tri_circ = MatrixNetwork([1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,7,7,8,8,9,9,10,10],
                         [2,3,4,5,1,3,4,5,1,2,4,5,1,2,3,5,1,2,3,4,6,5,7,10,6,8,7,9,8,10,9,6])
#println(louvain_motif(tri_circ; M=MatrixNetwork([1,2,3], [2,3,1])))
karate_edges = [[2,1],[3,1],[3,2],[4,1],[4,2],[4,3],[5,1],[6,1],[7,1],[7,5],[7,6],[8,1],[8,2],[8,3],[8,4],[9,1],
[9,3],[10,3],[11,1],[11,5],[1,6],[12,1],[13,1],[13,4],[14,1],[14,2],[14,3],[14,4],[17,6],[17,7],[18,1],[18,2],
[20,1],[20,2],[22,1],[22,2],[26,24],[26,25],[28,3],[28,24],[28,25],[29,3],[30,24],[30,27],[31,2],[31,9],[32,1],
[32,25],[32,26],[32,29],[33,3],[33,9],[33,15],[33,16],[33,19],[33,21],[33,23],[33,24],[33,30],[33,31],[33,32],[34,9],
[34,10],[34,14],[34,15],[34,16],[34,19],[34,20],[34,21],[34,23],[34,24],[34,27],[34,28],[34,29],[34,30],[34,31],[34,32],[34,33]]
karate_edges = reduce(vcat, [[i, reverse(i)] for i in karate_edges])
karate_network = MatrixNetwork(map(x->x[1],karate_edges),map(x->x[2],karate_edges))
#@time println(louvain_motif(karate_network))

# TODO: Community detection works for basic case where M is a single edge, but fails with larger motifs, usually having made no merges. Where do I implicitly assume edges in my implementation?
# TODO: Our algorithm is slightly different from the Louvain in that Louvain checks FOR EVERY COMMUNITY, if it would be best to leave its current community. Might implement this
# TODO: Don't count degenerated motifs!

# Doing a dyadic merge, at least initially, will only decrease modularity, as we increase intra-community nullcase without increase intra-community motifs.
# Merge hyperedges?
