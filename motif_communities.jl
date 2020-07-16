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

        if !isempty(motif_membership)
            motif_membership = sort(collect(motif_membership))
            if !haskey(shared_strengths, motif_membership)
                shared_strengths[motif_membership] = [0,0]
            end
            shared_strengths[motif_membership][1] += W_inc * weight
            shared_strengths[motif_membership][2] += p
        end

        inc_inds!(indices)
    end
    return (w/W) - (u/U), shared_strengths, W, U
end
# Somehow, we need to determine that, if communities A and B were combined, the motif would be fully contained within the union
# We could do something like checking that all the k_delt failures are from the same combination!
# But motifs that are shared between 3 communities will never be counted under this scheme...
# So it might be best to not completely memoize, to count the number of bi-shared motifs each time. This is still much better than total recomputation
# What we can memoize is all the edges from A to B
# Hmm, well we could represent motif sharing as weighted hyperedges! So we don't only track binary relationships
#function memo_modularity(A::MatrixNetwork, C::Vector{Vector{Int64}}; M::MatrixNetwork=MatrixNetwork([1],[2]))

#=
Okay, so here we are going to be memoizing the strength between communities for phi and the weight products between communities for omega.
We start out with a Dict wherein the keys are edges between communities and the values are the strength / weight products. So we only need manually compute modularity once.

=#

"""
`clump`
=======

Given a MatrixNetwork and a community partition of it, creates a new MatrixNetwork with each community as a node and edges between communities with weights equal to the previous weight between them (REFINE THIS!!!!).

Arguments
---------
    - `A::MatrixNetwork{Real}`: A weighted network
    - `C::Vector{Vector{Int64}}`: The communities of `A`
"""
function clump(A::MatrixNetwork, C::Vector{Vector{Int64}})
    sz = size(C,1)
    rp = collect(1:sz:sz*sz)
    push!(rp, sz*sz+1)
    cl = MatrixNetwork{Float64}(sz, rp, [i for j = 1:sz for i = 1:sz], zeros(sz*sz))
    #println(cl)

    weights = total_edge_weights(A)
    node_communities = zeros(Int64, A.n)
    for i = 1:size(C,1)
        for n in C[i]
            node_communities[n] = i
        end
    end

    for e in keys(weights)
        i,j = node_communities[e]
        cl.vals[cl.rp[i]+j] += weights[e]
    end

    return cl # Returns clumped graph
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
    C = [[i] for i = 1:A.n]
    C_last = deepcopy(C)

    mod, strengths, W, U = motif_modularity(A,C; M=M)
    #println(strengths)

    function merge_communities!(i, j)
        if i == j return end
        C[i] = [C[i]; C[j]]
        setdiff!(C, [C[j]])
        # First step: Find all groups that contain j
        comm_groups = deepcopy(collect(keys(strengths)))
        #j_groups = comm_groups[findall(x -> j in x, comm_groups)] # Find groups that contain j

        for group in comm_groups
            for c in group
                if c != j && haskey(strengths, sort([i,c]))
                    strengths[group] .+= strengths[sort([i,c])]
                end

                if c != i && haskey(strengths, sort([j,c]))
                    strengths[group] .+= strengths[sort([j,c])]
                end
            end
        end

        function relabel(group)
            labels =  setdiff(1:size(C,1), [j])
            return [findfirst(x -> x == i, labels) for i in group]
        end

        #println(comm_groups)
        for k = 1:size(comm_groups,1)
            g = comm_groups[k]
            g[g.==j] .= i # Replace j with i
            g = relabel(collect(Set(g))) # If i already included, remove duplicate
            comm_groups[k] = g
        end
        #println(setdiff(comm_groups, k))
        #map!(x -> relabel(collect(Set(replace(x,findfirst(y -> y == j, x),i)))), comm_groups) # Filter out j and relabel
        filter!(x -> size(x,1) > 1, comm_groups)
        #map!(sort, comm_groups)
        #println(comm_groups)

        new_strengths = Dict()
        for g = 1:length(comm_groups)
            println(comm_groups[g])
            if !isempty(comm_groups[g])
                new_strengths[sort(comm_groups[g])] = strengths[sort(collect(keys(strengths))[g])] #[1,1] is coming up...
            end
        end
        strengths = new_strengths
        #strengths = Dict(comm_groups[g] => strengths[keys(strengths)[g]] for g = 1:length(keys(strengths)))

        # My thoughts are feeling a bit jumbled about this...
        #=
        What am I trying to do here? I'm trying to update the strengths dict when we merge communities i and j. When i and j merge, any motivic strength between them disappears, as does any expected strength.
        This would be simpler if we were dealing with only pairs, but since we can have gruops of more than 2 communities, we need to think about merging within the lists as well. So when merging, we add all the



        Observation: If no communities within a group have pairwise strengths between them, then total motif strength is conserved with merges.
        Observation: Any pairwise strengths get added to the group strength upon the merging of the pair. BUT, the inter-community strength remains the same until they all become one.
        Observation: Only pairwise merges impact the modularity (trivially).
        =#
        println(C)
    end

    function louvain_iter!()
        for i = 1:size(C,1)
            if i > size(C,1) break end
            best = i
            for j = 1:size(C,1)
                if i == j continue end
                new_C = deepcopy(C)
                ordered = sort([i,j])
                #=push!(new_C, [C[i]; C[j]])
                setdiff!(new_C, [C[i], C[j]])
                post_mod = motif_modularity(A, new_C; M=M) # TODO: This is super gross; find formula?=#
                post_mod = mod + (strengths[ordered][1]/W) + (strengths[ordered][2]/U)
                if post_mod > mod
                    mod = post_mod
                    best = j
                end
            end
            merge_communities!(i, best)
        end
    end
    louvain_iter!()
    while C_last != C
        C_last = deepcopy(C)
        louvain_iter!()
    end
    return map(sort, C)
end

A = MatrixNetwork([(1,2), (1,3), (2,1), (2,3), (2,4), (3,1), (3,2), (4,2), (4,5), (5,4), (5,6), (5,7), (6,5), (6,7), (7,5), (7,6), (4,8), (8,4), (8,9), (8,10), (9,8), (9,10), (10,9), (10,8)], 10)
#println(motif_modularity(A, [[1,2,3,4], [5,6,7], [8,9,10]]))
@time println(louvain_motif(A))
#println(A.rp)

karate_edges = [[2,1],[3,1],[3,2],[4,1],[4,2],[4,3],[5,1],[6,1],[7,1],[7,5],[7,6],[8,1],[8,2],[8,3],[8,4],[9,1],
[9,3],[10,3],[11,1],[11,5],[1,6],[12,1],[13,1],[13,4],[14,1],[14,2],[14,3],[14,4],[17,6],[17,7],[18,1],[18,2],
[20,1],[20,2],[22,1],[22,2],[26,24],[26,25],[28,3],[28,24],[28,25],[29,3],[30,24],[30,27],[31,2],[31,9],[32,1],
[32,25],[32,26],[32,29],[33,3],[33,9],[33,15],[33,16],[33,19],[33,21],[33,23],[33,24],[33,30],[33,31],[33,32],[34,9],
[34,10],[34,14],[34,15],[34,16],[34,19],[34,20],[34,21],[34,23],[34,24],[34,27],[34,28],[34,29],[34,30],[34,31],[34,32],[34,33]]
karate_edges = reduce(vcat, [[i, reverse(i)] for i in karate_edges])
karate_network = MatrixNetwork(map(x->x[1],karate_edges),map(x->x[2],karate_edges))
#@time l = louvain_motif(karate_network)
#println(l)
