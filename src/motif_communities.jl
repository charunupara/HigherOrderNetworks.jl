using MatrixNetworks
using Dates
include("Utilities.jl")
include("clique_m.jl")
include("clique_j.jl")

"""
Common motif generators
    - `path(k::Int64)`: Returns a path of length `k`
    - `cycle(k::Int64)`: Returns a cycle of length `k`
    - `clique(k::Int64)`: Returns a clique of size `k`
"""
path(k::Int64) = MatrixNetwork([i for i = 1:k], [j+1 for j = 1:k])
cycle(k::Int64) = MatrixNetwork([i for i = 1:k], [[j+1 for j = 1:k-1]; [1]])
clique(k::Int64) = MatrixNetwork([i for i = 1:k for j = 1:k-1], [j for i = 1:k for j = 1:k if j != i])

function is_clique(M::MatrixNetwork)
    c = clique(M.n)
    return M.ci == c.ci && M.rp == c.rp # Checks whether a motif is a clique
end

"""
`motif_match`

"""
function motif_match(A::MatrixNetwork, M::MatrixNetwork)
    if is_clique(M)
        return find_clique(A, M.n), collect(1:M.n)
    end

    a_e = get_edges(A)
    m_e = get_edges(M)

    a_d = deg_distr(A)
    m_d = deg_distr(M)

    undirected = false

    ### Orders motif vertices such that non-matching branches will be cut off ASAP ###
    function greatest_constraint_first()
        highest = nodes_by_deg(M; descending=true)[1]
        vertices = setdiff(1:M.n, [highest])
        order = [highest]
        parents = [-1]
        next = -1

        are_connected(u,v) = return ([u,v] in m_e || [v,u] in m_e)

        while !isempty(vertices)
            m = length(order)
            next_rank = (-Inf, -Inf, -Inf)

            for u in setdiff(vertices, order)
                V1 = [v for v in order if are_connected(u,v)] # Members of path u is connected to
                V2 = [v for v in order if !isempty(setdiff(intersect(M.ci[M.rp[u]:M.rp[u+1]-1], M.ci[M.rp[v]:M.rp[v+1]-1]), V1))] # Members of path that u shares neighbors with
                #println("$(u-1), $(V2 .- 1)")
                V3 = [] # TODO: Get this

                if next_rank <= (length(V1), length(V2), length(V3)) || length(vertices) == 1
                    #println("Previous: $(next-1), $next_rank | Next: $(u-1), $((length(V1), length(V2), length(V3)))")
                    next = u
                    next_rank = (length(V1), length(V2), length(V3))
                end
            end
            pt = 0
            for v in order
                if are_connected(next,v) pt = v end
            end

            push!(order, next)
            push!(parents, pt)
            setdiff!(vertices, [next])
            #println(vertices)
        end

        return order, parents
    end

    ### Tests whether u ∈ V(M) and x ∈ V(A) are a valid match for graph isomorphism ###
    function meets_isomorphism(u, x, p1, p2)
        if u in p1
            #println("$u already in path")
            return false
        end # Condition 1
        if x in p2
            #println("$x already in path")
            return false
         end

        if a_d[x] < m_d[u]
            #println("Degree mismatch between $u and $x")
            return false
        end # Condition 3
        # Condition 4
        for i = 1:length(p1)
            if [u,p1[i]] in m_e && !([x,p2[i]] in a_e)
                #println("Edge $([x,p2[i]]) not in graph")
                return false
            end
            if [p1[i],u] in m_e && !([p2[i],x] in a_e)
                #println("Edge $([p2[i],x]) not in graph")
                return false
            end
        end
        #println("Match between $u and $x")
        return true # We made it!
    end

    mv_order, mv_parents = greatest_constraint_first()
    matches = Set()

    function match_path(p1, p2, new_v, available_vertices)
        l = length(p1) + 1
        if !isempty(p1) pt = p2[findfirst(x->x==mv_parents[l], p1)] end

        if isempty(p1) || (meets_isomorphism(mv_order[l], new_v, p1, p2) && ([pt, new_v] in a_e || [new_v, pt] in a_e))
            if l == M.n # Base case
                push!(matches, undirected ? sort([p2; [new_v]]) : [p2; [new_v]]) # If the graph is undirected, don't include "duplicate" motifs (i.e. don't include all vertex permutations)
                return
            end
            for x in available_vertices # Recursive case, try all next vertices
                match_path([p1; [mv_order[l]]], [p2; [new_v]], x, setdiff(available_vertices, [x]))
            end
        else
            return
        end
    end

    for v = 1:A.n
        if meets_isomorphism(mv_order[1], v, [], [])
            match_path([], [], v, setdiff(1:A.n, [v]))
        end
    end
    return collect(matches), mv_order
end

paper_example = MatrixNetwork([i+1 for i in [0,0,0,0,1,1,1,1,2,2,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,7,7,7,7,8,8,9,9,9,10]],
                              [i+1 for i in [9,3,4,1,0,4,5,2,1,5,0,6,0,6,1,7,5,2,1,4,7,8,9,3,4,7,6,4,5,8,5,7,0,6,10,9]])
my_example = MatrixNetwork([1,2,2,2,3,4,4,5,5,5], [2,1,4,5,5,2,5,2,3,4])
other_example = MatrixNetwork([(1,2), (2,1), (3,1), (3,2)], 3)
#println(motif_match(other_example, MatrixNetwork([(1,2), (3,2)], 3))) # Should find [2,1,3], [1,2,3]

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
        if length(inds) != length(Set(inds)) # Requires unique indices
            inc_inds!(inds)
        end
    end # TODO: Add undirected option?

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

        motif_membership = Set() # Communities that the motif is part of

        for e in m_e
            if !(indices[e] in a_e) # Motif edge not in graph
                w_inc = 0 # When we have efficient motif finding, get rid of this
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

        if !isempty(motif_membership) && reduce(+, length.(C[collect(motif_membership)])) >= M.n
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

function fast_modularity(A::MatrixNetwork, C::Vector{Vector{Int64}}; M::MatrixNetwork=MatrixNetwork([1],[2]))
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

    w,W,u,U = 0,0,0,0 # w: strength of motifs fully contained within communities, W: total strength of motifs
                      # u: w for a random, strength-preserving graph, U: W for the same

    motifs, order = motif_match(A, M)

    for m in motifs
        weight = 1
        p = 1

        motif_membership = Set(node_communities[m])
        d = Dict(order[i] => m[i] for i = 1:length(order))
        for e in m_e
            weight *= weights[[d[e[1]], d[e[2]]]]
            p *= out_w[d[e[1]]] * in_w[d[e[2]]]
        end

        if length(motif_membership) == 1
            w += weight
            u += p
        end
        W += weight
        U += p
        if length(motif_membership) > 1 && length(reduce(vcat, C[collect(motif_membership)])) >= M.n
            motif_membership = sort(collect(motif_membership))
            if !haskey(shared_strengths, motif_membership)
                shared_strengths[motif_membership] = [0,0]
            end
            shared_strengths[motif_membership][1] += weight
            shared_strengths[motif_membership][2] += p
        end
    end

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
    ### Collapse all communities in the group into one. Updates C and strengths ###
    function merge_community_group!(group)
        new_groups = collect(deepcopy(keys(strengths)))
        cpy = deepcopy(new_groups) # Remember old keys so we can get rid of ones that no longer exist at the end

        new_strengths = Dict()

        i = group[1]

        for k = 1:size(new_groups,1)
            g = new_groups[k]
            for c in group
                g[g.==c] .= i # Replace all other group members with i
            end

            g = sort(collect(Set(g))) # If i already included, remove duplicates and sort for consistency

            if length(g) < 2 continue end # Leave out degenerate groups

            if !haskey(new_strengths, g)
                new_strengths[g] = strengths[cpy[k]]
            else
                new_strengths[g] += strengths[cpy[k]]
            end
        end

        for g in setdiff(cpy, new_groups) # Erase entries that no longer exist
            delete!(new_strengths,g)
        end

        old_keys = collect(deepcopy(keys(new_strengths))) # THIS IS THE PART SLOWING THINGS DOWN!
        new_keys = []
        vals = []

        new_labels = collect(1:length(C))
        for c in group[2:length(group)]
            new_labels[c:length(new_labels)] .-= 1
        end

        for k = 1:length(old_keys) # Relabel community relationships to account for merge
            key = old_keys[k]
            push!(vals, new_strengths[key])
            key = [new_labels[k] for k in key]
            push!(new_keys, key)
        end

        # Update partition
        C[i] = sort(reduce(vcat, C[group]))
        for c in group[2:length(group)] setdiff!(C, [C_last[c]]) end

        strengths = Dict(new_keys[i] => vals[i] for i = 1:length(new_keys)) # Ta da!
    end

    ### Perform one iteration of the modified Louvain algorithm. Finds the group of communities that, if collapsed into one, would most increase modularity and performs it. ###
    function louvain_iter!()
        max_dQ = 0.0 # The merge must increase modularity
        best = -1
        str_k = filter(x -> strengths[x][1] > 0, collect(keys(strengths)))
        for k = 1:length(str_k)
            dQ = (strengths[str_k[k]][1]/W) - (strengths[str_k[k]][2]/U)
            if dQ > max_dQ
                max_dQ = dQ
                best = k
            end
        end
        if best != -1 merge_community_group!(str_k[best]) end
        mod += max_dQ
    end
    louvain_iter!() # Get the ball rolling, make C different from C_last

    while C_last != C # Stop merging when no merge increases modularity
        C_last = deepcopy(C)
        louvain_iter!()
    end

    return sort.(C)
end

A = MatrixNetwork([(1,2), (1,3), (2,1), (2,3), (2,4), (3,1), (3,2), (4,2), (4,5), (4,8), (5,4), (5,6), (5,7), (6,5), (6,7), (7,5), (7,6), (8,4), (8,9), (8,10), (9,8), (9,10), (10,9), (10,8)], 10)
#@time println(motif_modularity(A, [[1,2,3,4],[5,6,7],[8,9,10]]; M=clique(3)))
#@time println(fast_modularity(A,[[1,2,3,4],[5,6,7],[8,9,10]]))

star = MatrixNetwork([1,1,1,1,1,1,1,1,2,3,4,5,6,7,8,9], [2,3,4,5,6,7,8,9,1,1,1,1,1,1,1,1])

tri_circ = MatrixNetwork([1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,7,7,8,8,9,9,10,10],
                         [2,3,4,5,1,3,4,5,1,2,4,5,1,2,3,5,1,2,3,4,6,5,7,10,6,8,7,9,8,10,9,6])

karate_edges = [[2,1],[3,1],[3,2],[4,1],[4,2],[4,3],[5,1],[6,1],[7,1],[7,5],[7,6],[8,1],[8,2],[8,3],[8,4],[9,1],
[9,3],[10,3],[11,1],[11,5],[1,6],[12,1],[13,1],[13,4],[14,1],[14,2],[14,3],[14,4],[17,6],[17,7],[18,1],[18,2],
[20,1],[20,2],[22,1],[22,2],[26,24],[26,25],[28,3],[28,24],[28,25],[29,3],[30,24],[30,27],[31,2],[31,9],[32,1],
[32,25],[32,26],[32,29],[33,3],[33,9],[33,15],[33,16],[33,19],[33,21],[33,23],[33,24],[33,30],[33,31],[33,32],[34,9],
[34,10],[34,14],[34,15],[34,16],[34,19],[34,20],[34,21],[34,23],[34,24],[34,27],[34,28],[34,29],[34,30],[34,31],[34,32],[34,33]]

karate_edges = reduce(vcat, [[i, reverse(i)] for i in karate_edges])
K = MatrixNetwork(map(x->x[1],karate_edges),map(x->x[2],karate_edges))

# Todo ranked by priority
# TODO: Efficient motif detection: works for both directed and undirected, but we need to get U!
#   - It's actually slower than brute force! What the heck... I guess we ditch this?
# TODO: Try Charun's triangle detection

# Optional:
# TODO: Include option of specifying which edges k_delt should apply to
