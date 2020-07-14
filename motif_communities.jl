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

    a_e = get_edges(A)
    m_e = get_edges(M)
    out_w = outweight_distr(A)
    in_w = inweight_distr(A)
    weights = total_edge_weights(A)

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

    function motif_count() # Calculates the proportion of motif occurrences that do not cross communities
        indices = ones(Int64, M.n)
        z = zeros(Int64, M.n)
        n = 0
        N = 0

        while indices != z
            n_inc = 1
            N_inc = 1
            weight = 1
            for e in m_e
                if !(indices[e] in a_e)
                    n_inc = 0
                    N_inc = 0
                    break
                end
                if !k_delt(indices[e]...)
                    n_inc = 0
                end
                weight *= weights[[indices[e[1]], indices[e[2]]]]
            end
            n += n_inc * weight
            N += N_inc * weight
            inc_inds!(indices)
        end
        return n / N
    end

    function motif_expectation() # Calculates the proportion of intra-community motifs that would appear by chance
        indices = ones(Int64, M.n)
        z = zeros(Int64, M.n)
        n = 0
        N = 0

        while indices != z
            p = 1
            nf = 1
            for e in m_e
                if !k_delt(indices[e]...)
                    nf = 0
                end
                p *= out_w[indices[e[1]]] * in_w[indices[e[2]]]
            end

            n += p * nf
            N += p
            inc_inds!(indices)
        end

        return n / N
    end

    return motif_count() - motif_expectation()
end

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
#=function louvain_motifa(A::MatrixNetwork, C::Vector{Vector{Int64}}; M::MatrixNetwork=MatrixNetwork([1],[2]))
    node_c = collect(1:A.n)
    C_last = deepcopy(C)

    CL = clump(A,C)

    function merge_communities!()
        for i = 1:size(C,1)
            cur_mod = motif_modularity(A, C; M=M)
            best = i
            for c in CL.ci[CL.rp[i]:CL.rp[i+1]-1]
                new_C = deepcopy(C)
                setdiff!(new_C, [C[c], C[i]])
                push!(new_C, [C[c]; C[i]])
                post_mod = motif_modularity(A, new_C; M=M)
                if post_mod > cur_mod
                    cur_mod = post_mod
                    best = c
                end
            end
            setdiff!(C, [C[best], C[i]])
            push!(C, [C[best]; C[i]])
            CL = clump(A,C)
        end
    end

    function move_nodes!()
        for i = 1:A.n
            cur_mod = motif_modularity(A, C; M=M)
            best = node_c[i]
            for j in A.ci[A.rp[i]:A.rp[i+1]-1]
                new_C = deepcopy(C)
                #println(node_c[j])
                setdiff!(new_C[node_c[i]], [i])
                push!(new_C[node_c[j]], i)
                post_mod = motif_modularity(A, new_C; M=M)
                if post_mod > cur_mod
                    cur_mod = post_mod
                    best = j
                end
            end
            setdiff!(C[node_c[i]], [i])
            push!(C[node_c[best]], i)
            node_c[i] = node_c[best]
        end
        filter!(x -> !isempty(x), C)

        node_set = sort(collect(Set(node_c)))
        for i = 1:size(node_c,1)
            node_c[i] = findfirst(x -> x == node_c[i], node_set)
        end
        println(node_c)
    end
    #=move_nodes!()

    while C_last != C
        C_last = deepcopy(C)
        A = clump(A,C) # Node i represents the ith community
        move_nodes!()
    end
    return node_c=#

    #com_con = [CL.ci[[CL.rp[i]+j for j = 1:CL.rp[i+1]-CL.rp[i] if CL.vals[CL.rp[i]+j] > 0]] for i = 1:size(C,1)] # Get communities that each community is connected to
    com_con = [[] for i = 1:size(C,1)]
    for i = 1:size(C,1)
        for j = 1:CL.rp[i+1]-CL.rp[i]
            #if nonzero connection between communities i and j,
        #        add j to com_con[i]
            if i != size(C,1) && CL.vals[CL.rp[i]+j] > 0
                push!(com_con[i], j)
            end
        end
    end
    println(CL.ci)
    println(com_con)
    merge_communities!()
    while C_last != C
        C_last = deepcopy(C)
        merge_communities!()
    end
    return C
end=#

function louvain_motif(A::MatrixNetwork; M::MatrixNetwork=MatrixNetwork([1],[2]))
    C = [[i] for i = 1:A.n]
    C_last = deepcopy(C)
    function merge_communities!(i, j)
        if i == j return end
        C[i] = [C[i]; C[j]]
        setdiff!(C, [C[j]])
    end

    function louvain_iter!()
        for i = 1:size(C,1)
            if i > size(C,1) break end
            cur_mod = motif_modularity(A, C; M=M)
            best = i
            for j = 1:size(C,1)
                new_C = deepcopy(C)
                push!(new_C, [C[i]; C[j]])
                setdiff!(new_C, [C[i], C[j]])
                post_mod = motif_modularity(A, new_C; M=M) # TODO: This is super gross; find formula?
                if post_mod > cur_mod
                    cur_mod = post_mod
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
#println(louvain_motif(A))
#println(A.rp)

karate_edges = [[2,1],[3,1],[3,2],[4,1],[4,2],[4,3],[5,1],[6,1],[7,1],[7,5],[7,6],[8,1],[8,2],[8,3],[8,4],[9,1],
[9,3],[10,3],[11,1],[11,5],[1,6],[12,1],[13,1],[13,4],[14,1],[14,2],[14,3],[14,4],[17,6],[17,7],[18,1],[18,2],
[20,1],[20,2],[22,1],[22,2],[26,24],[26,25],[28,3],[28,24],[28,25],[29,3],[30,24],[30,27],[31,2],[31,9],[32,1],
[32,25],[32,26],[32,29],[33,3],[33,9],[33,15],[33,16],[33,19],[33,21],[33,23],[33,24],[33,30],[33,31],[33,32],[34,9],
[34,10],[34,14],[34,15],[34,16],[34,19],[34,20],[34,21],[34,23],[34,24],[34,27],[34,28],[34,29],[34,30],[34,31],[34,32],[34,33]]
karate_edges = reduce(vcat, [[i, reverse(i)] for i in karate_edges])
karate_network = MatrixNetwork(map(x->x[1],karate_edges),map(x->x[2],karate_edges))
@time l = louvain_motif(karate_network; M=MatrixNetwork([1,2],[2,3]))
println(l)
