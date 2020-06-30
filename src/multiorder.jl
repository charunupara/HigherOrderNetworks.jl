using GSL
using MatrixNetworks
"""
`k_model`
=========

Represents a graph model of order k as described in https://arxiv.org/pdf/1702.05499.pdf

Fields
------
    - `k::Int64`: The order of the model, which denotes the lengths of the paths
    - `P::Dict{Vector{Int64}, Float64}`: The probabilities that particular paths will be generated
"""
struct K_Model
    k::Int64
    P::Dict{Vector{Int64}, Float64}
end

"""
`k_model_k`
===========

`k_model` constructor that verifies that `k` is non-negative
"""
function k_model_k(k::Int64, P::Dict{Vector{Int64}, Float64})
    @assert k >= 0
    return K_Model(k, P)
end

"""
`optimal_k_model`
=================

Given a multiset of observed paths `S`, constructs the k_model that is most likely to generate `S`

Arguments
---------
    - `S::Vector{Vector{Int64}}`: Observed paths
    - `k::Int64`: Order of the desired model

Example
-------
~~~~
optimal_k_model([[1,2], [1,3], [1,2,3], [2,3], [3,2], [3,1,2]], 1)
~~~~
"""
function optimal_k_model(S::Vector{Vector{Int64}}, k::Int64)
    P = Dict{Vector{Int64}, Float64}() # Holds probabilities
    p_counts = Dict{Vector{Int64}, Int64}() # Number of times a path of length `k` appears in S
    w_counts = Dict{Vector{Int64}, Int64}() # Number of times a path of length `k-1` appears in S
    S_k = []
    for q in S
        for i = 1:size(q,1)-k # Get all subpaths of q
            p_k = q[i:i+k]
            w_k = p_k[1:k]

            push!(S_k, p_k)

            if !haskey(p_counts, p_k)
                p_counts[p_k] = 0
            end
            p_counts[p_k] += 1

            if !haskey(w_counts, w_k)
                w_counts[w_k] = 0
            end
            w_counts[w_k] += 1
        end
    end

    for p in Set(S_k)
        P[p] = p_counts[p] / w_counts[p[1:k]] # Pr = proportion of paths identical to p[1:k] that end in p[k]
    end
    return k_model_k(k, P)
end

"""
`generate_paths`
================

Recursively generate all the paths of a certain length given an adjacency list

Arguments
---------
    - `AL::Vector{Vector{Int64}}`: Adjacency list where `AL[i]` is all nodes adjacent to node `i`
    - `k::Int64`: The path length

Example
-------
~~~~
generate_paths([[2], [3,4], [1], [1]], 2)
~~~~
"""
function generate_paths(AL::Vector{Vector{Int64}}, k::Int64)
    if k == 1
        return [[i, j] for i = 1:size(AL,1) for j in AL[i]] # Return edges
    else
        k_paths = Vector{Vector{Int64}}()
        for p in copy(generate_paths(AL, k - 1)) # Recurse
            for n in AL[p[k]] # Get neighbors of last vertex in path
                p_c = copy(p)
                push!(p_c, n)
                push!(k_paths, p_c)
            end
        end
        return k_paths
    end
end

"""
`optimal_multiorder`
====================

Generates the optimal multi-order model for a set of observed paths

Arguments
---------
    - `S::Vector{Vector{Int64}}`: The observed paths
    - `K::Int64`: The maximum order of the model
"""
function optimal_multiorder(S::Vector{Vector{Int64}}, K::Int64)
    models = [optimal_k_model(S, k-1) for k = 1:K+1] # Generate each layer
    P = Dict{Vector{Int64}, Float64}()

    for p in S
        P[p] = 1.0
        for k = 1:min(K,size(p,1)) # Multiply by probability in each layer
            P[p] *= models[k].P[p[1:k]]
        end
        for i = K+1:size(p,1) # Multiply by M_K+1 probability to produce path
            P[p] *= models[K+1].P[p[i-K:i]]
        end
    end

    return k_model(K, P)
end

"""
`multiorder_likelihood`
=======================

Computes the likelihood that a multi-order model `M_k` produces the observed paths `S`

Precondition: S ⊆ M_k.P.keys
"""
function multiorder_likelihood(M_k::K_Model, S::Vector{Vector{Int64}})
    @assert issubset(Set(S), Set(key(M_k.P)))
    l = 1.0
    for p in S
        l *= M_k.P[p] # Product of likelihood for each path
    return l
end

"""
`K_opt`
=======

Compute the path length `k` by which the observed paths are most aptly modeled

Arguments
---------
    - `S::Vector{Vector{Int64}}`: The observed paths
    - `p::Float64 (=0.05)`: The significance threshold for `k`

Returns
-------
The optimal value for k and the corresponding multi-order model of order k
"""
function K_opt(S::Vector{Vector{Int64}}; p::Float64=0.05)
    K_opt = 1
    M_k = optimal_multiorder(S, K_opt)

    function p_value(M_k::K_Model, M_K::K_Model, S::Vector{Vector{Int64}})
        M_kl = multiorder_likelihood(M_k, S)
        M_Kl = multiorder_likelihood(M_K, S)

        function matrix_df(A::Vector{Float64,2}, k::Int64)
            Ak = A^k
            return sum([Ak[i,j] for i = 1:size(A,1) for j = 1:size(A,2)]) - sum([1 for i = 1:size(A,1) if !all(x -> x == 0, A[i,:]))
        end

        function df(M::K_Model)
            return (size(keys(M.P),1) - 1) + sum([matrix_df(MATRIX, k) for k = 1:M.K) #TODO: Figure out what MATRIX is
        end

        return 1 - sf_gamma_inc((df(M_K)-df(M_k))/2, -log10(M_kl/M_Kl))
    end

    while true
        M_K = optimal_multiorder(S, K_opt + 1)
        if p_value(M_k, M_K, S) < p
            break
        else
            M_k = M_K
            K_opt += 1
        end
    end
    return K_opt, M_k
end


AL = [[2], [4,3], [1], [1,2]]
paths = Vector{Vector{Int64}}()
for k = 1:3
    append!(paths, generate_paths(AL, k))
end
paths = generate_paths(AL, 4)
#=paths = [[2,4], [2,3], [4,1], [4,2], [1,2], [2,3,1], [1,2,4],
         [4,1,2], [2,4,2], [3,1,2], [4,2,4], [2,4,1], [1,2,3],
         [2,4,2,4], [4,1,2,4], [1,2,3,1], [1,2,4,2,4], [4,2,4,2,4],
         [3,1,2,4,2,4], [2,4,2,4,2,4]]=#
#println(optimal_k_model(paths,3).P)
#println(paths)
println(optimal_multiorder(paths, 3).P) # Numbers slightly off from paper...
