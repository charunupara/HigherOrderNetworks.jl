using MatrixNetworks
using SparseArrays
using StatsFuns
using Dates
"""
`K_Model`
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

Base.copy(M::K_Model) = K_Model(M.k, deepcopy(M.P))

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
        for i = 1:size(q,1)-k # Get all k-subpaths of q
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
        P[p] = log(p_counts[p] / w_counts[p[1:k]]) # Pr = proportion of paths identical to p[1:k] that end in p[k]
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
    if k == 0 # Return nodes
        return collect(1:size(AL,1))
    elseif k == 1
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
function optimal_multiorder(S::Vector{Vector{Int64}}, K::Int64; Log=true) # TODO: Log space?
    t1 = now()
    models = [optimal_k_model(S, k-1) for k = 1:K+1] # Generate each layer
    println("Layer models generated for K=$K. Took $(now() - t1).")
    P = Dict{Vector{Int64}, Float64}()

    for p in S
        P[p] = 0.0
        for k = 1:min(K,size(p,1)) # Multiply by probability in each layer. E.g., P(a,b,c) = P(a) * P(b|a) * P(c|a,b)
            P[p] += models[k].P[p[1:k]]
        end
        for i = K+1:size(p,1) # Multiply by M_K+1 probability to produce path. For paths longer than K
            P[p] += models[K+1].P[p[i-K:i]]
        end
    end

    return K_Model(K, P)
end

"""
`multiorder_likelihood`
=======================

Computes the likelihood that a multi-order model `M_k` produces the observed paths `S`

Precondition: S ⊆ M_k.P.keys
"""
function multiorder_likelihood(M_k::K_Model, S::Vector{Vector{Int64}}; Log=true)
    @assert issubset(Set(S), keys(M_k.P))

    if Log
        l = 0.0
        for p in S
            l += M_k.P[p]
        end
        return l
    end
    l = 1.0
    for p in S
        l *= M_k.P[p]
    end
    return l
end

"""
`K_opt`
=======

Compute the path length `k` by which the observed paths are most aptly modeled

Arguments
---------
    - `S::Vector{Vector{Int64}}`: The observed paths
    - `A::MatrixNetwork`: The graph underlying `S`
    - `p::Float64 (=0.05)`: The significance threshold for `k`

Returns
-------
The optimal value for k and the corresponding multi-order model of order k
"""
function K_opt(S::Vector{Vector{Int64}}, A::MatrixNetwork; p::Float64=0.01)
    K_opt = 1
    t1 = now()
    M_k = optimal_multiorder(S, K_opt)
    println("Initial k=1 model generated. Took $(now() - t1).")

    function p_value(M_k::K_Model, M_K::K_Model, S::Vector{Vector{Int64}})
        function matrix_df(A::Array{Bool,2}, k::Int64)
            Ak = A^k
            return sum(Ak) - sum([1 for i = 1:size(Ak,1) if sum(Ak[i,:]) != 0])
        end

        function df(M::K_Model)
            return (A.n - 1) + sum([matrix_df(Array(sparse(A)), k) for k = 1:M.k])
        end

        M_kl = multiorder_likelihood(M_k, S)
        M_Kl = multiorder_likelihood(M_K, S)
        println("k-Likelihood: $(M_kl)")
        println("K-Likelihood: $(M_Kl)")

        println(M_Kl - M_kl)

        dfMk = df(M_k)
        dfMK = df(M_K)
        println("$dfMK, $dfMk")
        return chisqcdf(dfMK-dfMk, 2(M_Kl-M_kl))
    end

    while true
        t2 = now()
        M_K = optimal_multiorder(S, K_opt + 1)
        println("Higher k=$(K_opt+1) model generated. Took $(now() - t1).")
        t3 = now()
        cp = p_value(M_k, M_K, S)
        println("Comparative p-value between k=$K_opt and k=$(K_opt+1): $cp. Took $(now() - t3).")
        if cp < p
            break
        else
            M_k = deepcopy(M_K)
            K_opt += 1
        end
    end
    println("Done. k=$K_opt. Total analysis took $(now() - t1).")
    return K_opt, M_k
end


AL = [[2], [4,3], [1], [1,2]]
paths = Vector{Vector{Int64}}()

#=for k = 1:20
    append!(paths, generate_paths(AL, k))
end=#
#=paths = [[2,4], [2,3], [4,1], [4,2], [1,2], [2,3,1], [1,2,4],
         [4,1,2], [2,4,2], [3,1,2], [4,2,4], [2,4,1], [1,2,3],
         [2,4,2,4], [4,1,2,4], [1,2,3,1], [1,2,4,2,4], [4,2,4,2,4],
         [3,1,2,4,2,4], [2,4,2,4,2,4]]=#

#println(optimal_k_model(paths,3).P)
#println(paths)
#println(optimal_multiorder(paths, 3).P) # Numbers slightly off from paper...
#println("wheef")
#println(K_opt(paths, MatrixNetwork([1,2,2,3,4,4], [2,4,3,1,1,2]))[1])

function read_paths(filepath::String; frequency=false)
    println("Begin reading $filepath.")
    nodes = Dict()
    paths = Array{Array{Int64,1},1}()
    edges = Set{Tuple{Int64,Int64}}()
    t1 = now()
    open(filepath) do file
        lines = [ln for ln in eachline(file)]
        for ln in lines
            separated = split(ln, ",")
            p = []
            for i = 1:(frequency ? length(separated)-1 : length(separated))
                n = separated[i]
                if !haskey(nodes,n)
                    nodes[n] = length(nodes) + 1
                end
                if i < length(separated)-1 && !haskey(nodes,separated[i+1])
                    nodes[separated[i+1]] = length(nodes) + 1
                end
                push!(p, nodes[n])
                if i < length(separated)-1
                    push!(edges, (nodes[n], nodes[separated[i+1]]))
                end
            end
            if frequency
                for i = 1:Int(parse(Float64, separated[length(separated)]))
                    push!(paths, p)
                end
            else
                push!(paths, p)
            end
        end
    end
    println("$(length(paths)) paths read, implying a network of $(length(nodes)) nodes and $(length(edges)) edges. Took $(now() - t1).")
    return paths, MatrixNetwork(collect(edges), length(nodes))
end

flights = "US_flights.ngram"
tube = "tube_paths_train.NGRAM"
wiki = "data\\wikipedia_clickstreams.NGRAM"
@time println(K_opt(read_paths("C:\\Users\\Elizabeth Turner\\Documents\\Josh\\MAP\\HO Networks\\Code\\$wiki")...)[1])

# TODO: Include layers in multi-order models
# TODO: Since multi model K is included in K+1, we can supply K to make generation faster
