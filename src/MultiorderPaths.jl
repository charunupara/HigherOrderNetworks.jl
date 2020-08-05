"""
A model for determining the memory of agents traversing a network given a list
of observed paths. For example, it can answer the question, "how many webpages
back do we have to look to predict where someone surfing the web will go next?"

Based on "When is a Network a Network? Multi-Order Graphical Model Selection in
Pathways and Temporal Networks" by Ingo Scholtes. The original paper may be
found at https://arxiv.org/pdf/1702.05499.pdf.

Also of note is Scholtes' argument that traditional network analytic methods are
only justified if the optimal memory for modeling paths is 1, i.e. the Markov
property holds.
"""

using MatrixNetworks # Used for realizing the network implied by observed paths
using SparseArrays # Used for going from MatrixNetwork to a matrix
using StatsFuns # Used for the χ² CDF
using Dates # Used for runtime calculations

"""
`HO_Markov`
=========

Represents a higher-order Markov model on a set of observed paths. (pp. 2)

Fields
------
    - `k::Int64`: The order of the model, which denotes the length of the paths
    - `P::Dict{Vector{Int64}, Float64}`: The log probabilities that particular paths
                                         will be generated. More specifically,
                                         given the first `k-1` steps of the path,
                                         the log probability that step `k` occurs

While traditional Markov models are memoryless, in that where they go next is
based solely on their current location, "higher-order" Markov models determine
the probability of each possible next location based on the `k` previous locations.
"""
struct HO_Markov
    k::Int64
    P::Dict{Vector{Int64}, Float64}
end

"""
`MultiOrderModel`
=================

Represents a multi-order model consisting of many layers of `HO_Markov` models.
(pp. 4)

Fields
------
    - `K::Int64`: The order of the highest `HO_Markov` model
    - `layers::Vector{HO_Markov}`: Array of `HO_Markov` models, going from orders
                                   0 to K
    - `P::Dict{Vector{Int64}, Float64}`: The cumulative log probabilities for the
                                         realization of each path based on the
                                         layers. For example, P(a,b,c) = P(a)P(b|a)P(c|a,b)
"""
struct MultiOrderModel
    K::Int64
    layers::Vector{HO_Markov}
    P::Dict{Vector{Int64}, Float64}
end

"""
`optimal_HO_Markov`
=================

Given a multiset of observed paths `S`, constructs the `HO_Markov` that is most
likely to generate `S`. (pp. 3)

Arguments
---------
    - `S::Vector{Vector{Int64}}`: Observed paths
    - `k::Int64`: Order of the desired model

Example
-------
~~~~
optimal_HO_Markov([[1,2], [1,3], [1,2,3], [2,3], [3,2], [3,1,2]], 1)
~~~~
"""
function optimal_HO_Markov(S::Vector{Vector{Int64}}, k::Int64)
    P = Dict{Vector{Int64}, Float64}() # Log probabilities
    p_counts = Dict{Vector{Int64}, Int64}() # Number of times a path of length `k` appears in S
    w_counts = Dict{Vector{Int64}, Int64}() # Number of times a path of length `k-1` appears in S
    S_k = []
    for q in S
        for i = 1:size(q,1)-k # Get all k-subpaths of q
            p_k = q[i:i+k] # k-subpath
            w_k = p_k[1:k] # p_k excluding the last node

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
                                                   # E.g., P(a,b,c) = |{x = (a,b,c) for x in S}| / |{x = (a,b,d) for x in S for all observed nodes d}|
    end
    return HO_Markov(k, P)
end

"""
`optimal_multiorder`
====================

Generates the optimal multi-order model of maximum order `K` for a set of
observed paths. (pp. 4)

Arguments
---------
    - `S::Vector{Vector{Int64}}`: The observed paths
    - `K::Int64`: The maximum order of the model
"""
function optimal_multiorder(S::Vector{Vector{Int64}}, K::Int64)
    t1 = now()
    layers = [optimal_HO_Markov(S, k-1) for k = 1:K+1] # Generate each layer
    println("Layer models generated for K=$K. Took $(now() - t1).")
    P = Dict{Vector{Int64}, Float64}()

    for p in S
        P[p] = 0.0
        for k = 1:min(K,size(p,1)) # Multiply by probability in each layer. E.g., P(a,b,c) = P(a) * P(b|a) * P(c|a,b)
            P[p] += layers[k].P[p[1:k]]
        end
        for i = K+1:size(p,1) # Multiply by M_K probability to produce path. For paths longer than K
            P[p] += layers[K+1].P[p[i-K:i]]
        end
    end

    return MultiOrderModel(K, layers, P)
end

"""
`increment_model`
=================

Generate a multi-order model of one order higher than the one given.

Arguments
---------
    - `model::MultiOrderModel`: The model to increment

This method takes advantage of the nestedness of multi-order models: that a model
of order K+1 contains all the layers of a model of order K.
"""
function increment_model(model::MultiOrderModel)
    S = collect(keys(model.P))
    K = model.K + 1
    layers = [model.layers; [optimal_HO_Markov(S, K)]]
    P = Dict{Vector{Int64}, Float64}()

    for p in S
        P[p] = 0.0
        for k = 1:min(K,size(p,1)) # Multiply by probability in each layer. E.g., P(a,b,c) = P(a)P(b|a)P(c|a,b)
            P[p] += layers[k].P[p[1:k]]
        end
        for i = K+1:size(p,1) # Multiply by M_K probability to produce path. For paths longer than K
            P[p] += layers[K+1].P[p[i-K:i]]
        end
    end

    return MultiOrderModel(model.K+1, layers, P)
end

"""
`multiorder_likelihood`
=======================

Compute the likelihood that a multi-order model `M_k` produces the observed paths `S`
by taking the product of all path probabilities. (pp. 4)

Precondition
------------
    - S ⊆ M_k.P.keys
"""
function multiorder_likelihood(M_k::MultiOrderModel, S::Vector{Vector{Int64}})
    @assert issubset(Set(S), keys(M_k.P))
    return sum([M_k.P[p] for p in S])
end

"""
`K_opt`
=======

Compute the memory length `K` by which the observed paths are optimally modeled.
(pp. 4,5)

Arguments
---------
    - `S::Vector{Vector{Int64}}`: The observed paths
    - `A::MatrixNetwork`: The graph underlying `S`
    - `p::Float64 (=0.01)`: The significance threshold that must be crossed to
                            prefer a simpler model to a more complex one

Returns
-------
    - The optimal value for `K`
    - The corresponding multi-order model of order `K`
"""
function K_opt(S::Vector{Vector{Int64}}, A::MatrixNetwork; p::Float64=0.01)
    K_opt = 0
    t1 = now()
    M_k = optimal_multiorder(S, K_opt)

    println("Initial K=0 model generated. Took $(now() - t1).")

    ### Compute the degrees of freedom in a binary adjacency matrix ###
    function matrix_df(A::Array{Bool,2}, k::Int64)
        Ak = A^k
        return sum(Ak) - sum([1 for i = 1:size(A,1) if any(Ak[i,:] .!= 0)])
    end

    ### Compute the degrees of freedom of a multi-order model ###
    function df(M::MultiOrderModel)
        return (A.n - 1) + sum([matrix_df(Array(sparse(A)), k) for k = 1:M.K])
    end

    M_kl = multiorder_likelihood(M_k, S)
    dfMk = df(M_k)

    ### Calculate the p-value of M_k over M_K###
    function p_value(M_K::MultiOrderModel) # pp. 5
        M_Kl = multiorder_likelihood(M_K, S)
        dfMK = df(M_K)

        return chisqcdf(dfMK-dfMk, 2(M_Kl-M_kl)), M_Kl, dfMK
    end

    while true # Keep increasing maximum order until the significance threshold
        t2 = now()
        M_K = increment_model(M_k)
        println("Higher K=$(K_opt+1) model generated. Took $(now() - t1).")
        t3 = now()
        cp, M_kl, dfMk = p_value(M_K)
        println("Comparative p-value between K=$K_opt and K=$(K_opt+1): $cp. Took $(now() - t3).\n")
        if cp < p
            break
        else
            M_k = deepcopy(M_K)
            K_opt += 1
        end
    end
    println("Done. K_opt=$K_opt. Total analysis took $(now() - t1).")
    return K_opt, M_k
end

"""
`read_paths`
============

Read a collection of paths from a file.

Each line must contain a single path with an optional frequency at the end to
indicate the number of times to count the path.

Arguments
---------
    - `filepath::String`: Where to read paths from
    - `separator::String (=",")`: The string that separates nodes on each line
    - `frequency::Bool (=false)`: Whether the file has frequencies at the end of
                                  each line


Returns
-------
    - An array of paths, where each path is an array of integers, in which each
      node from the file has been assigned a unique number
    - The `MatrixNetwork` implied by the paths
    - The dictionary relating node strings to node numbers

Example
-------
~~~~
read_paths("data.txt"; separator=" ", frequency=true)
~~~~
"""
function read_paths(filepath::String; separator::String=",", frequency::Bool=false)
    println("Begin reading $filepath.")
    nodes = Dict() # Maps strings in the file to positive integers
    paths = Array{Array{Int64,1},1}()
    edges = Set{Tuple{Int64,Int64}}() # Edges observed in paths
    t1 = now()
    open(filepath) do file
        lines = [ln for ln in eachline(file)]
        for ln in lines
            separated = split(ln, separator)
            p = []
            for i = 1:(frequency ? length(separated)-1 : length(separated)) # Parse strings in line
                n = separated[i]
                if !haskey(nodes,n)
                    nodes[n] = length(nodes) + 1
                end
                if i < length(separated)-1 && !haskey(nodes,separated[i+1]) # Get following string
                    nodes[separated[i+1]] = length(nodes) + 1
                end
                push!(p, nodes[n])
                if i < length(separated)-1
                    push!(edges, (nodes[n], nodes[separated[i+1]])) # Add edge between current and next
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
    return paths, MatrixNetwork(collect(edges), length(nodes)), nodes
end
