struct k_model
    k::Int64 # Memory of model
    P::Dict{Vector{Int64}, Float64} # Associates paths to probabilities
end

function k_model_k(k::Int64, P::Dict{Vector{Int64}, Float64})
    @assert k >= 0
    return k_model(k, P)
end

function optimal_k_model(S::Vector{Vector{Int64}}, k::Int64)
    P = Dict{Vector{Int64}, Float64}()
    p_counts = Dict{Vector{Int64}, Int64}()
    w_counts = Dict{Vector{Int64}, Int64}()
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
        P[p] = p_counts[p] / w_counts[p[1:k]]
    end
    return k_model_k(k, P)
end

function generate_paths(AL::Vector{Vector{Int64}}, k::Int64)
    if k == 1
        return [[i, j] for i = 1:size(AL,1) for j in AL[i]]
    else
        k_paths = Vector{Vector{Int64}}()
        for p in copy(generate_paths(AL, k - 1))
            for n in AL[p[k]]
                p_c = copy(p)
                push!(p_c, n)
                push!(k_paths, p_c)
            end
        end
        return k_paths
    end
end

function optimal_multiorder(S::Vector{Vector{Int64}}, K::Int64)
    models = [optimal_k_model(S, k) for k = 0:K] # models[i] = optimal_(i-1)_model
    P = Dict{Vector{Int64}, Float64}()

    S_K = [p[i:i+K] for p in S for i = 1:size(p,1)-K]

    for p in S_K
        pr = 1.0
        for k = 1:min(size(p,1), K-1)
            pr *= models[k].P[p[1:k]]
        end

        for i = K:size(p,1)
            pr *= models[K].P[p[i-K+1:i]]
        end
        P[p] = pr
    end

    return k_model(K, P)
end
AL = [[2], [4,3], [1], [1,2]]
paths = Vector{Vector{Int64}}()
for k = 1:11
    append!(paths, generate_paths(AL, k))
end
#=paths = [[2,4], [2,3], [4,1], [4,2], [1,2], [2,3,1], [1,2,4],
         [4,1,2], [2,4,2], [3,1,2], [4,2,4], [2,4,1], [1,2,3],
         [2,4,2,4], [4,1,2,4], [1,2,3,1], [1,2,4,2,4], [4,2,4,2,4],
         [3,1,2,4,2,4], [2,4,2,4,2,4]]=#
println(optimal_multiorder(paths, 1).P) # Numbers slightly off from paper...
