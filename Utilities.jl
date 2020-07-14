C(n,k) = factorial(big(n)) / (factorial(big(k)) * factorial(big(n-k))) # n choose k
vertex_floor(s::Te) where Te <: Real = Int(floor(s)) # Make sure a node (stub or vertex) is in vertex form, i.e. an integer

function deg_distr(A::MatrixNetwork)
    deg_dist = zeros(Int64, A.n)
    for i = 1:size(A.rp,1)-1
        deg_dist[i] += (A.rp[i+1] - A.rp[i])
    end

    for j in A.ci
        deg_dist[j] += 1
    end

    return deg_dist
end

function outweight_distr(A::MatrixNetwork) # For unweighted, equal to outdegree distribution
    deg_dist = zeros(Int64, A.n)
    if eltype(A.vals) <: Real
        for i = 1:size(A.rp,1)-1
            deg_dist[i] += sum(A.vals[A.rp[i]:A.rp[i+1]-1])
        end
    else
        for i = 1:size(A.rp,1)-1
            deg_dist[i] += (A.rp[i+1] - A.rp[i])
        end
    end
    return deg_dist
end

function inweight_distr(A::MatrixNetwork)
    deg_dist = zeros(Int64, A.n)
    if eltype(A.vals) <: Real
        for j = 1:size(A.ci,1)
            deg_dist[A.ci[j]] += A.vals[j]
        end
    else
        for j in A.ci
            deg_dist[j] += 1
        end
    end
    return deg_dist
end

function get_edges(A::MatrixNetwork)
    return [[i,j] for i = 1:size(A.rp,1)-1 for j in A.ci[A.rp[i]:A.rp[i+1]-1]]
end

function total_edge_weights(A::MatrixNetwork)
    weights = Dict{Vector{Int64}, eltype(A.vals)}()
    i = 1
    for j = 1:size(A.ci,1)
        if j >= A.rp[i+1]
            i += 1
         end

        if !haskey(weights, [i,A.ci[j]])
            weights[[i,A.ci[j]]] = 0
        end
        weights[[i,A.ci[j]]] += A.vals[j]
    end
    return weights
end
