C(n,k) = factorial(big(n)) / (factorial(big(k)) * factorial(big(n-k))) # n choose k
vertex_floor(s::Te) where Te <: Real = Int(floor(s)) # Make sure a node (stub or vertex) is in vertex form, i.e. an integer
function deg_distr(m::MatrixNetwork)
    deg_dist = zeros(Int64, m.n)
    for i = 1:size(m.rp,1)-1
        deg_dist[i] += (m.rp[i+1] - m.rp[i])
    end

    for j in m.ci
        deg_dist[j] += 1
    end

    return deg_dist ./= 2
end
