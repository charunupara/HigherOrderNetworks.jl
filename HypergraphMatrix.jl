"""
Represents a Hypergraph as an incidence matrix.

# Fields
- `ci::Vector{Int64}`: Column indices, correspond to nodes
- `ri::Vector{Int64}`: Row indices, correspond to hyperedges
- `vals::Vector{T}`: The nonzero values of the matrix
- `n::Int64`: The number of nodes
- `m::Int64`: The number of hyperedges
- `D::Vector{Int64}`: The degree sequence
- `K::Vector{Int64}`: The edge dimension sequence
"""
mutable struct HypergraphMatrix
    ci::Vector{Int64}
    ri::Vector{Int64}
    vals::Vector{Float64}
    n::Int64
    m::Int64
    D::Vector{Int64}
    K::Vector{Int64}
end

"""
Kernel hypergraph constructor
"""
function HypergraphMatrix(ci::Vector{Int64}, ri::Vector{Int64}, vals::Vector{Float64},
                          n::Int64, m::Int64)
    @assert size(ci) == size(ri) == size(vals)
    @assert maximum(ci) <= n && maximum(ri) <= m
    count(v, e) = sum([i == e ? 1 : 0 for i in v])

    D::Vector{Int64} = [count(ci, k) for k in 1:n]
    K::Vector{Int64} = [count(ri, k) for k in 1:m]
    return HypergraphMatrix(ci, ri, vals, n, m, D, K)
end

"""
Overloaded constructor, initializes `n` and `m` to the sizes of ci and ri
respectively. Requires fixing of last two arguments, as n and m are the number of _unique_ nodes and edges, where as nodes and edges
can appear multiple times in `ci` and `ri`.
"""
HypergraphMatrix(ci::Vector{Int64}, ri::Vector{Int64}, vals::Vector{Float64}) = HypergraphMatrix(ci, ri, vals, size(ci,1), size(ri,1))

"""
Overloaded constructor, initializes vals to ones
"""
HypergraphMatrix(ci::Vector{Int64}, ri::Vector{Int64}) = HypergraphMatrix(ci, ri, ones(size(ci)))

print(HypergraphMatrix([1,1,2,2,2,3,3,4,4], [1,3,1,2,3,1,2,1,3], ones(9), 4, 3))
