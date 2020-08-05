using SparseArrays
using LinearAlgebra
using Random


""" 
`hyperkron_graph`
==================
Code to generate a HyperKron Graph using the grass-hopping algorithm.

The original code for this method can be found at:
https://www.cs.purdue.edu/homes/dgleich/codes/hyperkron/

The model definition can be found in the following paper:
https://arxiv.org/abs/1809.03488


Input
-----
- `K`: The initiator tensor
- `k`: The power for the tensor
- `motif`: optional parameter. You can define your own motif to add to the model. Default is undirected triangle.

Functions
---------
-`hyperkron_graph(K,k)`
-`hyperkron_graph(K,k;triangle,False)`

Output
------
- A matrix network type for the hyperkron graph
- A list of hyperedges

Example
--------
A,hedges = hyperkron_graph(kron_params(0.99, 0.2, 0.3, 0.05), 5)
hyperkron_graph(K,k; motif=:triangle, onlytris::Bool=false)
 
The code calls "hypergraph_to_edges_with_" 'motif' in order
to generate the graph to return. Consequently, you can implement
your own motifs via this routine.
"""
function hyperkron_graph(K,k;
    motif_to_graph::Function = hyperedges_to_graph_with_triangles, kwargs...)

  hedges = fast_hyperkron_edges(K,k; kwargs...)
  return motif_to_graph(size(K,1)^k, hedges), hedges
end


"""
`_next_region!`
===============
Increment to the next region of a Kronecker graph

Input
-----
-`cur`: The current region
-`m`: The total number of elements in the initiated tensor (ie 8 for a 2x2x2 initiator tensor). 

Output
------
Outputs the next region
"""
@inline function _next_region!(cur::AbstractArray{T,1}, m::Integer) where T <: Integer
  lastk = length(cur)
  for k=length(cur):-1:1
    lastk = k              # save the last state
    cur[k] += 1            # increment the last element
    if cur[k] == m+1       # if there is spill
      if k > 1             # if there is an array
        continue # recur on prefix
      else
        cur[k] = 0
      end
    end
    break
  end
  for k=lastk+1:length(cur)
    cur[k] = cur[k-1]
  end
  cur
end

"""
`_first_region`
===============
Gives the first region.

Input
-----
-`a`: The vectorized initiator tensor
-`k`: The power of the HyperKron model
Output
------
An array of all ones, representing the first region.
"""
function _first_region(a::AbstractArray,k::Integer)
  return ones(UInt32,k)
end

"""
`num_multiset_permutations`
===========================
Returns the number of multiset permutations of elements in the given array.

Input
-----
-`ms`: A sorted array of numbers, which represent the values in the multiset.
-`start_pos`: Default is 1. 

Output
------
The number of multiset permutations.
"""
function num_multiset_permutations(ms::AbstractArray{T,1}, start_pos::Integer = 1) where T<:Integer
  nperm = 1
  count = 1
  countfac = 1
  @inbounds for i in start_pos+1:length(ms)
    nperm *= (i-start_pos+1)
    if ms[i] == ms[i-1]
      count += 1
      countfac *= count
    else
      nperm = div(nperm,countfac)
      # new element
      count = 1
      countfac = 1
    end
  end
  nperm = div(nperm,countfac)
end

"""
`_next_distinct_character`
==========================
Find the next distinct character after a[k]
Returns the index where it occurs or length(a)+1

Input
-----
-`a`: A region in the HyperKron tensor
-`k`: The current character (between 1 and m)
-`m`: Length of a. 

Output
------
The index where the next character occurs. 
"""
@inline function _next_distinct_character(a, k::Integer, m::Integer)
  # check for easy out
  if k>=m
    return m+1
  else
    ak = a[k]
    if ak != a[k+1] # many cases
      return k+1
    elseif ak == a[m] # some cases
      return m+1
    else
      # we have to search, so there are probably many dups
      # since we assume sorted, we can bisect
      lo = 1
      hi = m-k # we need to search [k+lo:k+hi]
      while hi - lo > 2
        mid = (hi + lo) >> 1
        if ak == a[k+mid]
          lo = mid
        else
          hi = mid
        end
      end
      ak != a[k+lo+1] && return k+lo+1
      return k+lo+2
    end
  end
end

"""
`unrank5!`
==========
This version doesn't use recursion, and unranks the sorted set ms in place
It also gets the multiplicity of each element to update nperm without
recomputing it from scratch, it also can take in the total number
of permutations and just update as we go.
This version tries to avoid integer division by storing nperm as a float

Input
-----
-`ms`: Multiset	
-`n`: Current rank
-`nperm`: The number of multiset permutations of the region
-`m`: Length of ms

Output
------
Nothing, it modifies ms. 
"""
@inline function unrank5!(ms::AbstractArray{T,1}, n::Integer, nperm, m) where T<:Integer
  multcur = 0
  imultcur = 1.0
  multnext = 0
  npermf = float(nperm)
  @inbounds for k=1:m # for each location in the output
    if n == 1
      break
    end

    if multcur == 0 # update
      multcur = _next_distinct_character(ms, k, m) - k # multiplicity of the current element
      imultcur = 1.0/multcur
    end

    # update nperm for removing multcur
    npermf *= multcur
    npermf /= (m-k+1) # shrink the length
    nperm = round(typeof(nperm),npermf)

    #=
    if nperm != num_multiset_permutations(ms, k+1)
      @show 0, ms, k, num_multiset_permutations(ms, k+1), nperm
    end
    @assert nperm == num_multiset_permutations(ms, k+1)
    =#

    if n <= nperm # this means we start with ms[k]
      multcur -= 1
      if multcur != 0
        imultcur = 1.0/multcur
      end
      continue
    else # we need to search
      n -= nperm

      if k+multcur <= m
        if ms[k] == ms[k+multcur]
          #@show ms, k, _next_distinct_character(ms,k,m), multcur
          @assert ms[k] != ms[k+multcur]
        end
      end

      pos = k+multcur
      while pos <= m
        # we are going to swap ms[k] and ms[pos]
        # but to update, we need to know the multiplicity of ms[pos]

        multnext = _next_distinct_character(ms,pos,m) - pos

        # so the element ms[pos] has multiplicity multnext
        # and element ms[k] has multiplicity multcur-1

        ms[pos], ms[k] = ms[k], ms[pos] # swap leading positions
        npermf *= multnext
        npermf *= imultcur
        nperm = round(typeof(nperm),npermf)

        #=
        if nperm != num_multiset_permutations(ms, k+1)
          @show 1, ms, k, num_multiset_permutations(ms, k+1), nperm
        end
        @assert nperm == num_multiset_permutations(ms, k+1)
        =#

        if n <= nperm
          break
        else
          n -= nperm
          pos += multnext
          multcur = multnext
          imultcur = 1.0/multcur
        end
      end

      multcur = 0 # reset the search
    end
  end
  return ms
end
unrank5!(ms, n) = unrank5!(ms, n, num_multiset_permutations(ms), length(ms))


"""
`_morton_decode3`
=================
Turns an index from a k-dimensional multiplication table into the corresponding index in the 3-dimensional tensor product.

Input
-----
-`mind`: Length k (where k is the number of kronecker products to take) array of indices in the k-dimensional multiplication table
-`n`: The side length of the initiator tensor
-`sz`: The full size of the initiator tensor

Output
------
Length 3 array, representing the corresponding indices in the tensor product. 

Example
-------
_morton_decode3([8,8,8,8],2,(2,2,2))
"""
@inline function _morton_decode3(mind::AbstractArray{T,1}, n, sz) where T<:Integer
  row = 0
  col = 0
  tub = 0
  base = 1
  for val in mind
    # CU: ind2sub() was deprecated in Julia 0.7, changed to a new way instead
    s1,s2,s3 = Base._ind2sub(sz, val) #basically backward map where base keeps track of the block size n^
    row += base*(s1-1)
    col += base*(s2-1)
    tub += base*(s3-1)
    base *= n
  end
  return (row+1, col+1, tub+1)
end


#_randgeo(p) = ceil(Int,-randexp() / log1p(-p)) # jsut a copy of Distributions.jl
#_randgeo(p) = ceil(Int,-randexp() / log1p(-p)) # jsut a copy of Distributions.jl
#_randgeo_true(p) = ceil(Int,-randexp() / log1p(-p)) # jsut a copy of Distributions.jl
#_randgeo(p) = _randgeo_true(max(0.01*p, p + 0.0*p*randn()))
_randgeo(p) = ceil(Int,-randexp() / log1p(-p))

"""
`_generate_region!`
==================
Given a region of the tensor, do the grass-hopping procedure to generate hyperedges.

Input
-----
-`ei`: First indices of hyperedges
-`ej`: Second indices of hyperedges
-`ek`: Third indices of hyperedges
-`r`: Region 
-`e`: Edge in Region
-`p`: The probability of edges within this particular region.
-`small_n`: The length of one dimension of the initiator matrix.
-`sz`: The full dimensions of the initiator tensor.
-`onlytris`: Boolean
-`onlysym`: Boolean

Output
------
Nothing, modifies input.
"""
function _generate_region!(
    ei::Array{Ta}, ej::Array{Ta}, ek::Array{Ta},
    r::AbstractArray{Tb,1},
    e::AbstractArray{Tb,1},
    p, small_n::Integer, sz,
    onlytris::Bool, onlysym::Bool) where {Ta, Tb <:Integer}

    N = num_multiset_permutations(r) # total area
    m = length(r)
    i = 0

    # ugh, do we need this? should make these tolerances.
    if N*p < 1.0e-12 || p < 1.0e-15
      return
    end

    gap = _randgeo(p)

    while i+gap <= N
      i += gap
      copy!(e, r)
      #@show e, i, N, m
      mult_edge = unrank5!(e, i, N, m) # generate the ith permutation
      #src, dst = _morton_decode(mult_edge, small_n)
      src,dst,tub = _morton_decode3(mult_edge, small_n, sz)
      # only add the edge if it's in the lower region of the matrix.
      # This means we are assuming a symmetric generation.
      if (onlysym && src <= dst <= tub) || !onlysym
        if (onlytris && (src != dst != tub != src)) || !onlytris
          push!(ei, src)
          push!(ej, dst)
          push!(ek, tub)
        end
      end
      gap = _randgeo(p)
    end
end

"""
`fast_hyperkron_edges`
======================
Returns the hyperedges of the Kronecker model

Input
-----
-`K`: Initial tensor for hyperkron model
-`k`: Power for the hyperkron model

Output
------
ei,ej,ek: Three arrays which hold the indices for hyperedges
"""
function fast_hyperkron_edges(K,k;
  noise::Float64=0.0,onlytris::Bool=false,onlysym::Bool=true)
  n = size(K,1)
  sz = size(K)
  v = vec(K) # vectorized

  vn = v/sum(v)
  V = collect(v for _ in 1:k) # create a set of parameters for each level
  #V = collect(max.(0.0,min.(2.0,sum(v).*normalize!(vec(kron_params(vn[1],vn[2]+noise*(2*rand()-1.0),vn[4]+noise*(2*rand()-1.0),vn[8])),1))) for _ in 1:k )
  #@show V

  # Dense = adjacency matrix, every spot is filled in, even if there's a 0. Very large graph = expensive to store
  # Sparse = Store a list of non-zeroes edges. For i and j, there will be two one-dimension array

  # hold edges
  ei = zeros(Int,0)
  ej = similar(ei)
  ek = similar(ei)

  nregions = 0
  # v is the vectoriezd tensor
  region = _first_region(v, k)

  edge_in_region = similar(region) # allocated once to avoid lots of allocs
  while region[1] != 0
    nregions += 1
    # compute the probability of erdos-reyes regions
    p = 1.0
    level = 1
    for j in region
      p *= V[level][j]
      level += 1
    end
    _generate_region!(ei, ej, ek, region, edge_in_region, p, n, sz, onlytris, onlysym)
    _next_region!(region, length(v))
  end
  return ei, ej, ek
end

"""
`kron_params`
=============
Turn 4 parameters into a 2x2x2 symmetric tensor representation for the hyperkron_graph function

Input
-----
-`a`: probability between 0 and 1
-`b`: probability between 0 and 1
-`c`: probability between 0 and 1
-`d`: probability between 0 and 1

Output
Correctly formatted tensor parameters

"""
function kron_params(a,b,c,d)
  K = zeros(2,2,2)
  K[1,1,1] = a
  K[1,1,2] = b; K[1,2,1] = b;   K[2,1,1] = b
  K[1,2,2] = c; K[2,1,2] = c;   K[2,2,1] = c
  K[2,2,2] = d
  return K
end

"""
`hyperedges_to_graph_with_triangles`
====================================
Turns a set of hyperedges into a graph.

Input
-----
-`n`: The size of the graph
-`hedges`: The hyperedges

Output
------
A MatrixNetwork graph based on the given hyperedges.
"""
function hyperedges_to_graph_with_triangles(n::Int, hedges)
  ei = zeros(Int64,0)
  ej = zeros(Int64,0)
  for i=1:length(hedges[1])
    push!(ei, hedges[1][i]); push!(ej, hedges[2][i]) # (i -> j)
    push!(ei, hedges[1][i]); push!(ej, hedges[3][i] ) # (i -> k)
    push!(ei, hedges[2][i]); push!(ej, hedges[3][i] ) # (j -> k)
  end
  A = sparse(ei,ej,1.0,n,n)
  A = max.(A,A') # symmetrize
  A = A - spdiagm(0 => diag(A)) #remove diagonal entries
  fill!(A.nzval, 1.0) # set all non-zeros
  return A
end

