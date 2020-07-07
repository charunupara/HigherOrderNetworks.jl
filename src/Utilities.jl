C(n,k) = factorial(big(n)) / (factorial(big(k)) * factorial(big(n-k))) # n choose k
vertex_floor(s::Te) where Te <: Real = Int(floor(s)) # Make sure a node (stub or vertex) is in vertex form, i.e. an integer
