"""
'triangle generalized_preferential_attachment_graph'
===========================================
Generate an instance of a triangle generalized preferential attachment graph which
follows the Eikmeier, Gleich description (https://arxiv.org/pdf/1904.12989.pdf).
This is an undirected graph that is generated as follows:
- Start with an empty graph
- Add n vertices where at each time step one of three events occurs:
A new node is added with probability p, along with two edges to existing nodes;
Three new nodes with two edges between them is added with probability 1 - p
Every m steps, coalesce the nodes added in time step 1

Functions
---------
The following functions are synonyms
- 'triangle_generalized_preferential_attachment_graph_v2'
- 'tgpa_graph_v2'
and
- 'triangle_generalized_preferential_attachment_edges_v2!'
- 'tgpa_edges_v2!'
The computational functions are
- 'tgpa_graph_v2(n,p,k0,m)' Generate a TGPA graph with n total nodes.
    This returns a MatrixNetwork type. Note that this process allows self loops.
The edge functions are
-   'tgpa_edges_v2!(n,p,m,edges,n0)' Add new edges to an existing set, by taking
    n0 time steps. Edges are added in one of two ways: From a new node to
    an existing node with probability p, between three new nodes with probability 1-p

Input
-----
- 'n': the number of nodes in the final graph.
- 'p': The probability of a node event, p must be a constant.
- 'k0': The number of nodes in the starting clique.
- 'm': The frequency of combining nodes
- 'edges': A list of edges to be manipulated in the process of generating
          new edges.
Output
------
- A matrix network type for the triangle generalized preferential attachment graph.
- 'edges': An updated list of edges.
Example:
triangle_generalized_preferential_attachment_graph_v2(100,1/3,1/2)
"""
:triangle_generalized_preferential_attachment_graph_v2, :tgpa_graph_v2,
:triangle_generalized_preferential_attachment_edges_v2!, :tgpa_edges_v2!

using Printf
using MatrixNetworks
using Distributions

function triangle_generalized_preferential_attachment_graph_v2(
    n::Int,p::Float64,k0::Int,m::Int)
    k0 >= 0 || throw(ArgumentError(@sprintf("k0=%i must be >= 0",n)))
    n >= k0 || throw(ArgumentError(@sprintf("n=%i must be >= k0=%i",n,k0)))
    0<=p<=1 || throw(ArgumentError(@sprintf("p=%0.3f must be between 0 and 1",p)))
    m >= 1 || throw(ArgumentError(@sprintf("m=%i must be >= 1",m)))
    m <= max(n,1) || throw(ArgumentError(@sprintf("m=%i must be <= n=%i",m,n)))
    edges = Vector{Tuple{Int,Int}}()
    # add the clique
    for i=1:k0
        for j=1:i-1
            push!(edges, (i,j))
            push!(edges, (j,i))
        end
    end
    return MatrixNetwork(triangle_generalized_preferential_attachment_edges_v2!(n,p,m,edges,k0),n)
end

function _edges_to_neighbors(edges, k0)
    neighbors = Vector{Vector{Int}}(undef, k0)
    for i=1:k0
        neighbors[i] = Vector{Int}()
    end
    for e in edges
        push!(neighbors[e[1]], e[2])
    end
    return neighbors
end


function triangle_generalized_preferential_attachment_edges_v2!(
    n::Int,p::Float64,m::Int,edges::Vector{Tuple{Int,Int}},k0::Int)
    i = k0 #number of nodes added

    # for this code, we need to have random neighbor access, so we'll keep
    neighbors = _edges_to_neighbors(edges, k0)

    while i < n
        #determine how many component events we should have
        #this says for every m steps (every m new nodes) we should have this many new components
        num_component_events = rand(Binomial(m,1-p),1)[1]

        #new component events
        for j = 1:num_component_events
            if i+3 <= n #only allow this step if there is room for three more nodes
                push!(edges, (i+1, i+2)) #add wedge to edge list
                push!(edges, (i+2, i+1))
                push!(edges, (i+1, i+3))
                push!(edges, (i+3, i+1))
                push!(neighbors,Vector{Int}()) ##update neighbor list
                push!(neighbors,Vector{Int}())
                push!(neighbors,Vector{Int}())
                push!(neighbors[i+1],i+2)
                push!(neighbors[i+1],i+3)
                push!(neighbors[i+2],i+1)
                push!(neighbors[i+3],i+1)
                i = i+3
            end
        end
        #new node events
        if i+1 <=n
            push!(neighbors,Vector{Int}())
            for j = 1:m
                if length(edges) > 0
                    v = rand(edges)[1]
                    v_neighbors = neighbors[v]
                else
                    v = i+1
                    v_neighbors = [];
                end
                #edge between new node and v
                push!(edges, (i+1, v))
                push!(edges, (v, i+1))
                #update neighbor vector
                push!(neighbors[i+1],v)
                push!(neighbors[v],i+1)

                #add another edge between new node and a neighbor of v
                if length(v_neighbors) > 0
                    v2 = rand(v_neighbors)
                    push!(edges, (i+1, v2))
                    push!(edges, (v2, i+1))
                    push!(neighbors[i+1],v2)
                    push!(neighbors[v2],i+1)
                end
            end
            i = i+1;
        end

    end
    return edges
end

tgpa_graph_v2 = triangle_generalized_preferential_attachment_graph_v2
tgpa_edges_v2! = triangle_generalized_preferential_attachment_edges_v2!
