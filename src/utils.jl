using LightGraphs
using ScHoLP
using SimpleWeightedGraphs

# Copied from ScHoLP-Tutorial/common.jl
"""Reads dataset and returns a HONData object."""
function read_txt_data(dataset::String)
    function read(filename::String)
        ret = Int64[]
        open(filename) do f
            for line in eachline(f)
                push!(ret, parse(Int64, line))
            end
        end
        return ret
    end
    return HONData(read("../data/$(dataset)/$(dataset)-simplices.txt"),
                   read("../data/$(dataset)/$(dataset)-nverts.txt"),
                   read("../data/$(dataset)/$(dataset)-times.txt"),
                   dataset)
end

"""Gather the edges that exist in the projected graph."""
function get_edge_list(ex::HONData, inverse::Bool=false, weighted::Bool=true)
    edges, vertices, vertex_id = Dict(), Set(), Dict()
    let idx = 0
        for nvert in ex.nverts
            for i in range(idx+1, stop=idx+nvert)
                for j in range(i+1, stop=idx+nvert)
                    v_i, v_j = ex.simplices[i], ex.simplices[j]
                    if !haskey(vertex_id, v_i)
                        vertex_id[v_i] = length(vertex_id) + 1
                    end
                    if !haskey(vertex_id, v_j)
                        vertex_id[v_j] = length(vertex_id) + 1
                    end
                    u, v = -1, -1
                    if vertex_id[v_i] < vertex_id[v_j]
                        u, v = vertex_id[v_i], vertex_id[v_j]
                    else
                        u, v = vertex_id[v_j], vertex_id[v_i]
                    end
                    edge = (u, v)
                    if !haskey(edges, edge)
                        edges[edge] = 1
                    elseif weighted
                        edges[edge] += 1
                    end
                    push!(vertices, u)
                    push!(vertices, v)
                end
            end
            idx += nvert
        end
    end

    n = length(vertices)
    edge_list = collect(edges)
    if inverse
        edge_list = [(e,1.0/w) for (e,w) in edge_list]
    end

    return n, edge_list
end

function get_adj_list(n::Int64, edge_list::Array)
    adj_list = SimpleWeightedGraph(n)
    for ((u, v), w) in edge_list
        add_edge!(adj_list, u, v, w)
        add_edge!(adj_list, v, u, w)
    end

    return adj_list
end

function get_adj_list_higher_deg(n::Int64, edge_list::Array)
    # Compute the degree of each vertex.
    degree = zeros(n)
    for ((u, v), _) in edge_list
        degree[u] += 1
        degree[v] += 1
    end

    # Construct an adjacency list representation of the graph where a vertex v
    # has an edge to a neighbor w if w has a higher degree than v. Checking if
    # edge(v, w) exists is a log(degree) operation in the following 
    # representation instead of constant, but it should suffice for our 
    # purposes.
    adj_list = SimpleWeightedGraph(n)
    for ((u, v), w) in edge_list
        if degree[u] > degree[v]
            add_edge!(adj_list, v, u, w)
        elseif degree[u] < degree[v]
            add_edge!(adj_list, u, v, w)
        elseif u < v
            add_edge!(adj_list, u, v, w)
        else
            add_edge!(adj_list, v, u, w)
        end
    end

    return adj_list
end