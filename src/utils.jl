using DataStructures
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

function get_vertex_id(ex::HONData)
    vertex_id = Dict()
    let idx = 0
        for nvert in ex.nverts
            for i in range(idx+1, stop=idx+nvert)
                v = ex.simplices[i]
                if !haskey(vertex_id, v)
                    vertex_id[v] = length(vertex_id) + 1
                end
            end
            idx += nvert
        end
    end

    return vertex_id
end

"""Gather the edges that exist in the projected graph."""
function get_edge_list(ex::HONData, p::Float64=1, weighted::Bool=true)
    edges, vertices = Dict(), Set()
    vertex_id = get_vertex_id(ex)
    let idx = 0
        for nvert in ex.nverts
            for i in range(idx+1, stop=idx+nvert)
                for j in range(i+1, stop=idx+nvert)
                    v_i, v_j = ex.simplices[i], ex.simplices[j]
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

    n = length(vertex_id)
    edge_list = collect(edges)
    edge_list = [(e, w^p) for (e, w) in edge_list]

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

function postprocess_counters(k::Int64,
                              x::Dict,
                              compute_weight::Function,
                              kprime::Int64=-1)
    if kprime == -1
        k = min(k, length(x))
        arr = [(compute_weight(t), t) for (t, _) in zip(keys(x), values(x))]
        top_k = nlargest(k, arr)
    else
        kprime = min(kprime, length(x))
        k = min(k, kprime)
        arr = [(cnt, t) for (t, cnt) in zip(keys(x), values(x))]
        top_kprime = nlargest(kprime, arr)
        top_k = nlargest(k, [(compute_weight(t), t) for (_, t) in top_kprime])
    end

    return top_k

end

function get_edges_to_simplices(m::Int64, ex::HONData)
    appearances = SimpleGraph(m + length(ex.nverts))
    num_simplices = 0
    edge_id, rev_edge_id = Dict(), Dict()
    vertex_id = get_vertex_id(ex)
    let idx = 0
        for nvert in ex.nverts
            num_simplices += 1
            for i in range(idx+1, stop=idx+nvert)
                for j in range(i+1, stop=idx+nvert)
                    v_i, v_j = ex.simplices[i], ex.simplices[j]
                    u, v = -1, -1
                    if vertex_id[v_i] < vertex_id[v_j]
                        u, v = vertex_id[v_i], vertex_id[v_j]
                    else
                        u, v = vertex_id[v_j], vertex_id[v_i]
                    end
                    edge = (u, v)
                    if !haskey(edge_id, edge)
                        edge_id[edge] = length(edge_id) + 1
                        rev_edge_id[edge_id[edge]] = edge
                    end
                    add_edge!(appearances, edge_id[edge], m+num_simplices)
                end
            end
            idx += nvert
        end
    end

    return num_simplices, appearances, edge_id, rev_edge_id
end

function get_simplices_to_vertices(n::Int64,
                                   ex::HONData,
                                   vertex_id::Dict)
    g = SimpleGraph(n + length(ex.nverts))
    num_simplices = 0
    let idx = 0
        for nvert in ex.nverts
            num_simplices += 1
            for i in range(idx+1, stop=idx+nvert)
                v = vertex_id[ex.simplices[i]]
                add_edge!(g, num_simplices, length(ex.nverts)+v)
                add_edge!(g, length(ex.nverts)+v, num_simplices)
            end
            idx += nvert
        end
    end

    return g
end