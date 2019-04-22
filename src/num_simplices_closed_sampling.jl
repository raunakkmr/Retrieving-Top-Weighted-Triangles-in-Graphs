include("./utils.jl")

using DataStructures
using DelimitedFiles
using Distributions
using ScHoLP
using SimpleWeightedGraphs
using SparseArrays
using StatsBase

"""Construct samplers for sampling neighbors of vertices. Only constructs vertices for the first idx vertices of the graph."""
function nbr_samplers(graph::SimpleGraph, idx::Int64)
    degrees = [length(neighbors(graph, v)) for v in vertices(graph)[1:idx]]
    wts = [ones(degree) for degree in degrees]
    probs = []
    for (i, (w, d)) in enumerate(zip(wts, degrees))
        if d == 0
            normalized_prob = Array([1.0])
        else
            normalized_prob = Array(w) ./ d
        end
        push!(probs, Categorical(normalized_prob))
    end
    graph_samplers = [sampler(prob) for prob in probs]

    return graph_samplers
end

"""Construct samplers for sampling edges from the bipartite graph between
edges and simplices, and sampling neighbors of vertices from the bipartite
graphs between edges and simplices, and simplices and vertices."""
function construct_samplers(graphs::Array,
                            m::Int64,
                            num_simplices::Int64)
    ledges, lgraph, rgraph = graphs
    ledge_prob = Categorical(ones(length(ledges)) / length(ledges))
    ledge_sampler = sampler(ledge_prob)

    # Since we only sample neighbors of edges in the first graph and the
    # neighbors of simplices in the second graph, only construct samplers for
    # those vertices.
    lgraph_samplers = nbr_samplers(lgraph, m)
    rgraph_samplers = nbr_samplers(rgraph, num_simplices)

    return [ledge_sampler, lgraph_samplers, rgraph_samplers]
end

"""Sample triangles / diamonds and return the top k closed triangles based on the number of simplices they appear in."""
function compute_weighted_triangles(n::Int64,
                                    m::Int64,
                                    kprime::Int64,
                                    k::Int64,
                                    num_simplices::Int64,
                                    graphs::Array,
                                    samplers::Array,
                                    ids::Array,
                                    ex::HONData)
    s = 1000000  # Number of samples.
    x = Dict()  # Counters.
    num_diamonds = 0

    # ledges is the set of edges between edges in the original graph and
    # simplices. There is an edge between e and s if edge e appears in simplex
    # s. lgraph is the adjacency list representation of this graph.
    # rgraph is the adjacency list for the graph between simplices and
    # vertices. There is an edge between s and v if vertex v appears in simplex
    # s.
    ledges, lgraph, rgraph = graphs
    ledge_sampler, lgraph_samplers, rgraph_samplers = samplers
    edge_id, rev_edge_id = ids

    # Sample (edge, simplex). Sample (simplex, vertex). If this forms a
    # triangle, then sample (edge, simplex'). If (simplex', vertex) exists then
    # increment the counter for the triangles corresponding to (edge, vertex).
    for l in range(1, stop=s)

        # a is edge id, c is m + simplex number
        a, c = ledges[rand(ledge_sampler)]
        # Need 3 different vertices in simplex to have a valid triangle.
        if length(neighbors(rgraph, c-m)) < 3
            continue
        end
        # b is num_simplices + vertex id
        b = neighbors(rgraph, c-m)[rand(rgraph_samplers[c-m])]
        if b-num_simplices == rev_edge_id[a][1] || b-num_simplices == rev_edge_id[a][2]
            continue
        end
        # cprime is m + simplex number
        cprime = neighbors(lgraph, a)[rand(lgraph_samplers[a])]
        if !has_edge(rgraph, cprime-m, b)
            continue
        end
        num_diamonds += 1
        t = Tuple(sort([rev_edge_id[a][1], rev_edge_id[a][2], b-num_simplices]))
        if !haskey(x, t)
            x[t] = 0
        end
        x[t] += 1
    end

    # Postprocessing. Compute the weight of the triangles corresponding to the
    # top kprime counters, and return the top k triangles from these.
    function compute_weight((a, b, c))
        ar = neighbors(rgraph, a+num_simplices)
        br = neighbors(rgraph, b+num_simplices)
        cr = neighbors(rgraph, c+num_simplices)
        as = sparsevec(ar, ones(length(ar)), num_simplices)
        bs = sparsevec(br, ones(length(br)), num_simplices)
        cs = sparsevec(cr, ones(length(cr)), num_simplices)
        weight = sum(as .* bs .* cs)

        return weight
    end
    top_k = postprocess_counters(k, kprime, x, compute_weight)

    return top_k
end

"""Construct the graph, and compute and return the triangles in sorted order."""
function construct_and_compute(n::Int64,
                               kprime::Int64,
                               k::Int64,
                               edge_list::Array,
                               ex::HONData)
    vertex_id = get_vertex_id(ex)
    m = length(edge_list)
    num_simplices, lgraph, edge_id, rev_edge_id = get_edges_to_simplices(m, ex)
    ledges = [(src(e), dst(e)) for e in edges(lgraph)]
    rgraph = get_simplices_to_vertices(n, ex, vertex_id)
    graphs = [ledges, lgraph, rgraph]
    samplers = construct_samplers(graphs, m, num_simplices)
    ids = [edge_id, rev_edge_id]
    triangles = compute_weighted_triangles(n, m, kprime, k, num_simplices, graphs, samplers, ids, ex)
end

function main()
    dataset_name, kprime, k = ARGS[1], parse(Int64, ARGS[2]), parse(Int64, ARGS[3])
    output_file = "../output/num_simplices_closed_sampling_$dataset_name.txt"

    # Load data in the form of simplices from the ScHoLP package.
    ex = read_txt_data(dataset_name)

    n, edge_list = get_edge_list(ex, 1.0)
    time = @elapsed triangles = construct_and_compute(n, kprime, k, edge_list, ex)
    writedlm(output_file, triangles)
    open(output_file, "a") do f
        write(f, "Time: $time\n")
        write(f, "k: $k\n")
        write(f, "kprime: $kprime")
    end
end

main()