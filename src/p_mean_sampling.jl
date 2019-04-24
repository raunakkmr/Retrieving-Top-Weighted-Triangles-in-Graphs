include("./utils.jl")

using DataStructures
using DelimitedFiles
using Distributions
using ScHoLP
using SimpleWeightedGraphs
using StatsBase

"""Construct a sampler to sample edges."""
function construct_samplers(n::Int64, edge_list::Array)
    edge_weights = [w for (e,w) in edge_list]
    Z = sum(edge_weights)
    edge_prob = Categorical(edge_weights ./ Z)
    edge_sampler = sampler(edge_prob)

    return edge_sampler
end

"""Sample triangles and return the top k triangles based on p-mean of their
edge weights."""
function compute_weighted_triangles(n::Int64,
                                    kprime::Int64,
                                    k::Int64,
                                    edge_list::Array,
                                    adj_list::SimpleWeightedGraph,
                                    edge_sampler::Distributions.AliasTable)
    s = 1000000  # Number of samples.
    edge_x, x = Dict(), Dict()  # Counters.
    num_triangles = 0

    # Sample l edges and maintain a counter for each edge.
    for l in range(1, stop=s)

        (a,b), _ = edge_list[rand(edge_sampler)]
        u, v = min(a,b), max(a,b)
        if !haskey(edge_x, (u,v))
            edge_x[(u,v)] = 0
        end
        edge_x[(u,v)] += 1
    end

    # For each sampled edge, iterate through the neighbors of its endpoints.
    # For each triangle formed by the edge and the neighbor increment its
    # counter by the counter of the edge.
    for (e, cnt) in edge_x
        u, v = e
        for w in neighbors(adj_list, u)
            if w == v || !has_edge(adj_list, v, w)
                continue
            end
            t = Tuple(sort([u, v, w]))
            if !haskey(x, t)
                x[t] = 0
            end
            x[t] += edge_x[e]
        end
        for w in neighbors(adj_list, v)
            if w == u || !has_edge(adj_list, u, w)
                continue
            end
            t = Tuple(sort([u, v, w]))
            if !haskey(x, t)
                x[t] = 0
            end
            x[t] += edge_x[e]
        end
    end

    # Postprocessing. Compute the weight of the triangles corresponding to the top kprime counters, and return the top k triangles from these.
    function compute_weight((a,b,c))
        return adj_list.weights[a,b] + adj_list.weights[b,c] + adj_list.weights[a,c]
    end
    top_k = postprocess_counters(k, kprime, x, compute_weight)

    return top_k
end

"""Construct the graph, and compute and return the triangles."""
function construct_and_compute(n::Int64,
                               kprime::Int64,
                               k::Int64,
                               edge_list::Array)
    adj_list = get_adj_list(n, edge_list)
    edge_sampler = construct_samplers(n, edge_list)
    triangles = compute_weighted_triangles(n, kprime, k, edge_list, adj_list, edge_sampler)

    return triangles
end

function main()
    dataset_name, kprime, k, p = ARGS[1], parse(Int64, ARGS[2]), parse(Int64, 
        ARGS[3]), parse(Float64, ARGS[4])
    output_file = "../output/mean_$(p)_sampling_$dataset_name.txt"

    # Load data in the form of simplices from the ScHoLP package.
    ex = read_txt_data(dataset_name)

    n, edge_list = get_edge_list(ex, p)
    time = @elapsed triangles = construct_and_compute(n, kprime, k, edge_list)
    writedlm(output_file, triangles)
    open(output_file, "a") do f
        write(f, "Time: $time\n")
        write(f, "k: $k\n")
        write(f, "kprime: $kprime\n")
        write(f, "p: $p\n")
    end
end

main()