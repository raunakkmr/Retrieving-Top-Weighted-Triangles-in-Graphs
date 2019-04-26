include("./utils.jl")

using ArgParse
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
function compute_weighted_triangles(s::Int64,
                                    n::Int64,
                                    k::Int64,
                                    edge_list::Array,
                                    adj_list::SimpleWeightedGraph,
                                    edge_sampler::Distributions.AliasTable)
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

    # Postprocessing. Compute the weight of the sampled triangles and return
    # the top k from these.
    function compute_weight((a,b,c))
        return adj_list.weights[a,b] + adj_list.weights[b,c] + adj_list.weights[a,c]
    end
    top_k = postprocess_counters(k, x, compute_weight)

    return top_k
end

"""Construct the graph, and compute and return the triangles."""
function construct_and_compute(s::Int64,
                               n::Int64,
                               k::Int64,
                               edge_list::Array)
    adj_list = get_adj_list(n, edge_list)
    edge_sampler = construct_samplers(n, edge_list)
    triangles = compute_weighted_triangles(s, n, k, edge_list, adj_list, edge_sampler)

    return triangles
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--dataset", "-d"
            help = "dataset name"
            arg_type = String
            required = true
        "-k"
            help = "parameter k for returning top-k triangles, default: 25"
            arg_type = Int64
            default = 25
        "--samples", "-s"
            help = "number of samples, default: 100000"
            arg_type = Int64
            default = 100000
        "--uniform", "-u"
            help = "ignore weights and use uniform sampling"
            action = :store_true
        "-p"
            help = "parameter p for p mean for weights, default: 1.0"
            arg_type = Float64
            default = 1.0
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    dataset_name, k, p = parsed_args["dataset"], parsed_args["k"], parsed_args["p"]
    s = parsed_args["samples"]
    uniform = parsed_args["uniform"]
    if uniform
        output_file = "../output/mean_$(p)_sampling_uniform_$dataset_name.txt"
    else
        output_file = "../output/mean_$(p)_sampling_$dataset_name.txt"
    end

    # Load data in the form of simplices from the ScHoLP package.
    ex = read_txt_data(dataset_name)

    n, edge_list = get_edge_list(ex, p)
    if uniform
        for i in range(1, stop=length(edge_list))
            e, _ = edge_list[i]
            edge_list[i] = (e, 1.0)
        end
    end
    time = @elapsed triangles = construct_and_compute(s, n, k, edge_list)
    writedlm(output_file, triangles)
    open(output_file, "a") do f
        write(f, "Time: $time\n")
        write(f, "k: $k\n")
        write(f, "p: $p\n")
        write(f, "s: $s")
    end
end

main()