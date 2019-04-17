include("./utils.jl")

using DataStructures
using DelimitedFiles
using Distributions
using LightGraphs
using ScHoLP
using SimpleWeightedGraphs
using StatsBase

function construct_samplers(n::Int64, edge_list::Array, adj_list::SimpleWeightedGraph)
    # Construct a sampler for sampling edge weights.
    # edge_weights = [w for (_,w) in edge_list]
    edge_weights = []
    for ((a,b),w) in edge_list
        a_nbrs = neighbors(adj_list, a)
        a_nbrs_weights = adj_list.weights[a, filter(e->e!=b, a_nbrs)]
        b_nbrs = neighbors(adj_list, b)
        b_nbrs_weights = adj_list.weights[b, filter(e->e!=a, b_nbrs)]
        push!(edge_weights, w * sum(a_nbrs_weights) * sum(b_nbrs_weights))
    end
    Z = sum(edge_weights)
    edge_prob = Categorical(edge_weights ./ Z)
    edge_sampler = sampler(edge_prob)

    # Weights of the neighbors of every vertex.
    nbr_weights = [adj_list.weights[v, neighbors(adj_list, v)] for v in vertices(adj_list)]

    return edge_sampler, nbr_weights, Z
end

function compute_weighted_triangles(n::Int64,
                                    kprime::Int64,
                                    k::Int64,
                                    Z::Float64,
                                    edge_list::Array,
                                    adj_list::SimpleWeightedGraph,
                                    edge_sampler::Distributions.AliasTable,
                                    nbr_weights::Array)
    s = 1000000  # Number of samples.
    x = Dict()  # Counters.
    num_triangles = 0

    for l in range(1, stop=s)

        (a,b), _ = edge_list[rand(edge_sampler)]

        if length(neighbors(adj_list, a)) == 1 || length(neighbors(adj_list, b)) == 1
            continue
        end

        a_nbrs = neighbors(adj_list, a)
        a_nbrs_weights = pweights(adj_list.weights[a, filter(e->e!=b, a_nbrs)])
        c = sample(filter(e->e!=b, a_nbrs), a_nbrs_weights)

        b_nbrs = neighbors(adj_list, b)
        b_nbrs_weights = pweights(adj_list.weights[b, filter(e->e!=a, b_nbrs)])
        cprime = sample(filter(e->e!=a, b_nbrs), b_nbrs_weights)

        if c == cprime
            num_triangles += 1
            t = Tuple(sort([a, b, c]))
            if !haskey(x, t)
                x[t] = 0
            end
            x[t] += 1
        end
    end

    # Postprocessing.
    kprime = min(kprime, length(x))
    k = min(k, kprime)
    function compute_weight((a,b,c))
        return adj_list.weights[a,b] * adj_list.weights[b,c] * adj_list.weights[a,c]
    end
    x_array = [(count,triangle) for (triangle,count) in zip(keys(x), values(x))]
    sorted_x = sort!(x_array, by = x -> x[1], rev=true)
    top_kprime = [(compute_weight(triangle), triangle) for (_,triangle) in sorted_x[1:kprime]]
    top_k = nlargest(k, top_kprime)

    return top_k
end

function construct_and_compute(n::Int64,
                               kprime::Int64,
                               k::Int64,
                               edge_list::Array)
    adj_list = get_adj_list(n, edge_list)
    edge_sampler, nbr_weights, Z = construct_samplers(n, edge_list, adj_list)
    triangles = compute_weighted_triangles(n, kprime, k, Z, edge_list, adj_list, edge_sampler, nbr_weights)

    return triangles
end

function main()
    dataset_name, kprime, k = ARGS[1], parse(Int64, ARGS[2]), parse(Int64, ARGS[3])
    output_file = "../output/geom_mean_sampling_wr_$dataset_name.txt"

    # Load data in the form of simplices from the ScHoLP package.
    ex = read_txt_data(dataset_name)

    n, edge_list = get_edge_list(ex, 1.0)
    time = @elapsed triangles = construct_and_compute(n, kprime, k, edge_list)
    writedlm(output_file, triangles)
    open(output_file, "a") do f
        write(f, "Time: $time\n")
        write(f, "k: $k\n")
        write(f, "kprime: $kprime")
    end
end

main()