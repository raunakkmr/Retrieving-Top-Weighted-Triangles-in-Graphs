include("./utils.jl")

using DataStructures
using DelimitedFiles
using Distributions
using ScHoLP
using SimpleWeightedGraphs
using StatsBase

function construct_samplers(n::Int64, edge_list::Array, adj_list::SimpleWeightedGraph)
    # # Construct a sampler for sampling edges.
    # edge_weights = [w for (_,w) in edge_list]
    # Z = sum(edge_weights)
    # edge_prob = Categorical(edge_weights ./ Z)
    # edge_sampler = sampler(edge_prob)

    # Construct sampler for sampling neighbors of a vertex.
    nbr_weights = [adj_list.weights[v,neighbors(adj_list, v)] for v in vertices(adj_list)]
    nbr_z = [sum(weights) for weights in nbr_weights]
    nbr_prob = []
    for (i, (weights, z)) in enumerate(zip(nbr_weights, nbr_z))
        normalized_prob = Array(weights) ./ z
        if length(neighbors(adj_list, i)) == 0
            normalized_prob = Array([1.0])
        end
        push!(nbr_prob, Categorical(normalized_prob))
    end
    nbr_samplers = [sampler(prob) for prob in nbr_prob]

    # Construct a sampler for sampling edges.
    edge_weights = [w*nbr_z[a]*nbr_z[b] for ((a,b),w) in edge_list]
    Z = sum(edge_weights)
    edge_prob = Categorical(edge_weights ./ Z)
    edge_sampler = sampler(edge_prob)

    return edge_sampler, nbr_samplers, Z
end

function compute_weighted_triangles(n::Int64,
                                    kprime::Int64,
                                    k::Int64,
                                    Z::Float64,
                                    edge_list::Array,
                                    adj_list::SimpleWeightedGraph,
                                    edge_sampler::Distributions.AliasTable,
                                    nbr_samplers::Array)
    s = 1000000  # Number of samples.
    x = Dict()  # Counters.
    num_triangles = 0

    for l in range(1, stop=s)

        (a,b), _ = edge_list[rand(edge_sampler)]

        if length(neighbors(adj_list, a)) == 1 || length(neighbors(adj_list, b)) == 1
            continue
        end

        c = b
        while c == b
            c = neighbors(adj_list, a)[rand(nbr_samplers[a])]
        end

        cprime = a
        while cprime == a
            cprime = neighbors(adj_list, b)[rand(nbr_samplers[b])]
        end

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
    function compute_weight((a,b,c))
        return adj_list.weights[a,b] * adj_list.weights[b,c] * adj_list.weights[a,c]
    end
    top_k = postprocess_counters(k, kprime, x, compute_weight)

    return top_k
end

function construct_and_compute(n::Int64,
                               kprime::Int64,
                               k::Int64,
                               edge_list::Array)
    adj_list = get_adj_list(n, edge_list)
    edge_sampler, nbr_samplers, Z = construct_samplers(n, edge_list, adj_list)
    triangles = compute_weighted_triangles(n, kprime, k, Z, edge_list, adj_list, edge_sampler, nbr_samplers)

    return triangles
end

function main()
    dataset_name, kprime, k = ARGS[1], parse(Int64, ARGS[2]), parse(Int64, ARGS[3])
    output_file = "../output/geom_mean_sampling_$dataset_name.txt"
    
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
