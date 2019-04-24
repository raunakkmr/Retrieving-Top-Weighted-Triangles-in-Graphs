include("./utils.jl")

using DataStructures
using DelimitedFiles
using Distributions
using ScHoLP
using SimpleWeightedGraphs
using StatsBase

"""Construct samplers for sampling edges and for sampling neighbors of a
vertex."""
function construct_samplers(n::Int64, edge_list::Array, adj_list::SimpleWeightedGraph)
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

"""Sample triangles and return the top k triangles based on geometric mean of
edge weights."""
function compute_weighted_triangles(n::Int64,
                                    k::Int64,
                                    Z::Float64,
                                    edge_list::Array,
                                    adj_list::SimpleWeightedGraph,
                                    edge_sampler::Distributions.AliasTable,
                                    nbr_samplers::Array)
    s = 1000000  # Number of samples.
    x = Dict()  # Counters.
    num_triangles = 0

    # Sample an edge, and a neighbor of each of the endpoints. If this forms a triangle, increment its counter.
    for l in range(1, stop=s)

        (a,b), _ = edge_list[rand(edge_sampler)]

        if length(neighbors(adj_list, a)) == 1 || length(neighbors(adj_list, b)) == 1
            continue
        end

        # Rejection sampling to sample neighbors of the endpoints. There are
        # other efficient ways to perform this sampling without rejection, but
        # this way is simple and works well in practice. Modifying the sampling
        # code below to use other strategies is not too difficult.

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

    # Postprocessing. Compute the weight of the sampled triangles and return
    # the top k from these.
    function compute_weight((a,b,c))
        return adj_list.weights[a,b] * adj_list.weights[b,c] * adj_list.weights[a,c]
    end
    top_k = postprocess_counters(k, x, compute_weight)

    return top_k
end

"""Construct the graph, and compute and return the triangles."""
function construct_and_compute(n::Int64,
                               k::Int64,
                               edge_list::Array)
    adj_list = get_adj_list(n, edge_list)
    edge_sampler, nbr_samplers, Z = construct_samplers(n, edge_list, adj_list)
    triangles = compute_weighted_triangles(n, k, Z, edge_list, adj_list, edge_sampler, nbr_samplers)

    return triangles
end

function main()
    dataset_name, k = ARGS[1], parse(Int64, ARGS[2])
    output_file = "../output/geom_mean_sampling_$dataset_name.txt"
    
    # Load data in the form of simplices from the ScHoLP package.
    ex = read_txt_data(dataset_name)

    n, edge_list = get_edge_list(ex, 1.0)
    time = @elapsed triangles = construct_and_compute(n, k, edge_list)
    writedlm(output_file, triangles)
    open(output_file, "a") do f
        write(f, "Time: $time\n")
        write(f, "k: $k")
    end
end

main()
