include("./utils.jl")

using DataStructures
using DelimitedFiles
using Distributions
using LightGraphs
using ScHoLP
using SimpleWeightedGraphs
using StatsBase

function construct_samplers(n::Int64, edge_list::Array, adj_list::SimpleWeightedGraph)
    degrees = [length(neighbors(adj_list, v)) for v in vertices(adj_list)]

    # Weights of the neighbors of every vertex.
    nbr_weights = [adj_list.weights[v, neighbors(adj_list, v)] for v in vertices(adj_list)]
    nbr_prefix = [[sum(@views weight[1:i]) for i in range(1, stop=length(weight))] for weight in nbr_weights]
    nbr_z = [sum(weights) for weights in nbr_weights]
    nbr_prefix_z = [sum(prefix) for prefix in nbr_prefix]

    # Construct a sampler for sampling edge weights.
    # edge_weights = [w for (_,w) in edge_list]
    edge_weights = []
    for ((a,b),w) in edge_list
        a_nbrs = neighbors(adj_list, a)
        a_nbrs_weights = pweights(adj_list.weights[a, filter(e->e!=b, a_nbrs)])
        b_nbrs = neighbors(adj_list, b)
        b_nbrs_weights = pweights(adj_list.weights[b, filter(e->e!=a, b_nbrs)])
        wt = (degrees[a]-1)*(degrees[b]-1)*w
        wt += (degrees[b]-1)*sum(a_nbrs_weights)
        wt += (degrees[a]-1)*sum(b_nbrs_weights)
        push!(edge_weights, wt)
    end
    Z = sum(edge_weights)
    edge_prob = Categorical(edge_weights ./ Z)
    edge_sampler = sampler(edge_prob)

    return degrees, edge_sampler, nbr_weights, nbr_prefix, nbr_z, nbr_prefix_z, Z
end

function compute_weighted_triangles(n::Int64,
                                    kprime::Int64,
                                    k::Int64,
                                    Z::Float64,
                                    edge_list::Array,
                                    adj_list::SimpleWeightedGraph,
                                    edge_sampler::Distributions.AliasTable,
                                    nbr_weights::Array,
                                    nbr_prefix::Array,
                                    nbr_z::Array,
                                    nbr_prefix_z::Array,
                                    degrees::Array)
    s = 1000000  # Number of samples.
    x = Dict()  # Counters.
    num_triangles = 0

    for l in range(1, stop=s)

        (a,b), _ = edge_list[rand(edge_sampler)]

        if length(neighbors(adj_list, a)) == 1 || length(neighbors(adj_list, b)) == 1
            continue
        end

        a_nbrs = neighbors(adj_list, a)
        u = rand()
        modified_u = (u * (nbr_z[a] - adj_list.weights[a,b]))
        modified_u -= ((degrees[b]-1)*adj_list.weights[a,b])
        modified_u -= (nbr_z[b]-adj_list.weights[a,b])
        modified_u /= (degrees[b]-1)
        b_index = searchsorted(a_nbrs, b)[1]
        index = searchsortedfirst(nbr_prefix[a], modified_u)
        if a_nbrs[index] != b
            c = a_nbrs[index]
        else
            if b_index == 1
                index = searchsortedfirst(nbr_prefix[a][2:end], modified_u-nbr_prefix[a][1])
                c = a_nbrs[1+index]
                c = a_nbrs[rand(2:end)]
            elseif b_index == degrees[a]
                index = searchsortedfirst(nbr_prefix[a][1:end-1], modified_u-nbr_prefix[a][1])
                c = a_nbrs[index]
                c = a_nbrs[rand(1:end-1)]
            else
                p = nbr_prefix[a][b_index-1] / nbr_z[a]
                coin_flip = rand(Bernoulli(p))
                if coin_flip == 1
                    c = a_nbrs[rand(1:b_index-1)]
                else
                    c = a_nbrs[rand(b_index+1:end)]
                end
            end
        end

        b_nbrs = neighbors(adj_list, b)
        u = rand()
        modified_u = (u * (nbr_z[b] - adj_list.weights[a,b]))
        modified_u -= (adj_list.weights[a,c] + adj_list.weights[a,b])
        a_index = searchsorted(b_nbrs, a)[1]
        index = searchsortedfirst(nbr_prefix[b], modified_u)
        if b_nbrs[index] != a
            cprime = b_nbrs[index]
        else
            if a_index == 1
                index = searchsortedfirst(nbr_prefix[b][2:end], modified_u-nbr_prefix[b][1])
                cprime = b_nbrs[1+index]
            elseif a_index == degrees[b]
                index = searchsortedfirst(nbr_prefix[b][1:end-1], modified_u-nbr_prefix[b][1])
                cprime = b_nbrs[index]
            else
                p = nbr_prefix[b][a_index-1] / nbr_z[b]
                coin_flip = rand(Bernoulli(p))
                if coin_flip == 1
                    cprime = b_nbrs[rand(1:a_index-1)]
                else
                    cprime = b_nbrs[rand(a_index+1:end)]
                end
            end
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
    kprime = min(kprime, length(x))
    k = min(k, kprime)
    function compute_weight((a,b,c))
        return adj_list.weights[a,b] + adj_list.weights[b,c] + adj_list.weights[a,c]
    end
    x_array = [(count,triangle) for (triangle,count) in zip(keys(x), values(x))]
    sorted_x = sort!(x_array, by = x -> x[1])
    top_kprime = [(compute_weight(triangle), triangle) for (_,triangle) in sorted_x[1:kprime]]
    top_k = nlargest(k, top_kprime)

    return top_k
end

function construct_and_compute(n::Int64,
                               kprime::Int64,
                               k::Int64,
                               edge_list::Array)
    adj_list = get_adj_list(n, edge_list)
    degrees, edge_sampler, nbr_weights, nbr_prefix, nbr_z, nbr_prefix_z, Z = construct_samplers(n, edge_list, adj_list)
    triangles = compute_weighted_triangles(n, kprime, k, Z, edge_list, adj_list, edge_sampler, nbr_weights, nbr_prefix, nbr_z, nbr_prefix_z, degrees)

    return triangles
end

function main()
    dataset_name, kprime, k, mode = ARGS[1], parse(Int64, ARGS[2]), parse(Int64, ARGS[3]), ARGS[4]  # mode is 'a' or 'h'
    output_file = "../output/$(mode)_mean_sampling_fast_$dataset_name.txt"

    # Load data in the form of simplices from the ScHoLP package.
    ex = read_txt_data(dataset_name)

    n, edge_list = get_edge_list(ex)
    time = @elapsed triangles = construct_and_compute(n, kprime, k, edge_list)
    writedlm(output_file, triangles)
    open(output_file, "a") do f
        write(f, "Time: $time\n")
        write(f, "k: $k\n")
        write(f, "kprime: $kprime")
    end
end

main()