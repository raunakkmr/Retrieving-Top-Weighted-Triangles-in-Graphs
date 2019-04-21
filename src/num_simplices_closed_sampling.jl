include("./utils.jl")

using DataStructures
using DelimitedFiles
using Distributions
using ScHoLP
using SimpleWeightedGraphs
using SparseArrays
using StatsBase

using Random
Random.seed!(0)

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

function nbr_samplers(graph::SimpleGraph)
    degrees = [length(neighbors(graph, v)) for v in vertices(graph)]
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

function construct_samplers(ledges::Array,
                            lgraph::SimpleGraph,
                            rgraph::SimpleGraph)
    ledge_prob = Categorical(ones(length(ledges)) / length(ledges))
    ledge_sampler = sampler(ledge_prob)

    lgraph_samplers = nbr_samplers(lgraph)
    rgraph_samplers = nbr_samplers(rgraph)

    return [ledge_sampler, lgraph_samplers, rgraph_samplers]
end

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

    ledges, lgraph, rgraph = graphs
    ledge_sampler, lgraph_samplers, rgraph_samplers = samplers
    edge_id, rev_edge_id = ids

    for l in range(1, stop=s)

        # a is edge id, c is m + simplex number
        a, c = ledges[rand(ledge_sampler)]
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

    # Postprocessing.
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

    kprime = min(kprime, length(x))
    k = min(k, kprime)
    x_array = [(count,triangle) for (triangle,count) in zip(keys(x), values(x))]
    sorted_x = sort!(x_array, by = x -> x[1], rev=true)
    top_kprime = [(compute_weight(triangle), triangle) for (_,triangle) in sorted_x[1:kprime]]
    top_k = nlargest(k, top_kprime)

    return top_k
end

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
    samplers = construct_samplers(ledges, lgraph, rgraph)
    graphs = [ledges, lgraph, rgraph]
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