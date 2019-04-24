include("./utils.jl")

using ArgParse
using DataStructures
using DelimitedFiles
using Distributions
using ScHoLP
using SimpleWeightedGraphs
using SparseArrays
using StatsBase

"""Construct samplers for sampling edges and sampling neighbors of vertices
from the bipartite graph between edges and simplices."""
function construct_samplers(graph::SimpleGraph, gedges::Array)
    edge_prob = Categorical(ones(length(gedges)) / length(gedges))
    edge_sampler = sampler(edge_prob)

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

    return [edge_sampler, graph_samplers]
end

"""Sample triangles / structures and return the top k triangles based on the number of simplices that contain at least one edge of this triangle."""
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
    num_structures = 0

    gedges, graph, adj_list = graphs
    edge_sampler, graph_samplers = samplers
    edge_id, rev_edge_id = ids

    # Sample (edge1, simplex), (edge2, simplex) and (edge3, simplex). If these edges from a triangle, then sample (edge1, simplex'). If (simplex', edg2) and (simplex', edge3) exist then increment the counter for the triangle corresponding to edge1, edge2 and edge3.
    for l in range(1, stop=s)

        # a is edge id, c is m + simplex number.
        a, c = gedges[rand(edge_sampler)]
        # Need 3 different edges to be connected to simplex to have a valid
        # triangle.
        if length(neighbors(graph, c)) < 3
            continue
        end
        # b, d are edge id.
        b = neighbors(graph, c)[rand(graph_samplers[c])]
        d = neighbors(graph, c)[rand(graph_samplers[c])]
        # Check if triangle.
        ae, be, de = rev_edge_id[a], rev_edge_id[b], rev_edge_id[d]
        t = sort(collect(Set([ae[1],ae[2],be[1],be[2],de[1],de[2]])))
        if length(t) != 3
            continue
        end
        u, v, w = t
        if !has_edge(adj_list, u, v) || !has_edge(adj_list, u, w) || !has_edge(adj_list, v, w)
            continue
        end
        # cprime is m + simplex number.
        cprime = neighbors(graph, a)[rand(graph_samplers[a])]
        if !has_edge(graph, cprime, b) || !has_edge(graph, cprime, d)
            continue
        end
        num_structures += 1
        if !haskey(x, t)
            x[t] = 0
        end
        x[t] += 1
    end

    # Postprocessing. Compute the weight of the triangles corresponding to the
    # top kprime counters, and return the top k triangles.
    function compute_weight((a, b, c))
        e1, e2, e3 = (a, b), (a, c), (b, c)
        ar = neighbors(graph, edge_id[e1])
        br = neighbors(graph, edge_id[e2])
        cr = neighbors(graph, edge_id[e3])
        as = sparsevec(ar, ones(length(ar)), m+num_simplices)
        bs = sparsevec(br, ones(length(br)), m+num_simplices)
        cs = sparsevec(cr, ones(length(cr)), m+num_simplices)
        weight = nnz(as .+ bs .+ cs)

        return weight
    end
    top_k = postprocess_counters(k, x, compute_weight, kprime)

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
    adj_list = get_adj_list(n, edge_list)
    num_simplices, graph, edge_id, rev_edge_id = get_edges_to_simplices(m, ex)
    gedges = [(src(e), dst(e)) for e in edges(graph)]
    samplers = construct_samplers(graph, gedges)
    graphs = [gedges, graph, adj_list]
    ids = [edge_id, rev_edge_id]
    triangles = compute_weighted_triangles(n, m, kprime, k, num_simplices, graphs, samplers, ids, ex)

    return triangles
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--dataset", "-d"
            help = "dataset name"
            arg_type = String
            required = true
        "--kprime"
            help = "budget for computing weights, default: 500"
            arg_type = Int64
            default = 500
        "--samples", "-s"
            help = "number of samples, default: 100000"
            arg_type = Int64,
            default = 100000
        "-k"
            help = "parameter k for returning top-k triangles, default: 25"
            arg_type = Int64
            default = 25
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    dataset_name = parsed_args["dataset"]
    kprime, k = parsed_args["kprime"], parsed_args["k"]
    output_file = "../output/num_simplices_all_sampling_$dataset_name.txt"

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