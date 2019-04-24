include("./utils.jl")

using ArgParse
using DataStructures
using DelimitedFiles
using SparseArrays
using ScHoLP
using SimpleWeightedGraphs

"""Enumerate triangles and sort them based on the number of simplices
that contain at least one edge of this triangle."""
function compute_weighted_triangles(num_simplices::Int64,
                                    adj_list::SimpleWeightedGraph, 
                                    appearances::SimpleGraph,
                                    edge_id::Dict)
    m = ne(adj_list)
    triangles = []
    seen = Set()
    for v in vertices(adj_list)
        neighbors_v = neighbors(adj_list, v)
        len = length(neighbors_v)
        for i in range(1, stop=len)
            for j in range(i+1, stop=len)
                u, w = neighbors_v[i], neighbors_v[j]
                triangle = Tuple(sort([u,v,w]))
                if in(triangle, seen)
                    continue
                end
                # Edges e = (a, b) follow the property a < b.
                e_vu = v < u ? (v,u) : (u,v)
                e_vw = v < w ? (v,w) : (w,v)
                avu = neighbors(appearances, edge_id[e_vu])
                avw = neighbors(appearances, edge_id[e_vw])
                svu = sparsevec(avu, ones(length(avu)), m+num_simplices)
                svw = sparsevec(avw, ones(length(avw)), m+num_simplices)
                if has_edge(adj_list, u, w)
                    u, w = u < w ? (u,w) : (w,u)
                    auw = neighbors(appearances, edge_id[(u,w)])
                    suw = sparsevec(auw, ones(length(auw)), m+num_simplices)
                    weight_uvw = nnz(svu .+ svw .+ suw)
                    # weight_uvw = sum(svu .* svw .* suw)
                    push!(triangles, (weight_uvw, triangle))
                    push!(seen, triangle)
                elseif has_edge(adj_list, w,u)
                    u, w = u < w ? (u,w) : (w,u)
                    auw = neighbors(appearances, edge_id[(w,u)])
                    suw = sparsevec(auw, ones(length(auw)), num_simplices)
                    weight_uvw = nnz(svu .+ svw .+ suw)
                    # weight_uvw = sum(svu .* svw .* suw)
                    push!(triangles, (weight_uvw, triangle))
                    push!(seen, triangle)
                end
            end
        end
    end
    # Sort by weight, break ties based on the ids of the vertices of the
    # triangle.
    sort!(triangles, by = x->[x[1],x[2][1],x[2][2],x[2][3]], rev=true)
    return triangles

end

"""Construct the graph, and compute and return the triangles in sorted order."""
function construct_and_compute(ex::HONData)
    n, edge_list = get_edge_list(ex, 1.0)
    adj_list = get_adj_list_higher_deg(n, edge_list)
    num_simplices, appearances, edge_id, _ = get_edges_to_simplices(length(edge_list), ex)
    triangles = compute_weighted_triangles(num_simplices, adj_list, appearances, edge_id)
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--dataset", "-d"
            help = "dataset name"
            arg_type = String
            required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    dataset_name = parsed_args["dataset"]
    output_file = "../output/num_simplices_all_exact_$dataset_name.txt"

    # Load data in the form of simplices from the ScHoLP package.
    ex = read_txt_data(dataset_name)

    time = @elapsed triangles = construct_and_compute(ex)
    writedlm(output_file, triangles)
    write(open(output_file, "a"), "Time: $time")
end

main()