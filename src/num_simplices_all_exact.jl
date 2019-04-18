include("./utils.jl")

using DataStructures
using DelimitedFiles
using SparseArrays
using ScHoLP
using SimpleWeightedGraphs

function get_edges_to_simplices(edge_list::Array, ex::HONData)
    appearances = SimpleGraph(length(edge_list))
    num_simplices = 0
    edge_id, vertex_id = Dict(), Dict()
    let idx = 0
        for nvert in ex.nverts
            num_simplices += 1
            for i in range(idx+1, stop=idx+nvert)
                for j in range(i+1, stop=idx+nvert)
                    v_i, v_j = ex.simplices[i], ex.simplices[j]
                    if !haskey(vertex_id, v_i)
                        vertex_id[v_i] = length(vertex_id) + 1
                    end
                    if !haskey(vertex_id, v_j)
                        vertex_id[v_j] = length(vertex_id) + 1
                    end
                    u, v = -1, -1
                    if vertex_id[v_i] < vertex_id[v_j]
                        u, v = vertex_id[v_i], vertex_id[v_j]
                    else
                        u, v = vertex_id[v_j], vertex_id[v_i]
                    end
                    edge = (u, v)
                    if !haskey(edge_id, edge)
                        edge_id[edge] = length(edge_id) + 1
                    end
                    add_edge!(appearances, edge_id[edge], num_simplices)
                end
            end
            idx += nvert
        end
    end

    return num_simplices, appearances, edge_id 
end

function compute_weighted_triangles(num_simplices::Int64,
                                    adj_list::SimpleWeightedGraph, 
                                    appearances::SimpleGraph,
                                    edge_id::Dict)
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
                svu = sparsevec(avu, ones(length(avu)), num_simplices)
                svw = sparsevec(avw, ones(length(avw)), num_simplices)
                if has_edge(adj_list, u, w)
                    auw = neighbors(appearances, edge_id[(u,w)])
                    suw = sparsevec(auw, ones(length(auw)), num_simplices)
                    weight_uvw = sum(svu .* svw .* suw)
                    push!(triangles, (weight_uvw, triangle))
                    push!(seen, triangle)
                elseif has_edge(adj_list, w,u)
                    auw = neighbors(appearances, edge_id[(w,u)])
                    suw = sparsevec(auw, ones(length(auw)), num_simplices)
                    weight_uvw = sum(svu .* svw .* suw)
                    push!(triangles, (weight_uvw, triangle))
                    push!(seen, triangle)
                end
            end
        end
    end
    sort!(triangles, by = x->[x[1],x[2][1],x[2][2],x[2][3]], rev=true)
    return triangles

end

function construct_and_compute(ex::HONData)
    n, edge_list = get_edge_list(ex, 1.0)
    adj_list = get_adj_list_higher_deg(n, edge_list)
    num_simplices, appearances, edge_id = get_edges_to_simplices(edge_list, ex)
    triangles = compute_weighted_triangles(num_simplices, adj_list, appearances, edge_id)
end

function main()
    dataset_name = ARGS[1]
    output_file = "../output/num_simplices_all_exact_$dataset_name.txt"

    # Load data in the form of simplices from the ScHoLP package.
    ex = read_txt_data(dataset_name)

    time = @elapsed triangles = construct_and_compute(ex)
    writedlm(output_file, triangles)
    write(open(output_file, "a"), "Time: $time")
end

main()