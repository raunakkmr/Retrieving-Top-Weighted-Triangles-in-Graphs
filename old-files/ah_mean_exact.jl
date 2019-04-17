include("./utils.jl")

using DataStructures
using DelimitedFiles
using ScHoLP
using SimpleWeightedGraphs

function compute_weighted_triangles(adj_list::SimpleWeightedGraph)
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
                weight_uvw = adj_list.weights[v,u] + adj_list.weights[v,w]
                if has_edge(adj_list, u, w)
                    weight_uvw += adj_list.weights[u,w]
                    push!(triangles, (weight_uvw, triangle))
                    push!(seen, triangle)
                elseif has_edge(adj_list, w,u)
                    weight_uvw += adj_list.weights[w,u]
                    push!(triangles, (weight_uvw, triangle))
                    push!(seen, triangle)
                end
            end
        end
    end
    sort!(triangles, by = x->[x[1],x[2][1],x[2][2],x[2][3]], rev=true)

    return triangles
end

function construct_and_compute(n::Int64, edge_list::Array)
    adj_list = get_adj_list_higher_deg(n, edge_list)
    triangles = compute_weighted_triangles(adj_list)

    return triangles
end

function main()
    dataset_name = ARGS[1]
    mode = ARGS[2]  # 'a' or 'h'
    output_file = "../output/$(mode)_mean_exact_$dataset_name.txt"

    # Load data in the form of simplices from the ScHoLP package.
    ex = read_txt_data(dataset_name)

    n, edge_list = get_edge_list(ex, mode == "a" ? 1.0 : -1.0)
    time = @elapsed triangles = construct_and_compute(n, edge_list)
    writedlm(output_file, triangles)
    write(open(output_file, "a"), "Time: $time")
end

main()