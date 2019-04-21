include("./utils.jl")

using DataStructures
using DelimitedFiles
using ScHoLP
using SimpleWeightedGraphs

"""Enumerate closed triangles and sort them based on the number of simplices that they appear in."""
function compute_weighted_triangles(ex::HONData)
    triangles = Dict()
    vertex_id = get_vertex_id(ex)
    let idx = 0
        for nvert in ex.nverts
            for i in range(idx+1, stop=idx+nvert)
                for j in range(i+1, stop=idx+nvert)
                    for k in range(j+1, stop=idx+nvert)
                        v_i, v_j, v_k = ex.simplices[i], ex.simplices[j], ex.simplices[k]
                        t = Tuple(sort([vertex_id[v_i], vertex_id[v_j], vertex_id[v_k]]))
                        if !haskey(triangles, t)
                            triangles[t] = 0
                        end
                        triangles[t] += 1
                    end
                end
            end
            idx += nvert
        end
    end
    triangles = collect(triangles)
    triangles = [(w, t) for (t, w) in triangles]
    # Sort by weight, break ties based on the ids of the vertices of the
    # triangle.
    sort!(triangles, by = x->[x[1],x[2][1],x[2][2],x[2][3]], rev=true)

    return triangles
end

"""Construct the graph, and compute and return the triangles in sorted order."""
function construct_and_compute(ex::HONData)
    triangles = compute_weighted_triangles(ex)
end

function main()
    dataset_name = ARGS[1]
    output_file = "../output/num_simplices_closed_exact_$dataset_name.txt"

    # Load data in the form of simplices from the ScHoLP package.
    ex = read_txt_data(dataset_name)

    time = @elapsed triangles = construct_and_compute(ex)
    writedlm(output_file, triangles)
    write(open(output_file, "a"), "Time: $time")
end

main()