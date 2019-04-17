include("./utils.jl")

using DataStructures
using DelimitedFiles
using ScHoLP
using SimpleWeightedGraphs

function compute_weighted_triangles(ex::HONData)
    triangles, vertex_id = Dict(), Dict()
    let idx = 0
        for nvert in ex.nverts
            for i in range(idx+1, stop=idx+nvert)
                for j in range(i+1, stop=idx+nvert)
                    for k in range(j+1, stop=idx+nvert)
                        v_i, v_j, v_k = ex.simplices[i], ex.simplices[j], ex.simplices[k]
                        if !haskey(vertex_id, v_i)
                            vertex_id[v_i] = length(vertex_id) + 1
                        end
                        if !haskey(vertex_id, v_j)
                            vertex_id[v_j] = length(vertex_id) + 1
                        end
                        if !haskey(vertex_id, v_k)
                            vertex_id[v_k] = length(vertex_id) + 1
                        end
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
    sort!(triangles, by = x->[x[1],x[2][1],x[2][2],x[2][3]], rev=true)

    return triangles
end

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