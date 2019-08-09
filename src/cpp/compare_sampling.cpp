#include <bits/stdc++.h>

#include "graph.h"
#include "triangle_sampler.h"

using namespace std;
using namespace wsdm_2019_graph;

#define PRINT_ARGS 1
#define PRINT_STATISTICS 0

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(0);
    cin.tie(0);
    srand(0); 

    // auto G = read_graph(argv[1]);
    auto GS = read_graph(argv[1], true);
    string tri_file = argv[2];
    int CHECK_TRIANGLES = atoi(argv[3]);
    int K = atoi(argv[4]);
    double MAX_TIME = atof(argv[5]);
    double INC = atof(argv[6]);

#if PRINT_ARGS
    cerr << "=============================================" << endl;
    cerr << "ARGUMENTS" << endl;
    cerr << "=============================================" << endl;
    cerr << "Dataset: " << argv[1] << endl;
    cerr << "Triangle File: " << argv[2] << endl;
    cerr << "K: " << K << endl;
    cerr << "MAX_TIME: " << MAX_TIME << endl;
    cerr << "INC: " << INC << endl;
    cerr << endl;
#endif

    // auto res = edge_sampler_time(GS, MAX_TIME, INC);
    auto res = edge_time_version(GS, MAX_TIME, INC);
    auto edge_sampling_tri = res.first;
    auto edge_sampling_times = res.second;
    // res = wedge_sampler_time(GS, MAX_TIME, INC);
    res = wedge_time_version(GS, MAX_TIME, INC);
    auto wedge_sampling_tri = res.first;
    auto wedge_sampling_times = res.second;
    // res = path_sampler_time(GS, MAX_TIME, INC);
    res = path_time_version(GS, MAX_TIME, INC);
    auto path_sampling_tri = res.first;
    auto path_sampling_times = res.second;

    if (CHECK_TRIANGLES) {
        auto all_tris = brute_force_sampler(GS, K);
        // auto all_tris = read_all_triangles_set(tri_file);
        cerr << "*** Comparing edge sampling ***" << endl;
        compare_statistics_time(all_tris, edge_sampling_tri, edge_sampling_times, K, true);
        cerr << "*** Comparing wedge sampling ***" << endl;
        compare_statistics_time(all_tris, wedge_sampling_tri, wedge_sampling_times, K, true);
        cerr << "*** Comparing path sampling ***" << endl;
        compare_statistics_time(all_tris, path_sampling_tri, path_sampling_times, K, true);

        Graph &G = GS.G;

        // Write out triangles to a file
        bool write_out_stats = false;
        if (write_out_stats) {
            ofstream triangle_file("triangles.txt");
            map<int, int> num_tris;
            for (auto t : all_tris) {
                triangle_file << t.weight << '\n';
                num_tris[get<0>(t.vertices)]++;
                num_tris[get<1>(t.vertices)]++;
                num_tris[get<2>(t.vertices)]++;
            }
            triangle_file.close();

            ofstream tris_to_weight_file("ntris-to-weight.txt");
            for (int i = 0; i < (int) G[i].size(); i++) {
                double weight_sum = 0;
                for (auto e : G[i]) {
                    weight_sum += e.wt;
                }
                tris_to_weight_file << weight_sum << " " << num_tris[i] << '\n';
            }
            tris_to_weight_file.close();

            ofstream degree_file("degree-to-weight.txt");
            for (int i = 0; i < (int) G[i].size(); i++) {
                double weight_sum = 0;
                for (auto e : G[i]) {
                    weight_sum += e.wt;
                }
                degree_file << G[i].size() << " " << weight_sum << '\n';
            }
            degree_file.close();
        }
    } else {
        cerr << "Skipping brute force triangle computations since it seems to be infeasible..." << endl; 
    }

    return 0;
}
