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

    vector<vector<weighted_triangle>> edge_sampling_tris;
    vector<double> times;
    double cur_time = 0;
    while (cur_time < MAX_TIME) {
        cur_time += INC;
        cerr << "Calling with  " << cur_time<< endl;
        times.push_back(cur_time);
        auto res = edge_sampler_parallel_time(GS, thread::hardware_concurrency(), cur_time, -1, true);
        edge_sampling_tris.push_back(res);
    }

    if (CHECK_TRIANGLES) {
        auto all_tris = brute_force_sampler(GS, K);
        // auto all_tris = read_all_triangles_set(tri_file);
        for (int i = 0; i < (int) times.size(); i++) {
            cerr << "*** Comparing parallel edge sampling ***" << endl;
            cerr << "Time: " << times[i] << endl;
            compare_statistics(all_tris, edge_sampling_tris[i], K, true);
        }
    } else {
        cerr << "Skipping brute force triangle computations since it seems to be infeasible..." << endl; 
    }

    return 0;
}
