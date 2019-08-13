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

    auto GS = read_graph(argv[1], true);
    string tri_file = argv[2];
    int CHECK_TRIANGLES = atoi(argv[3]);
    int K = atoi(argv[4]);
    double START_TIME = atof(argv[5]);
    double MAX_TIME = atof(argv[6]);
    double INC = atof(argv[7]);
    string sampler = argv[8];  // E or W (edge or wedge)

#if PRINT_ARGS
    cerr << "=============================================" << endl;
    cerr << "ARGUMENTS" << endl;
    cerr << "=============================================" << endl;
    cerr << "Dataset: " << argv[1] << endl;
    cerr << "Triangle File: " << argv[2] << endl;
    cerr << "K: " << K << endl;
    cerr << "START_TIME: " << START_TIME << endl;
    cerr << "MAX_TIME: " << MAX_TIME << endl;
    cerr << "INC: " << INC << endl;
    cerr << "Sampler: " << sampler << endl;
    cerr << endl;
#endif

    double cur_time = START_TIME;
    vector<double> times;
    vector<vector<weighted_triangle>> sampling_tris;
    while (cur_time <= MAX_TIME) {
        times.push_back(cur_time);
        cerr << "Runnning for " << cur_time << " (s)" << endl;
        if (sampler == "E") {
            sampling_tris.push_back(edge_parallel_time_version(GS, thread::hardware_concurrency(), cur_time, -1, false));
        } else if (sampler == "W") {
            sampling_tris.push_back(wedge_parallel_time_version(GS, thread::hardware_concurrency(), cur_time, -1, false));
        }
        cur_time += INC;
    }

    if (CHECK_TRIANGLES) {
        auto all_tris = brute_force_sampler(GS, K);
        // auto all_tris = read_all_triangles_set(tri_file);
        for (int i = 0; i < (int) times.size(); i++) {
            cerr << "*** Comparing parallel sampling ***" << endl;
            compare_statistics(all_tris, sampling_tris[i], K, true);
        }
    } else {
        cerr << "Skipping brute force triangle computations since it seems to be infeasible..." << endl; 
    }

    return 0;
}
