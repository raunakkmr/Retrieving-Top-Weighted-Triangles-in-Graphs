#include <bits/stdc++.h>

#include "turan.h"
#include "graph.h"
#include "clique_sampler.h"
#include "triangle_sampler.h"

using namespace std;
using namespace wsdm_2019_graph;

#define PRINT_ARGS 1
#define PRINT_STATISTICS 0

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(0);
    cin.tie(0);
    srand(0); 

    auto G = read_graph(argv[1], true);
    int NUM_SAMPLES_EDGE = atoi(argv[2]);
    int CLIQUE_SIZE = atoi(argv[3]);
    int NUM_SAMPLES_CLIQUE = atoi(argv[4]);

#if PRINT_ARGS
    cerr << "=============================================" << endl;
    cerr << "ARGUMENTS" << endl;
    cerr << "=============================================" << endl;
    cerr << "Dataset: " << argv[1] << endl;
    cerr << "NUM_SAMPLES_EDGE: " << NUM_SAMPLES_EDGE << endl;
    cerr << "CLIQUE_SIZE, NSAMPS: " << CLIQUE_SIZE << " " << NUM_SAMPLES_CLIQUE << endl;
#endif

#if PRINT_STATISTICS
    cerr << "=============================================" << endl;
    cerr << "Computing graph statistics" << endl;
    cerr << "=============================================" << endl;
    cerr << "Computing degeneracy..." << endl;
    auto retval = compute_degeneracy(G);
    auto degenOrder = retval.degenOrder;
    auto degeneracy = (--retval.degenFreq.end())->first;

    cerr << "Degeneracy distribution:" << endl;
    double sum = 0, num = 0;
    for (const auto& kv : retval.degenFreq) {
        //cerr << kv.first << " " << kv.second << endl;
        sum += kv.first * kv.second;
        num += kv.second;
    }
    cerr << "Average degeneracy: " << sum / num << endl;
    cerr << "Degeneracy: " << degeneracy << endl;
    cerr << endl;
#endif

    int nthreads = thread::hardware_concurrency();

    auto edge_sampling_tri = edge_sampler(G, NUM_SAMPLES_EDGE);
    auto edge_sampling_tri_parallel = edge_sampler_parallel(G, NUM_SAMPLES_EDGE, nthreads);
    // auto wedge_sampling_tri = wedge_sampler(G, NUM_SAMPLES_WEDGE);
    // auto path_sampling_tri = path_sampler(G, NUM_SAMPLES_PATH);
    // auto heavy_light_sampling_tri = heavy_light_sampler(G, 0.05);
    // auto adaptive_heavy_light_tri = adaptive_heavy_light(G, K);
    // auto auto_thresholded_heavy_light_tri = auto_thresholded_heavy_light(G, K);

    if (CLIQUE_SIZE > 1) {
        auto sampled_cliques = clique_sampler(G, CLIQUE_SIZE, NUM_SAMPLES_CLIQUE);
        auto sampled_cliques_parallel = clique_sampler_parallel(G, CLIQUE_SIZE, NUM_SAMPLES_CLIQUE, nthreads);
    }
    return 0;
}
