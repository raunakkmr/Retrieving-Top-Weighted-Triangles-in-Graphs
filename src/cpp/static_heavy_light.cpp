#include <bits/stdc++.h>

#include "graph.h"
#include "triangle_sampler.h"

using namespace std;
using namespace wsdm_2019_graph;

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(0);
    cin.tie(0);
    srand(0); 

    vector<string> datasets({"tags-stack-overflow", "eth"});
    vector<int> kvals({25, 1000, 40000});
    vector<double> thresholds({0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.90, 0.95, 1.00});

    for (const auto &dataset : datasets) {
        auto GS = read_graph(argv[1], true);
        for (const auto &k : kvals) {
            auto all_tris = brute_force_sampler(GS, k);
            for (const auto &threshold : thresholds) {
                cerr << "Threshold: " << threshold << endl;
                auto heavy_light_tri = heavy_light_sampler(GS, k, threshold);
                compare_statistics(all_tris, heavy_light_tri, k, true);
                cerr << endl;
            }
        }
    }

}