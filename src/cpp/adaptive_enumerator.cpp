#include <bits/stdc++.h>

#include "turan.h"
#include "graph.h"
#include "clique_sampler.h"
#include "triangle_sampler.h"

using namespace std;
using namespace wsdm_2019_graph;

#define WRITE_STATISTICS 0

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(0);
    cin.tie(0);
    srand(0); 

    auto G = read_graph(argv[1]);
    int K = atoi(argv[2]);

    auto adaptive_heavy_light_tri = adaptive_heavy_light(G, K);
    auto auto_thresh_heavy_light_tri = auto_thresholded_heavy_light(G, K);

#if WRITE_STATISTICS
    auto all_tris = brute_force_sampler(G);

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
#endif

    return 0;
}
