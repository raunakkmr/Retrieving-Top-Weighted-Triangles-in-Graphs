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

    auto G = read_graph(argv[1], true);
    string tri_file = argv[2];

    auto bf_tris = brute_force_sampler(G);
    auto read_tris = read_all_triangles_set(tri_file);

    int num_bf_not_in_read = 0;
    for (const auto &T : bf_tris) {
        if (read_tris.count(T) == 0) {
            // cerr << "T in bf but not in read" << '\n';
            num_bf_not_in_read++;
        }
    }

    cerr << "====================" << '\n';

    int num_read_not_in_bf = 0;
    for (const auto &T : read_tris) {
        if (bf_tris.count(T) == 0) {
            // cerr << "T in read but not in bf" << '\n';
            num_read_not_in_bf++;
        }
    }

    cerr << bf_tris.size() << " " << read_tris.size() << '\n';
    cerr << num_bf_not_in_read << " " << num_read_not_in_bf << '\n';

}
