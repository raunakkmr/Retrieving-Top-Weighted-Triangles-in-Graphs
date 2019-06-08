#include <bits/stdc++.h>

#include "graph.h"
#include "triangle_sampler.h"

using namespace std;
using namespace wsdm_2019_graph;

int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    srand(0);

    auto G = read_graph(argv[1]);
    string dataset = argv[2];
    double P = atof(argv[3]);

    modify_weights(G, P);
    for (int i = 0; i < G.size(); i++) {
        for (auto e : G[i]) {
            if (e.wt < 0) {
                cerr << "Negative weight" << endl;
            }
        }
    }
    return -1;

    auto all_tris = brute_force_sampler(G);

    ofstream triangle_file(dataset+"-triangles.txt");
    map<long long, long long> num_tris;
    for (auto t : all_tris) {
      triangle_file << t.weight << '\n';
      num_tris[get<0>(t.vertices)]++;
      num_tris[get<1>(t.vertices)]++;
      num_tris[get<2>(t.vertices)]++;
    }
    triangle_file.close();

    ofstream tris_to_weight_file(dataset+"-ntris-to-weight.txt");
    for (int i = 0; i < (int) G[i].size(); i++) {
      double weight_sum = 0;
      for (auto e : G[i]) {
        weight_sum += e.wt;
      }
      tris_to_weight_file << weight_sum << " " << num_tris[i] << '\n';
    }
    tris_to_weight_file.close();

    ofstream degree_file(dataset+"-degree-to-weight.txt");
    for (int i = 0; i < (int) G[i].size(); i++) {
      double weight_sum = 0;
      for (auto e : G[i]) {
        weight_sum += e.wt;
      }
      degree_file << G[i].size() << " " << weight_sum << '\n';
    }
    degree_file.close();

    return 0;
}