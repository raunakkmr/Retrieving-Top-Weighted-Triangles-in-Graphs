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
  auto G = read_graph(argv[1], true);
  string tri_file = argv[2];
  int CHECK_TRIANGLES = atoi(argv[3]);
  int K = atoi(argv[4]);

#if PRINT_ARGS
  cerr << "=============================================" << endl;
  cerr << "ARGUMENTS" << endl;
  cerr << "=============================================" << endl;
  cerr << "Dataset: " << argv[1] << endl;
  cerr << "Triangle File: " << argv[2] << endl;
  cerr << "K: " << K << endl;
  cerr << endl;
#endif

  auto heavy_light_tri = heavy_light_sampler(G);
  auto adaptive_heavy_light_tri = adaptive_heavy_light(G, K);;
  auto auto_thresholded_tri = auto_thresholded_heavy_light(G, K);

  if (CHECK_TRIANGLES) {
    auto all_tris = brute_force_sampler(G);
    // auto all_tris = read_all_triangles_set(tri_file);
    cerr << "*** Comparing heavy light***" << endl;
    compare_statistics(all_tris, heavy_light_tri, K);
    cerr << "*** Comparing adaptive heavy light***" << endl;
    compare_statistics(all_tris, adaptive_heavy_light_tri, K);
    cerr << "*** Comparing auto thresholded heavy light***" << endl;
    compare_statistics(all_tris, auto_thresholded_tri, K);

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
