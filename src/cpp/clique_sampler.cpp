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

  auto G = read_graph(argv[1]);
  int NUM_SAMPLES_EDGE = atoi(argv[2]);
  int NUM_SAMPLES_PATH = atoi(argv[3]);
  int CHECK_TRIANGLES = atoi(argv[4]);
  int K = atoi(argv[5]);
  double P = atof(argv[6]);
  // Since we are focusing on triangles for now, I am setting these to 0
  // manually to make it easier to specify arguments on the command line. We
  // should make the command line interface nicer at some point.
  // int CLIQUE_SIZE = atoi(argv[7]);
  // int NUM_SAMPLES_CLIQUE = atoi(argv[8]);
  int CLIQUE_SIZE = 0;
  int NUM_SAMPLES_CLIQUE = 0;

#if PRINT_ARGS
  cerr << "=============================================" << endl;
  cerr << "ARGUMENTS" << endl;
  cerr << "=============================================" << endl;
  cerr << "Dataset: " << argv[1] << endl;
  cerr << "NUM_SAMPLES_EDGE: " << NUM_SAMPLES_EDGE << endl;
  cerr << "NUM_SAMPLES_PATH: " << NUM_SAMPLES_PATH << endl;
  cerr << "K: " << K << endl;
  cerr << "P: " << P << endl;
  cerr << endl;
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

  modify_weights(G, P);

  auto edge_sampling_tri = edge_sampler(G, NUM_SAMPLES_EDGE);
  auto path_sampling_tri = path_sampler(G, NUM_SAMPLES_PATH);
  auto heavy_light_sampling_tri = heavy_light_sampler(G, 0.05);
  auto adaptive_heavy_light_tri = adaptive_heavy_light(G, K);
  auto auto_thresholded_heavy_light_tri = auto_thresholded_heavy_light(G, K);
  
  if (CHECK_TRIANGLES) {
    auto all_tris = brute_force_sampler(G);
    if (NUM_SAMPLES_EDGE > 0) {
      cerr << "*** Comparing edge sampling ***" << endl;
      compare_statistics(all_tris, edge_sampling_tri, K);
    }
    if (NUM_SAMPLES_PATH > 0) {
      cerr << "*** Comparing path sampling ***" << endl;
      compare_statistics(all_tris, path_sampling_tri, K);
    }
    cerr << "*** Comparing heavy light sampling ***" << endl;
    compare_statistics(all_tris, heavy_light_sampling_tri, K);
    cerr << "*** Comparing adaptive heavy light ***" << endl;
    compare_statistics(all_tris, adaptive_heavy_light_tri, K);
    cerr << "*** Comparing auto thresholded heavy light ***" << endl;
    compare_statistics(all_tris, auto_thresholded_heavy_light_tri, K);

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
  
  if (CLIQUE_SIZE > 1) {
    auto sampled_cliques = clique_sampler(G, CLIQUE_SIZE, NUM_SAMPLES_CLIQUE);
  }
  return 0;
}