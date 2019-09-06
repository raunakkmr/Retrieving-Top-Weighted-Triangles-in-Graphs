// This file runs the various heavy-light algorithms for all datasets, and brute
// force for all datasets except spotify since it takes over 24 hours to
// terminate.

#include <bits/stdc++.h>

#include "gflags/gflags.h"

#include "graph.h"
#include "triangle_sampler.h"

using namespace std;
using namespace wsdm_2019_graph;

#define PRINT_ARGS 1
#define PRINT_STATISTICS 0

DEFINE_string(filename, "", "Path to graph file.");
DEFINE_bool(binary, true, "Flag for if graph is in our binary format.");
DEFINE_string(format, "", "If binary is false this is required and indicates format of graph file. One of weighted, temporal, simplicial.");
DEFINE_int32(k, 1000, "Parameter k for top-k.");

int main(int argc, char* argv[]) {
  string usage("Deterministic algorithms for triangles.\n"
      "Sample usage:\n"
      "\t./compare_deterministic -filename=[fill_this_in] "
      "-binary=true -k=[fill_this_in]\n"
      "Additionally, these flags can be loaded from a single file "
      "with the option -flagfile=[filename].");

  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_filename.empty()) {
    std::cerr << "No file specified! Type ./compare_deterministic --help for a description of the program parameters." << std::endl;
    return 0;
  }

  ios::sync_with_stdio(0);
  cin.tie(0);
  srand(0); 

  auto GS = read_graph(FLAGS_filename, FLAGS_binary, FLAGS_format);
  int K = FLAGS_k;

#if PRINT_ARGS
  cerr << "=============================================" << endl;
  cerr << "ARGUMENTS" << endl;
  cerr << "=============================================" << endl;
  cerr << "Dataset: " << FLAGS_filename << endl;
  cerr << "Binary: " << FLAGS_binary << endl;
  cerr << "Format: " << FLAGS_format << endl;
  cerr << "K: " << K << endl;
  cerr << endl;
#endif

  auto heavy_light_tri = heavy_light_sampler(GS, K);
  auto adaptive_heavy_light_tri = dynamic_heavy_light(GS, K, false);
  auto auto_thresholded_tri = dynamic_heavy_light(GS, K, true);

  // Compute accuracy by comparing with dynamic heavy-light on the
  // Spotify dataset. For all others run brute force.
  auto all_tris = adaptive_heavy_light_tri;
  string spotify = "spotify";
  if (spotify.find(FLAGS_filename) == string::npos) {
    all_tris = brute_force_sampler(GS, K);
  }
  cerr << "*** Comparing heavy light***" << endl;
  compare_statistics(adaptive_heavy_light_tri, heavy_light_tri, K, true);
  cerr << "*** Comparing dynamic heavy light***" << endl;
  compare_statistics(all_tris, adaptive_heavy_light_tri, K, true);
  cerr << "*** Comparing auto heavy light***" << endl;
  compare_statistics(all_tris, auto_thresholded_tri, K, true);

  /*
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
   */

  return 0;
}
