// This file runs the dynamic heavy-light and auto heavy-light algorithms, and
// optionally writes out all triangles to file.

#include <bits/stdc++.h>

#include "gflags/gflags.h"

#include "turan.h"
#include "graph.h"
#include "clique_sampler.h"
#include "triangle_sampler.h"

using namespace std;
using namespace wsdm_2019_graph;

DEFINE_string(filename, "", "Path to graph file.");
DEFINE_bool(binary, true, "Flag for if graph is in our binary format.");
DEFINE_string(format, "", "If binary is false this is required and indicates format of graph file. One of weighted, temporal, simplicial.");
DEFINE_int32(k, 1000, "Parameter k for top-k.");

#define WRITE_STATISTICS 0

int main(int argc, char* argv[]) {
  string usage("Dynamic and auto heavy-light algorithms for triangles.\n"
      "Sample usage:\n"
      "\t./adaptive_enumerator -filename=[fill_this_in] "
      "-binary=true -k=[fill_this_in]\n"
      "Additionally, these flags can be loaded from a single file "
      "with the option -flagfile=[filename].");

  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_filename.empty()) {
    std::cerr << "No file specified! Type ./adaptive_enumerator --help for a description of the program parameters." << std::endl;
    return 0;
  }

  ios::sync_with_stdio(0);
  cin.tie(0);
  srand(0); 

  auto GS = read_graph(FLAGS_filename, FLAGS_binary, FLAGS_format);
  int K = FLAGS_k;

  auto adaptive_heavy_light_tri = dynamic_heavy_light(GS, K, false);
  auto auto_thresh_heavy_light_tri = dynamic_heavy_light(GS, K, true);

#if WRITE_STATISTICS
  auto all_tris = brute_force_sampler(GS);

  // Write out triangles to a file.
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
