// This file runs the static heavy-light algorithm for different thresholds in
// order to illustrate the trade-off between accuracy and time as threshold
// varies.

#include <bits/stdc++.h>

#include "gflags/gflags.h"

#include "graph.h"
#include "triangle_sampler.h"

using namespace std;
using namespace wsdm_2019_graph;

DEFINE_string(filename, "", "Path to graph file.");
DEFINE_bool(binary, true, "Flag for if graph is in our binary format.");
DEFINE_string(format, "", "If binary is false this is required and indicates format of graph file. One of weighted, temporal, simplicial.");
DEFINE_int32(k, 1000, "Parameter k for top-k.");

int main(int argc, char* argv[]) {
  string usage("Illustrate trade-off in static heavy-light algorithm.\n"
      "Sample usage:\n"
      "\t./static_heavy_lightt -filename=[fill_this_in] "
      "-binary=true -k=[fill_this_in]\n"
      "Additionally, these flags can be loaded from a single file "
      "with the option -flagfile=[filename].");

  ios::sync_with_stdio(0);
  cin.tie(0);
  srand(0); 

  // vector<string> datasets({"tags-stack-overflow", "eth"});
  // vector<int> kvals({25, 1000, 40000});
  vector<double> thresholds({0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.90, 0.95, 1.00});

  auto GS = read_graph(FLAGS_filename, FLAGS_binary, FLAGS_format);
  auto all_tris = brute_force_sampler(GS, FLAGS_k);
  for (const auto &threshold : thresholds) {
    cerr << "Threshold: " << threshold << endl;
    auto heavy_light_tri = heavy_light_sampler(GS, FLAGS_k, threshold);
    compare_statistics(all_tris, heavy_light_tri, FLAGS_k, true);
    cerr << endl;
  }

}
