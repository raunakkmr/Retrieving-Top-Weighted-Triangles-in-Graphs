// This file runs edge, wedge and path sampling (sequential, not parallel) for
// specified amounts of time.

#include <bits/stdc++.h>

#include "gflags/gflags.h"

#include "graph.h"
#include "triangle_sampler.h"

using namespace std;
using namespace wsdm_2019_graph;

DEFINE_string(filename, "", "Path to graph file.");
DEFINE_bool(binary, true, "Flag for if graph is in our binary format.");
DEFINE_string(format, "", "If binary is false this indicates format of graph file. One of weighted, temporal, simplicial. Required if binary is false.");
DEFINE_int32(k, 1000, "Parameter k for top-k.");
DEFINE_double(end_time, 0, "When to terminate (in seconds).");
DEFINE_double(increment, 0, "Time increments (in seconds).");

#define PRINT_ARGS 1
#define PRINT_STATISTICS 0

int main(int argc, char* argv[]) {
  string usage("Edge, wedge and path sampling for triangles.\n"
      "Sample usage:\n"
      "\t./compare_sampling -filename=[fill_this_in] "
      "-binary=true -k=[fill_this_in] "
      "-start_time=[fill_this_in] -end_time=[fill_this_in] -increment=[fill_this_in]\n"
      "For example, if end_time and increment are 6 and 2 respectively "
      "then the algorithms are run 3 times - for 2, 4, and 6 seconds.\n"
      "Additionally, these flags can be loaded from a single file "
      "with the option -flagfile=[filename].");

  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_filename.empty()) {
    std::cerr << "No file specified! Type ./compare_sampling --help for a description of the program parameters." << std::endl;
    return 0;
  }

  ios::sync_with_stdio(0);
  cin.tie(0);
  srand(0); 

  auto GS = read_graph(FLAGS_filename, FLAGS_binary, FLAGS_format);
  int K = FLAGS_k;
  double end_time = FLAGS_end_time;
  double increment = FLAGS_increment;

#if PRINT_ARGS
  cerr << "=============================================" << endl;
  cerr << "ARGUMENTS" << endl;
  cerr << "=============================================" << endl;
  cerr << "Dataset: " << FLAGS_filename << endl;
  cerr << "Binary: " << FLAGS_binary << endl;
  cerr << "Format: " << FLAGS_format << endl;
  cerr << "K: " << K << endl;
  cerr << "End time: " << end_time << endl;
  cerr << "Increment: " << increment << endl;
  cerr << endl;
#endif

  double cur_time = increment;
  vector<double> times;
  set<weighted_triangle> edge_sampling_tri, wedge_sampling_tri, path_sampling_tri;
  auto all_tris = dynamic_heavy_light(GS, K);
  while (cur_time <= end_time) {
    times.push_back(cur_time);
    edge_sampling_tri = (edge_time_version(GS, cur_time, false));
    wedge_sampling_tri = (wedge_time_version(GS, cur_time, false));
    path_sampling_tri = (path_time_version(GS, cur_time, false));
    cerr << "*** Comparing edge sampling ***" << endl;
    compare_statistics(all_tris, edge_sampling_tri, K, true);
    cerr << "*** Comparing wedge sampling ***" << endl;
    compare_statistics(all_tris, wedge_sampling_tri, K, true);
    cerr << "*** Comparing path sampling ***" << endl;
    compare_statistics(all_tris, path_sampling_tri, K, true);
    cur_time += increment;
  }


  #if PRINT_STATISTICS
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
#endif

  return 0;
}
