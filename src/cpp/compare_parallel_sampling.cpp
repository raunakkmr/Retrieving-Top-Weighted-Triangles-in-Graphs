// This file runs edge, wedge or path sampling in parallel for specified amounts
// of time.

#include <bits/stdc++.h>

#include "gflags/gflags.h"

#include "graph.h"
#include "triangle_sampler.h"

using namespace std;
using namespace wsdm_2019_graph;

#define PRINT_ARGS 1

DEFINE_string(filename, "", "Path to graph file.");
DEFINE_bool(binary, true, "Flag for if graph is in our binary format.");
DEFINE_string(format, "", "If binary is false this indicates format of graph file. One of weighted, temporal, simplicial. Required if binary is false.");
DEFINE_int32(k, 1000, "Parameter k for top-k.");
DEFINE_string(sampler, "edge", "Which sampling algorithm to use. One of edge, wedge or path.");
DEFINE_double(start_time, 0, "How long to run the algorithm for (in seconds).");
DEFINE_double(end_time, 0, "When to terminate (in seconds).");
DEFINE_double(increment, 0, "Time increments (in seconds).");

int main(int argc, char* argv[]) {
  string usage("Parallel sampling for triangles.\n"
      "Sample usage:\n"
      "\t./compare_parallel_sampling -filename=[fill_this_in] "
      "-binary=true -k=[fill_this_in] -sampler=[fill_this_in] "
      "-start_time=[fill_this_in] -end_time=[fill_this_in] -increment=[fill_this_in]\n"
      "For example, if start_time, end_time, increment are 2, 6 and 2 respectively "
      "then the chosen algorithm is run 3 times - for 2, 4, and 6 seconds.\n"
      "Additionally, these flags can be loaded from a single file "
      "with the option -flagfile=[filename].");

  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_filename.empty()) {
    std::cerr << "No file specified! Type ./compare_parallel_sampling --help for a description of the program parameters." << std::endl;
    return 0;
  }

  ios::sync_with_stdio(0);
  cin.tie(0);
  srand(0); 

  auto GS = read_graph(FLAGS_filename, FLAGS_binary, FLAGS_format);
  int K = FLAGS_k;
  string sampler = FLAGS_sampler;
  double start_time = FLAGS_start_time;
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
  cerr << "Sampler: " << sampler << endl;
  cerr << "Start time: " << start_time << endl;
  cerr << "End time: " << end_time << endl;
  cerr << "Increment: " << increment << endl;
  cerr << endl;
#endif

  double cur_time = start_time;
  vector<double> times;
  vector<weighted_triangle> sampling_tri;
  auto all_tris = dynamic_heavy_light(GS, K);
  while (cur_time <= end_time) {
    times.push_back(cur_time);
    cerr << "Runnning for sampler for " << cur_time << " (s)" << endl;
    if (sampler.find("edge") != string::npos) {
      sampling_tri = (edge_parallel_time_version(GS, thread::hardware_concurrency(), cur_time, -1, false));
    } else if (sampler.find("wedge") != string::npos) {
      sampling_tri = (wedge_parallel_time_version(GS, thread::hardware_concurrency(), cur_time, -1, false));
    } else if (sampler.find("path") != string::npos) {
      // TODO: Implement path_parallel_time_version.
    } else {
      cerr << "Unrecognized sampler." << endl;
    }
    cerr << "*** Comparing parallel sampling ***" << endl;
    compare_statistics(all_tris, sampling_tri, K, true);
    cur_time += increment;
  }

  return 0;
}
