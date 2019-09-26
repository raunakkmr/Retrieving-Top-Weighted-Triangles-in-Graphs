#include <bits/stdc++.h>

#include "graph.h"
#include "clique_sampler.h"
#include "triangle_sampler.h"
#include "gflags/gflags.h"

using namespace std;
using namespace wsdm_2019_graph;

#define PRINT_STATISTICS 0

DEFINE_string(filename, "", "Path to graph file.");
DEFINE_bool(binary, true, "Flag for if graph is in our binary format.");
DEFINE_string(format, "", "If binary is false this indicates format of graph file. One of weighted, temporal, simplicial. Required if binary is false.");
DEFINE_int32(clique_size, 4, "The size of the clique.");
DEFINE_double(start_time, 1, "How long to run the algorithm for (in seconds).");
DEFINE_double(end_time, 1, "When to terminate (in seconds).");
DEFINE_double(increment, 0, "Time increments (in seconds).");
DEFINE_int32(k, 1000, "Parameter k for top-k.");
DEFINE_int32(runs, 10, "Number of runs.");
DEFINE_bool(print_statistics, false, "Prints extra debug statistics about the graph.");

int main(int argc, char* argv[]) {
  std::string usage("Edge based clique sampler.\n"
      "Sample usage:\n"
      "\t./clique_sampler -filename=[fill_this_in] "
      "-binary=true -clique_size=[fill_this_in] "
      "-start_time=[fill_this_in] -end_time=[fill_this_in] -increment=[fill_this_in]\n"
      "For example, if start_time, end_time, increment are 2, 6 and 2 respectively "
      "then the chosen algorithm is run 3 times - for 2, 4, and 6 seconds.\n"
      "Additionally, these flags can be loaded from a single file "
      "with the option -flagfile=[filename].");

  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_filename.empty()) {
    std::cerr << "No file specified! Type ./clique_sampler --help for a description of the program parameters." << std::endl;
    return 0;
  }

  ios::sync_with_stdio(0);
  cin.tie(0);
  srand(0); 

  auto GS = read_graph(FLAGS_filename, FLAGS_binary, FLAGS_format);
  int CLIQUE_SIZE = FLAGS_clique_size;
  double start_time = FLAGS_start_time;
  double end_time = FLAGS_end_time;
  double increment = FLAGS_increment;

  cerr << "=============================================" << endl;
  cerr << "ARGUMENTS" << endl;
  cerr << "=============================================" << endl;
  cerr << "Dataset: " << FLAGS_filename << endl;
  cerr << "Binary: " << FLAGS_binary << endl;
  cerr << "Format: " << FLAGS_format << endl;
  cerr << "CLIQUE_SIZE: " << CLIQUE_SIZE << endl;
  cerr << "Start time: " << start_time << endl;
  cerr << "End time: " << end_time << endl;
  cerr << "Increment: " << increment << endl;

  if (FLAGS_print_statistics) {
    cerr << "=============================================" << endl;
    cerr << "Computing graph statistics" << endl;
    cerr << "=============================================" << endl;
    cerr << "Computing degeneracy..." << endl;
    auto retval = compute_degeneracy(GS.G);
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
  }

  int nthreads = thread::hardware_concurrency();

  if (CLIQUE_SIZE > 1) {
    auto all_cliques = clique_brute_force(GS.G, CLIQUE_SIZE);
    for (int run = 0; run < FLAGS_runs; run++) {
      cerr << "Run: " << run << endl;
      double cur_time = start_time;
      vector<double> times;
      vector<weighted_clique> sampling_cliques;
      while (cur_time <= end_time) {
        times.push_back(cur_time);
        cerr << "Runnning for sampler for " << cur_time << " (s)" << endl;
        // sampling_cliques = clique_sampler_parallel(GS, CLIQUE_SIZE, 1000, nthreads);
        sampling_cliques = clique_sampler_tmp(GS, CLIQUE_SIZE, nthreads, 1000);
        cerr << "*** Comparing parallel sampling ***" << endl;
        compare_statistics(all_cliques, sampling_cliques, CLIQUE_SIZE);
        cur_time += increment;
      }
    }
  }

  return 0;
}
