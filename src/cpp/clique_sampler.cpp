#include <bits/stdc++.h>

#include "turan.h"
#include "graph.h"
#include "clique_sampler.h"
#include "triangle_sampler.h"
#include "gflags/gflags.h"

using namespace std;
using namespace wsdm_2019_graph;

#define PRINT_STATISTICS 0

DEFINE_string(filename, "", "Path to graph file.");
DEFINE_bool(binary, true, "Flag for if graph is in our binary format.");
DEFINE_int32(clique_size, 4, "The size of the clique.");
DEFINE_int32(nsamples, 0, "Number of samples for clique edge sampler.");
DEFINE_bool(print_statistics, false, "Prints extra debug statistics about the graph.");

int main(int argc, char* argv[]) {
  std::string usage("Edge based clique sampler.\n"
      "Sample usage:\n"
      "\t./clique_sampler -filename=[fill_this_in] "
      "-binary=true -clique_size=4 -nsamples=100 \n"
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

  auto GS = read_graph(FLAGS_filename, true);
  int CLIQUE_SIZE = FLAGS_clique_size;
  int NUM_SAMPLES_CLIQUE = FLAGS_nsamples;

  cerr << "=============================================" << endl;
  cerr << "ARGUMENTS" << endl;
  cerr << "=============================================" << endl;
  cerr << "Dataset: " << FLAGS_filename << endl;
  cerr << "CLIQUE_SIZE, NSAMPS: " << CLIQUE_SIZE << " " << NUM_SAMPLES_CLIQUE << endl;

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

  // auto edge_sampling_tri = edge_samples_version(GS, NUM_SAMPLES_EDGE);
  // auto edge_sampling_tri_parallel = edge_parallel_samples_version(GS, nthreads, NUM_SAMPLES_EDGE);
  // auto wedge_sampling_tri = wedge_sampler(GS, NUM_SAMPLES_WEDGE);
  // auto path_sampling_tri = path_sampler(GS, NUM_SAMPLES_PATH);
  // auto heavy_light_sampling_tri = heavy_light_sampler(GS, 0.05);

  if (CLIQUE_SIZE > 1) {
    // auto sampled_cliques = clique_sampler(GS, CLIQUE_SIZE, NUM_SAMPLES_CLIQUE);
    auto sampled_cliques_parallel = clique_sampler_parallel(GS, CLIQUE_SIZE, NUM_SAMPLES_CLIQUE, nthreads);
    auto all_cliques = clique_brute_force(GS.G, CLIQUE_SIZE);
  }
  return 0;
}
