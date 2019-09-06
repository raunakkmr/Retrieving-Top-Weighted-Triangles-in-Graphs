// This file writes out the edge weights and their counts to file.

#include <bits/stdc++.h>

#include "gflags/gflags.h"

#include "graph.h"

using namespace std;
using namespace wsdm_2019_graph;

DEFINE_string(filename, "", "Path to graph file.");
DEFINE_bool(binary, true, "Flag for if graph is in our binary format.");
DEFINE_string(format, "", "If binary is false this is required and indicates format of graph file. One of weighted, temporal, simplicial.");
DEFINE_string(out_path, "", "Path to output file.");

int main(int argc, char *argv[]) {
  string usage("Print edge weights of a graph.\n"
      "Sample usage:\n"
      "\t./print_edge_weights -filename=[fill_this_in] "
      "-binary=true \n"
      "Additionally, these flags can be loaded from a single file "
      "with the option -flagfile=[filename].");

  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_filename.empty()) {
    std::cerr << "No file specified! Type ./compare_deterministic --help for a description of the program parameters." << std::endl;
    return 0;
  }

  ios_base::sync_with_stdio(0);
  cin.tie(0);

  auto GS = read_graph(FLAGS_filename, FLAGS_binary, FLAGS_format);
  string out_path = FLAGS_out_path;

  const auto &edges = GS.edges;
  map<long long, int> count;
  for (const auto &e : edges) {
  count[e.wt]++;
  }
  ofstream ofile(out_path, ios::out);
  for (const auto &kv : count) {
    ofile << kv.first << " " << kv.second << endl;
  }
  ofile.close();
  /*
     auto mle = [&]() {
     double xmin = count.begin()->first;
     cerr << "xmin: " << xmin << endl;
     double n = edges.size();
     double s = 0;
     for (const auto &e : edges) {
     s += log((e.wt*1.0) / (xmin));
     }
     double alpha = 1 + n / s;
     return alpha;
     };
     double alpha = mle();
     cerr << "alpha: " << alpha << endl;
   */
}
