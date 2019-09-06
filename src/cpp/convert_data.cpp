// This file converts the graph file to a specialized compressed binary format
// which is significantly smaller and faster to read.

#include <bits/stdc++.h>

#include "gflags/gflags.h"

#include "graph.h"

using namespace std;
using namespace wsdm_2019_graph;

DEFINE_string(filename, "", "Path to graph file.");
DEFINE_string(format, "", "If binary is false this is required and indicates format of graph file. One of weighted, temporal, simplicial.");
DEFINE_string(binary_path, "", "Path to output file.")

int main(int argc, char* argv[]) {
  string usage("Convert graph file to binary format.\n"
      "Sample usage:\n"
      "\t./convert_data -filename=[fill_this_in] -format=[fill_this_in] -binary_path=[fill_this_in]\n"
      "Additionally, these flags can be loaded from a single file "
      "with the option -flagfile=[filename].");

  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

	if (FLAGS_filename.empty()) {
		std::cerr << "No file specified! Type ./convert_data --help for a description of the program parameters." << std::endl;
		return 0;
	}

  ios::sync_with_stdio(0);
  cin.tie(0);
  srand(0);

  auto GS = read_graph(FLAGS_filename, false, FLAGS_format);
  string binary_path = FLAGS_binary_path;

  Graph G = GS.G;

  // If a node u has high degree then it will appear in the edge list a large
  // number of times. So assign it a small integer label that is smaller than 4
  // bytes. Since the binary reader and writer do not simply use 4 bytes for
  // every integer value but use a smaller number of bytes if possible, this
  // leads to significant speedup.
  int n = G.size();
  long long m = 0;
  vector<pair<int, int>> sort_by_deg;
  for (int u = 0; u < (int) G.size(); u++) {
    if (G[u].empty()) continue;
    m += G[u].size();
    sort_by_deg.push_back(make_pair((int) G[u].size(), u));
  }
  m /= 2;
  sort(sort_by_deg.rbegin(), sort_by_deg.rend());
  map<int, int> label;
  for (int i = 0; i < (int) sort_by_deg.size(); i++) {
    label[sort_by_deg[i].second] = i;
  }

  // Prune out degree 0 nodes
  n = label.size();

  ofstream out_file(binary_path+".binary", ios::binary | ios::out);

  binary_write(out_file, n);
  int m_int = (int) m;
  binary_write(out_file, m_int);

  size_t est_bytes = 8;
  for (int u = 0; u < (int) G.size(); u++) {
    int x = label[u];
    for (auto &e : G[u]) {
      int y = label[e.dst];
      if (x < y) {
        int bytes = 0;
        // Use last 4 bits to indicate whether edge weight is small or not. This
        // might be prematrure optimization but it saves a lot since a vast
        // majority of weights are extremely small.
        bytes |= (e.wt <= 12 ? e.wt-1 : 12 + get_bytes(e.wt));
        bytes <<= 2;
        // How many bytes are needed to represent y and x.
        bytes |= get_bytes(y);
        bytes <<= 2;
        bytes |= get_bytes(x);

        // Write out the number of bytes needed to represent the nodes and edge
        // weight, and then the nodes and weight themselves.
        binary_compressed_write(out_file, bytes);
        binary_compressed_write(out_file, x);
        binary_compressed_write(out_file, y);
        est_bytes += 3 + get_bytes(x) + get_bytes(y);
        if (e.wt > 12) {
          binary_compressed_write(out_file, e.wt);
          est_bytes += 1 + get_bytes(e.wt);
        }
      }
    }
  }

  cerr << "Estimated " << est_bytes << " bytes written." << endl;
  out_file.close();

  return 0;
}
