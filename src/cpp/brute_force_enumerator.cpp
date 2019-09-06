// This file runs the brute force algorithm to enumerate triangles.

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
DEFINE_string(out_path, "", "Path to output file.");

int main(int argc, char* argv[]) {
  string usage("Brute force enumerator for triangles.\n"
      "Sample usage:\n"
      "\t./brute_force_enumerator -filename=[fill_this_in] "
      "-binary=true -k=[fill_this_in] -out_path=[fill_this_in]\n"
      "Additionally, these flags can be loaded from a single file "
      "with the option -flagfile=[filename].");

  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_filename.empty()) {
    std::cerr << "No file specified! Type ./brute_force_enumerator --help for a description of the program parameters." << std::endl;
    return 0;
  }

  ios_base::sync_with_stdio(0);
  cin.tie(0);
  srand(0);

  auto GS = read_graph(FLAGS_filename, FLAGS_binary, FLAGS_format);
  int K = FLAGS_k;
  string out_path = FLAGS_out_path;

  auto all_tris = brute_force_sampler(GS, K);
  long long num_all_tris = all_tris.size();

  ofstream triangle_file(out_path+"-triangles.binary", ios::binary | ios::out);
  binary_write(triangle_file, num_all_tris);
  int est_bytes = 8;
  map<long long, long long> num_tris;
  for (auto t : all_tris) {
    unsigned int w = t.weight;
    int a = get<0>(t.vertices);
    int b = get<1>(t.vertices);
    int c = get<2>(t.vertices);
    int bytes = 0;
    // Use 2 as a special code that the weight is 3.
    bytes |= (w == 3 ? 2 : get_bytes(w));
    bytes <<= 2;
    bytes |= get_bytes(a);
    bytes <<= 2;
    bytes |= get_bytes(b);
    bytes <<= 2;
    bytes |= get_bytes(c);

    binary_compressed_write(triangle_file, bytes);
    binary_compressed_write(triangle_file, a);
    binary_compressed_write(triangle_file, b);
    binary_compressed_write(triangle_file, c);
    est_bytes += 4 + get_bytes(a) + get_bytes(b) + get_bytes(c);
    if (w > 3) {
      binary_compressed_write(triangle_file, w);
      est_bytes += 1 + get_bytes(w);
    }

    num_tris[get<0>(t.vertices)]++;
    num_tris[get<1>(t.vertices)]++;
    num_tris[get<2>(t.vertices)]++;
  }
  triangle_file.close();

  Graph &G = GS.G;

  ofstream tris_to_weight_file(out_path+"-ntris-to-weight.binary", ios::binary);
  int g_sz = G.size();
  binary_write(tris_to_weight_file, g_sz);
  for (int i = 0; i < (int) G.size(); i++) {
    double weight_sum = 0;
    for (auto e : G[i]) {
      weight_sum += e.wt;
    }
    binary_write(tris_to_weight_file, weight_sum);
    binary_write(tris_to_weight_file, num_tris[i]);
  }
  tris_to_weight_file.close();

  ofstream degree_file(out_path+"-degree-to-weight.binary", ios::binary);
  binary_write(degree_file, g_sz);
  for (int i = 0; i < (int) G.size(); i++) {
    double weight_sum = 0;
    for (auto e : G[i]) {
      weight_sum += e.wt;
    }
    int gi_sz = G[i].size();
    binary_write(degree_file, gi_sz);
    binary_write(degree_file, weight_sum);
  }
  degree_file.close();

  return 0;
}
