#include <bits/stdc++.h>

#include "turan.h"
#include "graph.h"
#include "clique_sampler.h"

using namespace std;
using namespace wsdm_2019_graph;

#define PRINT_STATISTICS 1

int main(int argc, char* argv[]) {
  ios::sync_with_stdio(0);
  cin.tie(0);
  srand(time(0)); 

  auto G = read_graph(argv[1]);
  int CLIQUE_SIZE = atoi(argv[2]);
  int NUM_SAMPLES_EDGE = atoi(argv[3]);
  int NUM_SAMPLES_PATH = atoi(argv[4]);
  int NUM_SAMPLES_CLIQUE = atoi(argv[5]);

#if PRINT_STATISTICS
  cerr << "=============================================" << endl;
  cerr << "Computing graph statistics" << endl;
  cerr << "=============================================" << endl;
  cerr << "Computing degeneracy..." << endl;
  auto retval = compute_degeneracy(G);
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
#endif

  auto edge_sampling_tri = edge_sampler(G, NUM_SAMPLES_EDGE);
  auto path_sampling_tri = path_sampler(G, NUM_SAMPLES_PATH);
  
  if (G.size() < 10000) {
    auto all_tris = brute_force_sampler(G);
    cerr << "*** Comparing edge sampling ***" << endl;
    compare_statistics(all_tris, edge_sampling_tri);
    cerr << "*** Comparing path sampling ***" << endl;
    compare_statistics(all_tris, path_sampling_tri);
  } else {
    cerr << "Skipping brute force triangle computations since they are infeasible..." << endl; 
  }
  
  if (CLIQUE_SIZE > 1) {
    auto sampled_cliques = clique_sampler(G, CLIQUE_SIZE, NUM_SAMPLES_CLIQUE);
  }
  return 0;
}