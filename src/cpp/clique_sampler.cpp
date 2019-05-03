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
  int NUM_SAMPLES = atoi(argv[3]);

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

  auto sampled_tris = edge_sampler(G, NUM_SAMPLES);
  
  if (G.size() < 10000) {
    auto all_tris = brute_force_sampler(G);
    compare_statistics(all_tris, sampled_tris);
  } else {
    cerr << "Skipping brute force triangle computations since they are infeasible..." << endl; 
  }
  
  auto sampled_cliques = clique_sampler(G, CLIQUE_SIZE, NUM_SAMPLES);

  return 0;
}