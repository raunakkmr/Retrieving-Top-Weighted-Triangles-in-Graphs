#include <bits/stdc++.h>

#include "graph.h"

using namespace std;
using namespace wsdm_2019_graph;

int main(int argc, char* argv[]) {
  ios::sync_with_stdio(0);
  cin.tie(0);
  srand(0); 

  auto G = read_graph(argv[1]);
  string dataset_path = argv[2];

  int n = G.size();
  int m = 0;
  for (const auto &u : G) {
      m += u.size();
  }
  m /= 2;

  ofstream out_file(dataset_path+".binary", ios::binary);

  binary_write(out_file, n);
  binary_write(out_file, m);

  for (int u = 0; u < (int) G.size(); u++) {
      for (const auto &e : G[u]) {
          if (u < e.dst) {
            int a = u, b = e.dst;
            int c = e.wt;
            binary_write(out_file, a);
            binary_write(out_file, b);
            binary_write(out_file, c);
          }
      }
  }

  out_file.close();

  return 0;
}