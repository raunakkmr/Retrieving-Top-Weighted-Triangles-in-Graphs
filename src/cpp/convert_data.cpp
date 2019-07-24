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

  ofstream out_file(dataset_path+".binary", ios::binary | ios::out);

  binary_write(out_file, n);
  binary_write(out_file, m);

  auto get_bytes = [](int x) {
    if (x <= numeric_limits<char>::max()) {
      return 0;
    } else if (x <= numeric_limits<short>::max()) {
      return 1;
    } else if (x <= numeric_limits<int>::max()) {
      return 3;
    }
    return -1;
  };

  int est_bytes = 8;
  for (int u = 0; u < (int) G.size(); u++) {
    for (auto &e : G[u]) {
      if (u < e.dst) {
        int bytes = 0;
        // uses 2 as a special code that the weight is 1.
        // call it premature optimization but it saves a lot
        bytes |= (e.wt == 1 ? 2 : get_bytes(e.wt));
        bytes <<= 2;
        bytes |= get_bytes(e.dst);
        bytes <<= 2;
        bytes |= get_bytes(u);

        binary_compressed_write(out_file, bytes);
        binary_compressed_write(out_file, u);
        binary_compressed_write(out_file, e.dst);
        est_bytes += 3 + get_bytes(u) + get_bytes(e.dst);
        if (e.wt > 1) {
          binary_compressed_write(out_file, e.wt);
          est_bytes += 1 + get_bytes(e.wt);
        }
      }
    }
  }

  cerr << "estimated " << est_bytes << " bytes written" << endl;
  out_file.close();

  return 0;
}