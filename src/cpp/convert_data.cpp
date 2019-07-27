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

  ofstream out_file(dataset_path+".binary", ios::binary | ios::out);

  binary_write(out_file, n);
  binary_write(out_file, m);

  // /*
  auto get_bytes = [](int x) {
    if (x <= numeric_limits<unsigned char>::max()) {
      return 0;
    } else if (x <= numeric_limits<short>::max()) {
      return 1;
    } else if (x <= numeric_limits<unsigned int>::max()) {
      return 3;
    }
    return -1;
  };
  // */

  size_t est_bytes = 8;
  for (int u = 0; u < (int) G.size(); u++) {
    int x = label[u];
    for (auto &e : G[u]) {
      int y = label[e.dst];
      if (x < y) {
        int bytes = 0;
        // uses last 4 bits to indicate whether its a small weight edge or not.
        // call it premature optimization but it saves a lot
        bytes |= (e.wt <= 12 ? e.wt-1 : 12 + get_bytes(e.wt));
        bytes <<= 2;
        bytes |= get_bytes(y);
        bytes <<= 2;
        bytes |= get_bytes(x);

        binary_compressed_write(out_file, bytes);
        binary_compressed_write(out_file, x);
        binary_compressed_write(out_file, y);
        est_bytes += 3 + get_bytes(x) + get_bytes(y);
        if (e.wt > 12) {
          binary_compressed_write(out_file, e.wt);
          est_bytes += 1 + get_bytes(e.wt);
        }
        // cerr << x << " " << y << " " << e.wt << '\n';
      }
    }
  }

  cerr << "estimated " << est_bytes << " bytes written" << endl;
  out_file.close();

  return 0;
}
