#include <bits/stdc++.h>
#include "turan.h"
#include "graph.h"

using namespace std;
using namespace wsdm_2019_graph;

int main(int argc, char* argv[]) {
  ios::sync_with_stdio(0);
  cin.tie(0);
  srand(time(0)); 

  auto res = freopen(argv[1], "r", stdin);
  if (res == nullptr) return 1;

  int CLIQUE_SIZE = atoi(argv[2]);

  int n, m;
  cin >> n >> m;

  Graph G(n), Gu(n);
  vector<int> ordering;
  ordering.resize(n);
  for (int i = 0; i < n; i++) {
    ordering[i] = i;
  }
  random_shuffle(ordering.begin(), ordering.end());

  int u, v;
  while (cin >> u >> v) {
    if (ordering[u] < ordering[v]) {
      G[u].push_back({v, 1});
    } else {
      G[v].push_back({u, 1});
    }
    Gu[u].push_back({v, 1});
    Gu[v].push_back({u, 1});
  }

  cerr << "=============================================" << endl;
  cerr << "Computing graph statistics" << endl;
  cerr << "=============================================" << endl;
  cerr << "Computing degeneracy..." << endl;
  auto retval = compute_degeneracy(Gu);
  auto degenOrder = retval.degenOrder;
  auto degeneracy = (--retval.degenFreq.end())->first;

  cerr << "Degeneracy distribution:" << endl;
  double sum = 0, num = 0;
  for (auto kv : retval.degenFreq) {
    cerr << kv.first << " " << kv.second << endl;
    sum += kv.first * kv.second;
    num += kv.second;
  }
  cerr << "Average degeneracy: " << sum / num << endl;

  cerr << "Number of vertice: " << n << endl;
  cerr << "Number of edges: " << m << endl;
  cerr << "Average degree: " << 2.0 * m / n << endl;
  cerr << "Degeneracy: " << degeneracy << endl;

  //TuranShadow(G, CLIQUE_SIZE);
  //TuranBK(Gu, degenOrder, CLIQUE_SIZE);
  //TuranBK(Gu, ordering, CLIQUE_SIZE);
  reverse(degenOrder.begin(), degenOrder.end());
  TuranBK(Gu, degenOrder, CLIQUE_SIZE);
  //auto whatAmIDoing = degenOrder;
  //random_shuffle(whatAmIDoing.begin() + n / 2, whatAmIDoing.end());
  //TuranBK(Gu, whatAmIDoing, CLIQUE_SIZE);
  return 0;
}