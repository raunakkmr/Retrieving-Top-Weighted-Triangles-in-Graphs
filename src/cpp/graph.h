#ifndef GRAPH_H
#define GRAPH_H

#include <bits/stdc++.h>
using namespace std;

namespace wsdm_2019_graph {

struct full_edge {
  int src, dst, wt;
  bool operator<(const full_edge& o) const {
    if (o.wt != wt) return wt < o.wt;
    return make_pair(src, dst) < make_pair(o.src, o.dst);
  }
};

struct half_edge {
  int dst, wt;
};

typedef vector<vector<half_edge>> Graph;
typedef map<int, vector<half_edge>> MappedGraph;

struct Shadow {
  MappedGraph g;
  int l;
  bool operator<(const Shadow& o) const {
    return l < o.l;
  }
};

template<class InputGraph>
inline MappedGraph create_shadow(int c, InputGraph& g) {
  map<int, bool> inside;
  for (auto v : g[c]) {
    inside[v.dst] = true;
  }

  MappedGraph ret;
  for (auto v : g[c]) {
    for (auto w : g[v.dst]) {
      if (!inside[w.dst]) continue;
      ret[v.dst].push_back(w);
      ret[w.dst];
    }
  }
  return ret;
}

struct degeneracy_info {
  map<int, int> degenFreq;
  vector<int> degenOrder;
};

degeneracy_info compute_degeneracy(Graph& G) {
  const int n = G.size();

  priority_queue<pair<int, int>> pq;
  vector<bool> done(n);
  vector<int> degree(n);
  for (int i = 0; i < n; i++) {
    degree[i] = G[i].size();
    pq.push(make_pair(-degree[i], i));
  }

  map<int, int> dfreq;
  int degeneracy = 0, cnt = 0;
  vector<int> degenOrder;
  while (!pq.empty()) {
    auto p = pq.top();
    pq.pop();
    int u = p.second;
    if (done[u]) continue;

    degenOrder.push_back(u);
    cnt++;
    degeneracy = max(degeneracy, -p.first);
    dfreq[-p.first]++;
    done[u] = true;
    for (auto v : G[u]) {
      if (done[v.dst]) continue;
      degree[v.dst]--;
      pq.push(make_pair(-degree[v.dst], v.dst));
    }
  }

  degeneracy_info retval;
  retval.degenFreq = dfreq;
  retval.degenOrder = degenOrder;
  return retval;
}

}

#endif /* GRAPH_H */