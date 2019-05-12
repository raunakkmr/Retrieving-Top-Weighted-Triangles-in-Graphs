#ifndef GRAPH_H
#define GRAPH_H

#include <bits/stdc++.h>
using namespace std;

namespace wsdm_2019_graph {

// generic data structures for graphs
struct full_edge {
  int src, dst, wt;
  bool operator<(const full_edge& o) const {
    if (o.wt != wt) return wt < o.wt;
    return make_pair(src, dst) < make_pair(o.src, o.dst);
  }
};

struct half_edge {
  int dst, wt;
  // THIS CANNOT BE CHANGED SINCE IT IS USED BY PATH SAMPLER
  const bool operator<(const half_edge& o) const {
    if (dst == o.dst) return wt > o.wt;
    return dst < o.dst;
  }
};

struct weighted_triangle {
  tuple<int, int, int> vertices;
  int weight;
  weighted_triangle(int u, int v, int w, int wt) {
    int verts[] = {u, v, w};
    sort(verts, verts+3);
    vertices = make_tuple(verts[0], verts[1], verts[2]);
    weight = wt;
  }

  const bool operator<(const weighted_triangle& o) const {
    if (weight != o.weight) return weight > o.weight;
    return vertices < o.vertices;
  }

  friend ostream& operator<<(ostream& stream, const weighted_triangle& t) {
    stream << '(' << get<0>(t.vertices) << ", " << get<1>(t.vertices) << ", " << get<2>(t.vertices) << ", " << t.weight << ")";
    return stream;
  }
};

struct weighted_clique {
  vector<int> vertices;
  int weight;
  weighted_clique(vector<int> V, int wt) {
    sort(V.begin(), V.end());
    vertices = V;
    weight = wt;
  }
  weighted_clique() {
    weight = 0;
  };

  const bool operator<(const weighted_clique& o) const {
    if (weight != o.weight) return weight > o.weight;
    return vertices < o.vertices;
  }

  friend ostream& operator<<(ostream& stream, const weighted_clique& t) {
    stream << '(';
    for (const int u : t.vertices) {
      stream << u << ", ";
    }
    stream << t.weight << ")";
    return stream;
  }
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

// Computes the induced subgraph of neighbours of c
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

// computes the degeneracy ordering of a graph G
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

// filename contains the path and common prefix of the datafiles
Graph read_graph(string filename) {
  ifstream simplices(filename + "-simplices.txt", ifstream::in);
  ifstream nverts(filename + "-nverts.txt", ifstream::in);
  ifstream times(filename + "-times.txt", ifstream::in);

  cerr << "reading in graph " << filename << endl;

  // no use for the times right now
  map<int, map<int, int>> weight;
  map<int, int> label;
  int ns, nnodes = 0;
  while (nverts >> ns) {
    vector<int> simplex(ns);
    for (int i = 0; i < ns; i++) {
      simplices >> simplex[i];
      if (!label.count(simplex[i])) {
        label[simplex[i]] = nnodes++;
      }
    }

    for (int i = 0; i < ns; i++) {
      for (int j = i+1; j < ns; j++) {
        weight[simplex[i]][simplex[j]]++;
        weight[simplex[j]][simplex[i]]++;
      }
    }
  }

  int nedges = 0;
  Graph G(nnodes);
  for (auto& e0 : weight) {
    for (auto& e1 : e0.second) {
      int u = label[e0.first], v = label[e1.first], w = e1.second;
      G[u].push_back({v, w});
      nedges++;
    }
  }

  cerr << "read in a graph with " << nnodes << " " << nedges / 2 << " edges" << endl;
  cerr << "Average degree: " << nedges / nnodes << endl;
  return G;
}

}

#endif /* GRAPH_H */