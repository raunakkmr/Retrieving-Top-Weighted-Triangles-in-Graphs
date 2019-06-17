#ifndef GRAPH_H
#define GRAPH_H

#include <bits/stdc++.h>
using namespace std;

namespace wsdm_2019_graph {

long long rand64() {
  return rand() * (1LL << 32) + rand();
}

// generic data structures for graphs
struct full_edge {
  int src, dst;
  long long wt;
  bool operator<(const full_edge& o) const {
    if (o.wt != wt) return wt < o.wt;
    return make_pair(src, dst) < make_pair(o.src, o.dst);
  }
};

struct half_edge {
  int dst;
  long long wt;
  // THIS CANNOT BE CHANGED SINCE IT IS USED BY PATH SAMPLER
  const bool operator<(const half_edge& o) const {
    if (dst == o.dst) return wt > o.wt;
    return dst < o.dst;
  }
};

struct weighted_triangle {
  tuple<int, int, int> vertices;
  long long weight;
  weighted_triangle(int u, int v, int w, long long wt) {
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
  long long weight;
  weighted_clique(vector<int> V, long long wt) {
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
  map<int, map<int, long long>> weight;
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
      int u = label[e0.first], v = label[e1.first];
      long long w = e1.second;
      G[u].push_back({v, w});
      nedges++;
    }
  }

  cerr << "read in a graph with " << nnodes << " nodes and " << nedges / 2 << " edges" << endl;
  cerr << "Average degree: " << nedges / nnodes << endl;
  return G;
}

// Modify the weights of the graph for p-means.
// If p=0, then w -> log(w). Otherwise, w -> w^p.
// However, a lot of the existing code assumes integer edge weights and this may
// no longer be true for certain values of p, such as 0, -1, 1.5, etc. So we
// scale the weights by multiplying it with a large constant and rounding them
// to an integer. The distribution of the weights does not change. If the
// maximum edge weight was too large, then prints a warning.
void modify_weights(Graph &G, double p=1.0) {
  // if ((long long) p == 1) {
  //   return;
  // }
  long long max_weight_edge = -1;
  long long num_edges = 0;
  for (int u = 0; u < (int) G.size(); u++) {
    for (auto &e : G[u]) {
      max_weight_edge = max(max_weight_edge, e.wt);
      num_edges++;
    }
  }
  num_edges /= 2;

  long long factor = 100000000;

  if (max_weight_edge > numeric_limits<long long>::max() / factor) {
    cerr << "WARNING: max_weight_edge is too large to truncate!!!" << endl;
  }

  long long large_constant = numeric_limits<long long>::max() / (num_edges*factor*max_weight_edge);

  cerr << max_weight_edge << " " << numeric_limits<long long>::max() << " " << (num_edges*factor*max_weight_edge) << endl;
  cerr << large_constant << endl;

  if (abs(p) < 1e-8) {
    for (int u = 0; u < (int) G.size(); u++) {
      for (auto &e : G[u]) {
        e.wt = (long long) (log(e.wt) * large_constant);
      }
    }
  } else {
    for (int u = 0; u < (int) G.size(); u++) {
      for (auto &e : G[u]) {
        e.wt = (long long) (pow(e.wt, p) * large_constant);
      }
    }
  }
}

}

#endif /* GRAPH_H */