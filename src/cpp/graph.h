#ifndef GRAPH_H
#define GRAPH_H

#include <bits/stdc++.h>
using namespace std;

namespace wsdm_2019_graph {

// From https://baptiste-wicht.com/posts/2011/06/write-and-read-binary-files-in-c.html
inline int binary_read(std::istream& stream, int nbytes){
  if (nbytes == 0) return 1;

  int res = 0;
  if (nbytes == 1) {
    char c;
    stream.read(reinterpret_cast<char*>(&c), nbytes);
    res = c;
  } else if (nbytes == 2) {
    short c;
    stream.read(reinterpret_cast<char*>(&c), nbytes);
    res = c;
  } else if (nbytes == 4) {
    int c;
    stream.read(reinterpret_cast<char*>(&c), nbytes);
    res = c;
  } else {
    throw "not yet implemented";
  }
  return res;
}

template<typename T>
void binary_write(std::ostream& stream, T value){
  stream.write(reinterpret_cast<char*>(&value), sizeof(T));
}

inline void binary_compressed_write(std::ostream& stream, int x) {
  if (x <= numeric_limits<char>::max()) {
    char value = x;
    stream.write(reinterpret_cast<char*>(&value), 1);
  } else if (x <= numeric_limits<short>::max()) {
    short value = x;
    stream.write(reinterpret_cast<char*>(&value), 2);
  } else if (x <= numeric_limits<int>::max()) {
    int value = x;
    stream.write(reinterpret_cast<char*>(&value), 4);
  } else {
    throw "not yet implemented";
  }
}

long long rand64() {
  return rand() * (1LL << 32) + rand();
}

// generic data structures for graphs
struct full_edge {
  int src, dst;
  long long wt;
  inline bool operator<(const full_edge& o) const {
    if (o.wt != wt) return wt < o.wt;
    return make_pair(src, dst) < make_pair(o.src, o.dst);
  }
};

struct half_edge {
  int dst;
  long long wt;
  // THIS CANNOT BE CHANGED SINCE IT IS USED BY PATH SAMPLER
  inline const bool operator<(const half_edge& o) const {
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

  inline const bool operator<(const weighted_triangle& o) const {
    if (weight != o.weight) return weight > o.weight;
    return vertices < o.vertices;
  }

  inline const bool operator==(const weighted_triangle& o) const {
    return (weight == o.weight && vertices == o.vertices);
  }

  inline const bool operator!=(const weighted_triangle& o) const {
    return (weight != o.weight || vertices != o.vertices);
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

// If binary is true then filename is the path to the .binary file which
// contains the graph in the following form: one line with the number of
// vertices n, one line with the number of edges m, and 3m lines representing
// edges (u, v, w).
// Otherwise, for the simplicial datasets filename contains the path and common
// prefix of the datafiles. For the temporal-reddit-reply dataset filename is
// the path to the temporal-reddit-reply.txt file.
Graph read_graph(string filename, bool binary=false) {

  // no use for the times right now
  map<int, map<int, long long>> weight;
  map<int, int> label;
  int nnodes = 0;
  int nedges = 0;
  int m = 0;

  cerr << "reading in graph " << filename << endl;

  if (binary) {
    const unsigned int blength = 1024 * 1024;
    char buffer[blength];

    ifstream data_file(filename, ios::binary | ios::in);
    data_file.rdbuf()->pubsetbuf(buffer, blength);
    data_file.sync_with_stdio(0);
    data_file.tie(0);

    nnodes = binary_read(data_file, 4);
    m = binary_read(data_file, 4);

    int bytes_read = 8;
    cerr << "nodes and edges: " << nnodes << " " << m << endl;

    //data_file.seekg (0, data_file.end);
    //int length = data_file.tellg();
    //data_file.seekg (0, data_file.beg);
    //cerr << "number of bytes in file: " << length << endl;
    for (int i = 0; i < m; i++) {
      char type = binary_read(data_file, 1);
      bytes_read += 1;

      int bytes[] = {type&3, (type>>2)&3, (type>>4)&3};
      if (bytes[2] == 2) { // special case when weight is 1 (majority of weights)
        bytes[2] = -1; //signals to binary read that nothing shall be read
      }
      //cerr << "byte types: " << bytes[0] << " " << bytes[1] << " " << bytes[2] << endl;
      int u = binary_read(data_file, 1+bytes[0]), 
          v = binary_read(data_file, 1+bytes[1]),
          w = binary_read(data_file, 1+bytes[2]);

      //cerr << u << " " << v << " " << w << '\n';
      weight[u][v] = w;
      weight[v][u] = w;

      bytes_read += 3 + bytes[0] + bytes[1] + bytes[2];
      //cerr << bytes_read << " " << data_file.tellg() << '\n';
    }
    cerr << bytes_read << " bytes read" << endl;
  } else {
    if (filename.find("reddit") != string::npos) {
      ifstream edges(filename, ifstream::in);

      int u, v, t;

      while (edges >> u) {
        edges >> v; edges >> t;
        if (!label.count(u)) {
          label[u] = nnodes++;
        }
        if (!label.count(v)) {
          label[v] = nnodes++;
        }

        weight[u][v]++;
        weight[v][u]++;
      }
    } else {
      ifstream simplices(filename + "-simplices.txt", ifstream::in);
      ifstream nverts(filename + "-nverts.txt", ifstream::in);
      ifstream times(filename + "-times.txt", ifstream::in);

      int ns = 0;

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
    }
  }

  double largest_weight = 0, sum_weight = 0, median_weight = 0;
  vector<double> all_weights;

  Graph G(nnodes);
  cerr << "constructing graph with " << nnodes << " nodes" << endl;
  for (auto& e0 : weight) {
    for (auto& e1 : e0.second) {
      int u, v;
      if (!binary) {
        u = label[e0.first], v = label[e1.first];
      } else {
        u = e0.first, v = e1.first;
      }
      long long w = e1.second;
      G[u].push_back({v, w});
      nedges++;
      if (u < v) {
        largest_weight = max(largest_weight, (double) w);
        sum_weight += w;
        all_weights.push_back(w);
      }
    }
  }

  cerr << "done constructing graph" << endl;

  nth_element(all_weights.begin(), all_weights.begin() + all_weights.size() / 2, all_weights.end());
  median_weight = all_weights[all_weights.size() / 2];
  cerr << "Largest edge weight is " << largest_weight << endl;
  cerr << "Mean edge weight is " << (2 * sum_weight) / (nedges) << endl;
  cerr << "Median edge weight is " << median_weight << endl;

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
  if ((long long) p == 1) {
    return;
  }
  long long max_weight_edge = -1;
  long long num_edges = 0;
  for (int u = 0; u < (int) G.size(); u++) {
    for (auto &e : G[u]) {
      max_weight_edge = max(max_weight_edge, e.wt);
      num_edges++;
    }
  }
  num_edges /= 2;

  long long factor = 1000000000;

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