// This file implements generic graph data structures, and graph reading /
// writing.

#ifndef GRAPH_H
#define GRAPH_H

#include <bits/stdc++.h>
#include <thread>
#include <parallel/algorithm>
#include <time.h>
#include <sys/sysinfo.h>

using namespace std;

namespace wsdm_2019_graph {

  // Returns the number of bytes - 1 needed to represented x for our specialized
  // binary files.
  inline int get_bytes(unsigned int x) {
    if (x <= numeric_limits<unsigned char>::max()) {
      return 0;
    } else if (x <= numeric_limits<short>::max()) {
      return 1;
    } else if (x <= numeric_limits<unsigned int>::max()) {
      return 3;
    }
    return -1;
  }

  // Class to read and write in our specialized binary format.
  class BinaryReader {
    const static int BLENGTH = 1024 * 1024;
    char buf[2 * BLENGTH + 4];
    int curr = 0, next_to_read = 0;

    istream& stream;

    public:
    BinaryReader(istream& stream) : stream(stream) {
      stream.read(buf, 2 * BLENGTH);
    }

    inline long long read(int nbytes) {
      long long res = 0;
      if (nbytes == 1) {
        res = *reinterpret_cast<unsigned char*>(buf + curr);
      } else if (nbytes == 2) {
        res = *reinterpret_cast<short*>(buf + curr);
      } else if (nbytes == 4) {
        res = *reinterpret_cast<unsigned int*>(buf + curr);
      } else if (nbytes == 8) {
        res = *reinterpret_cast<long long*>(buf + curr);
      } else {
        throw "not yet implemented";
      }

      curr += nbytes;
      if (curr >= BLENGTH && !next_to_read) {
        if (!stream.eof()) {
          stream.read(buf, BLENGTH);
        }
        // When reading into this circular buffer, an integer might get split at
        // the edges, so we assume the maximum contiguous chunk representing an
        // int is at most length 4, and copy over 3 bytes from the beginning to
        // make sure this never happens.
        buf[2 * BLENGTH] = buf[0];
        buf[2 * BLENGTH + 1] = buf[1];
        buf[2 * BLENGTH + 2] = buf[2];
        next_to_read = 1;
      } else if (curr >= 2 * BLENGTH && next_to_read) {
        curr %= 2 * BLENGTH;
        if (!stream.eof()) {
          stream.read(buf+BLENGTH, BLENGTH);
        }
        next_to_read = 0;
      }

      return res;
    }
  };

  // Write out value in binary format without any compression.
  template<typename T>
    void binary_write(std::ostream& stream, T value){
      stream.write(reinterpret_cast<char*>(&value), sizeof(T));
    }

  // Write out x in binary format with compression.
  inline void binary_compressed_write(std::ostream& stream, unsigned int x) {
    if (x <= numeric_limits<unsigned char>::max()) {
      unsigned char value = x;
      stream.write(reinterpret_cast<char*>(&value), 1);
    } else if (x <= numeric_limits<short>::max()) {
      short value = x;
      stream.write(reinterpret_cast<char*>(&value), 2);
    } else if (x <= numeric_limits<unsigned int>::max()) {
      unsigned int value = x;
      stream.write(reinterpret_cast<char*>(&value), 4);
    } else {
      throw "not yet implemented";
    }
  }

  // Generate a 64 bit random number.
  long long rand64() {
    return rand() * (1LL << 32) + rand();
  }

  // Generic data structures for graphs.
  struct full_edge {
    int src, dst;
    long long wt;
    full_edge() {
      wt = 0;
    }
    full_edge(int s, int d, long long w) : src(s), dst(d), wt(w) {}
    inline bool operator<(const full_edge& o) const {
      if (o.wt != wt) return wt < o.wt;
      return make_pair(src, dst) < make_pair(o.src, o.dst);
    }
    inline bool operator==(const full_edge& o) const {
      return (o.wt == wt) && (o.src == src) && (o.dst == dst);
    }
    inline bool operator>(const full_edge& o) const {
      if (o.wt != wt) return wt > o.wt;
      return make_pair(src, dst) > make_pair(o.src, o.dst);
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
    inline bool operator==(const half_edge& o) const {
      return (dst == o.dst) && (wt == o.wt);
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
    weighted_triangle() {
      weight = 0;
    }
    weighted_triangle(const weighted_triangle& t) :
      vertices(t.vertices), weight(t.weight) {
      }
    weighted_triangle(weighted_triangle&& t) noexcept : 
      vertices(move(t.vertices)), weight(t.weight) { 
      }
    weighted_triangle& operator=(weighted_triangle&& t) {
      vertices = move(t.vertices);
      weight = t.weight;
      return *this;
    }
    weighted_triangle& operator=(const weighted_triangle&) = default;

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
    weighted_clique(const weighted_clique& wc) : 
      vertices(wc.vertices), weight(wc.weight) {
      }
    weighted_clique(weighted_clique&& wc) noexcept : 
      vertices(move(wc.vertices)), weight(wc.weight) { 
      }
    weighted_clique& operator=(weighted_clique&& wc) {
      vertices = move(wc.vertices);
      weight = wc.weight;
      return *this;
    }
    weighted_clique& operator=(const weighted_clique&) = default;

    inline const bool operator<(const weighted_clique& o) const {
      if (weight != o.weight) return weight > o.weight;
      return vertices < o.vertices;
    }

    inline const bool operator==(const weighted_clique& o) const {
      return (weight == o.weight && vertices == o.vertices);
    }

    inline const bool operator!=(const weighted_clique& o) const {
      return (weight != o.weight || vertices != o.vertices);
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

  struct GraphStruct {
    int n, m;
    Graph G;
    vector<full_edge> edges;

    GraphStruct() {}
    GraphStruct(int n_, int m_, Graph G_, vector<full_edge> edges_) {
      n = n_;
      m = m_;
      G = G_;
      edges = edges_;
    }
  };

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

  // Reads the graph file and returns a graph data structure.
  // Arguments:
  //   filename: Path to graph file.
  //   binary: Reading our specialized binary files or not.
  //   format: If binary is false then this indicates the format of the graph
  //   file. One of simplicial, temporal or weighted.
  // If binary is true then filename is the path to the .binary file which
  // contains the graph in the following form: one line with the number of
  // vertices n, one line with the number of edges m, and 3m lines representing
  // edges (u, v, w). format is ignored in this case.
  // Otherwise, for the simplicial datasets (format = simplicial) filename
  // contains the path and common prefix of the datafiles. For the
  // temporal-reddit-reply dataset (format = temporal) filename is the path to
  // the temporal-reddit-reply.txt file. For weighted graphs (format = weighted)
  // filename is the path to the file which contains the graph in the following
  // form: one line # n m, and m lines representing edges u v w.
  GraphStruct read_graph(string filename, bool binary=false, string format="") {

    // No use for the times right now.
    unordered_map<int, unordered_map<int, long long>> weight;
    map<int, int> label;
    GraphStruct GS;
    int nnodes = 0;
    int nedges = 0;
    int m = 0;

    cerr << "Reading in graph: " << filename << endl;

    if (binary) {
      cerr << "Binary reading mode." << endl;
      ifstream data_file(filename, ios::binary | ios::in);
      BinaryReader reader(data_file);

      nnodes = reader.read(4);
      m = reader.read(4);

      size_t bytes_read = 8;
      cerr<< "Nodes: " << nnodes << ", edges: " << m << endl;

      GS.G.resize(nnodes);

      for (int i = 0; i < m; i++) {
        unsigned char type = reader.read(1);
        bytes_read += 1;

        int bytes[] = {type&3, (type>>2)&3, (type>>4)};

        int u = reader.read(1+bytes[0]), 
            v = reader.read(1+bytes[1]);

        unsigned int w;
        // Special case when weight is small, which is true for a majority of
        // the weights.
        if (bytes[2] < 12) {
          w = 1 + bytes[2];
          bytes[2] = -1;
        } else {
          bytes[2] -= 12;
          w = reader.read(1+bytes[2]);
        }

        GS.G[u].push_back({v, w});
        GS.G[v].push_back({u, w});

        bytes_read += 3 + bytes[0] + bytes[1] + bytes[2];
      }

      cerr << "Bytes read: " << bytes_read << endl;
    } else {
      cerr << "Non-binary reading mode." << endl;
      if (format.find("temporal") != string::npos) {
        ifstream edges(filename, ifstream::in);

        long long u, v, t;

        while (edges >> u) {
          edges >> v; edges >> t;
          if (!label.count(u)) {
            label[u] = nnodes++;
          }
          if (!label.count(v)) {
            label[v] = nnodes++;
          }

          u = label[u], v = label[v];
          // Weight of an edge is equal to the number of temporal edges between
          // the endpoints.
          if (u < v) {
            weight[u][v]++;
          } else if (u > v) {
            weight[v][u]++;
          }
        }
      } else if (format.find("weighted") != string::npos) {
        ifstream edge_file(filename, ifstream::in);

        string hash = "";
        edge_file >> hash >> nnodes >> m;
        cerr<< "Nodes: " << nnodes << ", edges: " << m << endl;

        GS.G.resize(nnodes);

        long long u, v, w;
        long long node_id = 0;
        for (int i = 0; i < m; i++) {
          edge_file >> u >> v >> w;
          if (!label.count(u)) {
            label[u] = node_id++;
          }
          if (!label.count(v)) {
            label[v] = node_id++;
          }

          u = label[u], v = label[v];
          GS.G[u].push_back({(int) v, w});
          GS.G[v].push_back({(int) u, w});
        }
      } else if (format.find("simplicial") != string::npos) {
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
            simplex[i] = label[simplex[i]];
          }

          // Weight of an edge is the number of times the endpoints appeared
          // together in a simplex.
          for (int i = 0; i < ns; i++) {
            for (int j = i+1; j < ns; j++) {
              // We filter out self loops in this step, but if a simplex has two
              // of the same node, we double the weight contribution to that
              // node
              if (simplex[i] < simplex[j]) {
                weight[simplex[i]][simplex[j]]++;
              } else if (simplex[i] > simplex[j]) {
                weight[simplex[j]][simplex[i]]++;
              }
            }
          }
        }
      } else {
        cerr << "ERROR! Unrecognized file format." << endl;
      }
    }

    double largest_weight = 0, sum_weight = 0, median_weight = 0;
    vector<double> all_weights;

    cerr << "Constructing graph with " << nnodes << " nodes." << endl;
    if (!binary && (format.find("weighted") == string::npos)) {
      GS.G.resize(nnodes);
      for (auto& e0 : weight) {
        for (auto& e1 : e0.second) {
          int u = e0.first, v = e1.first;
          long long w = e1.second;
          GS.G[u].push_back({v, w});
          GS.G[v].push_back({u, w});
          if (u < v) {
            GS.edges.push_back({u, v, w});
          } else {
            GS.edges.push_back({v, u, w});
          }
          nedges++;

          largest_weight = max(largest_weight, (double) w);
          sum_weight += w;
          all_weights.push_back(w);
        }
      }
    } else {
      for (int u = 0; u < (int) GS.G.size(); u++) {
        for (const auto& e : GS.G[u]) {
          int v = e.dst, w = e.wt;
          if (u < v) {
            GS.edges.push_back({u, v, w});
            nedges++;

            largest_weight = max(largest_weight, (double) w);
            sum_weight += w;
            all_weights.push_back(w);
          }
        }
      }
    }

    GS.n = nnodes;
    GS.m = nedges;

    cerr << "Done constructing graph." << endl;
    cerr << endl;

    omp_set_num_threads(thread::hardware_concurrency());
    omp_set_nested(1);
    cerr << "Starting parallel sort." << endl;
    struct timespec sort_st, sort_end;
    clock_gettime(CLOCK_MONOTONIC, &sort_st);
    __gnu_parallel::sort(GS.edges.rbegin(), GS.edges.rend());
    clock_gettime(CLOCK_MONOTONIC, &sort_end);
    double sort_time = (sort_end.tv_sec - sort_st.tv_sec);
    sort_time += (sort_end.tv_nsec - sort_st.tv_nsec) / 1000000000.0;
    cerr << "Sort time (s): " << sort_time << endl;

    nth_element(all_weights.begin(), all_weights.begin() + all_weights.size() / 2, all_weights.end());
    median_weight = all_weights[all_weights.size() / 2];
    cerr << "Largest edge weight is " << largest_weight << endl;
    cerr << "Mean edge weight is " << sum_weight / nedges << endl;
    cerr << "Median edge weight is " << median_weight << endl;

    cerr << endl;

    cerr << "Read in a graph with " << nnodes << " nodes and " << nedges << " edges." << endl;
    cerr << "Average degree: " << 2.0 * nedges / nnodes << endl;

    return GS;
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
