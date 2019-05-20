#include <bits/stdc++.h>

#include "turan.h"
#include "graph.h"
#include "clique_sampler.h"
#include "triangle_sampler.h"

using namespace std;
using namespace wsdm_2019_graph;

#define WRITE_STATISTICS 0

vector<int> greedy_indep_set(Graph& G) {
  cerr << "=============================================" << endl;
  cerr << "Running greedy independent set algorithm" << endl;
  cerr << "=============================================" << endl;
  double st = clock();

  // Compute a maximal indep set. Should be fairly big for these graphs
#if 0
  // This is dynamically ordering by degree, which seems to be quite slow in practice.
  set<pair<int, int>> q;
  vector<int> degree(G.size());
  int totedges = 0;
  for (int i = 0; i < (int) G.size(); i++) {
    q.insert(make_pair(int(G[i].size()), i));
    degree[i] = G[i].size();
    totedges += degree[i];
  }

  vector<bool> done(G.size());
  vector<int> indep_set;
  int nedges = 0;
  while (!q.empty()) {
    auto p = *q.begin();
    q.erase(p);
    int u = p.second;
    
    indep_set.push_back(u);
    done[u] = true;
    nedges += G[u].size();
    for (auto e : G[u]) {
      if (done[e.dst]) continue;
      done[e.dst] = true;
      q.erase(make_pair(degree[e.dst], e.dst));
      for (auto e_ : G[e.dst]) {
        if (!done[e_.dst]) {
          q.erase(make_pair(degree[e_.dst], e_.dst));
          degree[e_.dst]--;
          q.insert(make_pair(degree[e_.dst], e_.dst));
        }
      }
    }
  }
#else
  int totedges = 0;
  vector<pair<int, int>> verts;
  verts.reserve(G.size());
  for (int i = 0; i < (int) G.size(); i++) {
    verts.push_back(make_pair(int(G[i].size()), i));
    totedges += G[i].size();
  }
  sort(verts.begin(), verts.end());

  vector<int> indep_set;
  vector<bool> done(G.size());
  int nedges = 0;
  for (int i = 0; i < (int) G.size(); i++) {
    int u = verts[i].second;
    if (done[u]) continue;

    done[u] = true;
    indep_set.push_back(u);
    nedges += G[u].size();
    for (auto e : G[u]) {
      done[e.dst] = true;
    }
  }

#endif

  cerr << "maximal greedy indep. set size: " << indep_set.size() << " " << nedges << endl;
  cerr << "percentage of total nodes: " << 1.0 * indep_set.size() / G.size() << endl;
  cerr << "percentage of total edges: " << 1.0 * nedges / totedges << endl;

  double tot_time = (clock() - st) / CLOCKS_PER_SEC;
  cerr << "Total Time (s): " << tot_time << endl;
  cerr << endl;

  return indep_set;
}

vector<int> random_indep_set(Graph& G) {
  cerr << "=============================================" << endl;
  cerr << "Running random independent set algorithm" << endl;
  cerr << "=============================================" << endl;
  double st = clock();

  // Compute a random indep set
  set<int> picked;
  int totedges = 0;
  for (int i = 0; i < (int) G.size(); i++) {
    int degree = G[i].size();
    if (degree == 0 || rand()%degree == 0) {
      picked.insert(i);
    }
    totedges += degree;
  }

  int nedges = 0;
  set<int> indep_set = picked;
  for (int u : picked) {
    bool erased = false;
    for (auto e : G[u]) {
      if (indep_set.count(e.dst)) {
        indep_set.erase(u);
        erased = true;
        break;
      }
    }
    if (!erased) {
      nedges += G[u].size();
    }
  }

  cerr << "maximal greedy indep. set size: " << indep_set.size() << " " << nedges << endl;
  cerr << "percentage of total nodes: " << 1.0 * indep_set.size() / G.size() << endl;
  cerr << "percentage of total edges: " << 1.0 * nedges / totedges << endl;

  double tot_time = (clock() - st) / CLOCKS_PER_SEC;
  cerr << "Total Time (s): " << tot_time << endl;
  cerr << endl;

  return vector<int>(indep_set.begin(), indep_set.end());
}

void layered_random_indep_set_ordering(Graph& G) {
  cerr << "=============================================" << endl;
  cerr << "Running layered IS ordering algorithm" << endl;
  cerr << "=============================================" << endl;
  double st = clock();

  // Compute a random indep set
  vector<int> num_removed(G.size());
  set<int> active_vertices;
  for (int i = 0; i < (int) G.size(); i++) {
    active_vertices.insert(i);
  }

  int num_rounds = 0;
  while (active_vertices.size() > 0.01 * G.size()) {
    num_rounds++;

    set<int> picked;
    int totedges = 0;
    for (int u : active_vertices) {
      int degree = G[u].size() - num_removed[u];
      if (degree == 0 || rand()%degree == 0) {
        picked.insert(u);
      }
      totedges += degree;
    }

    int nedges = 0;
    set<int> indep_set = picked;
    for (int u : picked) {
      bool erased = false;
      for (auto e : G[u]) {
        if (num_removed[e.dst] == (int) G[e.dst].size()) {
          continue;
        }

        if (indep_set.count(e.dst)) {
          indep_set.erase(u);
          erased = true;
          break;
        }
      }

      if (!erased) {
        nedges += G[u].size() - num_removed[u];
        for (auto e : G[u]) {
          if (num_removed[e.dst] == (int) G[e.dst].size()) {
            continue;
          }
          num_removed[e.dst]++;
        }
        active_vertices.erase(u);
        num_removed[u] = G[u].size();
      }
    }

    cerr << "maximal greedy indep. set size: " << indep_set.size() << " " << nedges << endl;
    cerr << "percentage of total nodes: " << 1.0 * indep_set.size() / G.size() << endl;
    cerr << "percentage of total edges: " << 1.0 * nedges / totedges << endl;
  }
  double tot_time = (clock() - st) / CLOCKS_PER_SEC;
  cerr << "Total Time (s): " << tot_time << endl;
  cerr << "Total number of rounds taken: " << num_rounds << endl;
  cerr << endl;
}

void layered_greedy_indep_set_ordering(Graph& G) {
  cerr << "=============================================" << endl;
  cerr << "Running layered IS ordering algorithm" << endl;
  cerr << "=============================================" << endl;
  double st = clock();

  // Compute a random indep set
  vector<int> num_removed(G.size());
  set<int> active_vertices;
  for (int i = 0; i < (int) G.size(); i++) {
    active_vertices.insert(i);
  }

  vector<pair<int, int>> verts;
  verts.reserve(G.size());
  int num_rounds = 0;
  while (active_vertices.size() > 0.01 * G.size()) {
    num_rounds++;

    int totedges = 0;
    verts.clear();
    for (int u : active_vertices) {
      int degree = G[u].size() - num_removed[u];
      verts.push_back(make_pair(degree, u));
      totedges += degree;
    }
    sort(verts.begin(), verts.end());

    vector<int> indep_set;
    set<int> done;
    int nedges = 0;
    for (int i = 0; i < (int) verts.size(); i++) {
      int u = verts[i].second;
      if (done.count(u)) continue;
      done.insert(u);
      active_vertices.erase(u);
      nedges += G[u].size() - num_removed[u];
      num_removed[u] = G[u].size();
      
      indep_set.push_back(u);
      nedges += G[u].size() - num_removed[u];
      for (auto e : G[u]) {
        if (num_removed[e.dst] == (int) G[e.dst].size()) {
          continue;
        }
        num_removed[e.dst]++;
        done.insert(e.dst);
      }
    }

    cerr << "maximal greedy indep. set size: " << indep_set.size() << " " << nedges << endl;
    cerr << "percentage of total nodes: " << 1.0 * indep_set.size() / G.size() << endl;
    cerr << "percentage of total edges: " << 1.0 * nedges / totedges << endl;
  }
  double tot_time = (clock() - st) / CLOCKS_PER_SEC;
  cerr << "Total Time (s): " << tot_time << endl;
  cerr << "Total number of rounds taken: " << num_rounds << endl;
  cerr << endl;
}

int main(int argc, char* argv[]) {
  ios::sync_with_stdio(0);
  cin.tie(0);
  srand(0); 

  auto G = read_graph(argv[1]);

  greedy_indep_set(G);
  random_indep_set(G);
  layered_random_indep_set_ordering(G);
  layered_greedy_indep_set_ordering(G);

  return 0;
}