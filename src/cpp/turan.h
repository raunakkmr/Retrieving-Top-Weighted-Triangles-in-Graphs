#ifndef TURAN_H
#define TURAN_H

#include <bits/stdc++.h>
#include "graph.h"

using namespace std;

namespace wsdm_2019_graph {

void TuranShadow(Graph& G, int CLIQUE_SIZE) {
  cerr << "=============================================" << endl;
  cerr << "Running Turan combined with sampling" << endl;
  cerr << "=============================================" << endl;
  double st = clock();

  int n = G.size();

  priority_queue<tuple<int, int, Shadow>> q;
  for (int i = 0; i < n; i++) {
    if ((int) G[i].size() < CLIQUE_SIZE - 1) continue;
    Shadow s = {create_shadow(i, G), CLIQUE_SIZE - 1};
    q.push(make_tuple(-int(s.g.size()), i, s));
  }
  cerr << "Done initial pushes." << endl;

  int total_shadow_size = 0;
  int samples_needed = 0;
  while (!q.empty()) {
    if (q.size() % 33333 == 0) {
      cerr << "Queue size: " << q.size() << endl;
      cerr << "Top graph size is: " << get<2>(q.top()).g.size() << " " << -get<0>(q.top()) << endl;
      cerr << "Looking for cliques of size: " << get<2>(q.top()).l << endl;
    }

    auto subprob = q.top();
    int deg = -get<0>(subprob);
    //int v = get<1>(subprob);
    Shadow& S = get<2>(subprob);
    q.pop();

    if (deg < S.l) {
      continue;
    } else {
      int totedges = 0;
      for (auto keyval : S.g) {
        totedges += keyval.second.size();
      }

      int totverts = S.g.size();
      if ((int) S.g.size() == S.l || S.l <= 2) {
        // Do some counting here.
        samples_needed++;
        //cerr << "Edge case" << endl;
        continue;
      } else {
        if (totedges < S.l * (S.l - 1) / 2) {
          continue;
        } else if (totedges/(0.5 * totverts * (totverts - 1)) > 1 - 1.0/(S.l - 1)) {
          // Do some counting here.
          samples_needed++;
          //cerr << "Density case " << totedges << " " << totverts << endl;
          continue;
        } else {
          // Here we recurse.
          for (auto keyval : S.g) {
            if ((int) keyval.second.size() < S.l - 1) continue;

            // Create the graph and add the node to the PQ
            Shadow ns;
            ns.g = create_shadow(keyval.first, S.g);
            ns.l = S.l - 1;

            q.push(make_tuple(-int(ns.g.size()), keyval.first, ns));
            total_shadow_size++;
          }
        }
      }
    }
  }
  cerr << "Shadow size: " << total_shadow_size << " " << samples_needed << endl;
  cerr << "Time (s): " << (clock() - st) / CLOCKS_PER_SEC << endl;
}

int nedgecase = 0, nsamplingcase = 0;
int BKPivot(vector<int> R, vector<int> P, vector<int> X, Graph& G, int k) {
  // Get the edge density and such
  if ((int) P.size() < k) {
    return 0;
  }
  if (k <= 2) {
    //cerr << "edge or node case" << endl;
    //cerr << P.size() << " " << k << endl;
    nedgecase++;
    return 0;
  }

  // Compute edge density of P.
  set<int> sP(P.begin(), P.end());
  /*
  bool advanced = true;
  do {
    advanced = false;
    for (int u : P) {
      int d = 0;
      for (auto e : G[u]) {
        int v = e.dst;
        if (sP.count(v)) d++;
      }
      if (d < k - 1) {
        X.push_back(u);
        sP.erase(u);
        advanced = true;
      }
    }
    P = vector<int>(sP.begin(), sP.end());
  } while (advanced);
  //*/
  set<int> sX(X.begin(), X.end());

  int edges = 0;
  int max_deg = 0, max_v = -1;
  for (int u : P) {
    int d = 0;
    for (auto e : G[u]) {
      int v = e.dst;
      if (sP.count(v)) d++;
    }
    if (d > max_deg) {
      max_deg = d;
      max_v = u;
    }
    edges += d;
  }
  assert(edges%2 == 0);
  edges /= 2;

  if (edges < k * (k-1) / 2) {
    return 0;
  }
  if (edges == 0) return 0;

  int verts = P.size();
  if (edges/(0.5 * verts * (verts - 1)) > 1 - 1.0/(k - 1)) {
    //cerr << "sampling case" << endl;
    //cerr << edges << " " << verts << " " << k << endl;
    nsamplingcase++;
    //return 0;
  }

  for (int u : X) {
    int d = 0;
    for (auto e : G[u]) {
      int v = e.dst;
      if (sP.count(v)) d++;
    }
    if (d > max_deg) {
      max_deg = d;
      max_v = u;
    }
  }

  set<int> adjMax;
  for (auto e : G[max_v]) {
    int v = e.dst;
    if (sP.count(v)) {
      adjMax.insert(v);
    }
  }

  int shadow_size = 0;
  for (int v : P) {
    if (adjMax.count(v)) continue;
    vector<int> nR = R, nP, nX;
    for (auto e : G[v]) {
      int u = e.dst;
      if (sP.count(u)) nP.push_back(u);
      if (sX.count(u)) nX.push_back(u);
    }
    nR.push_back(v);
    shadow_size += BKPivot(nR, nP, nX, G, k-1);
    shadow_size++;

    sP.erase(v);
    sX.insert(v);
  }

  return shadow_size;
}

void TuranBK(Graph& G, const vector<int>& degenOrder, int CLIQUE_SIZE) {
  cerr << "=============================================" << endl;
  cerr << "Running BronKerbosh combined with sampling" << endl;
  cerr << "=============================================" << endl;
  double st = clock();

  int n = G.size();
  cerr << "Graph size: " << n << endl;
  vector<int> degenIndex(n);
  for (int i = 0; i < n; i++) {
    degenIndex[degenOrder[i]] = i;
  }

  int shadow_size = 0;
  nedgecase = 0;
  nsamplingcase = 0;
  for (int i = 0; i < n; i++) {
    vector<int> P;
    vector<int> X;
    vector<int> R;

    int u = degenOrder[i];
    for (auto e : G[u]) {
      int v = e.dst;
      if (degenIndex[v] < i) {
        X.push_back(v);
      } else {
        P.push_back(v);
      }
    }
    R.push_back(u);

    if (i % 33333 == 0) {
      cerr << "Running " << i << "th BKPivot call" << endl;
      cerr << "Current shadow size is: " << shadow_size << endl;
      cerr << R.size() << " " << P.size() << " " << X.size() << endl;
    }    
    shadow_size += BKPivot(R, P, X, G, CLIQUE_SIZE - 1);
  }

  cerr << "Total BK shadow size: " << shadow_size << endl;
  cerr << "Edge and sampling case: " << nedgecase << " " << nsamplingcase << endl; 
  cerr << "Time (s): " << (clock() - st) / CLOCKS_PER_SEC << endl;
}

}

#endif /* TURAN_H */