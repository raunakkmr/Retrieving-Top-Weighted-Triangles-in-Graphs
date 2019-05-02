#include <bits/stdc++.h>

using namespace std;

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

vector<int> ordering;

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

void TuranShadow(Graph& G, int CLIQUE_SIZE) {
  cerr << "=============================================" << endl;
  cerr << "Running Turan combined with sampling" << endl;
  cerr << "=============================================" << endl;
  double st = clock();

  int n = G.size();

  priority_queue<tuple<int, int, Shadow>> q;
  for (int i = 0; i < n; i++) {
    if (G[i].size() < CLIQUE_SIZE - 1) continue;
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
    int v = get<1>(subprob);
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
      if (S.g.size() == S.l || S.l <= 2) {
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
            if (keyval.second.size() < S.l - 1) continue;

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
  if (P.size() < k) {
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
  priority_queue<pair<int, int>> pq;
  vector<bool> done(n);
  vector<int> degree(n);
  for (int i = 0; i < n; i++) {
    degree[i] = Gu[i].size();
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
    for (auto v : Gu[u]) {
      if (done[v.dst]) continue;
      degree[v.dst]--;
      pq.push(make_pair(-degree[v.dst], v.dst));
    }
  }

  cerr << "Degeneracy distribution:" << endl;
  double sum = 0, num = 0;
  for (auto kv : dfreq) {
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