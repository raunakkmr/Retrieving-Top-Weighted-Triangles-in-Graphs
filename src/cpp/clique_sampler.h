#ifndef CLIQUE_SAMPLER_H
#define CLIQUE_SAMPLER_H

#include <bits/stdc++.h>

#include "graph.h"
#include "clique_enumerator.h"

using namespace std;

namespace wsdm_2019_graph {

    inline void prune_edges(Graph& G, vector<int>& degree, vector<bool>& removed, int thresh) {
        degree.clear();
        removed.clear();

        degree.resize(G.size());
        removed.resize(G.size());
        queue<int> q;
        for (int i = 0; i < (int) G.size(); i++) {
            removed[i] = false;
            degree[i] = G[i].size();

            if (degree[i] < thresh) {
                q.push(i);
                removed[i] = true;
            }
        }
        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (auto e : G[u]) {
                degree[e.dst]--;
                if (!removed[e.dst] && degree[e.dst] < thresh) {
                    removed[e.dst] = true;
                    q.push(e.dst);
                }
            }
        }
    }

    vector<weighted_clique> enumerate_cliques(Graph& G, int k) {
        vector<weighted_clique> retval;
        if (k == 1) {
            for (int i = 0; i < (int) G.size(); i++) {
                retval.push_back(weighted_clique(vector<int>({i}), 0));
            }
            return retval;
        } else if (k == 2) {
            for (int u = 0; u < (int) G.size(); u++) {
                for (auto e : G[u]) {
                    int v = e.dst;
                    if (u >= v) continue;
                    long long w = e.wt;
                    retval.push_back(weighted_clique(vector<int>({u, v}), w));
                }
            }
            return retval;
        }

        // prune the graph based on degeneracy
        vector<int> degree;
        vector<bool> removed;
        prune_edges(G, degree, removed, k-1);

        for (int i = 0; i < (int) G.size(); i++) {
            if (removed[i] || degree[i] < k-1) continue;

            map<int, int> label, unlabel;
            int cidx = 0;
            for (int j = 0; j < (int) G[i].size(); j++) {
                int u = G[i][j].dst;
                if (u < i || removed[u]) continue;
                label[u] = cidx;
                unlabel[cidx] = u;
                cidx++;
            }
            if ((int) label.size() < k-1) {
                removed[i] = true;
                continue;
            }

            int nedges = 0;
            Graph subgraph(label.size());
            for (const auto& kv : label) {
                int u = kv.first;
                for (const auto& e : G[u]) {
                    if (e.dst <= u || !label.count(e.dst)) continue;
                    int v = e.dst, w = e.wt;
                    subgraph[label[u]].push_back({label[v], w});
                    subgraph[label[v]].push_back({label[u], w});
                    nedges++;
                }
            }
            if (nedges < (k-1) * (k-2) / 2) continue;

            auto cliques = enumerate_cliques(subgraph, k-1);
            for (auto& clique : cliques) {
                // unlabelling phase
                set<int> seen;
                for (int& v : clique.vertices) {
                    v = unlabel[v];
                    seen.insert(v);
                }
                clique.vertices.push_back(i);

                for (const auto& e : G[i]) {
                    if (seen.count(e.dst)) {
                        clique.weight += e.wt;
                    }
                }
                sort(clique.vertices.begin(), clique.vertices.end());
                retval.push_back(move(clique));
            }

            removed[i] = true;

            for (auto e : G[i]) {
                degree[e.dst]--;
            }
        }
        return retval;
    }

    set<weighted_clique> clique_brute_force(Graph& G, int k) {
        cerr << "=============================================" << endl;
        cerr << "Running brute force for k-cliques" << endl;
        cerr << "=============================================" << endl;
        double st = clock();

        vector<weighted_clique> cliques = enumerate_cliques(G, k);
        // Min clique is actually max weight
        auto max_clique = min_element(cliques.begin(), cliques.end());
        double tot_time = (clock() - st) / CLOCKS_PER_SEC;

        cerr << "Found " << cliques.size() << " " << k << "-cliques." << endl;
        if (cliques.size()) cerr << "The maximum weight " << k << "-clique was " << *max_clique << endl;
        cerr << "Total Time (s): " << tot_time << endl;

        return set<weighted_clique>(cliques.begin(), cliques.end());
    }

    set<weighted_clique> clique_sampler(GraphStruct& GS, int k, int nsamples) {
        cerr << "=============================================" << endl;
        cerr << "Running edge sampling for k-cliques" << endl;
        cerr << "=============================================" << endl;
        double pre_st = clock();

        Graph &G = GS.G;
        const vector<full_edge>& edges = GS.edges;
        long long total_edge_weight = 0;
        vector<int> weight_index;
        int cur = 0;
        while (cur < (int) edges.size()) {
          weight_index.push_back(cur);
          long long cur_wt = edges[cur].wt;
          int nsteps = 5, found = 0;
          while (cur < (int) edges.size() && nsteps--) {
            cur++;
            if (edges[cur].wt < cur_wt) {
              found = 1;
              break;
            }
          }

          if (!found) {
            cur = lower_bound(edges.begin() + cur, edges.end(), full_edge(0, 0, cur_wt)) - edges.begin();
          }
          total_edge_weight += (cur - weight_index.back()) * cur_wt;
        }
        weight_index.push_back(cur);

        // Prune the graph first
        vector<int> degree;
        vector<bool> removed;
        prune_edges(G, degree, removed, k-1);

        cerr << "Precompute time (s): " << 1.0 * (clock() - pre_st)/CLOCKS_PER_SEC << endl;
        cerr << "Edge weight classes: " << int(weight_index.size())-1 << endl;
        cerr << "Total edge weight: " << total_edge_weight << endl;

        double st = clock();

        auto batched_sample_edges = [&](int num_samples){
          vector<int> sample_index(num_samples);
          vector<long long> sample_numbers(num_samples);
          for (int i = 0; i < num_samples; i++) {
            long long s = rand64() % total_edge_weight;
            sample_numbers[i] = s;
          }

          long long cur_weight = 0;
          int j = 0;
          for (int i = 0; i < int(weight_index.size()) - 1; i++) {
            cur_weight += (weight_index[i+1] - weight_index[i]) * edges[weight_index[i]].wt;
            while (j < num_samples && sample_numbers[j] <= cur_weight) {
              int e = rand() % (weight_index[i+1] - weight_index[i]);
              sample_index[j] = e + weight_index[i];
              j++;
            }
            if (j == num_samples) break;
          }
          return sample_index;
        };

        set<weighted_clique> counter;
        set<pair<int, int>> history;
        auto sample_index = batched_sample_edges(nsamples);
        for (int samp = 0; samp < nsamples; samp++) {
            auto sample = edges[sample_index[samp]];
            int u = sample.src, v = sample.dst;
            if (removed[u] || removed[v]) continue;

            long long w = sample.wt;
            if (history.count(make_pair(u, v))) {
                //cerr << "RESAMPLED!!" << endl;
                continue;
            }
            history.insert(make_pair(u, v));
            map<int, long long> vert_to_wt_u, vert_to_wt_v;
            for (auto e : G[u]) {
                if (removed[e.dst]) continue;
                vert_to_wt_u[e.dst] = e.wt;
            }

            vector<int> common_nbrs;
            for (auto e : G[v]) {
                if (removed[e.dst]) continue;
                vert_to_wt_v[e.dst] = e.wt;
                if (vert_to_wt_u.count(e.dst)) {
                    // todo: replace with p means
                    if ((int) G[e.dst].size() >= k-1) {
                        common_nbrs.push_back(e.dst);
                    }
                }
            }

            if ((int) common_nbrs.size() < k-2) continue;
            int nedges = 0;

            Graph subgraph(common_nbrs.size());
            map<int, int> label, unlabel;
            int cidx = 0;
            for (int u_ : common_nbrs) {
                label[u_] = cidx;
                unlabel[cidx] = u_;
                cidx++;
            }

            set<int> cn_set(common_nbrs.begin(), common_nbrs.end());
            for (int u_ : common_nbrs) {
                for (const auto& e : G[u_]) {
                    if (e.dst > u_ && cn_set.count(e.dst)) {
                        int v_ = e.dst, w_ = e.wt;
                        subgraph[label[u_]].push_back({label[v_], w_});
                        subgraph[label[v_]].push_back({label[u_], w_});
                        nedges++;
                    }
                }
            }

            if (nedges < (k-2) * (k-3) / 2) continue;
            vector<weighted_clique> cliques;
            if (k < 5) {
                cliques = enumerate_cliques(subgraph, k-2);
            } else {
                cliques = find_cliques(subgraph, k-2);
            }

            for (auto& clique : cliques) {
                // unlabelling phase
                set<int> seen;
                for (int& vert : clique.vertices) {
                    vert = unlabel[vert];
                    seen.insert(vert);
                }
                clique.vertices.push_back(u);
                clique.vertices.push_back(v);

                for (const auto& e : G[u]) {
                    if (seen.count(e.dst)) {
                        clique.weight += e.wt;
                    }
                }
                for (const auto& e : G[v]) {
                    if (seen.count(e.dst)) {
                        clique.weight += e.wt;
                    }
                }
                clique.weight += w;
                sort(clique.vertices.begin(), clique.vertices.end());
                counter.insert(move(clique));
            }
        }
        cerr << "Found " << counter.size() << " " << k << "-cliques." << endl;
        if (counter.size()) cerr << "The maximum weight " << k << "-clique was " << *counter.begin() << endl;

        double tot_time = (clock() - st) / CLOCKS_PER_SEC;
        cerr << "Total Time (s): " << tot_time << endl;
        cerr << "Time per sample (s): " << tot_time / nsamples << endl;

        return counter;
    }

    vector<weighted_clique> clique_sampler_parallel(GraphStruct& GS, int k, int nsamples, int nthreads) {
        cerr << "=============================================" << endl;
        cerr << "Running parallel clique sampling (k=" << k << ", " << nthreads << " threads)" << endl;
        cerr << "=============================================" << endl;

        struct timespec pre_start, pre_finish;
        double pre_elapsed;
        clock_gettime(CLOCK_MONOTONIC, &pre_start);

        Graph &G = GS.G;
        const vector<full_edge>& edges = GS.edges;
        long long total_edge_weight = 0;
        vector<int> weight_index;
        int cur = 0;
        while (cur < (int) edges.size()) {
          weight_index.push_back(cur);
          long long cur_wt = edges[cur].wt;
          int nsteps = 5, found = 0;
          while (cur < (int) edges.size() && nsteps--) {
            cur++;
            if (edges[cur].wt < cur_wt) {
              found = 1;
              break;
            }
          }

          if (!found) {
            cur = lower_bound(edges.begin() + cur, edges.end(), full_edge(0, 0, cur_wt)) - edges.begin();
          }
          total_edge_weight += (cur - weight_index.back()) * cur_wt;
        }
        weight_index.push_back(cur);

        // Prune the graph first
        vector<int> degree;
        vector<bool> removed;
        prune_edges(G, degree, removed, k-1);

        clock_gettime(CLOCK_MONOTONIC, &pre_finish);

        pre_elapsed = (pre_finish.tv_sec - pre_start.tv_sec);
        pre_elapsed += (pre_finish.tv_nsec - pre_start.tv_nsec) / 1000000000.0;
        cerr << "Pre-processing time: " << pre_elapsed << endl;
        cerr << "Edge weight classes: " << int(weight_index.size())-1 << endl;
        cerr << "Total edge weight: " << total_edge_weight << endl;

        struct timespec start, finish;
        double tot_time;
        clock_gettime(CLOCK_MONOTONIC, &start);

        vector<thread> threads(nthreads);
        vector<vector<weighted_clique>> counters(nthreads);
        vector<set<pair<int, int>>> histories(nthreads);
        int nsamples_per_thread = ceil(nsamples / nthreads);

        auto batched_sample_edges = [&](int num_samples){
          vector<int> sample_index(num_samples);
          vector<long long> sample_numbers(num_samples);
          for (int i = 0; i < num_samples; i++) {
            long long s = rand64() % total_edge_weight;
            sample_numbers[i] = s;
          }

          long long cur_weight = 0;
          int j = 0;
          for (int i = 0; i < int(weight_index.size()) - 1; i++) {
            cur_weight += (weight_index[i+1] - weight_index[i]) * edges[weight_index[i]].wt;
            while (j < num_samples && sample_numbers[j] <= cur_weight) {
              int e = rand() % (weight_index[i+1] - weight_index[i]);
              sample_index[j] = e + weight_index[i];
              j++;
            }
            if (j == num_samples) break;
          }
          return sample_index;
        };

        auto parallel_sampler = [&](int i) {
            auto sample_index = batched_sample_edges(nsamples);
            for (int samp = 0; samp < nsamples_per_thread; samp++) {
                auto sample = edges[sample_index[samp]];
                int u = sample.src, v = sample.dst;
                if (removed[u] || removed[v]) continue;

                long long w = sample.wt;
                bool cont = false;
                for (int j = 0; j < nthreads; j++) {
                    if (histories[j].count(make_pair(u, v))) {
                        cont = true;
                        break;
                    }
                }
                if (cont) continue;
                histories[i].insert(make_pair(u, v));
                map<int, long long> vert_to_wt_u, vert_to_wt_v;
                for (auto e : G[u]) {
                    if (removed[e.dst]) continue;
                    vert_to_wt_u[e.dst] = e.wt;
                }

                vector<int> common_nbrs;
                for (auto e : G[v]) {
                    if (removed[e.dst]) continue;
                    vert_to_wt_v[e.dst] = e.wt;
                    if (vert_to_wt_u.count(e.dst)) {
                        // todo: replace with p means
                        if ((int) G[e.dst].size() >= k-1) {
                            common_nbrs.push_back(e.dst);
                        }
                    }
                }

                if ((int) common_nbrs.size() < k-2) continue;
                int nedges = 0;

                Graph subgraph(common_nbrs.size());
                map<int, int> label, unlabel;
                int cidx = 0;
                for (int u_ : common_nbrs) {
                    label[u_] = cidx;
                    unlabel[cidx] = u_;
                    cidx++;
                }

                set<int> cn_set(common_nbrs.begin(), common_nbrs.end());
                for (int u_ : common_nbrs) {
                    for (const auto& e : G[u_]) {
                        if (e.dst > u_ && cn_set.count(e.dst)) {
                            int v_ = e.dst, w_ = e.wt;
                            subgraph[label[u_]].push_back({label[v_], w_});
                            subgraph[label[v_]].push_back({label[u_], w_});
                            nedges++;
                        }
                    }
                }

                if (nedges < (k-2) * (k-3) / 2) continue;
                vector<weighted_clique> cliques;
                if (k < 5) {
                    cliques = enumerate_cliques(subgraph, k-2);
                } else {
                    cliques = find_cliques(subgraph, k-2);
                }

                for (auto& clique : cliques) {
                    // unlabelling phase
                    set<int> seen;
                    for (int& vert : clique.vertices) {
                        vert = unlabel[vert];
                        seen.insert(vert);
                    }
                    clique.vertices.push_back(u);
                    clique.vertices.push_back(v);

                    for (const auto& e : G[u]) {
                        if (seen.count(e.dst)) {
                            clique.weight += e.wt;
                        }
                    }
                    for (const auto& e : G[v]) {
                        if (seen.count(e.dst)) {
                            clique.weight += e.wt;
                        }
                    }
                    clique.weight += w;
                    sort(clique.vertices.begin(), clique.vertices.end());
                    counters[i].push_back(move(clique));
                }
            }
            sort(counters[i].begin(), counters[i].end());
        };

        auto parallel_merger = [&](int i, int j) {
            vector<weighted_clique> W;
            W.reserve(counters[i].size() + counters[j].size());
            int a = 0, b = 0, Li = counters[i].size(), Lj = counters[j].size();
            while (a < Li && b < Lj) {
                if (counters[i][a] < counters[j][b]) {
                    if (counters[i][a] != W.back()) {
                        W.push_back(move(counters[i][a]));
                    }
                    a++;
                } else if (counters[i][a] == counters[j][b]) {
                    if (counters[i][a] != W.back()) {
                        W.push_back(move(counters[i][a]));
                    }
                    a++;
                    b++;
                } else {
                    if (counters[j][b] != W.back()) {
                        W.push_back(move(counters[j][b]));
                    }
                    b++;
                }
            }
            while (a < Li) {
                if (counters[i][a] != W.back()) {
                    W.push_back(move(counters[i][a]));
                }
                a++;
            }
            while (b < Lj) {
                if (counters[j][b] != W.back()) {
                    W.push_back(move(counters[j][b]));
                }
                b++;
            }
            counters[i].swap(W);
        };

        for (int i = 0; i < nthreads; i++) {
            thread th(parallel_sampler, i);
            threads[i] = move(th);
        }
        for (int i = 0; i < nthreads; i++) {
            threads[i].join();
        }

        struct timespec merge_start, merge_finish;
        double merge_elapsed;
        clock_gettime(CLOCK_MONOTONIC, &merge_start);

        // Parallel merging.

        int pow2_sz = 1, log2_sz = 0;
        while (pow2_sz < nthreads) {
            pow2_sz *= 2;
            log2_sz++;
        }

        for (int i = counters.size(); i < pow2_sz; i++) {
            counters.push_back(vector<weighted_clique>());
        }

        vector<thread> merge_threads(pow2_sz);
        int val = 1;
        for (int level = 1; level < log2_sz+1; level++) {
            val *= 2;
            for (int i = 0; i < (int) counters.size()/val; i++) {
                thread merge_th(parallel_merger, i*val, i*val+(int)val/2);
                merge_threads[i*val] = move(merge_th);
            }
            for (int i = 0; i < (int) counters.size()/val; i++) {
                merge_threads[i*val].join();
            }
        }

        clock_gettime(CLOCK_MONOTONIC, &merge_finish);
        merge_elapsed = (merge_finish.tv_sec - merge_start.tv_sec);
        merge_elapsed += (merge_finish.tv_nsec - merge_start.tv_nsec) / 1000000000.0;
        cerr << "Merge time: " << merge_elapsed << endl;

        cerr << "Found " << counters[0].size() << " cliques." << endl;
        if (counters[0].size()) cerr << "The maximum weight clique was " << *counters[0].begin() << endl;

        clock_gettime(CLOCK_MONOTONIC, &finish);
        tot_time = (finish.tv_sec - start.tv_sec);
        tot_time += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        cerr << "Total Time (s): " << tot_time << endl;
        cerr << "Time per sample (s): " << tot_time / nsamples << endl;
        cerr << endl;

        return counters[0];
    }

}

#endif /* CLIQUE_SAMPLER_H */
