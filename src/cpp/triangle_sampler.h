#ifndef TRIANGLE_SAMPLER_H
#define TRIANGLE_SAMPLER_H

#include <bits/stdc++.h>
#include <parallel/algorithm>

#include "graph.h"

using namespace std;

namespace wsdm_2019_graph {

  // Reads triangles from a given file and returns it as a vector.
  vector<weighted_triangle> read_all_triangles_vector(string filename) {

    ifstream tri_file(filename, ios::binary | ios::in);
    BinaryReader reader(tri_file);

    long long num_tris = reader.read(8);
    size_t bytes_read = 8;

    cerr << "Number of triangles: " << num_tris << '\n';

    vector<weighted_triangle> v;
    v.reserve(num_tris);

    for (int i = 0; i < num_tris; i++) {
      unsigned char type = reader.read(1);
      bytes_read += 1;

      int bytes[] = {type&3, (type>>2)&3, (type>>4)&3, (type>>6)};

      int a = reader.read(1+bytes[2]),
          b = reader.read(1+bytes[1]),
          c = reader.read(1+bytes[0]);

      int w;
      if (bytes[3] == 2) {
        // Special case when the triangle weight is 3.
        w = 3;
        bytes[3] = -1;
      } else {
        w = reader.read(1+bytes[3]);
      }

      bytes_read += 4 + bytes[0] + bytes[1] + bytes[2] + bytes[3];

      v.push_back(weighted_triangle(a, b, c, w));
    }

    return v;
  }

  // Reads triangles from a given file and returns it as a set.
  set<weighted_triangle> read_all_triangles_set(string filename) {
    vector<weighted_triangle> v = read_all_triangles_vector(filename);
    set<weighted_triangle> s(v.begin(), v.end());

    return s;
  }

  // Minbucket brute force algorithm.
  // Arguments:
  //   GS: Adjacency list and edge list.
  //   k: Parameter k for top-k. If k = -1, returns all triangles. Otherwise
  //   returns the top-k triangles.
  //   diagnostic: Whether to print diagnostic information.
  set<weighted_triangle> brute_force_sampler_minbucket(GraphStruct &GS, int k=-1, bool diagnostic=true) {
    if (diagnostic) {
      cerr << "=============================================" << endl;
      cerr << "Running brute force detection (minbucket) for triangles" << endl;
      cerr << "=============================================" << endl;
    }

    double st = clock();

    // Can use degeneracy ordering to speed up brute force
    // However, on certain tests it didnt have great performance.
    Graph &G = GS.G;
    set<weighted_triangle> counter;
    long long num_tris = 0;

    vector<long long> vert_to_wt(G.size());
    for (int u = 0; u < (int) G.size(); u++) {
      for (const auto& e : G[u]) {
        if (G[e.dst].size() > G[u].size() || 
            (G[e.dst].size() == G[u].size() && e.dst < u)) continue;
        vert_to_wt[e.dst] = e.wt;
      }

      for (const auto& e : G[u]) {
        int v = e.dst;
        long long w = e.wt;
        if (G[e.dst].size() > G[u].size() || 
            (G[e.dst].size() == G[u].size() && e.dst < u)) continue;

        for (const auto& ev : G[v]) {
          if (G[ev.dst].size() > G[v].size() || 
              (G[ev.dst].size() == G[v].size() && ev.dst < v)) continue;
          if (vert_to_wt[ev.dst]) {
            // todo: replace with p means
            long long val = ev.wt + vert_to_wt[ev.dst] + w;
            weighted_triangle tri = weighted_triangle(u, v, ev.dst, val);
            num_tris++;
            if (k < 0 || (int) counter.size() < k) {
              counter.insert(tri);
            } else {
              // Compare weights, not triangles.
              auto it = --(counter.end());
              if (val > it->weight) {
                counter.erase(it);
                counter.insert(tri);
              }
            }
          }
        }
      }

      for (const auto& e : G[u]) {
        if (G[e.dst].size() > G[u].size() || 
            (G[e.dst].size() == G[u].size() && e.dst < u)) continue;
        vert_to_wt[e.dst] = 0;
      }
    }

    cerr << "Found " << num_tris << " triangles." << endl;
    if (diagnostic) {
      if (counter.size()) cerr << "The maximum weight triangle was " << *counter.begin() << endl;

      double tot_time = (clock() - st) / CLOCKS_PER_SEC;
      cerr << "Total Time (s): " << tot_time << endl;
      cerr << endl;
    }
    return counter;
  }

  // Brute force algorithm that iterates over vertices in decreasing order of
  // degree, and for each vertex only enumerates triangles with lower degree
  // than itself.
  // Arguments:
  //   GS: Adjacency list and edge list.
  //   k: Parameter k for top-k. If k = -1, returns all triangles. Otherwise
  //   returns the top-k triangles.
  //   diagnostic: Whether to print diagnostic information.
  set<weighted_triangle> brute_force_sampler(GraphStruct &GS, int k=-1, bool diagnostic=true) {
    if (diagnostic) {
      cerr << "=============================================" << endl;
      cerr << "Running brute force detection for triangles" << endl;
      cerr << "=============================================" << endl;
    }

    double st = clock();

    // Can use degeneracy ordering to speed up brute force
    // However, on certain tests it didnt have great performance.
    Graph &G = GS.G;
    set<weighted_triangle> counter;
    long long num_tris = 0;

    vector<long long> vert_to_wt(G.size());
    for (int u = 0; u < (int) G.size(); u++) {
      for (auto e : G[u]) {
        if (e.dst < u) continue;
        vert_to_wt[e.dst] = e.wt;
      }

      for (auto e : G[u]) {
        int v = e.dst;
        long long w = e.wt;
        if (v < u) continue;

        for (auto ev : G[v]) {
          if (ev.dst < v) continue;
          if (vert_to_wt[ev.dst]) {
            // todo: replace with p means
            long long val = ev.wt + vert_to_wt[ev.dst] + w;
            weighted_triangle tri = weighted_triangle(u, v, ev.dst, val);
            num_tris++;
            if (k < 0 || (int) counter.size() < k) {
              counter.insert(tri);
            } else {
              // Compare weights, not triangles.
              auto it = --(counter.end());
              if (val > it->weight) {
                counter.erase(it);
                counter.insert(tri);
              }
            }
          }
        }
      }

      for (auto e : G[u]) {
        if (e.dst < u) continue;
        vert_to_wt[e.dst] = 0;
      }
    }

    cerr << "Found " << num_tris << " triangles." << endl;
    if (diagnostic) {
      if (counter.size()) cerr << "The maximum weight triangle was " << *counter.begin() << endl;

      double tot_time = (clock() - st) / CLOCKS_PER_SEC;
      cerr << "Total Time (s): " << tot_time << endl;
      cerr << endl;
    }
    return counter;
  }

  vector<weighted_triangle> edge_sampler(GraphStruct &GS, int nthreads=1, int max_samples=-1, double max_time=-1, bool include_setup=true) {
    cerr << "=============================================" << endl;
    cerr << "Running edge sampling for triangles (" << nthreads << " threads)" << endl;
    cerr << "=============================================" << endl;

    struct timespec pre_start, pre_finish;
    double pre_elapsed;
    clock_gettime(CLOCK_MONOTONIC, &pre_start);
  
    Graph &G = GS.G;
    const vector<full_edge>& edges = GS.edges;
    long long total_edge_weight = 0;
    vector<int> weight_index;
    vector<long long> weight_value;
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
        cur = lower_bound(edges.begin() + cur, edges.end(), full_edge(0, 0, cur_wt), greater<full_edge>()) - edges.begin();
      }
      total_edge_weight += (cur - weight_index.back()) * cur_wt;
      weight_value.push_back(total_edge_weight);
    }
    weight_index.push_back(cur);

    vector<thread> threads(nthreads);
    vector<vector<weighted_triangle>> counters(nthreads);
    vector<set<pair<int, int>>> histories(nthreads);
    int nsamples_per_thread = ceil(max_samples / nthreads);
  
    clock_gettime(CLOCK_MONOTONIC, &pre_finish);
    pre_elapsed = (pre_finish.tv_sec - pre_start.tv_sec);
    pre_elapsed += (pre_finish.tv_nsec - pre_start.tv_nsec) / 1000000000.0;
  
    cerr << "Pre-processing time: " << pre_elapsed << endl;
    cerr << "Edge weight classes: " << int(weight_index.size())-1 << endl;
    cerr << "Total edge weight: " << total_edge_weight << endl;
  
    struct timespec start, finish;
    double tot_time;
    clock_gettime(CLOCK_MONOTONIC, &start);
  
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
  
    auto sample_single_edge = [&](){
      long long s = rand64() % total_edge_weight;
      int index = lower_bound(weight_value.begin(), weight_value.end(), s) - weight_value.begin();
      int e = rand() % (weight_index[index+1] - weight_index[index]);
      return edges[e + weight_index[index]];
    };

    auto terminate = [&](int nsamples_) {
      if (max_samples != -1) {
        return nsamples_ >= nsamples_per_thread;
      } else {
        struct timespec cur;
        double tot_time;
        clock_gettime(CLOCK_MONOTONIC, &cur);
        if (include_setup) {
          tot_time = (cur.tv_sec - pre_start.tv_sec);
          tot_time += (cur.tv_nsec - pre_start.tv_nsec) / 1000000000.0;
        } else {
          tot_time = (cur.tv_sec - start.tv_sec);
          tot_time += (cur.tv_nsec - start.tv_nsec) / 1000000000.0;
        }
        return tot_time >= max_time;
      }
    };

    auto sampler = [&](int i) {
      vector<int> sample_index;
      if (max_samples != -1) {
        sample_index = batched_sample_edges(max_samples);
      }
      int nsamples_ = 0;
      while (!terminate(nsamples_)) {
        full_edge e;
        if (max_samples != -1) {
          e = edges[sample_index[nsamples_]];
        } else {
          e = sample_single_edge();
        }
        nsamples_++;
        int u = e.src, v = e.dst;
        long long w = e.wt;
        bool cont = false;
        for (int j = 0; j < nthreads; j++) {
          if (histories[j].count(make_pair(u, v))) {
            cont = true;
            break;
          }
        }
        if (cont) continue;
        histories[i].insert(make_pair(u, v));
        unordered_map<int, long long> vert_to_wt;
        for (const auto &eu : G[u]) {
          vert_to_wt[eu.dst] = eu.wt;
        }
  
        for (const auto &ev : G[v]) {
          if (vert_to_wt.count(ev.dst)) {
            auto tri = weighted_triangle(u, v, ev.dst, ev.wt + vert_to_wt[ev.dst] + w);
            counters[i].push_back(tri);
          }
        }
      }
      sort(counters[i].begin(), counters[i].end());
    };

    auto parallel_merger = [&](int i, int j) {
      vector<weighted_triangle> W;
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

    if (nthreads > 1) {
      for (int i = 0; i < nthreads; i++) {
        thread th(sampler, i);
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
        counters.push_back(vector<weighted_triangle>());
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
    } else {
      sampler(0);
    }

    cerr << "Found " << counters[0].size() << " triangles." << endl;
    if (counters[0].size()) cerr << "The maximum weight triangle was " << *counters[0].begin() << endl;

    clock_gettime(CLOCK_MONOTONIC, &finish);
    tot_time = (finish.tv_sec - start.tv_sec);
    tot_time += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    cerr << "Total Time (s): " << tot_time << endl;
    // cerr << "Time per sample (s): " << tot_time / nsamples << endl;
    cerr << endl;

    return counters[0];

  }

  vector<weighted_triangle> wedge_sampler(GraphStruct &GS, int nthreads, int max_samples=-1, double max_time=-1, bool include_setup=true) {
    cerr << "=============================================" << endl;
    cerr << "Running wedge sampling for triangles (" << nthreads << " threads)" << endl;
    cerr << "=============================================" << endl;

    struct timespec pre_start, pre_finish;
    double pre_elapsed;
    clock_gettime(CLOCK_MONOTONIC, &pre_start);
  
    Graph &G = GS.G;
  
    // Build sampling distribution over vertices.
    vector<long long> cumulative_weights(G.size());
    vector<vector<long long>> vertex_cumulative_weights_1(G.size());
    vector<vector<long long>> vertex_cumulative_weights_2(G.size());
    long long prev = 0;
    // todo: replace this with p means
    if (nthreads > 1) {
      omp_set_num_threads(thread::hardware_concurrency());
    } else {
      omp_set_nested(1);
    }
  
  #pragma omp parallel for
    for (int i = 0; i < (int) G.size(); i++) {
      long long vertex_weight = 0;
      vertex_cumulative_weights_1[i].reserve(G[i].size());
      vertex_cumulative_weights_2[i].reserve(G[i].size());
      for (const auto &e : G[i]) {
        vertex_weight += e.wt;
        vertex_cumulative_weights_2[i].push_back(vertex_weight);
      }
  
      long long total_weight = 0;
      for (const auto &e : G[i]) {
        total_weight += G[i].size() * e.wt + vertex_weight;
        vertex_cumulative_weights_1[i].push_back(total_weight);
      }
    }
  
    for (int i = 0; i < (int) G.size(); i++) {
      cumulative_weights[i] = vertex_cumulative_weights_2[i].back() + prev;
      prev = cumulative_weights[i];
    }
  
    // build an adjacency matrix where a(i, j) = weight of edge (i, j)
    // TODO: maybe we should lift this out of the functions and make a more general graph structure
    // vector<unordered_map<int, long long>> weight(G.size());
    // vector<google::dense_hash_map<int, long long>> weight(G.size());
    // for (int u = 0; u < (int) G.size(); u++) {
    //   weight[u].set_empty_key(-1);
    //   for (const auto &e : G[u]) {
    //     int v = e.dst;
    //     long long w = e.wt;
    //     weight[u][v] = w;
    //   }
    // }
  
    vector<thread> threads(nthreads);
    vector<vector<weighted_triangle>> counters(nthreads);
    int nsamples_per_thread = ceil(max_samples / nthreads);
  
    clock_gettime(CLOCK_MONOTONIC, &pre_finish);
    pre_elapsed = (pre_finish.tv_sec - pre_start.tv_sec);
    pre_elapsed += (pre_finish.tv_nsec - pre_start.tv_nsec) / 1000000000.0;
  
    cerr << "Pre-processing time: " << pre_elapsed << endl;
  
    struct timespec start, finish;
    double tot_time;
    clock_gettime(CLOCK_MONOTONIC, &start);
  
    auto sample_vertex = [&]() {
      long long s = rand64() % cumulative_weights.back();
      int idx = lower_bound(cumulative_weights.begin(), cumulative_weights.end(), s) - cumulative_weights.begin();
      return idx;
    };
  
    auto sample_neighbour_1 = [&](int v) {
      long long s = rand64() % vertex_cumulative_weights_1[v].back();
      int idx = lower_bound(vertex_cumulative_weights_1[v].begin(), vertex_cumulative_weights_1[v].end(), s) - vertex_cumulative_weights_1[v].begin();
      return G[v][idx];
    };
  
    auto sample_neighbour_2 = [&](int v, long long shift) {
      long long s = rand64() % (vertex_cumulative_weights_2[v].back() + shift * G[v].size());
      if (s >= vertex_cumulative_weights_2[v].back()) {
        return G[v][rand() % G[v].size()];
      } else {
        s = rand64() % vertex_cumulative_weights_2[v].back();
        int idx = lower_bound(vertex_cumulative_weights_2[v].begin(), vertex_cumulative_weights_2[v].end(), s) - vertex_cumulative_weights_2[v].begin();
        return G[v][idx];
      }
    };
  
    auto terminate = [&](int nsamples_) {
      if (max_samples != -1) {
        return nsamples_ >= nsamples_per_thread;
      } else {
        struct timespec cur;
        double tot_time;
        clock_gettime(CLOCK_MONOTONIC, &cur);
        if (include_setup) {
          tot_time = (cur.tv_sec - pre_start.tv_sec);
          tot_time += (cur.tv_nsec - pre_start.tv_nsec) / 1000000000.0;
        } else {
          tot_time = (cur.tv_sec - start.tv_sec);
          tot_time += (cur.tv_nsec - start.tv_nsec) / 1000000000.0;
        }
        return tot_time >= max_time;
      }
    };
  
    auto sampler = [&](int i) {
      int nsamples_ = 0;
      while (!terminate(nsamples_)) {
        int u = sample_vertex();
        auto ev = sample_neighbour_1(u);
        auto ew = sample_neighbour_2(u, ev.wt);
        if (ev.dst == ew.dst) continue;
        nsamples_++;
  
        if (G[ev.dst].size() > G[ew.dst].size()) {
          std::swap(ev, ew);
        }
  
        for (const auto &e : G[ev.dst]) {
          if (e.dst == ew.dst) {
            auto tri = weighted_triangle(u, ev.dst, ew.dst, ev.wt + ew.wt + e.wt);
            counters[i].push_back(tri);
            break;
          }
        }
        /*
           if (weight[ev.dst].count(ew.dst)) {
        // todo: replace with p means
        auto tri = weighted_triangle(u, ev.dst, ew.dst, ev.wt + ew.wt + weight[ev.dst][ew.dst]);
  
        bool cont = false;
        for (int j = 0; j < nthreads; j++) {
        if (histories[j].count(tri)) {
        cont = true;
        break;
        }
        }
        if (cont) continue;
        histories[i].insert(tri);
        counters[i].push_back(tri);
        }
         */
      }
      sort(counters[i].begin(), counters[i].end());
    };
  
    auto parallel_merger = [&](int i, int j) {
      vector<weighted_triangle> W;
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

    if (nthreads > 1) {
      for (int i = 0; i < nthreads; i++) {
        thread th(sampler, i);
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
        counters.push_back(vector<weighted_triangle>());
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
    } else {
      sampler(0);
    }
  
    cerr << "Found " << counters[0].size() << " triangles." << endl;
    if (counters[0].size()) cerr << "The maximum weight triangle was " << *counters[0].begin() << endl;
  
    clock_gettime(CLOCK_MONOTONIC, &finish);
    tot_time = (finish.tv_sec - start.tv_sec);
    tot_time += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    cerr << "Total Time (s): " << tot_time << endl;
    // cerr << "Time per sample (s): " << tot_time / nsamples << endl;
    cerr << endl;
  
    return counters[0];
  }

  vector<weighted_triangle> path_sampler(GraphStruct &GS, int nthreads=1, int max_samples=-1, double max_time=-1, bool include_setup=true) {
    cerr << "=============================================" << endl;
    cerr << "Running path sampling for triangles (" << nthreads << " threads)" << endl;
    cerr << "=============================================" << endl;

    struct timespec pre_start, pre_finish;
    double pre_elapsed;
    clock_gettime(CLOCK_MONOTONIC, &pre_start);
  
    Graph &G = GS.G;
    vector<full_edge> &edges = GS.edges;
  
    vector<double> weight_sum(G.size());
    vector<vector<long long>> node_sums(G.size());

    if (nthreads > 1) {
      omp_set_num_threads(thread::hardware_concurrency());
    } else {
      omp_set_nested(1);
    }

  #pragma omp parallel for
    for (int u = 0; u < (int) G.size(); u++) {
      sort(G[u].begin(), G[u].end());
  
      long long prev = 0;
      for (auto e : G[u]) {
        weight_sum[u] += e.wt;
        node_sums[u].push_back(prev + e.wt);
        prev = node_sums[u].back();
      }
    }
  
    vector<double> sum_edge_weight(GS.m);
    int count = 0;
    double prev = 0;
    for (auto e : edges) {
      double weight = e.wt * (weight_sum[e.src] - e.wt) * (weight_sum[e.dst] - e.wt);
      sum_edge_weight[count] = weight + prev;
      prev = sum_edge_weight[count++];
    }
  
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0, sum_edge_weight.back());
    auto sample_edge = [&]() {
      double s = distribution(generator);
      int idx = lower_bound(sum_edge_weight.begin(), sum_edge_weight.end(), s) - sum_edge_weight.begin();
      return edges[idx];
    };
  
    clock_gettime(CLOCK_MONOTONIC, &pre_finish);
    pre_elapsed = (pre_finish.tv_sec - pre_start.tv_sec);
    pre_elapsed += (pre_finish.tv_nsec - pre_start.tv_nsec) / 1000000000.0;
  
    cerr << "Pre-processing time: " << pre_elapsed << endl;
  
    struct timespec start, finish;
    double tot_time;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // set<int> sorted_nodes;
    auto sample_neighbour = [&](int node, int exclude) {
      int idx;
      // if (sorted_nodes.count(node) == 0) {
      //   sort(G[node].begin(), G[node].end());
      //   sorted_nodes.insert(node);
      // }
      // int exclude_idx = lower_bound(G[node].begin(), G[node].end(), half_edge{exclude, 0}) - G[node].begin();
      int exclude_idx = lower_bound(G[node].begin(), G[node].end(), half_edge{exclude, numeric_limits<long long>::max()}) - G[node].begin();

      // Decide whether to sample left side or right side.
      int s = rand64() % (node_sums[node].back() - G[node][exclude_idx].wt);
      if (exclude_idx == 0 || s >= node_sums[node][exclude_idx-1]) {
        // right side
        s = rand64() % (node_sums[node].back() - node_sums[node][exclude_idx]);
        idx = lower_bound(node_sums[node].begin() + exclude_idx, node_sums[node].end(), node_sums[node][exclude_idx] + s) - node_sums[node].begin();
      } else {
        // left side
        s = rand64() % node_sums[node][exclude_idx-1];
        idx = lower_bound(node_sums[node].begin(), node_sums[node].begin() + exclude_idx, s) - node_sums[node].begin();
      }
      return G[node][idx];
    };
  
    vector<set<weighted_triangle>> histories;
    vector<thread> threads(nthreads);
    vector<vector<weighted_triangle>> counters(nthreads);
    int nsamples_per_thread = ceil(max_samples / nthreads);

    auto terminate = [&](int nsamples_) {
      if (max_samples != -1) {
        return nsamples_ >= nsamples_per_thread;
      } else {
        struct timespec cur;
        double tot_time;
        clock_gettime(CLOCK_MONOTONIC, &cur);
        if (include_setup) {
          tot_time = (cur.tv_sec - pre_start.tv_sec);
          tot_time += (cur.tv_nsec - pre_start.tv_nsec) / 1000000000.0;
        } else {
          tot_time = (cur.tv_sec - start.tv_sec);
          tot_time += (cur.tv_nsec - start.tv_nsec) / 1000000000.0;
        }
        return tot_time >= max_time;
      }
    };

    auto sampler = [&](int i) {
      int nsamples_ = 0;
      while (!(terminate(nsamples_))) {
        auto edge = sample_edge();
        int u = edge.src, v = edge.dst;
        long long w = edge.wt;

        if (G[u].size() == 1 || G[v].size() == 1) {
          continue;
        } 
    
        auto c0 = sample_neighbour(u, v);
        auto c1 = sample_neighbour(v, u);
        if (c0.dst == c1.dst) {
          nsamples_++;
          auto tri = weighted_triangle(u, v, c0.dst, c0.wt + c1.wt + w);
          counters[i].push_back(tri);
          // if (history.count(tri) == 0) {
          //   history.insert(tri);
          //   counter.insert(tri);
          // }
        }
      }
      sort(counters[i].begin(), counters[i].end());
    };

    auto parallel_merger = [&](int i, int j) {
      vector<weighted_triangle> W;
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

    if (nthreads > 1) {
      for (int i = 0; i < nthreads; i++) {
        thread th(sampler, i);
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
        counters.push_back(vector<weighted_triangle>());
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
    } else {
      sampler(0);
    }

    cerr << "Found " << counters[0].size() << " triangles." << endl;
    if (counters[0].size()) cerr << "The maximum weight triangle was " << *counters[0].begin() << endl;
  
    clock_gettime(CLOCK_MONOTONIC, &finish);
    tot_time = (finish.tv_sec - start.tv_sec);
    tot_time += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    cerr << "Total Time (s): " << tot_time << endl;
    // cerr << "Time per sample (s): " << tot_time / nsamples << endl;
    cerr << endl;
  
    return counters[0];
  }
  
  vector<weighted_triangle> edge_time_version(GraphStruct &GS, double max_time, double inc, bool include_setup=true) {
    return edge_sampler(GS, 1, -1, max_time, include_setup);
  }
  
  vector<weighted_triangle> edge_samples_version(GraphStruct &GS, int nsamples) {
    return edge_sampler(GS, 1, nsamples, -1, false);
  }
  
  vector<weighted_triangle> edge_parallel_time_version(GraphStruct &GS, int nthreads, double max_time, double inc=-1, bool include_setup=true) {
    return edge_sampler(GS, nthreads, -1, max_time, include_setup);
  }
  
  vector<weighted_triangle> edge_parallel_samples_version(GraphStruct &GS, int nthreads, int nsamples) {
    return edge_sampler(GS, nthreads, nsamples, -1);
  }
  
  vector<weighted_triangle> wedge_time_version(GraphStruct &GS, double max_time, double inc, bool include_setup=true) {
    return wedge_sampler(GS, 1, -1, max_time, include_setup);
  }
  
  vector<weighted_triangle> wedge_samples_version(GraphStruct &GS, int nsamples) {
    return wedge_sampler(GS, 1, nsamples, -1, false);
  }
  
  vector<weighted_triangle> wedge_parallel_time_version(GraphStruct &GS, int nthreads, double max_time, double inc=-1, bool include_setup=true) {
    return wedge_sampler(GS, nthreads, -1, max_time, include_setup);
  }
  
  vector<weighted_triangle> wedge_parallel_samples_version(GraphStruct &GS, int nthreads, int nsamples) {
    return wedge_sampler(GS, nthreads, nsamples, -1);
  }
  
  vector<weighted_triangle> path_time_version(GraphStruct &GS, double max_time, double inc, bool include_setup=true) {
    return path_sampler(GS, 1, -1, max_time, include_setup);
  }
  
  vector<weighted_triangle> path_samples_version(GraphStruct &GS, int nsamples) {
    return path_sampler(GS, 1, nsamples, -1, false);
  }
  
  vector<weighted_triangle> path_parallel_time_version(GraphStruct &GS, int nthreads, double max_time, double inc=-1, bool include_setup=true) {
    return path_sampler(GS, nthreads, -1, max_time, include_setup);
  }
  
  vector<weighted_triangle> path_parallel_samples_version(GraphStruct &GS, int nthreads, int nsamples) {
    return path_sampler(GS, nthreads, nsamples, -1);
  }
  
  // Static heavy-light algorithm.
  // Arguments:
  //   GS: Adjacency list and edge list.
  //   k: Parameter k for top-k. If k = -1, returns all enumerated triangles.
  //   Otherwise returns the top-k triangles.
  //   p: What fraction of edges are considered heavy.
  set<weighted_triangle> heavy_light_sampler(GraphStruct &GS, int k=-1, double p = 0.1) {
    cerr << "=============================================" << endl;
    cerr << "Running heavy light sampling for triangles" << endl;
    cerr << "=============================================" << endl;
  
    vector<full_edge> &edges = GS.edges;
  
    double st = clock();
    // Since brute_force_sampler does not need to know the number of edges, GSh.m is not updated here.
    GraphStruct GSh;
    for (int i = 0; i < (int) (p * edges.size()); i++) {
      GSh.G.resize(max(edges[i].dst+1, (int) GSh.G.size()));
      GSh.G[edges[i].src].push_back({edges[i].dst, edges[i].wt});
      GSh.G[edges[i].dst].push_back({edges[i].src, edges[i].wt});
    }
    auto counter = brute_force_sampler(GSh, k, false);
  
    cerr << "Found " << counter.size() << " triangles." << endl;
    if (counter.size()) cerr << "The maximum weight triangle was " << *counter.begin() << endl;
  
    double tot_time = (clock() - st) / CLOCKS_PER_SEC;
    cerr << "Total Time (s): " << tot_time << endl;
    cerr << endl;
  
    return counter;
  
  }
  
  // Dynamic heavy-light and auto heavy-light algorithms.
  // Arguments:
  //   GS: Adjacency list and edge list.
  //   k: Parameter k for top-k.
  //   auto_threshold: If true then runs auto heavy-light, otherwise dynamic.
  //   alpha: Parameter for dynamic heavy-light. Ignored if running auto
  //   heavy-light.
  //   keep_all: If true then returns all triangles found, otherwise only top-k.
  set<weighted_triangle> dynamic_heavy_light(GraphStruct &GS, int k = 100, bool auto_threshold = false, double alpha = 1.25, bool keep_all = false) {
    string description = "dynamic";
    if (auto_threshold) {
      description = "auto";
    }
    cerr << "=============================================" << endl;
    cerr << "Running " << description << " heavy light for triangles" << endl;
    cerr << "=============================================" << endl;
  
    double pre_st = clock();
  
    Graph &G = GS.G;
    const vector<full_edge> &edges = GS.edges;
    vector<unordered_map<int, long long>> exists(G.size());
    vector<long long> vert_to_wt(G.size());
    vector<bool> computed(G.size());
    vector<set<int>> deleted(G.size());
  
    cerr << "Precompute time (s): " << 1.0 * (clock() - pre_st)/CLOCKS_PER_SEC << endl;
  
    double st = clock();
    set<weighted_triangle> counter, topk;
    long long num_tris = 0;
  
    Graph Gh;
    int hi = 0, hj = 0;
    // Initialize threshold to infinity and insert a dummy triangle with this
    // weight.
    long long threshold = numeric_limits<long long>::max();
    counter.insert(weighted_triangle(0, 0, 0, threshold));
    auto curr = counter.begin();
  
    // Use memoization to compute exists.
    auto compute_exists_per_node = [&](int u) {
      if (!computed[u]) {
        computed[u] = true;
        for (const auto& e : G[u]) {
          if (!deleted[u].count(e.dst)) {
            exists[u][e.dst] = e.wt;
            exists[e.dst][u] = e.wt;
          }
        }
      }
    };
  
    // Return edge weight between u and v.
    auto search = [&](int u, int v) {
      if (G[u].size() > G[v].size()) std::swap(u, v);
      if (deleted[u].count(v)) return 0LL;
      if (exists[u].count(v)) return exists[u][v];
      for (const auto& e : G[u]) {
        if (e.dst == v) return e.wt;
      }
      return 0LL;
    };
  
    // If triangle is heavier than threshold insert into topk. Otherwise insert
    // into counter and update curr pointer.
    auto insert_triangle = [&](long long weight, long long threshold, weighted_triangle T) {
      if (weight >= threshold) {
        topk.insert(T);
      } else {
        auto it = counter.insert(T).first;
        if (*it < *curr) {
          curr = it;
        }
      }
    };
  
    // Used by auto to choose which pointers to move.
    int edge_i_left = 0, edge_i_right = 0;
    int edge_j_left = 0, edge_j_right = 0;
    auto move_ptrs = [&](int& l, int& r, int i) {
      if (edges[i].wt <= edges[r].wt) l = i;
      while (r < (int) edges.size() && edges[r].wt == edges[l].wt) r++;
    };
    double delta_ei = 0, delta_ej = 0, ei_cost = 0, ej_cost = 0;
  
    while ((int) topk.size() < k+1 && hj < (int) edges.size()) {
      if (auto_threshold) {
        move_ptrs(edge_i_left, edge_i_right, hi);
        move_ptrs(edge_j_left, edge_j_right, hj);
      }
  
      auto ei = edges[hi], ej = edges[hj];
      Gh.resize(max(ej.dst+1, (int) Gh.size()));
      threshold = 2 * ej.wt + ei.wt;
  
      if (auto_threshold) {
        delta_ei = double(ei.wt - edges[edge_i_right].wt) / (edge_i_right - edge_i_left);
        delta_ej = double(ej.wt - edges[edge_j_right].wt) / (edge_j_right - edge_j_left);
        ei_cost = exists[ei.src].size() + exists[ei.dst].size() + 1 
          + G[ej.src].size() + G[ej.dst].size()
          - exists[ej.src].size() - exists[ej.dst].size();
        ej_cost = Gh[ej.src].size() + Gh[ej.dst].size();
      }
  
      // Auto heavy-light moves pointers based on computing the ratio of the
      // average decrease in weight and the cost of moving the corresponding
      // pointer, and moves the pointer with higher "bang-per-buck".
      // Dynamic heavy-light assumes a power law distribution and chooses which
      // pointer to move based on a parameter alpha, which is a function of the
      // power law parameter.
      bool advance_j = false;
      if (auto_threshold) {
        advance_j = delta_ej  * ei_cost >= delta_ei * ej_cost;
      } else {
        advance_j = pow(ej.wt, alpha) >= ei.wt;
      }
  
      if (hj == hi || advance_j) {
  
        // Enumerate triangles with 2 heavy edges.
        for (const auto& e : Gh[ej.src]) {
          long long wt = search(e.dst, ej.dst);
          if (wt) {
            long long weight = ej.wt + e.wt + wt;
            weighted_triangle T(e.dst, ej.dst, ej.src, weight);
            insert_triangle(weight, threshold, T);
            num_tris++;
          }
        }
  
        // Enumerate triangles with 2 heavy edges.
        for (const auto& e : Gh[ej.dst]) {
          long long wt = search(e.dst, ej.src);
          if (wt) {
            long long weight = ej.wt + e.wt + wt;
            weighted_triangle T(e.dst, ej.dst, ej.src, weight);
            insert_triangle(weight, threshold, T);
            num_tris++;
          }
        }
  
        // Enumerate triangles with 3 heavy edges.
        for (const auto& e : Gh[ej.src]) {
          vert_to_wt[e.dst] = e.wt;
        }
        for (const auto& e : Gh[ej.dst]) {
          if (vert_to_wt[e.dst]) {
            long long weight = ej.wt + e.wt + vert_to_wt[e.dst];
            weighted_triangle T(e.dst, ej.src, ej.dst, weight);
            insert_triangle(weight, threshold, T);
            num_tris++;
          }
        }
        for (const auto& e : Gh[ej.src]) {
          vert_to_wt[e.dst] = 0;
        }
  
        hj++;
  
        // Remove ej from light edges.
        deleted[ej.src].insert(ej.dst);
        deleted[ej.dst].insert(ej.src);
        Gh[ej.src].push_back({ej.dst, ej.wt});
        Gh[ej.dst].push_back({ej.src, ej.wt});
      } else {
  
        // Enumerate triangles with 1 heavy edge.
        compute_exists_per_node(ei.src);
        compute_exists_per_node(ei.dst);
        for (const auto& kv : exists[ei.src]) {
          vert_to_wt[kv.first] = kv.second;
        }
        for (const auto& kv : exists[ei.dst]) {
          if (vert_to_wt[kv.first]) {
            long long weight = ei.wt + kv.second + vert_to_wt[kv.first];
            weighted_triangle T(kv.first, ei.src, ei.dst, weight);
            insert_triangle(weight, threshold, T);
            num_tris++;
          }
        }
        for (const auto& kv : exists[ei.src]) {
          vert_to_wt[kv.first] = 0;
        }
        hi++;
      }
  
      // Add enumerated triangles heavier than threshold to topk.
      while (curr != counter.end() && curr->weight >= threshold) {
        topk.insert(*curr);
        auto prev = curr;
        curr++;
        if (prev != counter.begin()) {
          counter.erase(prev);
        }
      }
      if (curr == counter.end()) {
        curr--;
      }
    }
  
    // Removing the dummy triangle of weight INF. There should be one
    // in topk as well. So topk actually has one fewer triangle than it reports.
    counter.erase(counter.begin());
    topk.erase(topk.begin());
  
    cerr << "Found " << num_tris << " triangles." << endl;
    cerr << "Out of these, the top " << int(topk.size()) << " are found for sure." << endl;
    if (topk.size()) cerr << "The maximum weight triangle was " << *topk.begin() << endl;
  
    double tot_time = (clock() - st) / CLOCKS_PER_SEC;
    cerr << "Total Time (s): " << tot_time << endl;
    cerr << endl;
  
    if (keep_all) {
      for (const auto& T : counter) {
        topk.insert(T);
      }
    }
    return topk;
  
  }
  
  bool custom_find(const vector<weighted_triangle> &v,
      const weighted_triangle &T) {
    return binary_search(v.begin(), v.end(), T);
  }
  
  bool custom_find(const set<weighted_triangle> &s,
      const weighted_triangle &T) {
    return s.count(T);
  }
  
  // U should be vector or set of weighted triangle.
  template<class U>
  void compare_statistics(set<weighted_triangle> &all_triangles,
      U &sampled_triangles,
      int K,
      bool check_k=false) {
    cerr << "=============================================" << endl;
    cerr << "Comparing statistics" << endl;
    cerr << "=============================================" << endl;
  
    int num_found = 0;
    int curr_tri = 0;
    bool first_break = 0;
    int bidx = 0;
    vector<double> breakpoints({0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0 - 1e-6});
  
    int k = min(K, (int)sampled_triangles.size());
    vector<long long> ranks(k);
    // vector<long long> top_sampled_weights(k), top_true_weights(k);
    // vector<long double> percentiles(k);
  
    set<long long> unique_weights;
    for (const auto &T : all_triangles) {
      unique_weights.insert(T.weight);
    }
    vector<long long> weights(unique_weights.begin(), unique_weights.end());
  
    omp_set_num_threads(thread::hardware_concurrency());
    omp_set_nested(1);
    __gnu_parallel::sort(weights.rbegin(), weights.rend());
  
    if (check_k) {
      for (const auto &T : sampled_triangles) {
        num_found++;
        if (num_found < k+1) {
          ranks[num_found-1] = lower_bound(weights.begin(), weights.end(), T.weight, greater<long long>()) - weights.begin() + 1;
        } else {
          break;
        }
      }
    } else {
      for (const auto &T : all_triangles) {
        if (custom_find(sampled_triangles, T)) {
          num_found++;
          if (num_found < k+1) {
            ranks[num_found-1] = lower_bound(weights.begin(), weights.end(), T.weight, greater<long long>()) - weights.begin() + 1;
            /*
               percentiles[num_found-1] = 1.0 - (long double) (curr_tri+1.0)/all_triangles.size();
               top_sampled_weights[num_found-1] = T.weight;
             */
          }
        }
        curr_tri++;
        /*
           if (curr_tri < k+1) {
           top_true_weights[curr_tri-1] = T.weight;
           }
         */
  
        if (num_found != curr_tri && !first_break) {
          first_break = true;
          cerr << "Found top " << 100.0 * num_found / all_triangles.size() << " (" << num_found << ") percent of weighted triangles." << endl;
        }
  
        if (bidx < (int) breakpoints.size() && curr_tri == int(breakpoints[bidx] * all_triangles.size())) {
          cerr << "Found " << 100.0 * num_found / curr_tri << " percent of weighted triangles top " << int(breakpoints[bidx] * 100 + 1e-3) <<"%." << endl;
          bidx++;
        }
      }
    }
    /*
       cerr << "=============================================" << endl;
       cerr << "Ranks of top " << k << " triangles" << endl;
       cerr << "=============================================" << endl;
       for (auto rank : ranks) {
       cerr << rank << " ";
       }
       cerr << endl;
     */
  
    /*
       cerr << "=============================================" << endl;
       cerr << "Percentiles of top " << k << " triangles" << endl;
       cerr << "=============================================" << endl;
       for (auto percentile: percentiles) {
       cerr << percentile << " ";
       }
       cerr << endl;
     */
  
    long double acc = 0.0;
    vector<long long> all_weights;
    for (const auto &T : all_triangles) {
      all_weights.push_back(T.weight);
      if ((int) all_weights.size() == k) break;
    }
    // all_weights.reserve(all_triangles.size());
    // for (const auto &T : all_triangles) all_weights.push_back(T.weight);
    set<long long> weights_set(all_weights.begin(), all_weights.end());
    map<long long, long long> cnt_all, cnt_sampled;
    vector<long long> top_weights(all_weights.begin(), all_weights.begin()+k);
    for (const auto &w : all_weights) cnt_all[w]++;
    for (const auto &T : sampled_triangles) cnt_sampled[T.weight]++;
    for (const auto &w : weights_set) acc += min(cnt_all[w], cnt_sampled[w]);
    acc /= k;
    cerr << "=============================================" << endl;
    cerr << "Accuracy: " << acc << endl;
    cerr << "=============================================" << endl;
  
    // for (const auto &w : top_weights) cerr << w << " "; cerr << endl;
    // for (const auto &T : sampled_triangles) cerr << T.weight << " "; cerr << endl;
  
    cerr << endl;
  }
  
  void compare_statistics_time(set<weighted_triangle> &all_triangles,
      vector<set<weighted_triangle>> &vec_sampled_triangles, 
      vector<double> times, int K, bool check_k=false) {
    cerr << "=============================================" << endl;
    cerr << "Comparing sampling statistics" << endl;
    cerr << "=============================================" << endl;
  
    vector<double> breakpoints({0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0 - 1e-6});
  
    set<long long> unique_weights;
    for (const auto &T : all_triangles) {
      unique_weights.insert(T.weight);
    }
    vector<long long> weights(unique_weights.begin(), unique_weights.end());
  
    omp_set_num_threads(thread::hardware_concurrency());
    omp_set_nested(1);
    __gnu_parallel::sort(weights.rbegin(), weights.rend());
  
    set<weighted_triangle> sampled_triangles;
    for (int i = 0; i < (int) times.size(); i++) {
      cerr << "=============================================" << endl;
      cerr << "Time: " << times[i] << endl;
      for (const auto &T : vec_sampled_triangles[i]) {
        sampled_triangles.insert(T);
      }
      int num_found = 0;
      int curr_tri = 0;
      bool first_break = 0;
      int bidx = 0;
      int k = min(K, (int)sampled_triangles.size());
      vector<long long> ranks(k);
  
      if (!check_k) {
        for (auto T : all_triangles) {
          if (sampled_triangles.count(T)) {
            num_found++;
            if (num_found < k+1) {
              ranks[num_found-1] = lower_bound(weights.begin(), weights.end(), T.weight, greater<long long>()) - weights.begin() + 1;
            }
          }
          curr_tri++;
  
          // // to speed up finding time and inc for datasets, remove later
          // if (num_found == K || sampled_triangles.size() == 0) break;
  
          if (num_found != curr_tri && !first_break) {
            first_break = true;
            cerr << "Found top " << 100.0 * num_found / all_triangles.size() << " (" << num_found << ") percent of weighted triangles." << endl;
          }
  
          if (bidx < (int) breakpoints.size() && curr_tri == int(breakpoints[bidx] * all_triangles.size())) {
            cerr << "Found " << 100.0 * num_found / curr_tri << " percent of weighted triangles top " << int(breakpoints[bidx] * 100 + 1e-3) <<"%." << endl;
            bidx++;
          }
        }
      } else {
      }
  
      long double acc = 0.0;
      vector<long long> all_weights;
      for (const auto &T : all_triangles) {
        all_weights.push_back(T.weight);
        if ((int) all_weights.size() == k) break;
      }
      // all_weights.reserve(all_triangles.size());
      // for (const auto &T : all_triangles) all_weights.push_back(T.weight);
      set<long long> weights_set(all_weights.begin(), all_weights.end());
      map<long long, long long> cnt_all, cnt_sampled;
      vector<long long> top_weights(all_weights.begin(), all_weights.begin()+k);
      for (const auto &w : all_weights) cnt_all[w]++;
      for (const auto &T : sampled_triangles) cnt_sampled[T.weight]++;
      for (const auto &w : weights_set) acc += min(cnt_all[w], cnt_sampled[w]);
      acc /= k;
      cerr << "=============================================" << endl;
      cerr << "Accuracy: " << acc << endl;
      cerr << "=============================================" << endl;
  
      // for (const auto &w : top_weights) cerr << w << " "; cerr << endl;
      // for (const auto &T : sampled_triangles) cerr << T.weight << " "; cerr << endl;
    }
  
    cerr << endl;
  }

}

#endif

/*
  vector<weighted_triangle> edge_sampler_parallel_everything(GraphStruct &GS, int nthreads, int max_samples=-1, double max_time=-1, double inc=-1, bool include_setup=true) {
    cerr << "=============================================" << endl;
    cerr << "Running parallel edge sampling for triangles (" << nthreads << " threads)" << endl;
    cerr << "=============================================" << endl;
    struct timespec pre_start, pre_finish;
    double pre_elapsed;
    clock_gettime(CLOCK_MONOTONIC, &pre_start);
  
    Graph &G = GS.G;
    const vector<full_edge>& edges = GS.edges;
    long long total_edge_weight = 0;
    vector<int> weight_index;
    vector<long long> weight_value;
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
        cur = lower_bound(edges.begin() + cur, edges.end(), full_edge(0, 0, cur_wt), greater<full_edge>()) - edges.begin();
      }
      total_edge_weight += (cur - weight_index.back()) * cur_wt;
      weight_value.push_back(total_edge_weight);
    }
    weight_index.push_back(cur);
  
    vector<thread> threads(nthreads);
    vector<vector<weighted_triangle>> counters(nthreads);
    vector<set<pair<int, int>>> histories(nthreads);
    int nsamples_per_thread = ceil(max_samples / nthreads);
  
    clock_gettime(CLOCK_MONOTONIC, &pre_finish);
    pre_elapsed = (pre_finish.tv_sec - pre_start.tv_sec);
    pre_elapsed += (pre_finish.tv_nsec - pre_start.tv_nsec) / 1000000000.0;
  
    cerr << "Pre-processing time: " << pre_elapsed << endl;
    cerr << "Edge weight classes: " << int(weight_index.size())-1 << endl;
    cerr << "Total edge weight: " << total_edge_weight << endl;
  
    struct timespec start, finish;
    double tot_time;
    clock_gettime(CLOCK_MONOTONIC, &start);
  
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
  
    auto sample_single_edge = [&](){
      long long s = rand64() % total_edge_weight;
      int index = lower_bound(weight_value.begin(), weight_value.end(), s) - weight_value.begin();
      int e = rand() % (weight_index[index+1] - weight_index[index]);
      return edges[e + weight_index[index]];
    };
  
    auto terminate = [&](int nsamples_) {
      if (max_samples != -1) {
        return nsamples_ >= nsamples_per_thread;
      } else {
        struct timespec cur;
        double tot_time;
        clock_gettime(CLOCK_MONOTONIC, &cur);
        if (include_setup) {
          tot_time = (cur.tv_sec - pre_start.tv_sec);
          tot_time += (cur.tv_nsec - pre_start.tv_nsec) / 1000000000.0;
        } else {
          tot_time = (cur.tv_sec - start.tv_sec);
          tot_time += (cur.tv_nsec - start.tv_nsec) / 1000000000.0;
        }
        return tot_time >= max_time;
      }
    };
  
    auto parallel_sampler = [&](int i) {
      vector<int> sample_index;
      if (max_samples != -1) {
        sample_index = batched_sample_edges(max_samples);
      }
      int nsamples_ = 0;
      while (!terminate(nsamples_)) {
        full_edge e;
        if (max_samples != -1) {
          e = edges[sample_index[nsamples_]];
        } else {
          e = sample_single_edge();
        }
        nsamples_++;
        int u = e.src, v = e.dst;
        long long w = e.wt;
        bool cont = false;
        for (int j = 0; j < nthreads; j++) {
          if (histories[j].count(make_pair(u, v))) {
            cont = true;
            break;
          }
        }
        if (cont) continue;
        histories[i].insert(make_pair(u, v));
        unordered_map<int, long long> vert_to_wt;
        for (const auto &eu : G[u]) {
          vert_to_wt[eu.dst] = eu.wt;
        }
  
        for (const auto &ev : G[v]) {
          if (vert_to_wt.count(ev.dst)) {
            counters[i].push_back(weighted_triangle(u, v, ev.dst, ev.wt + vert_to_wt[ev.dst] + w));
          }
        }
      }
      sort(counters[i].begin(), counters[i].end());
    };
  
    auto parallel_merger = [&](int i, int j) {
      vector<weighted_triangle> W;
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
      counters.push_back(vector<weighted_triangle>());
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
  
    cerr << "Found " << counters[0].size() << " triangles." << endl;
    if (counters[0].size()) cerr << "The maximum weight triangle was " << *counters[0].begin() << endl;
  
    clock_gettime(CLOCK_MONOTONIC, &finish);
    tot_time = (finish.tv_sec - start.tv_sec);
    tot_time += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    cerr << "Total Time (s): " << tot_time << endl;
    // cerr << "Time per sample (s): " << tot_time / nsamples << endl;
    cerr << endl;
  
    return counters[0];
  }

  set<weighted_triangle> edge_sampler_everything(GraphStruct &GS, int max_samples=-1, double max_time=-1, double inc=-1, bool include_setup=true) {
    cerr << "=============================================" << endl;
    cerr << "Running edge sampling for triangles" << endl;
    cerr << "=============================================" << endl;

    double pre_st = clock();

    Graph &G = GS.G;
    const vector<full_edge>& edges = GS.edges;
    long long total_edge_weight = 0;
    vector<int> weight_index;
    vector<long long> weight_value;
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
        cur = lower_bound(edges.begin() + cur, edges.end(), full_edge(0, 0, cur_wt), greater<full_edge>()) - edges.begin();
      }
      total_edge_weight += (cur - weight_index.back()) * cur_wt;
      weight_value.push_back(total_edge_weight);
    }
    weight_index.push_back(cur);

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

    auto sample_single_edge = [&](){
      long long s = rand64() % total_edge_weight;
      int index = lower_bound(weight_value.begin(), weight_value.end(), s) - weight_value.begin();
      int e = rand() % (weight_index[index+1] - weight_index[index]);
      return edges[e + weight_index[index]];
    };

    set<weighted_triangle> counter;
    set<pair<int, int>> history;
    int nsamples = 0;
    double init_time = include_setup? pre_st : st;

    auto terminate = [&]() {
      if (max_samples != -1) {
        return nsamples >= max_samples;
      } else {
        double tot_time = (clock() - init_time) / CLOCKS_PER_SEC;
        return tot_time >= max_time;
      }
    };

    vector<int> sample_index;
    if (max_samples != -1) {
      sample_index = batched_sample_edges(max_samples);
    }
    while (!(terminate())) {
      full_edge e;
      if (max_samples != -1) {
        e = edges[sample_index[nsamples]];
      } else {
        e = sample_single_edge();
      }
      nsamples++;
      int u = e.src, v = e.dst;
      long long w = e.wt;
      // resampling isnt an issue from experimentation
      if (history.count(make_pair(u, v))) {
        //cerr << "RESAMPLED!!" << endl;
        continue;
      }
      history.insert(make_pair(u, v));
      map<int, long long> vert_to_wt;
      for (const auto &eu : G[u]) {
        vert_to_wt[eu.dst] = eu.wt;
      }

      for (const auto &ev : G[v]) {
        if (vert_to_wt.count(ev.dst)) {
          // todo: replace with p means
          counter.insert(weighted_triangle(u, v, ev.dst, ev.wt + vert_to_wt[ev.dst] + w));
        }
      }
    }

    cerr << "Found " << counter.size() << " triangles in counter." << endl;
    if (counter.size()) cerr << "The maximum weight triangle was " << *counter.begin() << endl;
    cerr << "Total Time (s): " << 1.0 * (clock() - pre_st) / CLOCKS_PER_SEC << endl;
    return counter;

  }

  set<weighted_triangle> path_sampler_everything(GraphStruct &GS, int max_samples=-1, double max_time=-1, double inc=-1, bool include_setup=true) {
    cerr << "=============================================" << endl;
    cerr << "Running path sampling for triangles" << endl;
    cerr << "=============================================" << endl;
  
    double pre_st = clock();
  
    Graph &G = GS.G;
    vector<full_edge> &edges = GS.edges;
  
    vector<double> weight_sum(G.size());
    vector<vector<long long>> node_sums(G.size());
    for (int u = 0; u < (int) G.size(); u++) {
      sort(G[u].begin(), G[u].end());
  
      long long prev = 0;
      for (auto e : G[u]) {
        weight_sum[u] += e.wt;
        node_sums[u].push_back(prev + e.wt);
        prev = node_sums[u].back();
      }
    }
  
    vector<double> sum_edge_weight(GS.m);
    int count = 0;
    double prev = 0;
    for (auto e : edges) {
      double weight = e.wt * (weight_sum[e.src] - e.wt) * (weight_sum[e.dst] - e.wt);
      // sum_edge_weight.push_back(weight + prev);
      // prev = sum_edge_weight.back();
      sum_edge_weight[count] = weight + prev;
      prev = sum_edge_weight[count++];
    }
  
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0, sum_edge_weight.back());
    auto sample_edge = [&]() {
      double s = distribution(generator);
      int idx = lower_bound(sum_edge_weight.begin(), sum_edge_weight.end(), s) - sum_edge_weight.begin();
      return edges[idx];
    };
    cerr << "Precompute time (s): " << 1.0 * (clock() - pre_st)/CLOCKS_PER_SEC << endl;
  
    double st = clock();
  
    auto sample_neighbour = [&](int node, int exclude) {
      int idx;
      int exclude_idx = lower_bound(G[node].begin(), G[node].end(), half_edge{exclude, 0}) - G[node].begin();
  
      // decide whether to sample left side or right side.
      int s = rand64() % (node_sums[node].back() - G[node][exclude_idx].wt);
      if (exclude_idx == 0 || s >= node_sums[node][exclude_idx-1]) {
        // right side
        s = rand64() % (node_sums[node].back() - node_sums[node][exclude_idx]);
        idx = lower_bound(node_sums[node].begin() + exclude_idx, node_sums[node].end(), node_sums[node][exclude_idx] + s) - node_sums[node].begin();
      } else {
        // left side
        s = rand64() % node_sums[node][exclude_idx-1];
        idx = lower_bound(node_sums[node].begin(), node_sums[node].begin() + exclude_idx, s) - node_sums[node].begin();
      }
      return G[node][idx];
    };
  
    set<weighted_triangle> counter, history;
    int nsamples = 0;
    double init_time = include_setup? pre_st : st;
  
    auto terminate = [&]() {
      if (max_samples != -1) {
        return nsamples >= max_samples;
      } else {
        double tot_time = (clock() - init_time) / CLOCKS_PER_SEC;
        return tot_time >= max_time;
      }
    };
  
    while (!(terminate())) {
      auto edge = sample_edge();
      int u = edge.src, v = edge.dst;
      long long w = edge.wt;
  
      auto c0 = sample_neighbour(u, v);
      auto c1 = sample_neighbour(v, u);
      if (c0.dst == c1.dst) {
        nsamples++;
        auto tri = weighted_triangle(u, v, c0.dst, c0.wt + c1.wt + w);
        if (history.count(tri) == 0) {
          history.insert(tri);
          counter.insert(tri);
        }
      }
    }
  
    cerr << "Found " << counter.size() << " triangles in counter." << endl;
    if (counter.size()) cerr << "The maximum weight triangle was " << *counter.begin() << endl;
    cerr << "Total Time (s): " << 1.0 * (clock() - pre_st) / CLOCKS_PER_SEC << endl;
    return counter;
  
  }

  set<weighted_triangle> wedge_sampler_everything(GraphStruct &GS, int max_samples=-1, double max_time=-1, double inc=-1, bool include_setup=true) {
    cerr << "=============================================" << endl;
    cerr << "Running wedge sampling for triangles" << endl;
    cerr << "=============================================" << endl;

    double pre_st = clock();

    Graph &G = GS.G;

    // build sampling distribution over vertices
    vector<long long> cumulative_weights(G.size());
    vector<vector<long long>> vertex_cumulative_weights_1(G.size());
    vector<vector<long long>> vertex_cumulative_weights_2(G.size());
    long long prev = 0;
    // todo: replace this with p means
    for (int i = 0; i < (int) G.size(); i++) {
      long long vertex_weight = 0;
      vertex_cumulative_weights_1[i].reserve(G[i].size());
      vertex_cumulative_weights_2[i].reserve(G[i].size());
      for (const auto &e : G[i]) {
        vertex_weight += e.wt;
        vertex_cumulative_weights_2[i].push_back(vertex_weight);
      }

      long long total_weight = 0;
      for (const auto &e : G[i]) {
        total_weight += G[i].size() * e.wt + vertex_weight;
        vertex_cumulative_weights_1[i].push_back(total_weight);
      }
      // cumulative_weights.push_back(total_weight);
      // cumulative_weights[cumulative_weights.size() - 1] += prev;
      // prev = cumulative_weights.back();
      cumulative_weights[i] = total_weight + prev;
      prev = cumulative_weights[i];
    }

    // build an adjacency matrix where a(i, j) = weight of edge (i, j)
    // TODO: maybe we should lift this out of the functions and make a more general graph structure
    // vector<unordered_map<int, long long>> weight(G.size());
    // vector<google::dense_hash_map<int, long long>> weight(G.size());
    // for (int u = 0; u < (int) G.size(); u++) {
    //   weight[u].set_empty_key(-1);
    //   for (const auto &e : G[u]) {
    //     int v = e.dst;
    //     long long w = e.wt;
    //     weight[u][v] = w;
    //   }
    // }
    cerr << "Precompute time (s): " << 1.0 * (clock() - pre_st)/CLOCKS_PER_SEC << endl;

    double st = clock();

    auto sample_vertex = [&]() {
      long long s = rand64() % cumulative_weights.back();
      int idx = lower_bound(cumulative_weights.begin(), cumulative_weights.end(), s) - cumulative_weights.begin();
      return idx;
    };

    auto sample_neighbour_1 = [&](int v) {
      long long s = rand64() % vertex_cumulative_weights_1[v].back();
      int idx = lower_bound(vertex_cumulative_weights_1[v].begin(), vertex_cumulative_weights_1[v].end(), s) - vertex_cumulative_weights_1[v].begin();
      return G[v][idx];
    };

    auto sample_neighbour_2 = [&](int v, long long shift) {
      long long s = rand64() % (vertex_cumulative_weights_2[v].back() + shift * G[v].size());
      if (s >= vertex_cumulative_weights_2[v].back()) {
        return G[v][rand() % G[v].size()];
      } else {
        s = rand64() % vertex_cumulative_weights_2[v].back();
        int idx = lower_bound(vertex_cumulative_weights_2[v].begin(), vertex_cumulative_weights_2[v].end(), s) - vertex_cumulative_weights_2[v].begin();
        return G[v][idx];
      }
    };

    set<weighted_triangle> counter, history;
    int nsamples = 0;
    double init_time = include_setup? pre_st : st;

    auto terminate = [&]() {
      if (max_samples != -1) {
        return nsamples >= max_samples;
      } else {
        double tot_time = (clock() - init_time) / CLOCKS_PER_SEC;
        return tot_time >= max_time;
      }
    };

    while (!(terminate())) {
      int u = sample_vertex();
      auto ev = sample_neighbour_1(u);
      auto ew = sample_neighbour_2(u, ev.wt);
      if (ev.dst == ew.dst) continue;
      nsamples++;

      if (G[ev.dst].size() > G[ew.dst].size()) {
        std::swap(ev, ew);
      }

      for (const auto &e : G[ev.dst]) {
        if (e.dst == ew.dst) {
          auto tri = weighted_triangle(u, ev.dst, ew.dst, ev.wt + ew.wt + e.wt);
          if (history.count(tri) == 0) {
            history.insert(tri);
            counter.insert(tri);
          }
        }
      }
//          if (weight[ev.dst].count(ew.dst)) {
// // todo: replace with p means
// auto tri = weighted_triangle(u, ev.dst, ew.dst, ev.wt + ew.wt + weight[ev.dst][ew.dst]);
// if (history.count(tri) == 0) {
// history.insert(tri);
// counter.insert(tri);
// }
// }
    }


    cerr << "Found " << counter.size() << " triangles in counter." << endl;
    if (counter.size()) cerr << "The maximum weight triangle was " << *counter.begin() << endl;
    cerr << "Total Time (s): " << 1.0 * (clock() - pre_st) / CLOCKS_PER_SEC << endl;
    return counter;

  }

  vector<weighted_triangle> wedge_sampler_parallel_everything(GraphStruct &GS, int nthreads, int max_samples=-1, double max_time=-1, double inc=-1, bool include_setup=true) {
    cerr << "=============================================" << endl;
    cerr << "Running parallel wedge sampling for triangles (" << nthreads << " threads)" << endl;
    cerr << "=============================================" << endl;
    struct timespec pre_start, pre_finish;
    double pre_elapsed;
    clock_gettime(CLOCK_MONOTONIC, &pre_start);
  
    Graph &G = GS.G;
  
    // build sampling distribution over vertices
    vector<long long> cumulative_weights(G.size());
    vector<vector<long long>> vertex_cumulative_weights_1(G.size());
    vector<vector<long long>> vertex_cumulative_weights_2(G.size());
    long long prev = 0;
    // todo: replace this with p means
    omp_set_num_threads(thread::hardware_concurrency());
    //omp_set_nested(1);
  
  #pragma omp parallel for
    for (int i = 0; i < (int) G.size(); i++) {
      long long vertex_weight = 0;
      vertex_cumulative_weights_1[i].reserve(G[i].size());
      vertex_cumulative_weights_2[i].reserve(G[i].size());
      for (const auto &e : G[i]) {
        vertex_weight += e.wt;
        vertex_cumulative_weights_2[i].push_back(vertex_weight);
      }
  
      long long total_weight = 0;
      for (const auto &e : G[i]) {
        total_weight += G[i].size() * e.wt + vertex_weight;
        vertex_cumulative_weights_1[i].push_back(total_weight);
      }
    }
  
    for (int i = 0; i < (int) G.size(); i++) {
      cumulative_weights[i] = vertex_cumulative_weights_2[i].back() + prev;
      prev = cumulative_weights[i];
    }
  
    // build an adjacency matrix where a(i, j) = weight of edge (i, j)
    // TODO: maybe we should lift this out of the functions and make a more general graph structure
    // vector<unordered_map<int, long long>> weight(G.size());
    // vector<google::dense_hash_map<int, long long>> weight(G.size());
    // for (int u = 0; u < (int) G.size(); u++) {
    //   weight[u].set_empty_key(-1);
    //   for (const auto &e : G[u]) {
    //     int v = e.dst;
    //     long long w = e.wt;
    //     weight[u][v] = w;
    //   }
    // }
  
    vector<thread> threads(nthreads);
    vector<vector<weighted_triangle>> counters(nthreads);
    int nsamples_per_thread = ceil(max_samples / nthreads);
  
    clock_gettime(CLOCK_MONOTONIC, &pre_finish);
    pre_elapsed = (pre_finish.tv_sec - pre_start.tv_sec);
    pre_elapsed += (pre_finish.tv_nsec - pre_start.tv_nsec) / 1000000000.0;
  
    cerr << "Pre-processing time: " << pre_elapsed << endl;
  
    struct timespec start, finish;
    double tot_time;
    clock_gettime(CLOCK_MONOTONIC, &start);
  
    auto sample_vertex = [&]() {
      long long s = rand64() % cumulative_weights.back();
      int idx = lower_bound(cumulative_weights.begin(), cumulative_weights.end(), s) - cumulative_weights.begin();
      return idx;
    };
  
    auto sample_neighbour_1 = [&](int v) {
      long long s = rand64() % vertex_cumulative_weights_1[v].back();
      int idx = lower_bound(vertex_cumulative_weights_1[v].begin(), vertex_cumulative_weights_1[v].end(), s) - vertex_cumulative_weights_1[v].begin();
      return G[v][idx];
    };
  
    auto sample_neighbour_2 = [&](int v, long long shift) {
      long long s = rand64() % (vertex_cumulative_weights_2[v].back() + shift * G[v].size());
      if (s >= vertex_cumulative_weights_2[v].back()) {
        return G[v][rand() % G[v].size()];
      } else {
        s = rand64() % vertex_cumulative_weights_2[v].back();
        int idx = lower_bound(vertex_cumulative_weights_2[v].begin(), vertex_cumulative_weights_2[v].end(), s) - vertex_cumulative_weights_2[v].begin();
        return G[v][idx];
      }
    };
  
    auto terminate = [&](int nsamples_) {
      if (max_samples != -1) {
        return nsamples_ >= nsamples_per_thread;
      } else {
        struct timespec cur;
        double tot_time;
        clock_gettime(CLOCK_MONOTONIC, &cur);
        if (include_setup) {
          tot_time = (cur.tv_sec - pre_start.tv_sec);
          tot_time += (cur.tv_nsec - pre_start.tv_nsec) / 1000000000.0;
        } else {
          tot_time = (cur.tv_sec - start.tv_sec);
          tot_time += (cur.tv_nsec - start.tv_nsec) / 1000000000.0;
        }
        return tot_time >= max_time;
      }
    };
  
    auto parallel_sampler = [&](int i) {
      int nsamples_ = 0;
      while (!terminate(nsamples_)) {
        int u = sample_vertex();
        auto ev = sample_neighbour_1(u);
        auto ew = sample_neighbour_2(u, ev.wt);
        if (ev.dst == ew.dst) continue;
        nsamples_++;
  
        if (G[ev.dst].size() > G[ew.dst].size()) {
          std::swap(ev, ew);
        }
  
        for (const auto &e : G[ev.dst]) {
          if (e.dst == ew.dst) {
            auto tri = weighted_triangle(u, ev.dst, ew.dst, ev.wt + ew.wt + e.wt);
            counters[i].push_back(tri);
            break;
          }
        }
        //    if (weight[ev.dst].count(ew.dst)) {
        // // todo: replace with p means
        // auto tri = weighted_triangle(u, ev.dst, ew.dst, ev.wt + ew.wt + weight[ev.dst][ew.dst]);
  
        // bool cont = false;
        // for (int j = 0; j < nthreads; j++) {
        // if (histories[j].count(tri)) {
        // cont = true;
        // break;
        // }
        // }
        // if (cont) continue;
        // histories[i].insert(tri);
        // counters[i].push_back(tri);
        // }
        //
      }
      sort(counters[i].begin(), counters[i].end());
    };
  
    auto parallel_merger = [&](int i, int j) {
      vector<weighted_triangle> W;
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
      counters.push_back(vector<weighted_triangle>());
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
  
    cerr << "Found " << counters[0].size() << " triangles." << endl;
    if (counters[0].size()) cerr << "The maximum weight triangle was " << *counters[0].begin() << endl;
  
    clock_gettime(CLOCK_MONOTONIC, &finish);
    tot_time = (finish.tv_sec - start.tv_sec);
    tot_time += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    cerr << "Total Time (s): " << tot_time << endl;
    // cerr << "Time per sample (s): " << tot_time / nsamples << endl;
    cerr << endl;
  
    return counters[0];
  }

*/  