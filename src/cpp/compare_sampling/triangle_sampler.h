#ifndef TRIANGLE_SAMPLER_H
#define TRIANGLE_SAMPLER_H

#include <bits/stdc++.h>

#include "graph.h"

using namespace std;

namespace wsdm_2019_graph {

// TODO: implement sample without replacement

set<weighted_triangle> brute_force_sampler(Graph& G, bool diagnostic = true) {
	if (diagnostic) {
		cerr << "=============================================" << endl;
		cerr << "Running brute force detection for triangles" << endl;
		cerr << "=============================================" << endl;
	}
	double st = clock();

	// Can use degeneracy ordering to speed up brute force
	// However, on certain tests it didnt have great performance.
	set<weighted_triangle> counter;

	// TODO: Can replace maps with normal arrays in most of these procedures.
	// Will only be useful in single threaded case.
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
					counter.insert(weighted_triangle(u, v, ev.dst, val));
				}
			}
		}

		for (auto e : G[u]) {
			if (e.dst < u) continue;
			vert_to_wt[e.dst] = 0;
		}
	}
	if (diagnostic) {
		cerr << "Found " << counter.size() << " triangles." << endl;
		if (counter.size()) cerr << "The maximum weight triangle was " << *counter.begin() << endl;

		double tot_time = (clock() - st) / CLOCKS_PER_SEC;
		cerr << "Total Time (s): " << tot_time << endl;
		cerr << endl;
	}
	return counter;
}

pair<vector<set<weighted_triangle>>, vector<double>> edge_sampler(Graph &G, double max_time, double inc) {
	cerr << "=============================================" << endl;
	cerr << "Running edge sampling for triangles" << endl;
	cerr << "=============================================" << endl;
	double st = clock();

	// build distribution over edges
	map<int, vector<full_edge>> edge_distribution;
	for (int u = 0; u < (int) G.size(); u++) {
		for (auto e : G[u]) {
			int v = e.dst;
			long long w = e.wt;
			if (u > v) continue;
			edge_distribution[e.wt].push_back({u, v, w});
		}
	}

	vector<long long> cumulative_weights;
	map<int, long long> index_to_weight;
	int count = 0;
	long long prev = 0;
	// todo: replace this with p means
	for (const auto& kv : edge_distribution) {
		//cerr << kv.first << " " << kv.second.size() << endl;
		cumulative_weights.push_back(kv.second.size() * kv.first);
		cumulative_weights[cumulative_weights.size() - 1] += prev;
		index_to_weight[count++] = kv.first;
		prev = cumulative_weights.back();
	}
	cerr << "Edge weight classes: " << edge_distribution.size() << endl;
	cerr << "Total edge weight: " << cumulative_weights.back() << endl;

	auto sample_edge = [&](){
		long long s = rand64() % cumulative_weights.back();
		int idx = lower_bound(cumulative_weights.begin(), cumulative_weights.end(), s) - cumulative_weights.begin();
		//cerr << "sampled weight: " << s << endl;
		//cerr << "sampled index: " << idx << " " << index_to_weight[idx] << endl;

		long long weight = index_to_weight[idx];
		auto& edges = edge_distribution[weight];
		return edges[rand() % edges.size()];
	};

	set<weighted_triangle> counter;
	set<pair<int, int>> history;
	vector<set<weighted_triangle>> counters;
	vector<double> times;
	int nsamples = 0;
	double last_time = 0;

	while (1) {
		double tot_time = (clock() - st) / CLOCKS_PER_SEC;
		if (tot_time >= max_time) {
			break;
		}
		if (tot_time - last_time >= inc) {
			last_time = tot_time;
			counters.push_back(counter);
			counter.clear();
			times.push_back(tot_time);
		}
		auto e = sample_edge();
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
		for (auto eu : G[u]) {
			vert_to_wt[eu.dst] = eu.wt;
		}

		for (auto ev : G[v]) {
			if (vert_to_wt.count(ev.dst)) {
				// todo: replace with p means
				counter.insert(weighted_triangle(u, v, ev.dst, ev.wt + vert_to_wt[ev.dst] + w));
			}
		}

	}

	return make_pair(counters, times);

}

pair<vector<set<weighted_triangle>>, vector<double>> wedge_sampler(Graph &G, double max_time, double inc) {
	cerr << "=============================================" << endl;
	cerr << "Running wedge sampling for triangles" << endl;
	cerr << "=============================================" << endl;
	double st = clock();

	// build sampling distribution over vertices
	vector<long long> cumulative_weights;
	vector<vector<long long>> vertex_cumulative_weights_1(G.size());
	vector<vector<long long>> vertex_cumulative_weights_2(G.size());
	long long prev = 0;
	// todo: replace this with p means
	for (int i = 0; i < (int) G.size(); i++) {
		long long vertex_weight = 0;
		for (auto e : G[i]) {
			vertex_weight += e.wt;
			vertex_cumulative_weights_2[i].push_back(vertex_weight);
		}

		long long total_weight = 0;
		for (auto e : G[i]) {
			total_weight += G[i].size() * e.wt + vertex_weight;
			vertex_cumulative_weights_1[i].push_back(total_weight);
		}
		cumulative_weights.push_back(total_weight);
		cumulative_weights[cumulative_weights.size() - 1] += prev;
		prev = cumulative_weights.back();
	}

	// build an adjacency matrix where a(i, j) = weight of edge (i, j)
	// TODO: maybe we should lift this out of the functions and make a more general graph structure
	map<int, map<int, long long>> weight;
	for (int u = 0; u < (int) G.size(); u++) {
		for (const auto &e : G[u]) {
			int v = e.dst;
			long long w = e.wt;
			weight[u][v] = w;
		}
	}

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
	vector<set<weighted_triangle>> counters;
	vector<double> times;
	int nsamples = 0;
	double last_time = 0;

	while(1) {
		double tot_time = (clock() - st) / CLOCKS_PER_SEC;
		if (tot_time >= max_time) {
			break;
		}
		if (tot_time - last_time >= inc) {
			last_time = tot_time;
			counters.push_back(counter);
			counter.clear();
			times.push_back(tot_time);
		}
		int u = sample_vertex();
		auto ev = sample_neighbour_1(u);
		auto ew = sample_neighbour_2(u, ev.wt);
		if (ev.dst == ew.dst) continue;
		nsamples++;

		if (weight[ev.dst].count(ew.dst)) {
			// todo: replace with p means
			auto tri = weighted_triangle(u, ev.dst, ew.dst, ev.wt + ew.wt + weight[ev.dst][ew.dst]);
			if (history.count(tri) == 0) {
				history.insert(tri);
				counter.insert(tri);
			}
		}
	}

	return make_pair(counters, times);

}

pair<vector<set<weighted_triangle>>, vector<double>> path_sampler(Graph &G, double max_time, double inc) {
	cerr << "=============================================" << endl;
	cerr << "Running path sampling for triangles" << endl;
	cerr << "=============================================" << endl;
	double st = clock();

	vector<full_edge> edges;
	map<int, double> weight_sum;
	vector<vector<long long>> node_sums(G.size());
	for (int u = 0; u < (int) G.size(); u++) {
		sort(G[u].begin(), G[u].end());

		long long prev = 0;
		for (auto e : G[u]) {
			weight_sum[u] += e.wt;
			if (u < e.dst) {
				edges.push_back({u, e.dst, e.wt});
			}
			node_sums[u].push_back(prev + e.wt);
			prev = node_sums[u].back();
		}
	}

	vector<double> sum_edge_weight;
	double prev = 0;
	for (auto e : edges) {
		double weight = e.wt * (weight_sum[e.src] - e.wt) * (weight_sum[e.dst] - e.wt);
		sum_edge_weight.push_back(weight + prev);
		prev = sum_edge_weight.back();
	}

	default_random_engine generator;
	uniform_real_distribution<double> distribution(0.0, sum_edge_weight.back());
	auto sample_edge = [&]() {
		double s = distribution(generator);
		int idx = lower_bound(sum_edge_weight.begin(), sum_edge_weight.end(), s) - sum_edge_weight.begin();
		return edges[idx];
	};

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
	vector<set<weighted_triangle>> counters;
	vector<double> times;
	int nsamples = 0;
	double last_time = 0;

	while (1) {
		double tot_time = (clock() - st) / CLOCKS_PER_SEC;
		if (tot_time >= max_time) {
			break;
		}
		if (tot_time - last_time >= inc) {
			last_time = tot_time;
			counters.push_back(counter);
			counter.clear();
			times.push_back(tot_time);
		}
	}

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

	return make_pair(counters, times);

}

void compare_statistics_old(set<weighted_triangle> &all_triangles,
						vector<set<weighted_triangle>> &vec_sampled_triangles,
						vector<double> times, int K) {
	cerr << "=============================================" << endl;

	set<long long> unique_weights;
	for (const auto &T : all_triangles) {
		unique_weights.insert(T.weight);
	}
	vector<long long> weights(unique_weights.begin(), unique_weights.end());
	sort(weights.rbegin(), weights.rend());

	set<weighted_triangle> sampled_triangles;
	vector<int> ranks;
	for (int i = 0; i < times.size(); i++) {
		// sampled_triangles.merge(vec_sampled_triangles[i]);
		for (const auto &T : vec_sampled_triangles[i]) {
			sampled_triangles.insert(T);
		}
		int num_found = 0;
		int k = min(K, (int)sampled_triangles.size());
		long double recall = 0.0;
		if (k == 0) {
			recall = -1;
		} else {
			ranks.clear(); ranks.resize(k);
			for (const auto &T : sampled_triangles) {
				ranks[num_found] = lower_bound(weights.begin(), weights.end(), T.weight, greater<long long>()) - weights.begin() + 1;
				num_found++;
				if (num_found == k) {
					break;
				}
			}
			for (int i = 0; i < k; i++) {
				recall += (ranks[i] <= K);
			}
			recall /= K;
		}

		cerr << "Time: " << times[i] << ", recall: " << recall << endl;

	}

	cerr << "=============================================" << endl;

}

void compare_statistics(set<weighted_triangle> &all_triangles,
						vector<set<weighted_triangle>> &vec_sampled_triangles, 
						vector<double> times, int K) {
	cerr << "=============================================" << endl;
	cerr << "Comparing sampling statistics" << endl;
	cerr << "=============================================" << endl;

	vector<double> breakpoints({0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0 - 1e-6});

	set<long long> unique_weights;
	for (const auto &T : all_triangles) {
		unique_weights.insert(T.weight);
	}
	vector<long long> weights(unique_weights.begin(), unique_weights.end());
	sort(weights.rbegin(), weights.rend());

	set<weighted_triangle> sampled_triangles;
	for (int i = 0; i < times.size(); i++) {
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

		for (auto T : all_triangles) {
			if (sampled_triangles.count(T)) {
				num_found++;
				if (num_found < k+1) {
					ranks[num_found-1] = lower_bound(weights.begin(), weights.end(), T.weight, greater<long long>()) - weights.begin() + 1;
				}
			}
			curr_tri++;

			if (num_found != curr_tri && !first_break) {
				first_break = true;
				cerr << "Found top " << 100.0 * num_found / all_triangles.size() << " (" << num_found << ") percent of weighted triangles." << endl;
			}

			if (bidx < (int) breakpoints.size() && curr_tri == int(breakpoints[bidx] * all_triangles.size())) {
				cerr << "Found " << 100.0 * num_found / curr_tri << " percent of weighted triangles top " << int(breakpoints[bidx] * 100 + 1e-3) <<"%." << endl;
				bidx++;
			}
		}

		long double recall = 0.0;
		if (sampled_triangles.size() == 0) {
			recall = -1;
		} else {
			for (int i = 0; i < k; i++) {
				recall += (ranks[i] <= k);
			}
		}

		recall /= K;
		cerr << "Recall: " << recall << endl;
		cerr << "=============================================" << endl;
	}

	cerr << endl;
}

}

#endif
