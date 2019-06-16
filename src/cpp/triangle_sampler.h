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

set<weighted_triangle> edge_sampler(Graph& G, int nsamples) {
	cerr << "=============================================" << endl;
	cerr << "Running edge sampling for triangles" << endl;
	cerr << "=============================================" << endl;
	// double st = clock();
	time_t st;
	time(&st);

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
		long long s = rand() % cumulative_weights.back();
		int idx = lower_bound(cumulative_weights.begin(), cumulative_weights.end(), s) - cumulative_weights.begin();
		//cerr << "sampled weight: " << s << endl;
		//cerr << "sampled index: " << idx << " " << index_to_weight[idx] << endl;

		long long weight = index_to_weight[idx];
		auto& edges = edge_distribution[weight];
		return edges[rand() % edges.size()];
	};

	set<weighted_triangle> counter;
	set<pair<int, int>> history;
	for (int samp = 0; samp < nsamples; samp++) {
		auto e = sample_edge();
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
	cerr << "Found " << counter.size() << " triangles." << endl;
	if (counter.size()) cerr << "The maximum weight triangle was " << *counter.begin() << endl;

	// double tot_time = (clock() - st) / CLOCKS_PER_SEC;
	time_t en;
	time(&en);
	double tot_time = difftime(en, st);
	cerr << "Total Time (s): " << tot_time << endl;
	cerr << "Time per sample (s): " << tot_time / nsamples << endl;
	cerr << endl;

	return counter;
}

set<weighted_triangle> edge_sampler_parallel(Graph &G, int nsamples, int nthreads) {
	cerr << "=============================================" << endl;
	cerr << "Running parallel edge sampling for triangles" << endl;
	cerr << "=============================================" << endl;
	// double st = clock();
	time_t st;
	time(&st);

	// build distribution over edges
	map<int, vector<full_edge>> edge_distribution;
	for (int u = 0; u < (int) G.size(); u++) {
		for (const auto &e : G[u]) {
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
	for (const auto &kv : edge_distribution) {
		cumulative_weights.push_back(kv.second.size() * kv.first);
		cumulative_weights[cumulative_weights.size() - 1] += prev;
		index_to_weight[count++] = kv.first;
		prev = cumulative_weights.back();
	}
	cerr << "Edge weight classes: " << edge_distribution.size() << endl;
	cerr << "Total edge weight: " << cumulative_weights.back() << endl;

	vector<thread> threads(nthreads);
	vector<set<weighted_triangle>> counters(nthreads);
	vector<set<pair<int, int>>> histories(nthreads);
	int nsamples_per_thread = ceil(nsamples / nthreads);

	auto sample_edge = [&](){
		long long s = rand() % cumulative_weights.back();
		int idx = lower_bound(cumulative_weights.begin(), cumulative_weights.end(), s) - cumulative_weights.begin();

		long long weight = index_to_weight[idx];
		auto &edges = edge_distribution[weight];
		return edges[rand() % edges.size()];
	};

	auto parallel_sampler = [&](int i){
		for (int samp = 0; samp < nsamples_per_thread; samp++) {
			auto e = sample_edge();
			int u = e.src, v = e.dst;
			long long w = e.wt;

			if (histories[i].count(make_pair(u, v))) {
				continue;
			}
			histories[i].insert(make_pair(u, v));
			map<int, long long> vert_to_wt;
			for (auto eu : G[u]) {
				vert_to_wt[eu.dst] = eu.wt;
			}

			for (auto ev : G[v]) {
				if (vert_to_wt.count(ev.dst)) {
					counters[i].insert(weighted_triangle(u, v, ev.dst, ev.wt + vert_to_wt[ev.dst] + w));
				}
			}
		}
	};

	for (int i = 0; i < nthreads; i++) {
		thread th(parallel_sampler, i);
		threads[i] = move(th);
	}
	for (int i = 0; i < nthreads; i++) {
		threads[i].join();
	}

	set<weighted_triangle> counter;
	for (const auto &c : counters) {
		for (const auto &t : c) {
			counter.insert(t);
		}
	}

	cerr << "Found " << counter.size() << " triangles." << endl;
	if (counter.size()) cerr << "The maximum weight triangle was " << *counter.begin() << endl;

	// double tot_time = (clock() - st) / CLOCKS_PER_SEC;
	time_t en;
	time(&en);
	double tot_time = difftime(en, st);
	cerr << "Total Time (s): " << tot_time << endl;
	cerr << "Time per sample (s): " << tot_time / nsamples << endl;
	cerr << endl;

	return counter;
}

// 4x speedup possible with an optimal random number implementation
// this isnt enough to get asymptotically past
set<weighted_triangle> path_sampler(Graph& G, int nsamples) {
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
		int s = rand() % (node_sums[node].back() - G[node][exclude_idx].wt);
		if (exclude_idx == 0 || s >= node_sums[node][exclude_idx-1]) {
			// right side
			s = rand() % (node_sums[node].back() - node_sums[node][exclude_idx]);
			idx = lower_bound(node_sums[node].begin() + exclude_idx, node_sums[node].end(), node_sums[node][exclude_idx] + s) - node_sums[node].begin();
		} else {
			// left side
			s = rand() % node_sums[node][exclude_idx-1];
			idx = lower_bound(node_sums[node].begin(), node_sums[node].begin() + exclude_idx, s) - node_sums[node].begin();
		}
		return G[node][idx];
	};

	set<weighted_triangle> counter;
	for (int samps = 0; samps < nsamples; samps++) {
		auto edge = sample_edge();
		int u = edge.src, v = edge.dst;
		long long w = edge.wt;

		auto c0 = sample_neighbour(u, v);
		auto c1 = sample_neighbour(v, u);
		if (c0.dst == c1.dst) {
			counter.insert(weighted_triangle(u, v, c0.dst, c0.wt + c1.wt + w));
		}

		/*
		map<int, long long> seen;
		for (int s = 0; s < G[v].size(); s++) {
			auto c0 = sample_neighbour(u, v);
			seen[c0.dst] = c0.wt;
		}

		for (int s = 0; s < G[v].size(); s++) {
			auto c1 = G[v][s];//sample_neighbour(v, u);
			if (seen.count(c1.dst)) {
				// todo: use product weight here
				counter.insert(weighted_triangle(u, v, c1.dst, seen[c1.dst] + c1.wt + w));
			}
		}
		*/
	}
	cerr << "Found " << counter.size() << " triangles." << endl;
	if (counter.size()) cerr << "The maximum weight triangle was " << *counter.begin() << endl;

	double tot_time = (clock() - st) / CLOCKS_PER_SEC;
	cerr << "Total Time (s): " << tot_time << endl;
	cerr << "Time per sample (s): " << tot_time / nsamples << endl;
	cerr << endl;

	return counter;
}

set<weighted_triangle> heavy_light_sampler(Graph& G, double p = 0.1) {
	cerr << "=============================================" << endl;
	cerr << "Running heavy light sampling for triangles" << endl;
	cerr << "=============================================" << endl;
	double st = clock();

	vector<full_edge> edges;
	for (int i = 0; i < (int) G.size(); i++) {
		for (auto e : G[i]) {
			if (i < e.dst) {
				edges.push_back({i, e.dst, e.wt});
			}
		}
	}
	sort(edges.rbegin(), edges.rend());

	Graph Gh, Gl;
	for (int i = 0; i < (int) (p * edges.size()); i++) {
		Gh.resize(max(edges[i].dst+1, (int) Gh.size()));
		Gh[edges[i].src].push_back({edges[i].dst, edges[i].wt});
		Gh[edges[i].dst].push_back({edges[i].src, edges[i].wt});
	}
	auto counter = brute_force_sampler(Gh, false);

/*
	for (int i = (int) (p * edges.size()); i < (int) edges.size(); i++) {
		Gl.resize(max(edges[i].dst+1, (int) Gl.size()));
		Gl[edges[i].src].push_back({edges[i].dst, edges[i].wt});
		Gl[edges[i].dst].push_back({edges[i].src, edges[i].wt});
	}
	auto counter2 = brute_force_sampler(Gl, false);
	for (auto t : counter2) counter.insert(t);
//*/

///*
	for (int i = (int) (p * edges.size()); i < (int) edges.size(); i++) {
		int u = edges[i].src, v = edges[i].dst;
		long long w = edges[i].wt;
		if (u >= (int) Gh.size() || v >= (int) Gh.size()) continue;
		map<int, long long> seen;
		for (auto e : Gh[u]) {
			seen[e.dst] = e.wt;
		}

		for (auto e : Gh[v]) {
			if (seen.count(e.dst)) {
				counter.insert(weighted_triangle(u, v, e.dst, e.wt + seen[e.dst] + w));
			}
		}
	}
//*/

/*
	for (int i = 0; i < (int) (p * edges.size()); i++) {
		int u = edges[i].src, v = edges[i].dst;
		long long w = edges[i].wt;
		if (u >= (int) Gl.size() || v >= (int) Gl.size()) continue;
		map<int, long long> seen;
		for (auto e : Gl[u]) {
			seen[e.dst] = e.wt;
		}

		for (auto e : Gl[v]) {
			if (seen.count(e.dst)) {
				counter.insert(weighted_triangle(u, v, e.dst, e.wt + seen[e.dst] + w));
			}
		}
	}
//*/
	cerr << "Found " << counter.size() << " triangles." << endl;
	if (counter.size()) cerr << "The maximum weight triangle was " << *counter.begin() << endl;

	double tot_time = (clock() - st) / CLOCKS_PER_SEC;
	cerr << "Total Time (s): " << tot_time << endl;
	cerr << endl;

	return counter;

}

set<weighted_triangle> adaptive_heavy_light(Graph& G, int k = 100) {
	cerr << "=============================================" << endl;
	cerr << "Running adaptive heavy light for triangles" << endl;
	cerr << "=============================================" << endl;
	double st = clock();

	vector<full_edge> edges;
	vector<map<int, long long>> exists(G.size());
	for (int i = 0; i < (int) G.size(); i++) {
		for (auto e : G[i]) {
			if (i < e.dst) {
				edges.push_back({i, e.dst, e.wt});
			}
			exists[i][e.dst] = e.wt;
			exists[e.dst][i] = e.wt;
		}
	}
	sort(edges.rbegin(), edges.rend());

	set<weighted_triangle> counter, topk;

	Graph Gh;
	int hi = 0, hj = 0;
	// double threshold = numeric_limits<double>::max();
	long long threshold = numeric_limits<long long>::max();
	counter.insert(weighted_triangle(0, 0, 0, threshold));
	auto curr = counter.begin();
	while ((int) topk.size() < k+1 && hj < (int) edges.size()) {
		// Version where we use a threshold
		auto ei = edges[hi], ej = edges[hj];
		Gh.resize(max(ej.dst+1, (int) Gh.size()));

		// TODO: COMPUTE THE MAGIC EXPONENT THROUGH SOME THEORY
		// This exponent should be roughly 2 - O(poly(1/beta)), where beta is 
		// the exponent of the power law governing the edge weights.
		// Back of the envelope calculations indicate
		// 	magic = 2 - 2/(b + 1).
		// Experimental evidence suggests the exponent is around 1.7 to 2.
		// This gives a magic range of [1.25, 1.33] which fits with what 
		// we observed.
		double magic = 1.25;
		if (pow(ej.wt, magic) > ei.wt) {
			// Advance j, H2 and H3 cases (at least two heavy)
			// Check for P2s with incoming ej
			for (auto e : Gh[ej.src]) {
				if (exists[e.dst].count(ej.dst)) {
					long long weight = ej.wt + e.wt + exists[ej.dst][e.dst];
					weighted_triangle T(e.dst, ej.dst, ej.src, weight);
					if (weight >= threshold) {
						topk.insert(T);
					}
					counter.insert(T);
				}
			}
			for (auto e : Gh[ej.dst]) {
				if (exists[e.dst].count(ej.src)) {
					long long weight = ej.wt + e.wt + exists[ej.src][e.dst];
					weighted_triangle T(e.dst, ej.src, ej.dst, weight);
					if (weight >= threshold) {
						topk.insert(T);
					}
					counter.insert(T);
				}
			}

			// Check for all 3 heavy
			map<int, long long> vert_to_wt;
			for (auto e : Gh[ej.src]) {
				vert_to_wt[e.dst] = e.wt;
			}
			for (auto e : Gh[ej.dst]) {
				if (vert_to_wt.count(e.dst)) {
					long long weight = ej.wt + e.wt + vert_to_wt[e.dst];
					weighted_triangle T(e.dst, ej.src, ej.dst, weight);
					if (weight >= threshold) {
						topk.insert(T);
					}
					counter.insert(T);
				}
			}

			// Remove from light edges
			exists[ej.src].erase(ej.dst);
			exists[ej.dst].erase(ej.src);
			hj++;

			Gh[ej.src].push_back({ej.dst, ej.wt});
			Gh[ej.dst].push_back({ej.src, ej.wt});
		} else {
			// Advance i, H1 case (exactly 1 heavy)
			map<int, long long> vert_to_wt;
			for (auto kv : exists[ei.src]) {
				vert_to_wt[kv.first] = kv.second;
			}
			for (auto kv : exists[ei.dst]) {
				if (vert_to_wt.count(kv.first)) {
					long long weight = ei.wt + kv.second + vert_to_wt[kv.first];
					weighted_triangle T(kv.first, ei.src, ei.dst, weight);
					if (weight >= threshold) {
						topk.insert(T);
					}
					counter.insert(T);
				}
			}
			hi++;
		}
		
		threshold = 2 * ej.wt + ei.wt - 1;
		while (curr != counter.end() && curr->weight >= threshold) {
			topk.insert(*curr);
			curr++;
		}
		if (curr == counter.end()) {
			curr--;
		}
	}
	// Removing the dummy triangle of weight INF. There should be one
	// in topk as well. So topk actually has one fewer triangle than it reports.
	counter.erase(counter.begin());

	cerr << "Found " << counter.size() << " triangles." << endl;
	cerr << "Out of these, the top " << int(topk.size()) - 1 << " are found for sure." << endl;
	if (counter.size()) cerr << "The maximum weight triangle was " << *counter.begin() << endl;

	double tot_time = (clock() - st) / CLOCKS_PER_SEC;
	cerr << "Total Time (s): " << tot_time << endl;
	cerr << endl;

	return counter;

}

// Same as adaptive heavy light EXCEPT the magic 
// constant is computed on the fly automatically
set<weighted_triangle> auto_thresholded_heavy_light(Graph& G, int k = 100) {
	cerr << "=============================================" << endl;
	cerr << "Running auto thresholded heavy light for triangles" << endl;
	cerr << "=============================================" << endl;
	double st = clock();

	vector<full_edge> edges;
	vector<map<int, long long>> exists(G.size());
	map<long long, int> edge_distribution;
	edge_distribution[0] = 1; // Assuming positive weight edges
	for (int i = 0; i < (int) G.size(); i++) {
		for (auto e : G[i]) {
			if (i < e.dst) {
				edges.push_back({i, e.dst, e.wt});
				edge_distribution[e.wt]++;
			}
			exists[i][e.dst] = e.wt;
			exists[e.dst][i] = e.wt;
		}
	}
	sort(edges.rbegin(), edges.rend());

	set<weighted_triangle> counter, topk;

	Graph Gh;
	int hi = 0, hj = 0;
	// double threshold = numeric_limits<double>::max();
	long long threshold = numeric_limits<long long>::max();
	counter.insert(weighted_triangle(0, 0, 0, threshold));
	auto curr = counter.begin();
	while ((int) topk.size() < k+1 && hj < (int) edges.size()) {
		// Version where we use a threshold
		auto ei = edges[hi], ej = edges[hj];
		Gh.resize(max(ej.dst+1, (int) Gh.size()));

		double delta_ei = double(ei.wt - (--edge_distribution.lower_bound(ei.wt))->first) / edge_distribution[ei.wt];
		double delta_ej = double(ej.wt - (--edge_distribution.lower_bound(ej.wt))->first) / edge_distribution[ej.wt];
		double ei_cost = exists[ei.src].size() + exists[ei.dst].size();
		double ej_cost = 2 * (Gh[ej.src].size() + Gh[ej.dst].size());

		if (delta_ej / ej_cost > delta_ei / ei_cost) {
			// Advance j, H2 and H3 cases (at least two heavy)
			// Check for P2s with incoming ej
			for (auto e : Gh[ej.src]) {
				if (exists[e.dst].count(ej.dst)) {
					long long weight = ej.wt + e.wt + exists[ej.dst][e.dst];
					weighted_triangle T(e.dst, ej.dst, ej.src, weight);
					if (weight >= threshold) {
						topk.insert(T);
					}
					counter.insert(T);
				}
			}
			for (auto e : Gh[ej.dst]) {
				if (exists[e.dst].count(ej.src)) {
					long long weight = ej.wt + e.wt + exists[ej.src][e.dst];
					weighted_triangle T(e.dst, ej.src, ej.dst, weight);
					if (weight >= threshold) {
						topk.insert(T);
					}
					counter.insert(T);
				}
			}

			// Check for all 3 heavy
			map<int, long long> vert_to_wt;
			for (auto e : Gh[ej.src]) {
				vert_to_wt[e.dst] = e.wt;
			}
			for (auto e : Gh[ej.dst]) {
				if (vert_to_wt.count(e.dst)) {
					long long weight = ej.wt + e.wt + vert_to_wt[e.dst];
					weighted_triangle T(e.dst, ej.src, ej.dst, weight);
					if (weight >= threshold) {
						topk.insert(T);
					}
					counter.insert(T);
				}
			}

			// Remove from light edges
			exists[ej.src].erase(ej.dst);
			exists[ej.dst].erase(ej.src);
			hj++;

			Gh[ej.src].push_back({ej.dst, ej.wt});
			Gh[ej.dst].push_back({ej.src, ej.wt});
		} else {
			// Advance i, H1 case (exactly 1 heavy)
			map<int, long long> vert_to_wt;
			for (auto kv : exists[ei.src]) {
				vert_to_wt[kv.first] = kv.second;
			}
			for (auto kv : exists[ei.dst]) {
				if (vert_to_wt.count(kv.first)) {
					long long weight = ei.wt + kv.second + vert_to_wt[kv.first];
					weighted_triangle T(kv.first, ei.src, ei.dst, weight);
					if (weight >= threshold) {
						topk.insert(T);
					}
					counter.insert(T);
				}
			}
			hi++;
		}
		
		threshold = 2 * ej.wt + ei.wt - 1;
		while (curr != counter.end() && curr->weight >= threshold) {
			topk.insert(*curr);
			curr++;
		}
		if (curr == counter.end()) {
			curr--;
		}
	}
	// Removing the dummy triangle of weight INF. There should be one
	// in topk as well. So topk actually has one fewer triangle than it reports.
	counter.erase(counter.begin());

	cerr << "Found " << counter.size() << " triangles." << endl;
	cerr << "Out of these, the top " << int(topk.size()) - 1 << " are found for sure." << endl;
	if (counter.size()) cerr << "The maximum weight triangle was " << *counter.begin() << endl;

	double tot_time = (clock() - st) / CLOCKS_PER_SEC;
	cerr << "Total Time (s): " << tot_time << endl;
	cerr << endl;

	return counter;

}

void compare_statistics(set<weighted_triangle>& all_triangles, set<weighted_triangle>& sampled_triangles, int K) {
	cerr << "=============================================" << endl;
	cerr << "Comparing sampling statistics" << endl;
	cerr << "=============================================" << endl;

	int num_found = 0;
	int curr_tri = 0;
	bool first_break = 0;
	int bidx = 0;
	vector<double> breakpoints({0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0 - 1e-6});
	//vector<double> breakpoints({0.1, 0.25, 0.5});

	int k = min(K, (int)sampled_triangles.size());
	vector<long long> ranks(k), top_sampled_weights(k), top_true_weights(k);
	vector<long double> percentiles(k);

	set<long long> unique_weights;
	for (const auto &T : all_triangles) {
		unique_weights.insert(T.weight);
	}
	vector<long long> weights(unique_weights.begin(), unique_weights.end());
	sort(weights.rbegin(), weights.rend());

	for (auto T : all_triangles) {
		if (sampled_triangles.count(T)) {
			num_found++;
			if (num_found < k+1) {
				// ranks[num_found-1] = curr_tri+1;
				ranks[num_found-1] = lower_bound(weights.begin(), weights.end(), T.weight, greater<long long>()) - weights.begin() + 1;
				percentiles[num_found-1] = 1.0 - (long double) (curr_tri+1.0)/all_triangles.size();
				top_sampled_weights[num_found-1] = T.weight;
			}
		}
		curr_tri++;
		if (curr_tri < k+1) {
			top_true_weights[curr_tri-1] = T.weight;
		}

		if (num_found != curr_tri && !first_break) {
			first_break = true;
			cerr << "Found top " << 100.0 * num_found / all_triangles.size() << " (" << num_found << ") percent of weighted triangles." << endl;
		}

		if (bidx < (int) breakpoints.size() && curr_tri == int(breakpoints[bidx] * all_triangles.size())) {
			cerr << "Found " << 100.0 * num_found / curr_tri << " percent of weighted triangles top " << int(breakpoints[bidx] * 100 + 1e-3) <<"%." << endl;
			bidx++;
		}
	}

	cerr << "=============================================" << endl;
	cerr << "Ranks of top " << k << " triangles" << endl;
	cerr << "=============================================" << endl;
	for (auto rank : ranks) {
		cerr << rank << " ";
	}
	cerr << endl;

	cerr << "=============================================" << endl;
	cerr << "Percentiles of top " << k << " triangles" << endl;
	cerr << "=============================================" << endl;
	for (auto percentile: percentiles) {
		cerr << percentile << " ";
	}
	cerr << endl;

	long double accuracy = 0.0;
	for (int i = 0; i < k; i++) {
		accuracy += (ranks[i] <= k);
	}
	accuracy /= k;
	cerr << "=============================================" << endl;
	cerr << "Accuracy: " << accuracy << endl;
	cerr << "=============================================" << endl;

	cerr << "=============================================" << endl;
	cerr << "Weights of true top " << k << " triangles" << endl;
	cerr << "=============================================" << endl;
	for (const auto &wt : top_true_weights) {
		cerr << wt << " ";
	}
	cerr << endl;
	cerr << "=============================================" << endl;
	cerr << "Weights of sampled top " << k << " triangles" << endl;
	cerr << "=============================================" << endl;
	for (const auto &wt : top_sampled_weights) {
		cerr << wt << " ";
	}
	cerr << endl;

	cerr << endl;
}

}

#endif
