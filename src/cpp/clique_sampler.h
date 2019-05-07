#ifndef CLIQUE_SAMPLER_H
#define CLIQUE_SAMPLER_H

#include <bits/stdc++.h>
#include "clique_enumerator.h"

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

	set<weighted_triangle> counter;
	for (int u = 0; u < (int) G.size(); u++) {
		for (auto e : G[u]) {
			int v = e.dst, w = e.wt;
			if (u > v) continue;
			map<int, int> vert_to_wt;
			for (auto e : G[u]) {
				vert_to_wt[e.dst] = e.wt;
			}

			for (auto e : G[v]) {
				if (vert_to_wt.count(e.dst)) {
					// todo: replace with p means
					counter.insert(weighted_triangle(u, v, e.dst, e.wt + vert_to_wt[e.dst] + w));
				}
			}
		}
	}
	if (diagnostic) {
		cerr << "Found " << counter.size() << " triangles." << endl;
		if (counter.size()) cerr << "The maximum weight triangle was " << *counter.begin() << endl;

		double tot_time = (clock() - st) / CLOCKS_PER_SEC;
		cerr << "Total Time (s): " << tot_time << endl;
	}
	return counter;
}

set<weighted_triangle> edge_sampler(Graph& G, int nsamples) {
	cerr << "=============================================" << endl;
	cerr << "Running edge sampling for triangles" << endl;
	cerr << "=============================================" << endl;
	double st = clock();

	// build distribution over edges
	map<int, vector<full_edge>> edge_distribution;
	for (int u = 0; u < (int) G.size(); u++) {
		for (auto e : G[u]) {
			int v = e.dst, w = e.wt;
			if (u > v) continue;
			edge_distribution[e.wt].push_back({u, v, w});
		}
	}

	vector<int> cumulative_weights;
	map<int, int> index_to_weight;
	int count = 0, prev = 0;
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
		int s = rand() % cumulative_weights.back();
		int idx = lower_bound(cumulative_weights.begin(), cumulative_weights.end(), s) - cumulative_weights.begin();
		//cerr << "sampled weight: " << s << endl;
		//cerr << "sampled index: " << idx << " " << index_to_weight[idx] << endl;

		int weight = index_to_weight[idx];
		auto& edges = edge_distribution[weight];
		return edges[rand() % edges.size()];
	};

	set<weighted_triangle> counter;
	set<pair<int, int>> history;
	for (int samp = 0; samp < nsamples; samp++) {
		auto e = sample_edge();
		int u = e.src, v = e.dst, w = e.wt;
		// resampling isnt an issue from experimentation
		if (history.count(make_pair(u, v))) {
			//cerr << "RESAMPLED!!" << endl;
			continue;
		}
		history.insert(make_pair(u, v));
		map<int, int> vert_to_wt;
		for (auto e : G[u]) {
			vert_to_wt[e.dst] = e.wt;
		}

		for (auto e : G[v]) {
			if (vert_to_wt.count(e.dst)) {
				// todo: replace with p means
				counter.insert(weighted_triangle(u, v, e.dst, e.wt + vert_to_wt[e.dst] + w));
			}
		}
	}
	cerr << "Found " << counter.size() << " triangles." << endl;
	if (counter.size()) cerr << "The maximum weight triangle was " << *counter.begin() << endl;

	double tot_time = (clock() - st) / CLOCKS_PER_SEC;
	cerr << "Total Time (s): " << tot_time << endl;
	cerr << "Time per sample (s): " << tot_time / nsamples << endl;

	return counter;
}

set<weighted_triangle> path_sampler(Graph& G, int nsamples) {
	cerr << "=============================================" << endl;
	cerr << "Running path sampling for triangles" << endl;
	cerr << "=============================================" << endl;
	double st = clock();
	vector<full_edge> edges;
	map<int, double> weight_sum;
	vector<vector<int>> node_sums(G.size());
	for (int u = 0; u < (int) G.size(); u++) {
		sort(G[u].begin(), G[u].end());

		int prev = 0;
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

		// decide weather to sample left side or right side.
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
		int u = edge.src, v = edge.dst, w = edge.wt;

		auto c0 = sample_neighbour(u, v);
		auto c1 = sample_neighbour(v, u);
		if (c0.dst == c1.dst) {
			counter.insert(weighted_triangle(u, v, c0.dst, c0.wt + c1.wt + w));
		}

		/*
		map<int, int> seen;
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

	return counter;
}

set<weighted_triangle> heavy_light_sampler(Graph& G, double p = 0.1) {
	cerr << "=============================================" << endl;
	cerr << "Running heavy light sampling for triangles" << endl;
	cerr << "=============================================" << endl;
	double st = clock();

	vector<full_edge> edges;
	for (int i = 0; i < (int) G[i].size(); i++) {
		for (auto e : G[i]) {
			if (i < e.dst) {
				edges.push_back({i, e.dst, e.wt});
			}
		}
	}
	sort(edges.rbegin(), edges.rend());

	Graph Gh;
	for (int i = 0; i < (int) (p * edges.size()); i++) {
		Gh.resize(max(edges[i].dst+1, (int) Gh.size()));
		Gh[edges[i].src].push_back({edges[i].dst, edges[i].wt});
		Gh[edges[i].dst].push_back({edges[i].src, edges[i].wt});
	}
	auto counter = brute_force_sampler(Gh, false);

	for (int i = (int) (p * edges.size()); i < (int) edges.size(); i++) {
		int u = edges[i].src, v = edges[i].dst, w = edges[i].wt;
		if (u >= (int) Gh.size() || v >= (int) Gh.size()) continue;
		map<int, int> seen;
		for (auto e : Gh[u]) {
			seen[e.dst] = e.wt;
		}

		for (auto e : Gh[v]) {
			if (seen.count(e.dst)) {
				counter.insert(weighted_triangle(u, v, e.dst, e.wt + seen[e.dst] + w));
			}
		}
	}

	cerr << "Found " << counter.size() << " triangles." << endl;
	if (counter.size()) cerr << "The maximum weight triangle was " << *counter.begin() << endl;

	double tot_time = (clock() - st) / CLOCKS_PER_SEC;
	cerr << "Total Time (s): " << tot_time << endl;

	return counter;

}

void compare_statistics(set<weighted_triangle>& all_triangles, set<weighted_triangle>& sampled_triangles) {
	cerr << "=============================================" << endl;
	cerr << "Comparing sampling statistics" << endl;
	cerr << "=============================================" << endl;

	int num_found = 0;
	int curr_tri = 0;
	bool first_break = 0;
	int bidx = 0;
	vector<double> breakpoints({0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0 - 1e-6});
	//vector<double> breakpoints({0.1, 0.25, 0.5});

	for (auto T : all_triangles) {
		if (sampled_triangles.count(T)) {
			num_found++;
		}
		curr_tri++;

		if (num_found != curr_tri && !first_break) {
			first_break = true;
			cerr << "Found top " << 100.0 * num_found / all_triangles.size() << " percent of weighted triangles." << endl;
		}

		if (bidx < (int) breakpoints.size() && curr_tri == int(breakpoints[bidx] * all_triangles.size())) {
			cerr << "Found " << 100.0 * num_found / curr_tri << " percent of weighted triangles top " << int(breakpoints[bidx] * 100 + 1e-6) <<"%." << endl;
			bidx++;
		}
	}
}


set<weighted_clique> enumerate_cliques(Graph& G, int k) {
	set<weighted_clique> retval;
	if (k == 1) {
		for (int i = 0; i < (int) G.size(); i++) {
			retval.insert(weighted_clique(vector<int>({i}), 0));
		}
		return retval;
	} else if (k == 2) {
		for (int u = 0; u < (int) G.size(); u++) {
			for (auto e : G[u]) {
				int v = e.dst, w = e.wt;
				retval.insert(weighted_clique(vector<int>({u, v}), w));
			}
		}
		return retval;
	}

	// todo: make this its own function
	// prune the graph based on degeneracy
	vector<int> degree(G.size());
	vector<bool> removed(G.size());
	queue<int> q;
	for (int i = 0; i < (int) G.size(); i++) {
		removed[i] = false;
		degree[i] = G[i].size();

		if (degree[i] < k-1) {
			q.push(i);
			removed[i] = true;
		}
	}
	while (!q.empty()) {
		int u = q.front();
		q.pop();

		for (auto e : G[u]) {
			degree[e.dst]--;
			if (!removed[e.dst] && degree[e.dst] < k-1) {
				removed[e.dst] = true;
				q.push(e.dst);
			}
		}
	}

	map<int, map<int, int>> adjmat;
	for (int i = 0; i < (int) G.size(); i++) {
		if (removed[i]) continue;
		for (auto e : G[i]) {
			if (removed[e.dst]) continue;
			adjmat[i][e.dst] = e.wt;
		}
	}

	for (int i = 0; i < (int) G.size(); i++) {
		if (removed[i] || degree[i] < k-1) continue;

		int nedges = 0;
		Graph subgraph;

		map<int, int> label, unlabel;
		int cidx = 0;
		for (auto e0 : G[i]) {
			int u = e0.dst;
			if (removed[u]) continue;
			if (!label.count(u)) {
				label[u] = cidx;
				unlabel[cidx] = u;
				cidx++;
			}

			for (auto e1 : G[i]) {
				int v = e1.dst;
				if (removed[v]) continue;
				if (v <= u) continue;

				if (!label.count(v)) {
					label[v] = cidx;
					unlabel[cidx] = v;
					cidx++;
				}
				
				if (label.size() > subgraph.size()) {
					subgraph.resize(label.size());
				}

				if (adjmat[u].count(v)) {
					int w = adjmat[u][v];
					subgraph[label[u]].push_back({label[v], w});
					subgraph[label[v]].push_back({label[u], w});
					nedges++;
				}
			}
		}
		subgraph.resize(label.size());

		if (nedges < (k-1) * (k-2) / 2) continue;

		auto cliques = enumerate_cliques(subgraph, k-1);
		for (auto clique : cliques) {
			// unlabelling phase
			set<int> seen;
			for (int& v : clique.vertices) {
				v = unlabel[v];
				seen.insert(v);
			}
			clique.vertices.push_back(i);

			for (auto e : G[i]) {
				if (seen.count(e.dst)) {
					clique.weight += e.wt;
				}
			}
			sort(clique.vertices.begin(), clique.vertices.end());
			retval.insert(clique);
		}

		removed[i] = true;

		for (auto e : G[i]) {
			degree[e.dst]--;
		}
	}
	return retval;
}

set<weighted_clique> clique_sampler(Graph& G, int k, int nsamples) {
	cerr << "=============================================" << endl;
	cerr << "Running edge sampling for k-cliques" << endl;
	cerr << "=============================================" << endl;
	double st = clock();

	// Prune the graph first
	vector<int> degree(G.size());
	vector<bool> removed(G.size());
	{
		queue<int> q;
		for (int i = 0; i < (int) G.size(); i++) {
			removed[i] = false;
			degree[i] = G[i].size();

			if (degree[i] < k-1) {
				q.push(i);
				removed[i] = true;
			}
		}
		while (!q.empty()) {
			int u = q.front();
			q.pop();

			for (auto e : G[u]) {
				degree[e.dst]--;
				if (!removed[e.dst] && degree[e.dst] < k-1) {
					removed[e.dst] = true;
					q.push(e.dst);
				}
			}
		}
	}

	// build distribution over edges
	map<int, vector<full_edge>> edge_distribution;
	map<int, map<int, int>> adjmat;
	for (int u = 0; u < (int) G.size(); u++) {
		if (removed[u]) continue;
		for (auto e : G[u]) {
			int v = e.dst, w = e.wt;
			if (u > v) continue;
			if (removed[v]) continue;
			edge_distribution[e.wt].push_back({u, v, w});
			adjmat[u][v] = adjmat[v][u] = w;
		}
	}

	vector<int> cumulative_weights;
	map<int, int> index_to_weight;
	int count = 0, prev = 0;
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
		int s = rand() % cumulative_weights.back();
		int idx = lower_bound(cumulative_weights.begin(), cumulative_weights.end(), s) - cumulative_weights.begin();
		//cerr << "sampled weight: " << s << endl;
		//cerr << "sampled index: " << idx << " " << index_to_weight[idx] << endl;

		int weight = index_to_weight[idx];
		auto& edges = edge_distribution[weight];
		return edges[rand() % edges.size()];
	};

	set<weighted_clique> counter;
	set<pair<int, int>> history;
	for (int samp = 0; samp < nsamples; samp++) {
		auto sample = sample_edge();
		int u = sample.src, v = sample.dst, w = sample.wt;
		if (history.count(make_pair(u, v))) {
			//cerr << "RESAMPLED!!" << endl;
			continue;
		}
		history.insert(make_pair(u, v));
		map<int, int> vert_to_wt_u, vert_to_wt_v;
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
			if (!label.count(u_)) {
				label[u_] = cidx;
				unlabel[cidx] = u_;
				cidx++;
			}
			for (int v_ : common_nbrs) {
				if (v_ <= u_) continue;
				if (!label.count(v_)) {
					label[v_] = cidx;
					unlabel[cidx] = v_;
					cidx++;
				}
				
				if (adjmat[u_].count(v_)) {
					int w_ = adjmat[u_][v_];
					subgraph[label[u_]].push_back({label[v_], w_});
					subgraph[label[v_]].push_back({label[u_], w_});
					nedges++;
				}
			}
		}

		if (nedges < (k-2) * (k-3) / 2) continue;
		set<weighted_clique> cliques;
		if (k < 5) {
			cliques = enumerate_cliques(subgraph, k-2);
		} else {
			cliques = find_cliques(subgraph, k-2);
		}
#if 0
		/*
		cerr << "Computing on subgraph with " << common_nbrs.size() << " nodes and " << nedges << " edges" << endl;
		int max_seen = 0;
		static map<int, int> freq;
		for (int u : common_nbrs) {
			freq[u]++;
			max_seen = max(freq[u], max_seen);
		}
		cerr << "Current max: " << max_seen << endl;
		*/

		// Lets do a two-way pruning:
		// After we find the (k-2)-cliques corresponding to this edge,
		// we can take each clique and try to complete it into a k-clique.
		// Then we can delete all of the edges here.
		for (auto clique : cliques) {
			// unlabelling phase
			map<int, int> seen, wsum;
			vector<int> nbrs;
			for (int& vert : clique.vertices) {
				vert = unlabel[vert];
				for (auto e : G[vert]) {
					seen[e.dst]++;
					wsum[e.dst] += e.wt;
					if (seen[e.dst] == k-2) {
						nbrs.push_back(e.dst);
					}
				}
			}

			cerr << "found " << nbrs.size() << " common neighbours" << endl;
			for (int i = 0; i < (int) nbrs.size(); i++) {
				for (int j = i + 1; j < (int) nbrs.size(); j++) {
					if (adjmat[nbrs[i]].count(nbrs[j])) {
						weighted_clique cnew = clique;
						cnew.vertices.push_back(nbrs[i]);
						cnew.vertices.push_back(nbrs[j]);
						cnew.weight += wsum[nbrs[i]] + wsum[nbrs[j]] + adjmat[nbrs[i]][nbrs[j]];
						sort(cnew.vertices.begin(), cnew.vertices.end());
						counter.insert(cnew);
					}
				}
			}
		}
#else
		for (auto clique : cliques) {
			// unlabelling phase
			set<int> seen;
			for (int& vert : clique.vertices) {
				vert = unlabel[vert];
				seen.insert(vert);
			}
			clique.vertices.push_back(u);
			clique.vertices.push_back(v);

			for (auto e : G[u]) {
				if (seen.count(e.dst)) {
					clique.weight += e.wt;
				}
			}
			for (auto e : G[v]) {
				if (seen.count(e.dst)) {
					clique.weight += e.wt;
				}
			}
			clique.weight += w;
			sort(clique.vertices.begin(), clique.vertices.end());
			counter.insert(clique);
		}
#endif

	}
	cerr << "Found " << counter.size() << " " << k << "-cliques." << endl;
	if (counter.size()) cerr << "The maximum weight " << k << "-clique was " << *counter.begin() << endl;

	double tot_time = (clock() - st) / CLOCKS_PER_SEC;
	cerr << "Total Time (s): " << tot_time << endl;
	cerr << "Time per sample (s): " << tot_time / nsamples << endl;

	return counter;
}

}

#endif /* CLIQUE_SAMPLER_H */