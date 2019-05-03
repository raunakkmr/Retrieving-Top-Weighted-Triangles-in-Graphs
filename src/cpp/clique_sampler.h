#ifndef CLIQUE_SAMPLER_H
#define CLIQUE_SAMPLER_H

#include <bits/stdc++.h>
using namespace std;

namespace wsdm_2019_graph {

// TODO: implement sample without replacement

set<weighted_triangle> brute_force_sampler(Graph& G) {
	cerr << "=============================================" << endl;
	cerr << "Running brute force detection for triangles" << endl;
	cerr << "=============================================" << endl;
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
	cerr << "Found " << counter.size() << " triangles. The maximum weight triangle was " << *counter.begin() << endl;

	double tot_time = (clock() - st) / CLOCKS_PER_SEC;
	cerr << "Total Time (s): " << tot_time << endl;

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
		if (history.count(make_pair(u, v))) {
			cerr << "RESAMPLED!!" << endl;
			history.insert(make_pair(u, v));
		}
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

set<weighted_triangle> path_sampler(Graph& G) {
	// TODO: implement raunak's path sampling algorithm
	return set<weighted_triangle>();
}

void compare_statistics(set<weighted_triangle>& all_triangles, set<weighted_triangle>& sampled_triangles) {
	cerr << "=============================================" << endl;
	cerr << "Comparing sampling statistics" << endl;
	cerr << "=============================================" << endl;

	int num_found = 0;
	int curr_tri = 0;
	bool first_break = 0;
	for (auto T : all_triangles) {
		if (sampled_triangles.count(T)) {
			num_found++;
		}
		curr_tri++;

		if (num_found != curr_tri && !first_break) {
			first_break = true;
			cerr << "Found top " << 100.0 * num_found / all_triangles.size() << " percent of weighted triangles." << endl;
		}

		if (curr_tri == int(0.1 * all_triangles.size())) {
			cerr << "Found " << 100.0 * num_found / curr_tri << " percent of weighted triangles top 10%." << endl;
		}

		if (curr_tri == int(0.25 * all_triangles.size())) {
			cerr << "Found " << 100.0 * num_found / curr_tri << " percent of weighted triangles top 25%." << endl;
		}

		if (curr_tri == int(0.50 * all_triangles.size())) {
			cerr << "Found " << 100.0 * num_found / curr_tri << " percent of weighted triangles top 50%." << endl;
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
		if (removed[i]) continue;

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
	}
	return retval;
}

set<weighted_clique> clique_sampler(Graph& G, int k, int nsamples) {
	cerr << "=============================================" << endl;
	cerr << "Running edge sampling for k-cliques" << endl;
	cerr << "=============================================" << endl;
	double st = clock();

	// build distribution over edges
	map<int, vector<full_edge>> edge_distribution;
	map<int, map<int, int>> adjmat;
	for (int u = 0; u < (int) G.size(); u++) {
		for (auto e : G[u]) {
			int v = e.dst, w = e.wt;
			if (u > v) continue;
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
			history.insert(make_pair(u, v));
			cerr << "RESAMPLED!!!" << endl;
		}
		map<int, int> vert_to_wt_u, vert_to_wt_v;
		for (auto e : G[u]) {
			vert_to_wt_u[e.dst] = e.wt;
		}

		vector<int> common_nbrs;
		for (auto e : G[v]) {
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
		set<weighted_clique> cliques = enumerate_cliques(subgraph, k-2);

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