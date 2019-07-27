#ifndef CLIQUE_SAMPLER_H
#define CLIQUE_SAMPLER_H

#include <bits/stdc++.h>

#include "graph.h"
#include "clique_enumerator.h"

using namespace std;

namespace wsdm_2019_graph {

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
				int v = e.dst;
				long long w = e.wt;
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

	map<int, map<int, long long>> adjmat;
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
					long long w = adjmat[u][v];
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
	vector<map<int, long long>> adjmat(G.size());
	int nedges = 0;
	for (int u = 0; u < (int) G.size(); u++) {
		if (removed[u]) continue;
		for (auto e : G[u]) {
			int v = e.dst;
			long long w = e.wt;
			if (u > v) continue;
			if (removed[v]) continue;
			edge_distribution[e.wt].push_back({u, v, w});
			adjmat[u][v] = adjmat[v][u] = w;
			nedges++;
		}
	}

	vector<long long> cumulative_weights;
	vector<long long> index_to_weight(nedges);
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

	set<weighted_clique> counter;
	set<pair<int, int>> history;
	for (int samp = 0; samp < nsamples; samp++) {
		auto sample = sample_edge();
		int u = sample.src, v = sample.dst;
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
			map<int, int> seen;
			long long wsum;
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
