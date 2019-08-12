#include <bits/stdc++.h>

#include "graph.h"

using namespace std;
using namespace wsdm_2019_graph;

int main(int argc, char *argv[]) {
    ios_base::sync_with_stdio(0);
    cin.tie(0);

    auto GS = read_graph(argv[1], true);
    string out_path = argv[2];
    const auto &edges = GS.edges;
    map<long long, int> count;
    for (const auto &e : edges) {
        count[e.wt]++;
    }
    ofstream ofile(out_path, ios::out);
    for (const auto &kv : count) {
        ofile << kv.first << " " << kv.second << endl;
    }
    ofile.close();
}