#include <bits/stdc++.h>

#include "graph.h"

using namespace std;
using namespace wsdm_2019_graph;

int main(int argc, char *argv[]) {
    ios::sync_with_stdio(0);
    cin.tie(0);
    srand(0);

    // struct timespec start1, finish1;
    // double tot_time1;
    // clock_gettime(CLOCK_MONOTONIC, &start1);

    // auto G = read_graph(argv[1]);

    // clock_gettime(CLOCK_MONOTONIC, &finish1);

    // tot_time1 = (finish1.tv_sec - start1.tv_sec);
    // tot_time1 += (finish1.tv_nsec - start1.tv_nsec) / 1000000000.0;
    // cerr << "Non binary read for " << argv[1] << " took " << tot_time1 << " seconds." << endl;

    cerr << endl << endl;

    struct timespec start2, finish2;
    double tot_time2;
    clock_gettime(CLOCK_MONOTONIC, &start2);

    // G = read_graph(argv[2], true);
    auto G = read_graph(argv[2], true);

    clock_gettime(CLOCK_MONOTONIC, &finish2);

    tot_time2 = (finish2.tv_sec - start2.tv_sec);
    tot_time2 += (finish2.tv_nsec - start2.tv_nsec) / 1000000000.0;
    cerr << "binary read for " << argv[2] << " took " << tot_time2 << " seconds." << endl;

    return 0;
}
