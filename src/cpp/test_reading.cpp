// This file compares time taken to read in binary vs non-binary format.

#include <bits/stdc++.h>

#include "gflags/gflags.h"

#include "graph.h"

using namespace std;
using namespace wsdm_2019_graph;

DEFINE_string(filename, "", "Path to graph file in non-binary format.");
DEFINE_string(format, "", "If binary is false this is required and indicates format of graph file. One of weighted, temporal, simplicial.");
DEFINE_string(binaryname, "", "Path to graph file in binary format.");

int main(int argc, char *argv[]) {
  string usage("Test reading time of graph file.\n"
      "Sample usage:\n"
      "\t./test_reading -filename=[fill_this_in] "
      "-format=[fill_this_in] -binaryname=[fill_this_in]\n"
      "Additionally, these flags can be loaded from a single file "
      "with the option -flagfile=[filename].");

  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_filename.empty()) {
    std::cerr << "No file specified! Type ./test_reading --help for a description of the program parameters." << std::endl;
    return 0;
  }

  ios::sync_with_stdio(0);
  cin.tie(0);
  srand(0);

  struct timespec start1, finish1;
  double tot_time1;
  clock_gettime(CLOCK_MONOTONIC, &start1);

  auto GS = read_graph(FLAGS_filename, false, FLAGS_format);

  clock_gettime(CLOCK_MONOTONIC, &finish1);
  tot_time1 = (finish1.tv_sec - start1.tv_sec);
  tot_time1 += (finish1.tv_nsec - start1.tv_nsec) / 1000000000.0;
  cerr << "Non binary read for " << FLAGS_filename << " took " << tot_time1 << " seconds." << endl;

  cerr << endl << endl;

  struct timespec start2, finish2;
  double tot_time2;
  clock_gettime(CLOCK_MONOTONIC, &start2);

  GS = read_graph(FLAGS_binaryname, true);

  clock_gettime(CLOCK_MONOTONIC, &finish2);
  tot_time2 = (finish2.tv_sec - start2.tv_sec);
  tot_time2 += (finish2.tv_nsec - start2.tv_nsec) / 1000000000.0;
  cerr << "Binary read for " << FLAGS_binaryname << " took " << tot_time2 << " seconds." << endl;

  return 0;
}
