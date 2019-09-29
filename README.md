# Retrieving Top Weighted Triangles in Graphs

This code accompanies the paper

* **TODO**: Cite paper here.

## Requirements

* C++ (for running the core algorithms)
* Python 3 (for reproducing plots)
* [gflags](https://github.com/gflags/gflags) (for command line arguments in C++)

## Datasets

Download the datasets from [this
webpage](http://www.cs.cornell.edu/~arb/data/index.html).
Our scripts assume that the binary files (described below) are in the
`data` folder, but the scripts can be edited to specify a custom path.

## Code

### Data formats

The datasets that we use in the paper come in one of the following 3 formats:

1. Timestamped sequence of simplices, where a simplex is a set of nodes from a vertex set. A weighted graph is constructed by forming the projected graph of the list of simplices.
2. List of temporal edges, i.e., a list of edges with timestamps. We ignore the timestamps and set the weight of each edge to be the number of times it appeared in the data.
3. List of weighted edges.

Due to the enormous size of these datasets, simply reading in the data file
takes a long time. So we first generated a specialized compressed binary
file. The first line of the file is an integer representing the number of
nodes (n) and the second line of the file is an integer representing the
number of edges (m). The remaining 3m lines describe the edges. Lines 3i + 2
and 3i + 3 are integers representing endpoints of the ith edge and line 3i +
4 is a long long representing the weight of the ith edge. The details of the
compression can be found in the code documentation.
Essentially, we reduce the number of bytes needed in 2 ways. First, we relabel
nodes so that high degree nodes are assigned a smaller label. If this label
can be represented in less than 4 bytes then we read / write only those number
of bytes. Additionally, since a majority of the weights are small we do not
always need 4 bytes to represent them either. This leads to significant benefit,
especially when running multiple experiments. For example, reading in the raw txt
file for Spotify took almost 1.5 hours whereas reading in our compressed binary
file only took around 5 minutes.

### Code structure
The following files in the `src/cpp` directory implement data reading,
conversion, proposed algorithms, and comparisons.

* `graph.h`: Implements generic graph data structures, and graph reading / writing.
* `triangle_sampler.h`: Implements all proposed algorithms.
* `convert_data.cpp`: Converts data to format as described above.
* `compare_deterministic.cpp`: Runs brute force, static heavy-light, dynamic heavy-light and auto heavy-light algorithms on a specified dataset and value of k. 
* `compare_parallel_sampling.cpp`: Runs brute force, parallel edge / wedge sampling on a specified dataset and value of k. It runs the sampling algorithms at intervals of time specified by the user.

## Reproduce results

In order to reproduce the results run the following scripts in the `scripts` directory (you may have to edit the scripts to specify the correct path to the datasets.):
* `convert_data.sh`: Converts the data as described above.
* `compare_deterministic.sh`: Runs the static heavy-light, dynamic heavy-light and auto heavy-light algorithms.
* `compare_parallel_edge.sh`: Runs parallel edge sampling.
* `compare_parallel_wedge.sh`: Runs parallel wedge sampling.

Run `python print_table.py` in the `src` directory to print the table of results.

### Complete Example
Suppose we want to run the deterministic and parallel edge sampling algorithms on the `email-Enron` dataset for `k = 1000`. Suppose the dataset is stored in the `../../data/` directory and we want to store the output in the `../../output/` directory. Then, run the following in the `src/cpp` directory:
* Compile: `make clean && make all`
* Convert simplicial data to binary format and store it as `../../data/binaries/email-Enron.binary`: `./convert_data -filename=../../data/email-Enron/email-Enron -format=simplicial -binary_path=../../data/binaries/email-Enron`
* Run deterministic algorithms and save the output in the `../../output/compare_deterministic_1000` directory: `./compare_deterministic -filename=../../data/binaries/email-Enron.binary -binary=true -k=1000 &> ../../output/compare_deterministic_1000/email-Enron`
* Run edge sampling algorithm from 0.4 to 2 seconds at an interval of 0.4 seconds and save the output in the `../../output/compare_parallel_edge_1000` directory: `./compare_parallel_sampling -filename=../../data/binaries/email-Enron.binary -binary=true -k=1000 -sampler=edge -start_time=0.4 -end_time=2.0 -increment=0.4 &> ../../output/compare_parallel_edge_1000/email-Enron`
