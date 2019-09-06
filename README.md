# Top-k Triangles

This code accompanies the paper

* **TODO**: Cite paper here.

## Requirements

* C++ (for running the core algorithms)
* Python 3 (for reproducing plots)
* [gflags](https://github.com/gflags/gflags) (for command line arguments in C++)

## Datasets

Download the datasets from [this webpage](http://www.cs.cornell.edu/~arb/data/index.html). **TODO**: Place them in a particular folder?

## Code

### Data formats

The datasets that we use in the paper come in one of the following 3 formats:

1. Timestamped sequence of simplices, where a simplex is a set of nodes from a vertex set. A weighted graph is constructed by forming the projected graph of the list of simplices.
2. List of temporal edges, i.e., a list of edges with timestamps. We ignore the timestamps and set the weight of each edge to be the number of times it appeared in the data.
3. List of weighted edges.

Due to the enormous size of these datasets, simply reading in the data file takes a long time. So we first generated a specialized compressed binary file. The first line of the file is an integer representing the number of nodes (n) and the second line of the file is an integer representing the number of edges (m). The remaining 3m lines describe the edges. Lines 3i + 2 and 3i + 3 are integers representing endpoints of the ith edge and line 3i + 4 is a long long representing the weight of the ith edge.
The details of the compression can be found in the code documentation. **TODO**: Do we want to talk about details of compression here?

This leads to significant benefit in running multiple experiments. For example, reading in the raw txt file for Spotify took almost 1.5 hours whereas reading in our compressed binary file only took around 5 minutes.

### Code structure
The following files in the `src/cpp` directory implement data reading, conversion, proposed algorithms, and comparisons.

* `graph.h`: Implements generic graph data structures, and graph reading / writing.
* `triangle_sampler.h`: Implements all proposed algorithms.
* `convert_data.cpp`: Converts data to format as described above.
    * Usage: **TODO**
* `compare_deterministic.cpp`: Runs brute force, static heavy-light, dynamic heavy-light and auto heavy-light algorithms on a specified dataset and value of k. 
    * Usage: **TODO**
* `compare_parallel_sampling.cpp`: Runs brute force, parallel edge / wedge sampling on a specified dataset and value of k. It runs the sampling algorithms at intervals of time specified by the user.
    * Usage: **TODO**

## Reproduce results

In order to reproduce the results run the following scripts in the `scripts` directory:
* `convert_data.sh`: Converts the data as described above.
    * Usage: **TODO**
* `compare_deterministic.sh`: Runs the static heavy-light, dynamic heavy-light and auto heavy-light algorithms for k = 1,000 and 100,000.
    * Usage: **TODO**
* `compare_parallel_edge.sh`: Runs parallel edge sampling for k = 1,000 and 100,000. 
    * Usage: **TODO**
* `compare_parallel_wedge.sh`: Runs parallel wedge sampling for k = 1,000 and 100,000. 
    * Usage: **TODO**

Run `python print_table.py` in the `src` directory to print the table of results.