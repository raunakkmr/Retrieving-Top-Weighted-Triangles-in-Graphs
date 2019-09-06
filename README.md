# Top-k Triangles

This code accompanies the paper

* Cite paper here.

## Requirements

* C++ (for running the core algorithms)
* Python 3 (for reproducing plots)
* [gflags](https://github.com/gflags/gflags) (for command line arguments in C++)

## Datasets

Download the datasets from [this webpage](http://www.cs.cornell.edu/~arb/data/index.html).

## Code

### Reproduce results

In order to reproduce the results run the following scripts in the `scripts` directory:
* `convert_data.sh`: Converts the data to binary format as described above.
* `compare_heavy_light.sh`: Runs the static heavy-light, dynamic heavy-light and auto heavy-light algorithms for k = 1,000 and 100,000.
* `compare_parallel_edge.sh`: Runs parallel edge sampling for k = 1,000 and 100,000.
* `compare_parallel_wedge.sh`: Runs parallel wedge sampling for k = 1,000 and 100,000.

Run `python print_table.py` in the `src` directory to print the table of results.