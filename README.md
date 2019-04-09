# Top-k Triangles

## Data
Follow the instructions at the end of the Data section in [this repository](https://github.com/arbenson/ScHoLP-Tutorial#data)
to download the data. Place it in a folder named *./data*.

## Code

### src
This folder contains julia scripts and jupyter notebooks.
* *analyze_dataset.ipynb* loads a dataset and plots the distribution of edge weights of the 
graph, shows the top 10 vertices with "heavy" neighboring edges, etc.
* *compare_results_\*.ipynb* compares the quality of the solutions obtained from the sampling 
algorithm and the exact solution. It displays the distribution of the weights of the top 
triangles (in the exact solution and the estimated solution), the distribution of 
percentile of the estimated triangles, etc.
* *utils.jl* contains functions to read the dataset from file, construct the projected graph,
get the adjacency list, etc.
* *geom_mean_\*.jl* implement orders triangles by geometric mean of edge weights.
*ah_mean_\*.jl* implement orders triangles by arithmetic or harmonic mean of edge weights.
* *\*exact.jl* implement an exact algorithm. They enumerate triangles and then extract the top 
k triangles from those.
* *geom_mean_sampling.jl* implements the sampling algorithm with rejection sampling. 
*geom_mean_sampling_wr.jl* implements the sampling algorithm without rejection sampling but the slow version.
* *ah_mean_sampling_fast.sh* implements the sampling algorithm without rejection sampling but
the fast version. *ah_mean_sampling_wr.sh* implements the sampling algorithm without rejection 
sampling but the slow vertion. 
* Usage: Suppose we choose the *email-Enron* dataset and set `kprime = 500` and 
`k = 25`.
    * Run the exact algorithms as `julia geom_mean_exact.sh email-Enron` or 
     `julia ah_mean_exact.sh email-Enron a` (replace `a` with `h` for harmonic mean).
    * Run the sampling algorithms as `julia geom_mean_sampling*.jl email-Enron 500 25` or
    `julia ah_mean_sampling*.jl email-Enron 500 25 a` (replace `a` with `h` for harmonic mean).
* Output is placed in a txt file in a folder named *./output*.

### scripts
This folder contains scripts to run a particular algorithm on multiple datasets.
Each of them uses `kprime = 500` and `k = 25`.
* For ordering triangles using the arithmetic mean or harmonic mean, use "a" or "h" 
respectively (without quotes) as a command line argument to the script. 
* *ah_mean_exact.sh* runs the exact algorithm. *ah_mean_sampling_wr.sh* runs the sampling 
algorithm without rejection sampling but the slow version. *ah_mean_sampling.sh* runs the 
sampling algorithm without rejection sampling but the fast version.
* For ordering triangles using the geometric mean, run *geom_mean_\*.sh*. *geom_mean_exact.sh*
runs the exact algorithm. *geom_mean_sampling_wr.sh* runs the sampling algorithm without rejection but the slow version. *geom_mean_sampling.sh* runs the sampling algorithm with rejection sampling.