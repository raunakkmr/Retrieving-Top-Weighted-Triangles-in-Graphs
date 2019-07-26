#!/bin/bash
declare -a datasets=(
    coauth-DBLP
    coauth-MAG-Geology
    coauth-MAG-History
    congress-bills
    tags-stack-overflow
    threads-math-sx
    threads-stack-overflow
)
declare -a times=(
    3
    3
    1
    5
    3
    3
    100
)
declare -a incs=(
    0.3
    0.3
    0.1
    0.5
    0.3
    0.3
    10
)
make -C ../src/cpp clean
make -C ../src/cpp compare_sampling
mkdir ../output/
mkdir ../output/compare_sampling_25
for i in "${!datasets[@]}"
do
    dataset=${datasets[$i]}
    time=${times[$i]}
    inc=${incs[$i]}
    ../src/cpp/compare_sampling ../data/binaries/${dataset}.binary ../output/brute_force_enumerator/${dataset}-triangles.binary 1 25 $time $inc &> ../output/compare_sampling_25/${dataset}
done
# ../src/cpp/compare_sampling/clique_sampler ../data/temporal-reddit-reply.txt 1 25 1.0 500 50 &> ../output/compare_sampling/temporal-reddit-reply
