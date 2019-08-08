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
make -C ../src/cpp clean
make -C ../src/cpp compare_deterministic
mkdir ../output/
mkdir ../output/compare_deterministic_25
for i in "${!datasets[@]}"
do
    dataset=${datasets[$i]}
    ../src/cpp/compare_deterministic ../data/binaries/${dataset}.binary ../output/brute_force_enumerator/${dataset}-triangles.binary 1 25 &> ../output/compare_deterministic_25/${dataset}
done
# ../src/cpp/compare_deterministic/clique_sampler ../data/temporal-reddit-reply.txt 1 25 1.0 0 0 &> ../output/compare_deterministic/temporal-reddit-reply
