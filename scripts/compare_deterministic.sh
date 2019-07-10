#!/bin/bash
declare -a datasets=(
    coauth-DBLP
    coauth-MAG-Geology
    coauth-MAG-History
    congress-bills
    DAWN
    email-Eu
    tags-math-sx
    tags-stack-overflow
    threads-ask-ubuntu
    threads-math-sx
    threads-stack-overflow
)
make -C ../src/cpp/compare_deterministic clean
make -C ../src/cpp/compare_deterministic
mkdir ../output/compare_deterministic
for i in "${!datasets[@]}"
do
    dataset=${datasets[$i]}
    ../src/cpp/compare_deterministic/clique_sampler ../data/$dataset/$dataset 1 25 1.0 0 0 &> ../output/compare_deterministic/${dataset}
done