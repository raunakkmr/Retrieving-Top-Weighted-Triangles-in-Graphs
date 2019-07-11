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
declare -a times=(
    3
    3
    1
    5
    0.1
    0.2
    0.2
    3
    0.05
    3
    100
)
declare -a incs=(
    0.3
    0.3
    0.1
    0.5
    0.01
    0.02
    0.02
    0.3
    0.5
    0.3
    10
)
make -C ../src/cpp/compare_sampling clean
make -C ../src/cpp/compare_sampling
mkdir ../output/
mkdir ../output/compare_sampling
for i in "${!datasets[@]}"
do
    dataset=${datasets[$i]}
    time=${times[$i]}
    inc=${incs[$i]}
    ../src/cpp/compare_sampling/clique_sampler ../data/$dataset/$dataset 1 25 1.0 $time $inc &> ../output/compare_sampling/${dataset}
done
../src/cpp/compare_sampling/clique_sampler ../data/temporal-reddit-reply.txt 1 25 1.0 500 50 &> ../output/compare_sampling/temporal-reddit-reply