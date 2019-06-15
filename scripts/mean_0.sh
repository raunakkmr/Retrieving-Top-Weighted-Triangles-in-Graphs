#!/bin/bash
declare -a datasets=(
    coauth-DBLP
    coauth-MAG-Geology
    coauth-MAG-History
    congress-bills
    contact-high-school
    contact-primary-school
    DAWN
    email-Enron
    email-Eu
    NDC-classes
    NDC-substances
    tags-ask-ubuntu
    tags-math-sx
    tags-stack-overflow
    threads-ask-ubuntu
    threads-math-sx
    threads-stack-overflow
)
declare -a samples=(
    1000000
    1000000
    1000000
    7500
    250
    500
    250
    100
    250
    250
    500
    250
    250
    2500
    5000
    5000
    1000000
)
declare -a kvals=(
    25
    50
    250
)
make clean
make
for i in "${!datasets[@]}"
do
    dataset=${datasets[$i]}
    sample=${samples[$i]}
    for p in "${pvals[@]}"
    do
        ../src/cpp/clique_sampler ~/Datasets/$dataset/$dataset- $sample 0 1 $k 0 &> ../../output/mean_0/${dataset}_${k}_${sample}
    done
done