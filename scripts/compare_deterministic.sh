#!/bin/bash

# Runs the compare_deterministic executable for all datasets for 10 runs and
# writes output to ../output/compare_deterministic_${k}_${run}/[dataset-name].

declare -a datasets=(
    tags-stack-overflow
    threads-math-sx
    threads-stack-overflow
    wikipedia
    eth
    aminer
    MAG
)

declare -a kvals=(
    1000
    100000
)

declare runs=10

make -C ../src/cpp clean
make -C ../src/cpp compare_deterministic

mkdir ../output/
for k in "${kvals[@]}"
do
    for run in {1..10..1}
    do
        mkdir ../output/compare_deterministic_${k}_${run}
    done
done

for run in {1..10..1}
do
    for k in "${kvals[@]}"
    do
        for dataset in "${datasets[@]}"
        do
            ../src/cpp/compare_deterministic -filename=../data/binaries/${dataset}.binary -binary=true -k=${k} &> ../output/compare_deterministic_${k}_${run}/${dataset}
        done
    done
done

for run in {1..10..1}
do
    for k in "${kvals[@]}"
    do
        ../src/cpp/compare_deterministic -filename=../data/binaries/temporal-reddit-reply.binary -binary=true -k=${k} &> ../output/compare_deterministic_${k}_${run}/temporal-reddit-reply

        ../src/cpp/compare_deterministic -filename=/mnt/disks/additional_ssd_dir/spotify.binary -binary=true -k=${k} &> ../output/compare_deterministic_${k}_${run}/spotify
    done
done