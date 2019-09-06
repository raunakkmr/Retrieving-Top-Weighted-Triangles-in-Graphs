#!/bin/bash

# Runs the compare_parallel_sampling (edge sampling) executable for all
# datasets for 10 runs and writes output to
# ../output/compare_parallel_edge_${k}_${run}/[dataset-name].

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

make -C ../src/cpp clean
make -C ../src/cpp compare_parallel_sampling

mkdir ../output/
for k in "${kvals[@]}"
do
    mkdir ../output/compare_parallel_edge_${k}_${run}
done

for run in {1..10..1}
do
    for k in "${kvals[@]}"
    do
        for dataset in "${datasets[@]}"
        do
            ../src/cpp/compare_parallel_sampling -filename=../data/binaries/${dataset}.binary -binary=true -k=${k} -sampler=edge -start_time=0.4 -end-time=2.0 -increment=0.4 &> ../output/compare_parallel_edge_${k}_${run}/${dataset}
        done
    done
done

for run in {1..10..1}
do
    ../src/cpp/compare_parallel_sampling -filename=../data/binaries/temporal-reddit-reply.binary -binary=true -k=1000 -sampler=edge -start_time=1.5 -end-time=3.5 -increment=0.5 &> ../output/compare_parallel_edge_1000_${run}/temporal-reddit-reply
done

for run in {1..10..1}
do
    ../src/cpp/compare_parallel_sampling -filename=../data/binaries/temporal-reddit-reply.binary -binary=true -k=100000 -sampler=edge -start_time=2.0 -end-time=6.0 -increment=1.0 &> ../output/compare_parallel_edge_100000_${run}/temporal-reddit-reply
done

for run in {1..10..1}
do
    ../src/cpp/compare_parallel_sampling -filename=../data/binaries/spotify.binary -binary=true -k=1000 -sampler=edge -start_time=30 -end-time=30 -increment=30 &> ../output/compare_parallel_edge_1000_${run}/spotify
done

for run in {1..10..1}
do
    ../src/cpp/compare_parallel_sampling -filename=../data/binaries/spotify.binary -binary=true -k=100000 -sampler=edge -start_time=20 -end-time=20 -increment=20 &> ../output/compare_parallel_edge_100000_${run}/spotify
done