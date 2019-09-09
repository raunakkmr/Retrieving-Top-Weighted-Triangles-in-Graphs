#!/bin/bash
declare -a simplicial_datasets=(
    tags-stack-overflow
    threads-math-sx
    threads-stack-overflow
)

declare -a weighted_datasets=(
    wikipedia
    eth
    aminer
)

make -C ../src/cpp/ clean
make -C ../src/cpp/ convert_data

for i in "${!simplicial_datasets[@]}"
do
    dataset=${datasets[$i]}
    ../src/cpp/convert_data -filename=../data/${dataset}/${dataset} -format=simplicial -binary_path../data/binaries/${dataset}
done

for i in "${!weighted_datasets[@]}"
do
    dataset=${datasets[$i]}
    ../src/cpp/convert_data -filename=../data/${dataset}-weighted.txt -format=weighted -binary_path../data/binaries/${dataset}
done

../src/cpp/convert_data -filename=../data/temporal-reddit-reply.txt -format=temporal -binary_path=../data/binaries/temporal-reddit-reply

../src/cpp/convert_data -filename=/mnt/disks/additional_ssd_dir/MAG-weighted.txt -format=weighted -binary_path../data/binaries/MAG

../src/cpp/convert_data -filename=/mnt/disks/additional_ssd_dir/spotify-weighted.txt -format=weighted -binary_path=/mnt/disks/additional_ssd_dir/spotify
