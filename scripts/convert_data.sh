#!/bin/bash
declare -a datasets=(
    congress-bills
    tags-stack-overflow
    threads-math-sx
    threads-stack-overflow
    wikipedia
    eth
    aminer
    temporal-reddit-reply
)
make -C ../src/cpp/ clean
make -C ../src/cpp/ convert_data
for i in "${!datasets[@]}"
do
    dataset=${datasets[$i]}
    ../src/cpp/convert_data ../data/${dataset}/${dataset} ../data/binaries/${dataset}
done
../src/cpp/convert_data /mnt/disks/additional_ssd_dir/MAG-weighted.txt ../data/binaries/MAG
../src/cpp/convert_data /mnt/disks/additional_ssd_dir/spotify-weighted.txt /mnt/disks/additional_ssd_dir/spotify