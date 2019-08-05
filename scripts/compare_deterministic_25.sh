#!/bin/bash
declare -a datasets=(
    congress-bills
    tags-stack-overflow
    threads-math-sx
    threads-stack-overflow
    wikipedia
    eth
    aminer
    MAG
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
../src/cpp/compare_deterministic /mnt/disks/additional_ssd_dir/spotify.binary notusingrightnow 1 25 &> ../output/compare_deterministic_25/spotify
