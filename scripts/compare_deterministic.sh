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

declare -a kvals=(
#     25
    1000
#     40000
    100000
)

make -C ../src/cpp clean
make -C ../src/cpp compare_deterministic

mkdir ../output/
for k in "${kvals[@]}"
do
    mkdir ../output/compare_deterministic_${k}
done

for i in "${!datasets[@]}"
do
    dataset=${datasets[$i]}
    for k in "${kvals[@]}"
    do
#       ../src/cpp/compare_deterministic ../data/binaries/${dataset}.binary ../output/brute_force_enumerator/${dataset}-triangles.binary 1 $k &> ../output/compare_deterministic_${k}/${dataset}
      ../src/cpp/compare_deterministic ../data/binaries/${dataset}.binary notusingrightnow 1 $k &> ../output/compare_deterministic_${k}/${dataset}
    done
done

for k in "${kvals[@]}"
do
    ../src/cpp/compare_deterministic ../data/binaries/temporal-reddit-reply.binary notusingrightnow 1 $k &> ../output/compare_deterministic_${k}/temporal-reddit-reply
done

for k in "${kvals[@]}"
do
    ../src/cpp/compare_deterministic /mnt/disks/additional_ssd_dir/spotify.binary notusingrightnow 0 $k &> ../output/compare_deterministic_${k}/spotify
done
