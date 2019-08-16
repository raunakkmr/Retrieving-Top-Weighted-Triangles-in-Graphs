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
    1000
    100000
)

make -C ../src/cpp clean
make -C ../src/cpp compare_parallel_sampling

mkdir ../output/
for k in "${kvals[@]}"
do
    mkdir ../output/compare_parallel_edge_${k}
done

for i in "${!datasets[@]}"
do
    dataset=${datasets[$i]}
    # for k in "${kvals[@]}"
    # do
    #   ../src/cpp/compare_parallel_sampling ../data/binaries/${dataset}.binary notusingrightnow 1 $k 0.4 1.6 0.4 E &> ../output/compare_parallel_edge_${k}/${dataset}
    # done
    ../src/cpp/compare_parallel_sampling ../data/binaries/${dataset}.binary notusingrightnow 1 1000 0.4 2.0 0.4 E &> ../output/compare_parallel_edge_1000/${dataset}
    ../src/cpp/compare_parallel_sampling ../data/binaries/${dataset}.binary notusingrightnow 1 100000 0.5 4.5 1.0 E &> ../output/compare_parallel_edge_100000/${dataset}
done

# for k in "${kvals[@]}"
# do
#     ../src/cpp/compare_parallel_sampling ../data/binaries/temporal-reddit-reply.binary notusingrightnow 1 $k 1.5 3.0 0.5 E &> ../output/compare_parallel_edge_${k}/temporal-reddit-reply
# done
../src/cpp/compare_parallel_sampling ../data/binaries/temporal-reddit-reply.binary notusingrightnow 1 1000 1.5 3.5 0.5 E &> ../output/compare_parallel_edge_1000/temporal-reddit-reply
../src/cpp/compare_parallel_sampling ../data/binaries/temporal-reddit-reply.binary notusingrightnow 1 100000 2.0 6.0 1.0 E &> ../output/compare_parallel_edge_100000/temporal-reddit-reply

for k in "${kvals[@]}"
do
    ../src/cpp/compare_parallel_sampling /mnt/disks/additional_ssd_dir/spotify.binary notusingrightnow 0 $k 0.4 2.0 0.4 E &> ../output/compare_parallel_edge_${k}/spotify
done
