#!/bin/bash
declare -a datasets=(
#     congress-bills
#     tags-stack-overflow
#     threads-math-sx
#     threads-stack-overflow
    wikipedia
    eth
    aminer
    MAG
)

declare -a times=(
#     5
#     3
#     3
#     100
    100
    150
    200
    200
)

declare -a incs=(
#     0.5
#     0.3
#     0.3
#     10
    10
    15
    20
    20
)

declare -a kvals=(
    25
    1000
    40000
)
# make -C ../src/cpp clean
# make -C ../src/cpp compare_sampling

mkdir ../output/
for k in "${kvals[@]}"
do
    mkdir ../output/compare_sampling_${k}
done

for i in "${!datasets[@]}"
do
    dataset=${datasets[$i]}
    tim=${times[$i]}
    inc=${incs[$i]}
    echo $dataset
    echo $tim
    echo $inc
    for k in "${kvals[@]}"
    do
      ../src/cpp/compare_sampling ../data/binaries/${dataset}.binary notusingrightnow 1 $k $tim $inc &> ../output/compare_sampling_${k}/${dataset}
    done
done

for k in "${kvals[@]}"
do
    ../src/cpp/compare_sampling ../data/binaries/temporal-reddit-reply.binary notusingrightnow 1 $k 200 20 &> ../output/compare_sampling_${k}/temporal-reddit-reply
done

for k in "${kvals[@]}"
do
    ../src/cpp/compare_sampling /mnt/disks/additional_ssd_dir/spotify.binary notusingrightnow 0 $k 500 50 &> ../output/compare_sampling_${k}/spotify
done
