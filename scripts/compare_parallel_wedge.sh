#!/bin/bash

# Runs the compare_parallel_sampling (wedge sampling) executable for all
# datasets for 10 runs and writes output to
# ../output/compare_parallel_edge_${k}_${run}/[dataset-name].

declare -a datasets=(
    tags-stack-overflow
    MAG
    temporal-reddit-reply
)

declare -a starts=(
  6
  25
  300
)

declare -a ends=(
  12
  37
  500
)

declare -a incs=(
  3
  6
  200
)

declare -a kvals=(
    1000
)

make -C ../src/cpp clean
make -C ../src/cpp compare_parallel_sampling

mkdir ../output/
for k in "${kvals[@]}"
do
    for run in {1..10..1}
    do
        mkdir ../output/compare_parallel_wedge_${k}_${run}
    done
done

for run in {1..10..1}
do
  for k in "${kvals[@]}"
  do
    for i in "!${datasets[@]}"
    do
      dataset=${datasets[$i]}
      st=${starts[$i]}
      en=${ends[$i]}
      inc=${incs[$i]}
      ../src/cpp/compare_parallel_sampling -filename=../data/binaries/${dataset}.binary -binary=true k=1000 -sampler=wedge -start_time=$st -end_time=$en -increment=$inc &> ../output/compare_parallel_wedge_${k}_${run}/${dataset}
    done
  done
done
