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
    wikipedia
    # eth
    # aminer
    # MAG
)
make -C ../src/cpp clean
make -C ../src/cpp compare_deterministic
mkdir -p ../output/
mkdir -p ../output/compare_deterministic_1000
for i in "${!datasets[@]}"
do
    dataset=${datasets[$i]}
    ../src/cpp/compare_deterministic ../data/binaries/${dataset}.binary ../output/brute_force_enumerator/${dataset}-triangles.binary 1 1000 &> ../output/compare_deterministic_1000/${dataset}
done
# ../src/cpp/compare_deterministic /mnt/disks/additional_ssd_dir/spotify.binary notusingrightnow 1 1000 &> ../output/compare_deterministic_1000/spotify
