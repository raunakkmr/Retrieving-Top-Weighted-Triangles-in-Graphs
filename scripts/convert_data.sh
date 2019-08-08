#!/bin/bash
declare -a datasets=(
    congress-bills
<<<<<<< HEAD
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
=======
    tags-stack-overflow
>>>>>>> c2225f06343169f420bd2b5ff7da22566691aee7
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
