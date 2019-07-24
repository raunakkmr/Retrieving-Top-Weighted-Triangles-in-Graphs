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
)
make -C ../src/cpp/ clean
make -C ../src/cpp/ convert_data
for i in "${!datasets[@]}"
do
    dataset=${datasets[$i]}
    ../src/cpp/convert_data ../data/${dataset}/${dataset} ../data/binaries/${dataset}
done
# ../src/cpp/convert_data ../data/temporal-reddit-reply.txt ../data/temporal-reddit-reply