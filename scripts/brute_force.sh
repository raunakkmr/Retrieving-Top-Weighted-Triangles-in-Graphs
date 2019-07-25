#!/bin/bash
declare -a datasets=(
#     coauth-DBLP
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
make -C ../src/cpp/ brute_force_enumerator
mkdir ../output/
mkdir ../output/brute_force_enumerator
for i in "${!datasets[@]}"
do
    dataset=${datasets[$i]}
    ../src/cpp/brute_force_enumerator ../data/binaries/${dataset}.binary ../output/brute_force_enumerator/${dataset} B &> ../output/brute_force_enumerator/${dataset}-timeinfo.txt
done
# ../src/cpp/brute_force_enumerator #
# ../data/binaries/temporal-reddit-reply.binary ../output/brute_force_enumerator/temporal-reddit-reply B &> ../output/brute_force_enumerator/temporal-reddit-reply-timeinfo.txt
