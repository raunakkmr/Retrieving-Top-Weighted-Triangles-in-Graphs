#!/bin/bash
declare -a arr=(
    email-Enron
    email-Eu
    contact-primary-school
    contact-high-school
    NDC-classes
    tags-math-sx
    congress-bills
)
for i in "${arr[@]}"
do
    julia ../src/ah_mean_sampling_fast.jl $i 500 25 $1
done
