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
    julia ../src/geom_mean_exact.jl $i
done
