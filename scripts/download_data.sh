#!/bin/bash

wget https://github.com/arbenson/ScHoLP-Data/archive/1.0.tar.gz
tar -xzvf 1.0.tar.gz
gunzip ScHoLP-Data-1.0/*/*.gz
mv ScHoLP-Data-1.0/* ../data/
rm -rf ScHoLP-Data-1.0
mv 1.0.tar.gz ../data/