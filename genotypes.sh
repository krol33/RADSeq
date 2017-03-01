#!/bin/sh

input=$1
bin_dir=$2
map=$3
cpu=$4
OPT=$5

echo "module load compiler/gcc-4.9.1 ; $bin_dir/genotypes -b 1 -P $input -t GEN -s -o joinmap $OPT"
module load compiler/gcc-4.9.1 ; $bin_dir/genotypes -b 1 -P $input -t GEN -s -o joinmap $OPT 

