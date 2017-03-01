#!/bin/sh

input=$1
bin_dir=$2
map=$3
cpu=$4
OPT=$5

echo "module load compiler/gcc-4.9.1 ; $bin_dir/populations -P $input -M $map -b 1 -k -t $cpu --fstats $OPT" 
module load compiler/gcc-4.9.1 ; $bin_dir/populations -P $input -M $map -b 1 -k -t $cpu --fstats $OPT
