#!/bin/sh

read1=$1
id=$2
out=$3
bin_dir=$4
cov=$5
mismatch=$6
mismatch2=$7
maxLocus=$8
opt=$9

# search stack for one individual
echo "$bin_dir/ustacks -t gzfastq -f $read1 -i $id -o $out -d -r -m $cov -M $mismatch -N $mismatch2 --max_locus_stacks $maxLocus $opt"
$bin_dir/ustacks -t gzfastq -f $read1 -i $id -o $out -d -r -m $cov -M $mismatch -N $mismatch2 --max_locus_stacks $maxLocus $opt
