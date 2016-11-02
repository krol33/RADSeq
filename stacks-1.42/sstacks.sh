#!/bin/sh

input=$1
cstack_dir=$2
out=$3
bin_dir=$4
opt=$5

echo "module load compiler/gcc-4.9.1 ; $bin_dir/sstacks -b 1 $opt -c $cstack_dir/batch_1 -s $input -o $out"
module load compiler/gcc-4.9.1 ; $bin_dir/sstacks -b 1 $opt -c $cstack_dir/batch_1 -s $input -o $out
