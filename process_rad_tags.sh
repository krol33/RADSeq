#!/bin/sh

#$ -S /bin/sh

#$ -q workq


#récupération des input
f1=$1	# fichier fastq 1
outdir=$2
enz=$3
stacks_dir=$4	# dossier d'executable stacks
max_len=$5

type=`echo $f1 | awk '{if(match($1,".gz") > 0 ){print "gzfastq"}else{print "fastq"}}'`

if [[ "$max_len" -gt 0 ]]
then
echo "CMD:		$stacks_dir/process_radtags -f $f1 -o $outdir -e $enz -i $type -y fastq --filter_illumina -r -q -E phred33 -D -s 20 -t $max_len --retain_header" 
$stacks_dir/process_radtags -f $f1 -o $outdir -e $enz -i $type -y fastq --filter_illumina -r -q -E phred33 -D -s 20 -t $max_len --retain_header
else
echo "CMD:		$stacks_dir/process_radtags -f $f1 -o $outdir -e $enz -i $type -y fastq --filter_illumina -r -q -E phred33 -D -s 20 --retain_header" 
$stacks_dir/process_radtags -f $f1 -o $outdir -e $enz -i $type -y fastq --filter_illumina -r -q -E phred33 -D -s 20 --retain_header
fi

if [[ "$type" != "gzfastq" && ! -h $f1 ]]
then
gzip $f1
fi
