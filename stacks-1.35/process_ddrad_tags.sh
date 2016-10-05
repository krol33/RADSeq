#!/bin/sh

#$ -S /bin/sh

#$ -q workq


#récupération des input
f1=$1	# fichier fastq 1
f2=$2
outdir=$3
enz1=$4
enz2=$5
stacks_dir=$6	# dossier d'executable stacks
max_len=$7

type=`echo $f1 | awk '{if(match($1,".gz") > 0 ){print "gzfastq"}else{print "fastq"}}'`

if [[ "$max_len" -gt 0 ]]
then
echo "CMD:		$stacks_dir/process_radtags -1 $f1 -2 $f2 -o $outdir --renz_1 $enz1  --renz_2 $enz2 -i $type -y fastq --filter_illumina -r -q -E phred33 -D -s 20 -t $max_len --retain_header" 
$stacks_dir/process_radtags -1 $f1 -2 $f2 -o $outdir --renz_1 $enz1 --renz_2 $enz2 -i $type -y fastq --filter_illumina -r -q -E phred33 -D -s 20 -t $max_len --retain_header
else
echo "CMD:		$stacks_dir/process_radtags -1 $f1 -2 $f2 -o $outdir --renz_1 $enz1  --renz_2 $enz2 -i $type -y fastq --filter_illumina -r -q -E phred33 -D -s 20 --retain_header" 
$stacks_dir/process_radtags -1 $f1 -2 $f2 -o $outdir --renz_1 $enz1  --renz_2 $enz2 -i $type -y fastq --filter_illumina -r -q -E phred33 -D -s 20 --retain_header
fi

if [[ "$type" != "gzfastq" ]]
then
if [[ ! -h $f1 ]]; then gzip $f1 ;fi
if [[ ! -h $f2 ]]; then gzip $f2 ;fi
fi

indiv=`basename $f1 | sed 's/_1.fq/ /' | awk '{print $1}'`
gzip $outdir/${indiv}_*
