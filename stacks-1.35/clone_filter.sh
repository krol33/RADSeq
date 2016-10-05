#!/bin/sh

#$ -S /bin/sh

#$ -q workq


#récupération des input
f1=$1	# fichier fastq 1
f2=$2	# fichier fastq 1
outdir=$3
stacks_dir=$4	# dossier d'executable stacks


type=`echo $f1 | awk '{if(match($1,".gz") > 0 ){print "gzfastq"}else{print "fastq"}}'`

if [[ "$type" == "gzfastq" ]]
then
gzip -dc $f1 > $outdir/`basename $f1`
gzip -dc $f2 > $outdir/`basename $f2`
f1=$outdir/`basename $f1`
f2=$outdir/`basename $f2`
fi

echo "CMD:		$stacks_dir/clone_filter -1 $f1 -2 $f2 -o $outdir -i fastq -y fastq" 
$stacks_dir/clone_filter -1 $f1 -2 $f2 -o $outdir -i fastq -y fastq

base=`basename $f1 | sed 's/_1.fq/ /' | awk '{print $1}'`

if [[ "$type" != "gzfastq" && ! -h $f1 ]]
then
	gzip $f1 $f2
else
	rm $f1 $f2
fi

mv $outdir/`basename $f1`.1.fq $outdir/${base}_cloneF_1.fq
mv $outdir/`basename $f2`.2.fq $outdir/${base}_cloneF_2.fq

gzip $outdir/${base}_cloneF*
