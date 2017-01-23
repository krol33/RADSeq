#!/bin/sh

#$ -S /bin/sh

#$ -q workq


#récupération des input
script_dir=$1
f1=$2	# fichier fastq 1
f2=$3	# fichier fastq 2
prefix=$4	# dossier de sortie
barcode=$5	# fichier de barcode
mismatch=$6
trim=$7

if [[ $trim = 1 ]]
then
echo  "CMD:		perl $script_dir/splitbc.pl $f1 $f2 --mismatches $mismatch --bcfile $barcode --bol --trim --prefix-r1 ${prefix}%_1.fq --prefix-r2 ${prefix}%_2.fq --no_adapt"
perl $script_dir/splitbc.pl $f1 $f2 --mismatches $mismatch --bcfile $barcode --bol --trim --prefix-r1 ${prefix}%_1.fq --prefix-r2 ${prefix}%_2.fq --no_adapt
else
# attention le barcode sur la lecture 2 doit être de la même taille que sur la lecture 1. La recherche se fait uniquement sur la lecture 1
echo  "CMD:		perl $script_dir/splitbc.pl $f1 $f2 --mismatches $mismatch --bcfile $barcode --bol --trim2 --prefix-r1 ${prefix}%_1.fq --prefix-r2 ${prefix}%_2.fq --no_adapt"
perl $script_dir/splitbc.pl $f1 $f2 --mismatches $mismatch --bcfile $barcode --bol --trim2 --prefix-r1 ${prefix}%_1.fq --prefix-r2 ${prefix}%_2.fq --no_adapt
fi
