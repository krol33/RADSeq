#!/bin/sh

#$ -o [Project_path_dir]/SGE_out/stacks_denovo_map.out
#$ -e [Project_path_dir]/SGE_out/stacks_denovo_map.err
#$ -N [Proj_Name]_stacks_denovo_map

#$ -S /bin/sh

#$ -q unlimitq

#$ -m bea

date
date >&2
Proj_Name=[Proj_Name]

# VARIABLES GENERALES
stacks_dir=/usr/local/bioinfo/src/Stacks/stacks-1.35/bin
RAD_DIR=[Project_path_dir]
SCRIPT_DIR=[Script_path_dir]

# input
INDIV_FILE=$RAD_DIR/indiv_barcode.txt
POP_FILE=$RAD_DIR/population.map
IN_DIR=$RAD_DIR/preprocessing/stacks_input
OUT_DIR=$RAD_DIR/denovo_map_populations
SGE=$RAD_DIR/SGE_out/denovo_map_populations
STAT_DIR=$RAD_DIR/stat/denovo_map_populations

# OPTIONS
# ustacks
# attention 
# si ustacks alors min coverage per stacks ~= allele
# si ustacks alors min coverage per location ~= locus!
MIN_DEPTH=3
MAX_PRIM_DIST=2
MAX_SEC_DIST=4
# --max_locus_stacks option pour individus diploïde. (Pour les individus happloïdes doublés cette option sera à 2)
MAC_LOCUS_STACKS=3
# autres options telles que :  --model_type --alpha --bound_low --bound_high --bc_err_freq. Ecrire la chaine de caractère entre "quotte"
STACKS_OPT=""

#cstacks
MISMATCH=1

################### checks ##################################
mkdir -p $OUT_DIR $STAT_DIR $SGE

dos2unix $INDIV_FILE
check=`cut -f 2 $INDIV_FILE | sort | uniq -c | awk '$1>1'`
if [[ $check == "" ]]
then echo "all sample names are uniq"
else
	echo "You have sample name(s) that are non uniq in your $INDIV_FILE file:"
	echo $check | awk '{for (i =1;i<= NF; i+=2){name=i+1; print "\t"$name" is present "$i" times"}}'
	exit 1
fi

check=`cut -f 2 $INDIV_FILE | grep "-" | wc -l`
if [[ $check -gt 0 ]]
then
	echo "You use \"-\" in your sample names, please search and replace those with \"_\" in you sample names (only)"
fi
################### Command line ############################

sample=`awk -v I=$IN_DIR '{s=s"-s "I"/"$2".fq.gz "}END{print s}' $INDIV_FILE`

# -t = -d -r dans ustacks
cmd_line="perl /usr/local/bioinfo/src/Stacks/stacks-1.35/scripts/denovo_map.pl $sample -o $OUT_DIR -t -m $MIN_DEPTH -M $MAX_PRIM_DIST -N $MAX_SEC_DIST -n 1 -T 10 -O $POP_FILE -b 1 -i 1 -e $stacks_dir -S "
if [[ $MAC_LOCUS_STACKS != "" ]]
then
  cmd_line=$cmd_line" -X \"ustacks:--max_locus_stacks 3\""
fi
echo $cmd_line  > $SGE/denovo_map.sh

echo "DENOVO_MAP"

id_denovo_job=`qsub -terse -N ${Proj_Name}_denovomap -l mem=8G -l h_vmem=16G -pe parallel_smp 10 -o $SGE -e $SGE $SGE/denovo_map.sh`
r=0;
while [ $r == 0 ] ; do qstat -j $id_denovo_job >& /dev/null ; r=$? ; sleep 10s ;done

# séparation locus poly/monomorphe
grep -E 'Catalog|consensus' $OUT_DIR/batch_1.haplotypes.tsv > $OUT_DIR/locus_monomorphes.txt
	# pour les locus polymorphes ajout des couvertures total et allelique à chaque genotype individuel
grep -v consensus $OUT_DIR/batch_1.haplotypes.tsv > $OUT_DIR/tmp_poly
datamash transpose < $OUT_DIR/OUT_DIR/tmp_poly_transposed
rm $OUT_DIR/tmp_poly
python $SCRIPT_DIR/add_cov_to_tab.py -i $OUT_DIR/tmp_poly_transposed -d $OUT_DIR -o $OUT_DIR/tmp_poly_cov_transposed
rm $OUT_DIR/tmp_poly_transposed 
datamash transpose < $OUT_DIR/tmp_poly_cov_transposed > $OUT_DIR/locus_polymorphes_cov.txt
rm $OUT_DIR/tmp_poly_cov_transposed

# génération des fichiers de comptage allèlique/genotypes
python $SCRIPT_DIR/genetics_stat.py -d -g -i $OUT_DIR/locus_polymorphes_cov.txt -p $POP_FILE -o $OUT_DIR -n all_locus_polymorphes

echo "###################"
echo "Stacks summary"
python $SCRIPT_DIR/stacks_summary.py --res-dir $OUT_DIR --logfile $OUT_DIR/denovo_map.log --pop-map $POP_FILE --summary $STAT_DIR/summary.html --stacks-prog denovo_map.pl
