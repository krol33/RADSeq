#$ -o [Project_path_dir]/SGE_out/stacks4.out
#$ -e [Project_path_dir]/SGE_out/stacks4.err
#$ -N RAD_Dev_stacks4
#$ -S /bin/sh

#$ -q unlimitq

#$ -m bea

date
date >&2
Proj_Name=RAD_Dev

# VARIABLES DEDIEES
# type d'analyse, croisement (ANALYSIS_TYPE="genotypes") ou populations (ANALYSIS_TYPE="populations")
ANALYSIS_TYPE="populations"
# nombre de locus minimum pour conserver un individus si vous ne voulez pas faire de filtre laissez à vide
seuil_locus=

# VARIABLES GENERALES
stacks_dir=/usr/local/bioinfo/src/Stacks/stacks-1.35/bin
RAD_DIR=[Project_path_dir]
SCRIPT_DIR=[Script_path_dir]
libR=[libR_path_dir]

# INPUT
STACKS1_DIR=$RAD_DIR/ustacks
CSTACKS_DIR=$RAD_DIR/cstacks
SSTACKS_DIR=$RAD_DIR/sstacks
INDIV_FILE=$RAD_DIR/indiv_barcode.txt
dos2unix $INDIV_FILE
POP_FILE=$RAD_DIR/population.map
dos2unix $POP_FILE
# si besoins de générer le tableau d'haplotypage sur une selection de population en particulier indiquez entre "" les numéros de population séparer par un espace.
EXTRACT_POP=""
# genotypes/populations options  : autre que -P -M -b -k -t --fstats (notamment les options de format de sortie  --genepop --fasta --structure --phase --beagle --plink --phylip --vcf ..., de  boostrap, de fusion et phasage, de filtre (attention à la version de stacks utilisé)
OPT="--vcf --fasta "

# output
OUT_DIR=$RAD_DIR/${ANALYSIS_TYPE}
SGE=$RAD_DIR/SGE_out/`basename $OUT_DIR`
STAT_DIR=$RAD_DIR/stat/`basename $OUT_DIR`

mkdir -p $SGE $OUT_DIR $STAT_DIR

###########################################################################################

if [[ -n $EXTRACT_POP ]]
then
	if [ -e $OUT_DIR/pop_extract.map ]
	then 
	rm $OUT_DIR/pop_extract.map
	fi
	for p in $EXTRACT_POP
	do
		awk -v P=$p '$2==P' $POP_FILE >> $OUT_DIR/pop_extract.map
	done
	POP_FILE=$OUT_DIR/pop_extract.map
	cat $POP_FILE | while read line
	do
		indiv=`echo $line | awk '{print $1}'`
		ln -s $STACKS1_DIR/${indiv}* $OUT_DIR/
		ln -s $SSTACKS_DIR/${indiv}* $OUT_DIR/
	done
	ln -s $CSTACKS_DIR/* $OUT_DIR/
else
	ln -s $STACKS1_DIR/* $OUT_DIR/
	ln -s $CSTACKS_DIR/* $OUT_DIR/
	ln -s $SSTACKS_DIR/* $OUT_DIR/
fi

# supression des individus sans correspondance
awk -v D=$OUT_DIR -v RM=$SSTACKS_DIR/indiv_to_be_removed.txt 'BEGIN{while(getline<RM>0){tab[$1]=1; system("rm "D"/*"$1"*")}}{if(tab[$1]!=1){print $0}}' $POP_FILE | sort -k 2,2n > $OUT_DIR/tmp.map
mv $OUT_DIR/tmp.map $OUT_DIR/pop_extract.map
POP_FILE=$OUT_DIR/pop_extract.map

echo qsub -N ${Proj_Name}_$ANALYSIS_TYPE -o $SGE/$ANALYSIS_TYPE.out -e $SGE/$ANALYSIS_TYPE.err -l h_vmem=32G -l mem=16G -pe parallel_smp 10 $SCRIPT_DIR/$ANALYSIS_TYPE.sh $OUT_DIR $stacks_dir $POP_FILE $SCRIPT_DIR $libR $STAT_DIR \"$OPT\"
id=`qsub -b y -terse -N ${Proj_Name}_$ANALYSIS_TYPE -o $SGE/$ANALYSIS_TYPE.out -e $SGE/$ANALYSIS_TYPE.err -l h_vmem=32G -l mem=16G -pe parallel_smp 10 $SCRIPT_DIR/$ANALYSIS_TYPE.sh $OUT_DIR $stacks_dir $POP_FILE $SCRIPT_DIR $libR $STAT_DIR \"$OPT\"| awk -F "." '{print $1}'`

r=0;
while [ $r == 0 ] ; do qstat -j $id >& /dev/null ; r=$? ; sleep 10s ; done 

if [[ -n $seuil_locus ]]
then
awk -v S=$seuil_locus '{if($3+$4+$5<S){print $1}}' $STAT_DIR/zygosity_rate.txt > $OUT_DIR/bad_indiv
echo `wc -l $OUT_DIR/bad_indiv | awk '{print $1}' ` "individus genotypés pour moins de "$seuil_locus" locus" >> $STAT_DIR/summary_statistics.txt
fi

