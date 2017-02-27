#!/bin/sh

#$ -o [Project_path_dir]/SGE_out/stacks1.out
#$ -e [Project_path_dir]/SGE_out/stacks1.err
#$ -N [Proj_Name]_stacks1
#$ -S /bin/sh

#$ -q unlimitq

#$ -m bea



date
date >&2
Proj_Name=[Proj_Name]

# VARIABLES GENERALES
# fonctionne aussi avec 1.44
stacks_dir=/usr/local/bioinfo/src/Stacks/stacks-1.42/bin
RAD_DIR=[Project_path_dir]
SCRIPT_DIR=[Script_path_dir]

# input
PREPROCESSING_DIR=$RAD_DIR/preprocessing/stacks_input
INDIV_FILE=$RAD_DIR/indiv_barcode.txt
dos2unix $INDIV_FILE
POP_FILE=$RAD_DIR/population.map
dos2unix $POP_FILE
POP=`cut -f 2 $POP_FILE | sort -u`

# path to genome fasta index (genome.fa, genome.fa.ann, genome.fa.amb, genome.fa.bwt, genome.fa.pac, genome.fa.sa). If not indexed please see bwa index command help
GENOME=""

# output
OUT_DIR=$RAD_DIR/stacks1
SGE=$RAD_DIR/SGE_out/`basename $OUT_DIR`
STAT_DIR=$RAD_DIR/stat/`basename $OUT_DIR`

mkdir -p $SGE $OUT_DIR $STAT_DIR

# OPTION
# ID of a control population double haploide
POP_HD=""

# projet DD Rad ? oui : DDRAD=1, non DDRAD=0
DDRAD=0

# attention 
# si ustacks alors min coverage per stacks ~= allele
# si pstacks alors min coverage per location ~= locus!
MIN_DEPTH=3

# autres options telles que :  --model_type --alpha --bound_low --bound_high --bc_err_freq. Ecrire la chaine de caractère entre "quotte"
STACKS_OPT=""

# SI USTACKS
MAX_PRIM_DIST=2
MAX_SEC_DIST=4
# --max_locus_stacks option pour individus diploïde. (Pour les individus happloïdes doublés cette option sera à 2)
MAC_LOCUS_STACKS=3


## TODO
#~ STAT sur qualité de clustering (dernière colonne de consensus)

####################################################################################################################
############################ U/PSTACKS ##############################
echo "############################ U/PSTACKS ##############################"
echo "############################ UP/STACKS ##############################" >&2

if [[ -e $SGE/stacks1.array ]]
then
rm $SGE/stacks1.array
fi

cat $INDIV_FILE | while read line
do
id=`echo $line | awk '{print $1}'`
indiv=`echo $line | awk '{print $2}'`
run=`echo $line | awk '{print $3}'`

# PSTACKS
if [[ $GENOME != "" ]]
then
cmd_line="sh $SCRIPT_DIR/pstacks.sh --genomeIndex=$GENOME --read1="`ls $PREPROCESSING_DIR/${indiv}_*1.fq.gz | grep -v "rem" `

rem1=""
if [[ -e `ls $PREPROCESSING_DIR/${indiv}_* | grep "rem.1.fq.gz" ` ]]
then cmd_line=$cmd_line" --rem1="`ls $PREPROCESSING_DIR/${indiv}_*rem.1.fq.gz`
fi

read2=""
if [[ -e `ls $PREPROCESSING_DIR/${indiv}_*2.fq.gz | grep -v "rem" ` ]]
then cmd_line=$cmd_line" --read2="`ls $PREPROCESSING_DIR/${indiv}_*2.fq.gz | grep -v "rem" `
fi

rem2=""
if [[ -e `ls $PREPROCESSING_DIR/${indiv}_* | grep "rem.2.fq.gz"` ]]
then cmd_line=$cmd_line" --rem2="`ls $PREPROCESSING_DIR/${indiv}_*rem.2.fq.gz`
fi

if [[ $DDRAD == 1 ]]
then  cmd_line=$cmd_line" --ddrad"
fi

echo $cmd_line" --id=$id --indiv=$indiv --out=$OUT_DIR --bin_dir=$stacks_dir --cov=$MIN_DEPTH --opt=\"$STACKS_OPT\"" >> $SGE/stacks1.array
#USTACKS
else
# configuration pour les individus happloïdes doublés ==> identifié par NUM_POP
is_HD=`awk -v I=$indiv -v P="$POP_HD" '{if($1==I && match(P,$2)){print "hd"}}' $POP_FILE`
if [[ "$is_HD" == "hd" ]]
then
echo "sh $SCRIPT_DIR/ustacks.sh $PREPROCESSING_DIR/${indiv}.fq.gz $id $OUT_DIR $stacks_dir $MIN_DEPTH $MAX_PRIM_DIST $MAX_SEC_DIST 2 \"$STACKS_OPT\"" >> $SGE/stacks1.array
else
echo "sh $SCRIPT_DIR/ustacks.sh $PREPROCESSING_DIR/${indiv}.fq.gz $id $OUT_DIR $stacks_dir $MIN_DEPTH $MAX_PRIM_DIST $MAX_SEC_DIST $MAC_LOCUS_STACKS \"$STACKS_OPT\"" >> $SGE/stacks1.array
fi
fi
done
 

id_stacks=`qarray -terse -N ${Proj_Name}_stacks1 -o $SGE -e $SGE $SGE/stacks1.array | awk -F "." '{print $1}'`
r=0;
while [ $r == 0 ] ; do qstat -j $id_stacks >& /dev/null ; r=$? ; sleep 10s ; done 

####################################################################################################################
############################ STAT ##############################
echo "############################ STAT ##############################"
echo "############################ STAT ##############################" >&2

# PSTACKS
if [[ $GENOME != "" ]]
then
python $SCRIPT_DIR/stacks_summary.py --stacks-prog pstacks --res-dir $OUT_DIR --pop-map $POP_FILE --summary $STAT_DIR/pstacks_summary.html
else
python $SCRIPT_DIR/stacks_summary.py --stacks-prog ustacks --res-dir $OUT_DIR --pop-map $POP_FILE --summary $STAT_DIR/ustacks_summary.html
fi
