#!/bin/sh

#$ -o [Project_path_dir]/SGE_out/stacks2.out
#$ -e [Project_path_dir]/SGE_out/stacks2.err
#$ -N [Proj_Name]_stacks2
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

# INPUT
STACKS1_DIR=$RAD_DIR/stacks1
INDIV_FILE=$RAD_DIR/indiv_barcode.txt
dos2unix $INDIV_FILE
POP_FILE=$RAD_DIR/population.map
dos2unix $POP_FILE

# output
CSTACKS_DIR=$RAD_DIR/cstacks
SGE=$RAD_DIR/SGE_out/`basename $CSTACKS_DIR`
STAT_DIR=$RAD_DIR/stat/`basename $CSTACKS_DIR`

mkdir -p $SGE $CSTACKS_DIR $STAT_DIR

# VARIABLE DEDIEE
MISMATCH=1
# voulez vous baser le catalogue sur un autre ? (attention au version) 
KNOWN_CAT=
# travaillez  vous sur des données alignées ou non. oui : align=1, non : align=0
align=0

# INDIVIDUALS SELECTION (1 option au choix, si aucune option tous les individus sont pris en compte)
POP_IN=""
POP_OUT=""
INDIV_IN=""
INDIV_OUT=""

# VARIABLE CLUSTER
MEM="16G"
VMEM="24G"
THREADS=10

###########################################################################################
### catalogue 

### SELECTION DES INDIV 

if [[ -n $POP_IN ]]
then
	if [[ -n $POP_OUT || -n $INDIV_IN || -n $INDIV_OUT ]]
	then
	echo "WRONG1, you must choose one selection option exclusively" >&2
	exit
	fi
else
	if [[ -n $POP_OUT ]]
	then
		if [[ -n $INDIV_IN || -n $INDIV_OUT ]]
		then
		echo "WRONG2, you must choose one selection option exclusively" >&2
		exit
		fi
	else
		if [[ -n $INDIV_IN && -n $INDIV_OUT ]]
		then
		echo "WRONG3, you must choose one selection option exclusively" >&2
		exit
		fi
	fi
fi

args=""

# tous ceux des pop POP_IN
if [[ -n $POP_IN ]]
then
args=`awk -v P="$POP_IN" -v U=$STACKS1_DIR '{if (match(P,$2)){s=s" -s "U"/"$1}}END{print s}' $POP_FILE `
fi

# tous sauf ceux des pop POP_OUT
if [[ -n $POP_IN ]]
then
args=`awk -v P="$POP_OUT" -v U=$STACKS1_DIR '{if (! match(P,$2)){s=s" -s "U"/"$1}}END{print s}' $POP_FILE `
fi

# tous les individus INDIV_IN
if [[ -n $INDIV_IN ]]
then
args=`awk -v P="$INDIV_IN" -v U=$STACKS1_DIR '{if (match(P,$1)){s=s" -s "U"/"$1}}END{print s}' $POP_FILE `
fi

# tous les individus sauf INDIV_OUT
if [[ -n $INDIV_OUT ]]
then
args=`awk -v P="$INDIV_OUT" -v U=$STACKS1_DIR '{if (! match(P,$1)){s=s" -s "U"/"$1}}END{print s}' $POP_FILE `
fi

if [[ ! -n "$args" ]]
then
args=`awk -v U=$STACKS1_DIR '{s=s" -s "U"/"$1}END{print s}' $POP_FILE`
fi

cmd="module load compiler/gcc-4.9.1 ; $stacks_dir/cstacks -b 1 $args -o $CSTACKS_DIR -p $THREADS -n $MISMATCH "

# à partir d'un catalogue 
if [[ -n $KNOWN_CAT ]]
then
cmd=`echo $cmd --catalog $KNOWN_CAT`
fi

# à partir de données alignées
if [[ $align == 1 ]]
then
cmd=`echo $cmd -g`
fi

echo $cmd > $SGE/cstacks.tmp.sh
id_cstacks=`qsub -terse -N ${Proj_Name}_cstacks -o $SGE/cstacks.out -e $SGE/cstacks.err -l mem=$MEM -l h_vmem=$VMEM -pe parallel_smp $THREADS $SGE/cstacks.tmp.sh`
r=0;
while [ $r == 0 ] ; do qstat -j $id_cstacks >& /dev/null ; r=$? ; sleep 10s ;done 

####################################################################################################################
###########################   CSTACK stat : $CSTACKS_DIR ######################################################
echo -e "\n###########################   CSTACK stat : $CSTACKS_DIR ######################################################"
echo -e "\n###########################   CSTACK stat : $CSTACKS_DIR ######################################################" >&2 

python $SCRIPT_DIR/stacks_summary.py --stacks-prog cstacks --res-dir $CSTACKS_DIR --pop-map $POP_FILE --summary $STAT_DIR/cstacks_summary.html --logfile $SGE/cstacks.err
