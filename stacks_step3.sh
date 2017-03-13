#!/bin/sh

#$ -o [Project_path_dir]/SGE_out/stacks3.out
#$ -e [Project_path_dir]/SGE_out/stacks3.err
#$ -N [Proj_Name]_stacks3
#$ -S /bin/sh

#$ -q unlimitq

#$ -m bea

date
date >&2
Proj_Name=[Proj_Name]

#VARIABLE DEDIEE
# travaillez  vous sur des données alignées ou non. oui : align=1, non : align=0
align=0
# other sstacks options ? available option : -x
OPT=""

# VARIABLES GENERALES
stacks_dir=/usr/local/bioinfo/src/Stacks/stacks-1.42/bin
RAD_DIR=[Project_path_dir]
SCRIPT_DIR=[Script_path_dir]

# INPUT
STACKS1_DIR=$RAD_DIR/stacks1
CSTACKS_DIR=$RAD_DIR/cstacks
INDIV_FILE=$RAD_DIR/indiv_barcode.txt
dos2unix $INDIV_FILE
POP_FILE=$RAD_DIR/population.map
dos2unix $POP_FILE

# output
SSTACKS_DIR=$RAD_DIR/sstacks
SGE=$RAD_DIR/SGE_out/`basename $SSTACKS_DIR`
STAT_DIR=$RAD_DIR/stat/`basename $SSTACKS_DIR`
mkdir -p $SGE $SSTACKS_DIR $STAT_DIR

# TO DO
# plot score de vraissemblance dernière colonne

###########################################################################################

if [[ -e $SGE/sstacks.qarray ]]
then
rm $SGE/sstacks.qarray
fi

if [[ $align == 1 ]]
then 
OPT=`echo $OPT" -g"`
fi

cat $INDIV_FILE | while read line
do
indiv=`echo $line | awk '{print $2}'`
echo "sh $SCRIPT_DIR/sstacks.sh $STACKS1_DIR/${indiv} $CSTACKS_DIR $SSTACKS_DIR $stacks_dir $OPT" >> $SGE/sstacks.qarray
done

id_sstacks=`qarray -terse -N ${Proj_Name}_sstacks -o $SGE -e $SGE $SGE/sstacks.qarray | awk -F "." '{print $1}'`
r=0;
while [ $r == 0 ] ; do qstat -j $id_sstacks >& /dev/null ; r=$? ; sleep 10s ; done 
#############################################################################################
###########################   SSTACK stat : $SSTACKS_DIR ######################################################
echo -e "\n###########################   SSTACK stat : $SSTACKS_DIR ######################################################"
echo -e "\n###########################   SSTACK stat : $SSTACKS_DIR ######################################################" >&2 

python $SCRIPT_DIR/stacks_summary.py --stacks-prog sstacks --res-dir $SSTACKS_DIR --pop-map $POP_FILE --summary $STAT_DIR/sstacks_summary.html
