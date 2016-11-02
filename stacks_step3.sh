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
# pour les stat nombre mini acceptable de locus avec correspondance
seuil_locus=0
# travaillez  vous sur des données alignées ou non. oui : align=1, non : align=0
align=1
# other sstacks options ? available option : -x
OPT=""

# VARIABLES GENERALES
stacks_dir=/usr/local/bioinfo/src/Stacks/stacks-1.42/bin
RAD_DIR=[Project_path_dir]
SCRIPT_DIR=[Script_path_dir]

# INPUT
STACKS1_DIR=$RAD_DIR/ustacks
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

if [[ $align == 1 ]]
then 
OPT=`echo $OPT" -g"`
fi

cat $INDIV_FILE | while read line
do
indiv=`echo $line | awk '{print $2}'`
echo "$SCRIPT_DIR/sstacks.sh $STACKS1_DIR/${indiv} $CSTACKS_DIR $SSTACKS_DIR $stacks_dir $OPT" >> $SGE/sstacks.qarray
done

id_sstacks=`qarray -terse -N ${Proj_Name}_sstacks -o $SGE -e $SGE $SGE/sstacks.qarray | awk -F "." '{print $1}'`
r=0;
while [ $r == 0 ] ; do qstat -j $id_sstacks >& /dev/null ; r=$? ; sleep 10s ; done 
#############################################################################################
###########################   SSTACK stat : $SSTACKS_DIR ######################################################
echo -e "\n###########################   SSTACK stat : $SSTACKS_DIR ######################################################"
echo -e "\n###########################   SSTACK stat : $SSTACKS_DIR ######################################################" >&2 
echo -e "###########################   SSTACK stat : $SSTACKS_DIR ######################################################"  > $STAT_DIR/summary_statistics.txt


echo -e "pop\tindiv\tmatching_locus\tmatching_haplotype" > $STAT_DIR/matches.count

cat $POP_FILE | while read line 
do
indiv=`echo $line | awk '{print $1}'`
pop=`echo $line | awk '{print $2}'`
stat=`zcat $SSTACKS_DIR/${indiv}.matches.tsv.gz | cut -f 3 | sort | uniq -c | colstat.sh `
loc=`echo $stat | awk '{print $2}'`
hap=`echo $stat | awk '{print $4}'`
echo -e "$pop\t$indiv\t$loc\t$hap" >> $STAT_DIR/matches.count
done

awk '{if($3==0){print $2}}' $STAT_DIR/matches.count > $SSTACKS_DIR/indiv_to_be_removed.txt

# stat sur le nombre restant de cluster par indiv
echo "stat sur le nombre correspondance de cluster par indiv" >> $STAT_DIR/summary_statistics.txt
grep -v pop $STAT_DIR/matches.count | cut -f 3 | colstat.sh >> $STAT_DIR/summary_statistics.txt
echo "détail du nombre de locus/haplotypes individuels correspondant à ceux du catalog: "$STAT_DIR/matches.count >> $STAT_DIR/summary_statistics.txt

echo "" >> $STAT_DIR/summary_statistics.txt
awk -v S=$seuil_locus 'BEGIN{n=0; m=0; o=0; p=0; q=0}{if ($3<10){n++} if ($3<100){m++} if($3<1000){o++} if($3<10000){p++} if($3<S){q++}}END{print "nb indiv < 10 correspondances : "n; print "nb indiv < 100 correspondances : "m; print "nb indiv < 1000 correspondances : "o; print "nb indiv < 10000 correspondances : "p; if (S>0){print "nb indiv < "S" correspondances : "q; }}' $STAT_DIR/matches.count >> $STAT_DIR/summary_statistics.txt

