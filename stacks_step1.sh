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
stacks_dir=/usr/local/bioinfo/src/Stacks/stacks-1.42/bin
libR=[libR_path_dir]
RAD_DIR=[Project_path_dir]
SCRIPT_DIR=[Script_path_dir]

# input
PREPROCESSING_DIR=$RAD_DIR/preprocessing

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
POP_HD=""
# projet DD Rad ? oui : DDRAD=1, non DDRAD=0
DDRAD=0
# attention 
# si ustacks alors min coverage per stacks ~= allele
# si ustacks alors min coverage per location ~= locus!
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
	cmd_line="$SCRIPT_DIR/pstacks.sh --genomeIndex=$GENOME --read1="`ls $PREPROCESSING_DIR/stacks_input/${indiv}_*1.fq.gz | grep -v "rem" `
	
	rem1=""
	if [[ -e `ls $PREPROCESSING_DIR/stacks_input/${indiv}_* | grep "rem.1.fq.gz" ` ]]
		then cmd_line=$cmd_line" --rem1="`ls $PREPROCESSING_DIR/stacks_input/${indiv}_*rem.1.fq.gz`
	fi
	
	read2=""
	if [[ -e `ls $PREPROCESSING_DIR/stacks_input/${indiv}_*2.fq.gz | grep -v "rem" ` ]]
		then cmd_line=$cmd_line" --read2="`ls $PREPROCESSING_DIR/stacks_input/${indiv}_*2.fq.gz | grep -v "rem" `
	fi
	
	rem2=""
	if [[ -e `ls $PREPROCESSING_DIR/stacks_input/${indiv}_* | grep "rem.2.fq.gz"` ]]
		then cmd_line=$cmd_line" --rem2="`ls $PREPROCESSING_DIR/stacks_input/${indiv}_*rem.2.fq.gz`
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
		echo "$SCRIPT_DIR/ustacks.sh $PREPROCESSING_DIR/stacks_input/${indiv}.fq.gz $id $OUT_DIR $stacks_dir $MIN_DEPTH $MAX_PRIM_DIST $MAX_SEC_DIST 2 \"$STACKS_OPT\"" >> $SGE/stacks1.array
	else
		echo "$SCRIPT_DIR/ustacks.sh $PREPROCESSING_DIR/stacks_input/${indiv}.fq.gz $id $OUT_DIR $stacks_dir $MIN_DEPTH $MAX_PRIM_DIST $MAX_SEC_DIST $MAC_LOCUS_STACKS \"$STACKS_OPT\"" >> $SGE/stacks1.array
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
if [[ -e $SGE/gzip-d-stacks.qarray ]]
then 
rm $SGE/gzip-d-stacks.qarray
fi

cat $INDIV_FILE | while read line
do
indiv=`echo $line | awk '{print $2}'`
echo "gzip -d $OUT_DIR/${indiv}.*tsv.gz " >> $SGE/gzip-d-stacks.qarray
done

id_gzip_d=`qarray -terse -N ${Proj_Name}_gzip-d -o /dev/null -e /dev/null $SGE/gzip-d-stacks.qarray`

cat $INDIV_FILE | awk -v S=$SGE -v P=$Proj_Name -v POP=$POP_FILE -v I=$id_stacks 'BEGIN{while(getline<POP>0){tab[$1]=$2}}{cpt++; pop=tab[$2]; system( "mv "S"/"P"_stacks1.o"I"."cpt" "S"/"P"_stacks1.o"I"."cpt"."pop) ; system( "mv "S"/"P"_stacks1.e"I"."cpt" "S"/"P"_stacks1.e"I"."cpt"."pop) }'



echo "###########################   U/PSTACKS stat ######################################################" > $STAT_DIR/summary_statistics.txt
echo "###########################   U/PSTACKS stat ######################################################" >&2
echo "#READS"  >> $STAT_DIR/summary_statistics.txt
echo "#READS"  >&2

echo -e "\nnombre de lectures utilisées" >&2
echo -e "\nnombre de lectures utilisées" >> $STAT_DIR/summary_statistics.txt
grep "Number of utilized reads" $SGE/${Proj_Name}_stacks1.e*  | awk '{print $NF}' | colstat.sh >> $STAT_DIR/summary_statistics.txt

for p in $POP
do
echo -e "\tread used for the population $p : "  >> $STAT_DIR/summary_statistics.txt
grep "Number of utilized reads" $SGE/${Proj_Name}_stacks1.e*.*.$p | awk '{print $NF}'  | colstat.sh 1  >> $STAT_DIR/summary_statistics.txt
done

r=0;
while [ $r == 0 ] ; do qstat -j $id_gzip_d >& /dev/null ; r=$? ; sleep 10s ; done 

echo -e "\n#CLUSTER" >> $STAT_DIR/summary_statistics.txt
echo -e "\n#CLUSTER" >&2
# nombre de cluster
echo -e "nombre de clusters créés" >&2
echo -e "nombre de clusters créés" >> $STAT_DIR/summary_statistics.txt
if [[ -e $STAT_DIR/cluster.count ]]
then
rm $STAT_DIR/cluster.count
fi

echo -e "#pop\tindiv\ttotal_cluster\twhitlisted_cluster" > $STAT_DIR/cluster.count
cat $POP_FILE | while read line
do
indiv=`echo $line | awk '{print $1}'`
pop=`echo $line | awk '{print $2}'`
tags=$OUT_DIR/$indiv.tags.tsv
grep consensus $tags | awk -v I=$indiv -v P=$pop '{n++; if($9!=1){m++}}END{print P,I,n,m}' >> $STAT_DIR/cluster.count
done

echo -e "\ttotal" >> $STAT_DIR/summary_statistics.txt
grep -v "#" $STAT_DIR/cluster.count | colstat.sh 3  >> $STAT_DIR/summary_statistics.txt
echo -e "\twhitelisted" >> $STAT_DIR/summary_statistics.txt
grep -v "#" $STAT_DIR/cluster.count | colstat.sh 4  >> $STAT_DIR/summary_statistics.txt

echo -e "\ttotal/whitelisted cluster per population : " >> $STAT_DIR/summary_statistics.txt
for p in $POP
do
echo -e "\tpopulation $p" >> $STAT_DIR/summary_statistics.txt
awk -v P=$p '$1==P' $STAT_DIR/cluster.count | colstat.sh 3 >> $STAT_DIR/summary_statistics.txt
awk -v P=$p '$1==P' $STAT_DIR/cluster.count | colstat.sh 4 >> $STAT_DIR/summary_statistics.txt
done

$libR/plot_hist.R "nb cluster per individuals" "individual" "nb cluster" $STAT_DIR/nbcluster_per_indiv $STAT_DIR/cluster.count 4

#~ barplot(t$V4, xlab="individual", ylab="number of whitelisted cluster", main="number of cluster par individual\nconf m5", col=ifelse(t$V1==1,"red",ifelse(t$V1==2,"green",ifelse(t$V1==3,"blue",ifelse(substr(t$V2,1,7)=="GenoRad","grey50",ifelse(substr(t$V2,1,3)=="DEX","grey70","black"))))), ylim=c(1,100000)) 
#~ legend("top",c("maigre","grasse","control","fondateur"),col=c("red","green","blue","black"),lty=c(1,1,1,1), ncol=4)



echo -e "\n\t(Note: La suite des analyses ne concerne que les cluster whitelistés)" >> $STAT_DIR/summary_statistics.txt
grep consensus $OUT_DIR/*tags.tsv  |  awk '{if($9==1){print $2"_"$3}}' > $STAT_DIR/blacklisted_cluster.txt

# couverture des clusters
echo -e "\ncouverture des clusters" >&2
echo -e "\ncouverture des clusters"  >> $STAT_DIR/summary_statistics.txt
echo -e "\tall"  >> $STAT_DIR/summary_statistics.txt
grep -vE 'consensus|model|#' $OUT_DIR/*tags.tsv  |awk -v IND=$INDIV_FILE -v POP=$POP_FILE -v BF=$STAT_DIR/blacklisted_cluster.txt 'BEGIN{while(getline<IND>0){tabI[$1]=$2} while(getline<POP>0){tabP[$1]=$2} while(getline<BF>0){tab[$1]=1}}{if(tab[$2"_"$3]!=1){p=tabP[tabI[$2]];  print p,$2,$3}}' | uniq -c > $STAT_DIR/tmp

cat $STAT_DIR/tmp | colstat.sh  >> $STAT_DIR/summary_statistics.txt


for p in $POP
do
echo -e "\tpopulation $p">> $STAT_DIR/summary_statistics.txt
awk -v P=$p '$2==P' $STAT_DIR/tmp | colstat.sh   >> $STAT_DIR/summary_statistics.txt
done

# nombre de lecture secondaires par cluster
echo -e "\nnombre de lecture secondaires" >&2
echo -e "\nnombre de lecture secondaires" >> $STAT_DIR/summary_statistics.txt
echo -e "\tall" >> $STAT_DIR/summary_statistics.txt
grep secondary $OUT_DIR/*tags.tsv  | awk -v IND=$INDIV_FILE -v POP=$POP_FILE -v BF=$STAT_DIR/blacklisted_cluster.txt  'BEGIN{while(getline<IND>0){tabI[$1]=$2} while(getline<POP>0){tabP[$1]=$2} while(getline<BF>0){tab[$1]=1}}{if(tab[$2"_"$3]!=1){p=tabP[tabI[$2]]; print p,$2,$3}}' | uniq -c > $STAT_DIR/tmp
cat $STAT_DIR/tmp | colstat.sh  >> $STAT_DIR/summary_statistics.txt

for p in $POP
do
echo -e "\tpopulation $p" >> $STAT_DIR/summary_statistics.txt
awk -v P=$p '$2==P' $STAT_DIR/tmp | colstat.sh   >> $STAT_DIR/summary_statistics.txt
done

# nombre de stacks
echo -e "\n#STACKS" >> $STAT_DIR/summary_statistics.txt
echo -e "\n#STACKS"  >&2
echo -e "nombre de stacks par cluster" >&2
echo -e "nombre de stacks par cluster" >> $STAT_DIR/summary_statistics.txt
echo -e "\tall" >> $STAT_DIR/summary_statistics.txt
grep primary $OUT_DIR/*tags.tsv  |  awk -v IND=$INDIV_FILE -v POP=$POP_FILE -v BF=$STAT_DIR/blacklisted_cluster.txt 'BEGIN{while(getline<IND>0){tabI[$1]=$2} while(getline<POP>0){tabP[$1]=$2} while(getline<BF>0){tab[$1]=1}}{if(tab[$2"_"$3]!=1){p=tabP[tabI[$2]]; print p,$2,$3,$5}}' | sort -u | awk '{print $1,$2,$3}' | sort | uniq -c > $STAT_DIR/tmp
cat $STAT_DIR/tmp | colstat.sh  >> $STAT_DIR/summary_statistics.txt

for p in $POP
do
echo -e "\tpopulation $p" >> $STAT_DIR/summary_statistics.txt
awk -v P=$p '$2==P' $STAT_DIR/tmp | colstat.sh   >> $STAT_DIR/summary_statistics.txt
done

# nombre de lectures par stacks
echo -e "\nnombre de lectures par stacks" >> $STAT_DIR/summary_statistics.txt
echo -e "\nnombre de lectures par stacks" >&2
echo -e "\tall" >> $STAT_DIR/summary_statistics.txt
grep primary $OUT_DIR/*tags.tsv  | awk -v IND=$INDIV_FILE -v POP=$POP_FILE -v BF=$STAT_DIR/blacklisted_cluster.txt 'BEGIN{while(getline<IND>0){tabI[$1]=$2} while(getline<POP>0){tabP[$1]=$2}  while(getline<BF>0){tab[$1]=1}}{if(tab[$2"_"$3]!=1){p=tabP[tabI[$2]]; print p,$2,$3,$5}}' | uniq -c > $STAT_DIR/tmp
cat $STAT_DIR/tmp | colstat.sh  >> $STAT_DIR/summary_statistics.txt

for p in $POP
do
echo -e "\tpopulation $p" >> $STAT_DIR/summary_statistics.txt
awk -v P=$p '$2==P' $STAT_DIR/tmp | colstat.sh   >> $STAT_DIR/summary_statistics.txt
done

# nombre de SNP détectés
echo -e "\nSNP" >> $STAT_DIR/summary_statistics.txt
echo -e "\nSNP" >&2
echo -e "nombre de SNPs détectés ">&2
echo -e "nombre de SNPs détectés " >> $STAT_DIR/summary_statistics.txt
echo -e "\tall:" >> $STAT_DIR/summary_statistics.txt

grep -v "#" $OUT_DIR/*.snps.tsv | awk '$5=="E"' | cut -f 2,3 | uniq -c | awk -v IND=$INDIV_FILE -v POP=$POP_FILE -v BF=$STAT_DIR/blacklisted_cluster.txt 'BEGIN{while(getline<IND>0){tabI[$1]=$2} while(getline<POP>0){tabP[$1]=$2}  while(getline<BF>0){tab[$1]=1}} {if(tab[$2"_"$3]!=1){p=tabP[tabI[$2]]; print p,$0}}' > $STAT_DIR/tmp

stat=`cat $STAT_DIR/tmp | colstat.sh 2`
echo $stat  >> $STAT_DIR/summary_statistics.txt
tot=`echo $stat | awk '{print $2}'`
cat $STAT_DIR/tmp | awk -v TOT=$tot '{if($2==1){n++} if($2==2){m++} if($2==3){p++}}END{print "1:"n*100/TOT"; 2: "m*100/TOT"; 3: "p*100/TOT}'   >> $STAT_DIR/summary_statistics.txt

for p in $POP
do
echo -e "\tpopulation $p" >> $STAT_DIR/summary_statistics.txt
stat=`awk -v P=$p '$1==P' $STAT_DIR/tmp | colstat.sh 2`
echo $stat  >> $STAT_DIR/summary_statistics.txt
tot=`echo $stat | awk '{print $2}'`
awk -v P=$p '$1==P' $STAT_DIR/tmp | awk -v TOT=$tot '{if($2==1){n++} if($2==2){m++} if($2==3){p++}}END{print "1:"n*100/TOT"; 2: "m*100/TOT"; 3: "p*100/TOT}'   >> $STAT_DIR/summary_statistics.txt
done
# couverture des clusters avec SNP
echo -e "\ncouverture des clusters avec SNP " >&2
echo -e "\ncouverture des clusters avec SNP " >> $STAT_DIR/summary_statistics.txt
echo -e "\tall:" >> $STAT_DIR/summary_statistics.txt
grep -v "#" $OUT_DIR/*.alleles.tsv  | awk -v IND=$INDIV_FILE -v POP=$POP_FILE -v BF=$STAT_DIR/blacklisted_cluster.txt 'BEGIN{while(getline<IND>0){tabI[$1]=$2} while(getline<POP>0){tabP[$1]=$2} while(getline<BF>0){vu[$1]=1}}{if(vu[$2"_"$3]!=1){p=tabP[tabI[$2]]; tab[p"\t"$2"_"$3]=tab[p"\t"$2"_"$3]+$6}}END{for (i in tab){print i,tab[i]}}' > $STAT_DIR/tmp
cat $STAT_DIR/tmp | colstat.sh 3  >> $STAT_DIR/summary_statistics.txt

for p in $POP
do
echo -e "\tpopulation $p" >> $STAT_DIR/summary_statistics.txt
awk -v P=$p '$1==P' $STAT_DIR/tmp | colstat.sh 3  >> $STAT_DIR/summary_statistics.txt
done
# nombre de cluster en fonction du nombre d'haplotype
echo -e "\nnombre de cluster en fonction du nombre d'haplotype " >&2
echo -e "\nnombre de cluster en fonction du nombre d'haplotype " >> $STAT_DIR/summary_statistics.txt
echo -e "\tall" >> $STAT_DIR/summary_statistics.txt
awk 'NF==6' $OUT_DIR/*.alleles.tsv | grep -v "#" |  awk -v IND=$INDIV_FILE -v POP=$POP_FILE -v BF=$STAT_DIR/blacklisted_cluster.txt 'BEGIN{while(getline<IND>0){tabI[$1]=$2} while(getline<POP>0){tabP[$1]=$2} while(getline<BF>0){vu[$1]=1}}{if(vu[$2"_"$3]!=1){p=tabP[tabI[$2]];  print p,$2,$3}}' | uniq -c > $STAT_DIR/tmp
awk '$1==1' $STAT_DIR/tmp | wc -l | awk '{print "\tnombre de cluster avec 1 happlotypes: "$1}'  >> $STAT_DIR/summary_statistics.txt
awk '$1==2' $STAT_DIR/tmp | wc -l | awk '{print "\tnombre de cluster avec 2 happlotypes: "$1}'  >> $STAT_DIR/summary_statistics.txt
awk '$1==3' $STAT_DIR/tmp | wc -l| awk '{print "\tnombre de cluster avec 3 happlotypes: "$1}'  >> $STAT_DIR/summary_statistics.txt
awk '$1==4' $STAT_DIR/tmp | wc -l| awk '{print "\tnombre de cluster avec 4 happlotypes: "$1}'  >> $STAT_DIR/summary_statistics.txt
awk '$1>=5' $STAT_DIR/tmp | wc -l| awk '{print "\tnombre de cluster avec >=5 happlotypes: "$1}'  >> $STAT_DIR/summary_statistics.txt
awk '$1>=10' $STAT_DIR/tmp | wc -l| awk '{print "\tnombre de cluster avec >=10 happlotypes: "$1}'  >> $STAT_DIR/summary_statistics.txt
cat $STAT_DIR/tmp | sort -k1,1n | tail -n 1 | awk '{print "\tnombre max d haplotypes: "$1}'  >> $STAT_DIR/summary_statistics.txt

for p in $POP
do
echo -e "\tpopulation $p" >> $STAT_DIR/summary_statistics.txt
awk -v P=$p '$1==1 && $2==P' $STAT_DIR/tmp | wc -l | awk '{print "\tnombre de cluster avec 1 happlotypes: "$1}'  >> $STAT_DIR/summary_statistics.txt
awk -v P=$p '$1==2 && $2==P' $STAT_DIR/tmp | wc -l | awk '{print "\tnombre de cluster avec 2 happlotypes: "$1}'  >> $STAT_DIR/summary_statistics.txt
awk -v P=$p '$1==3 && $2==P' $STAT_DIR/tmp | wc -l| awk '{print "\tnombre de cluster avec 3 happlotypes: "$1}'  >> $STAT_DIR/summary_statistics.txt
awk -v P=$p '$1==4 && $2==P' $STAT_DIR/tmp | wc -l| awk '{print "\tnombre de cluster avec 4 happlotypes: "$1}'  >> $STAT_DIR/summary_statistics.txt
awk -v P=$p '$1>=5 && $2==P' $STAT_DIR/tmp | wc -l| awk '{print "\tnombre de cluster avec >=5 happlotypes: "$1}'  >> $STAT_DIR/summary_statistics.txt
awk -v P=$p '$1>=10 && $2==P' $STAT_DIR/tmp | wc -l| awk '{print "\tnombre de cluster avec >=10 happlotypes: "$1}'  >> $STAT_DIR/summary_statistics.txt
awk -v P=$p '$2==P' $STAT_DIR/tmp | sort -k1,1n | tail -n 1 | awk '{print "\tnombre max d haplotypes: "$1}'  >> $STAT_DIR/summary_statistics.txt
done

rm $STAT_DIR/tmp

if [[ -e $SGE/gzip_ustacks.qarray ]]
then 
rm $SGE/gzip_ustacks.qarray
fi

cat $INDIV_FILE | while read line
do
indiv=`echo $line | awk '{print $2}'`
echo "gzip $OUT_DIR/${indiv}.*tsv " >> $SGE/gzip_ustacks.qarray
done
id_gzip=`qarray -terse -N ${Proj_Name}_gzip_ustacks  -o /dev/null -e /dev/null $SGE/gzip_ustacks.qarray`
r=0;
while [ $r == 0 ] ; do qstat -j $id_gzip >& /dev/null ; r=$? ; sleep 10s ; done 
