#!/bin/sh

#$ -o [Project_path_dir]/SGE_out/stacks0.out
#$ -e [Project_path_dir]/SGE_out/stacks0.err
#$ -N [Proj_Name]_stacks0
#$ -S /bin/sh

#$ -q unlimitq

#$ -m bea

# WARNING
# dans le fichier data les lectures doivent être stockées dans du _1.fq.gz et _2.fq.gz


date
date >&2
Proj_Name=[Proj_Name]

# VARIABLES GENERALES
stacks_dir=/usr/local/bioinfo/src/Stacks/stacks-1.35/bin
RAD_DIR=[Project_path_dir]
SCRIPT_DIR=[Script_path_dir]

# input
DATA_DIR=$RAD_DIR/data
INDIV_FILE=$RAD_DIR/indiv_barcode.txt
dos2unix $INDIV_FILE
POP_FILE=$RAD_DIR/population.map
dos2unix $POP_FILE

# output
OUT_DIR=$RAD_DIR/preprocessing
mkdir -p $OUT_DIR
SGE=$RAD_DIR/SGE_out/`basename $OUT_DIR`
mkdir -p $SGE
STAT_DIR=$RAD_DIR/stat/`basename $OUT_DIR`
mkdir -p $STAT_DIR
cd $OUT_DIR
Run=`cut -f 3 $INDIV_FILE  | sort -u | awk '{dir=dir $1" "}END{print dir}'`
mkdir -p $Run

# VARIABLE DEDIEES
# enzyme utilisée, séparée par un espace en cas de double digestion.
# enzymes possibles:
#	  'ageI', 'aluI', 'apeKI', 'apoI', 'bamHI', 'bgIII', 'bstYI', 'claI', 
#      'ddeI', 'dpnII', 'eaeI', 'ecoRI', 'ecoRV', 'ecoT22I', 'hindIII', 'kpnI', 
#      'mluCI', 'mseI', 'mspI', 'ndeI', 'nheI', 'nlaIII', 'notI', 'nsiI', 
#      'pstI', 'rsaI', 'sacI', 'sau3AI', 'sbfI', 'sexAI', 'sgrAI', 'speI', 
#      'sphI', 'taqI', 'xbaI', or 'xhoI'
Enzyme="sbfI"
double_digest=`echo $Enzyme | awk '{if(NF>1){print "1"}else{print "0"}}'`

#voulez vous travaillez sur des données alignées ou non. oui : align=1, non : align=0
align=0

#avez vous des données pairées ? oui : paired=1, non : paired=0
paired=1
if [[ $double_digest = 1 && $paired = 0 ]]
then
echo "Error. You can not have double digest and single end sequencing"
fi

# voulez vous supprimer les duplicats PCR oui : dereplication=1, non : dereplication=0. 
# !! attention vallable unique sur single digest pair end!!
dereplication=0
if [[ $double_digest = 1 && $dereplication = 1 ]]
then
echo "Error. You can not dereplicate double digest data"
exit 1
fi

if [[ $paired = 0 && $dereplication = 1 ]]
then
echo "Error. You can not dereplicate single end data"
exit 1
fi


# nombre de mismatch autorisé sur les barcodes BMM sur le tag de restriction qui suit TMM
BMM=0
# si le barcode n'est présent que sur la lecture 1 alors TRIM=1 sinon TRIM=2
TRIM=1
# Si vous avez des taille de barcodes différentes, alors après trimming de barcodes vos lectures seront de taille différente. Il est nécessaire pour la suite d'avoir des lectures de même taille. Indiquez ici la taille maximum que vous souhaitez conserver.
LEN=


############################ DEMULTIPLEXAGE ##############################
rm -rf $SGE/*

echo "DEMULTIPLEXAGE"
echo "DEMULTIPLEXAGE" >&2

if [[ $(ls */ | grep -c barcode ) -gt 0 ]]
then
	rm */barcode*
fi

# si les séquences de vos barcodes sont de différentes tailles alors il faut démultipléxer en plusieurs fois en commençant par les baroce de plus grande taille.
cat $INDIV_FILE | awk '{bl=length($4); print $2,$4 >> $3"/barcode_"bl}'

# Pour chaque run on va lancer le démultiplexage sur chaque fichier barcode en fonction de leur taille.
nj=0
for run in $Run
do
lastjobid=""
f1=`echo $DATA_DIR/${run}_1.fq.gz`
f2=`echo $DATA_DIR/${run}_2.fq.gz`
for bl in `ls $run/barcode_*  | sed -e s%$run/barcode_%% | sort -nr `
do
name=${Proj_Name}_preprocessing_${run}_${bl}
if [ "$lastjobid" = "" ]
then
echo "CMD:	qsub -N $name -o $SGE/${Proj_Name}_preprocessing_${run}_${bl}.out -e $SGE/${Proj_Name}_preprocessing_${run}_${bl}.err $SCRIPT_DIR/preprocessing.sh $SCRIPT_DIR  $f1 $f2 $OUT_DIR/$run/ $run/barcode_${bl} $BMM $TRIM"
qsub -N $name -o $SGE/${Proj_Name}_preprocessing_${run}_${bl}.out -e $SGE/${Proj_Name}_preprocessing_${run}_${bl}.err $SCRIPT_DIR/preprocessing.sh $SCRIPT_DIR  $f1 $f2 $OUT_DIR/$run/ $run/barcode_${bl} $BMM $TRIM
else
echo "CMD:	qsub -hold_jid $lastjobid -N $name -o $SGE/${Proj_Name}_preprocessing_${run}_${bl}.out -e $SGE/${Proj_Name}_preprocessing_${run}_${bl}.err $SCRIPT_DIR/preprocessing.sh $SCRIPT_DIR $f1 $f2 $OUT_DIR/$run/ $run/barcode_${bl} $BMM $TRIM"
qsub -hold_jid $lastjobid -N $name -o $SGE/${Proj_Name}_preprocessing_${run}_${bl}.out -e $SGE/${Proj_Name}_preprocessing_${run}_${bl}.err $SCRIPT_DIR/preprocessing.sh $SCRIPT_DIR $f1 $f2 $OUT_DIR/$run/ $run/barcode_${bl} $BMM $TRIM
fi
lastjobid=$name
f1=$OUT_DIR/$run/unmatched_1.fq
f2=$OUT_DIR/$run/unmatched_2.fq
let nj=$nj+1
done
done


# nombre de jobs finis
finished=0
jf=0
# tant que tous les jobs ne sont pas terminés
until [[ $(ls $SGE/ | grep -c preprocessing_) -gt 0 && $jf -eq $nj ]]
do 
# compter le nombre de jobs terminés
if [[ $(ls $SGE | grep preprocessing | wc -l) -gt 0 ]]
then
jf=`grep "Epilog" $SGE/${Proj_Name}_preprocessing*.out | wc -l`
else
jf=0
fi
# affiché le nombre de jobs terminés s'il est différent d'avant
if [[ "$finished" -ne "$jf" ]]
then
echo "    $jf paire filtered among $nj"
finished=$jf
fi
sleep 10s
done

jf=`grep "Epilog" $SGE/${Proj_Name}_preprocessing*.out | wc -l`

if [[ "$finished" -ne "$jf" ]]
then
jf=`grep "Epilog" $SGE/${Proj_Name}_preprocessing*.out | wc -l`
echo "    $jf paire filtered among $nj"
fi

for run in $Run
do
if [[ -e $OUT_DIR/$run/unmatched_1.fq ]]
then
qsub -N ${Proj_Name}_gzip_$run -o /dev/null -e /dev/null -b y "gzip $OUT_DIR/$run/ambiguous_* $OUT_DIR/$run/unmatched_* $OUT_DIR/$run/*2rad*"
fi

done
date
############################ FILTRE CLONE ##############################
if [[ $dereplication = 1 ]]
then
date
echo -e "\nFILTRE CLONE: clone_filter"
echo -e "\nFILTRE CLONE: clone_filter" >&2
if [[ $(ls $SGE/ | grep -c clone_filter) -gt 0  ]]
then
rm $SGE/*clone_filter*
rm $run/clone_filter/*
fi

cat $INDIV_FILE | while read line
do
run=`echo $line | awk '{print $3}'`
mkdir -p $run/clone_filter
indiv=`echo $line | awk '{print $2}'`
echo "$SCRIPT_DIR/clone_filter.sh `ls $run/${indiv}_1.fq*` `ls $run/${indiv}_2.fq*` $run/clone_filter $stacks_dir" >> $SGE/clone_filter_$run.qarray
done

job_id_clone_filter=""
for run in $Run
do
id=`qarray -terse -N ${Proj_Name}_clone_filter_$run -o $SGE -e $SGE $SGE/clone_filter_$run.qarray`
job_id_clone_filter=`echo $job_id_clone_filter $id`
done

for id in $job_id_clone_filter
do
r=0;
while [ $r == 0 ] ; do qstat -j $id >& /dev/null ; r=$? ; sleep 10s ; done 
done
fi

############################ FILTRE QUAL ##############################
date
echo -e "\nFILTRE QUALITE: process_rad_tags"
echo -e "\nFILTRE QUALITE: process_rad_tags" >&2
if [[ $(ls $SGE/ | grep -c process_rad_tags) -gt 0  ]]
then
rm $SGE/*process_rad_tags*
fi

cat $INDIV_FILE | while read line
do
run=`echo $line | awk '{print $3}'`
mkdir -p $run/checkRAD
indiv=`echo $line | awk '{print $2}'`
if [[ $double_digest = 0 ]]
then
if [[ $dereplication = 1 ]]
then
echo "sh $SCRIPT_DIR/process_rad_tags.sh `ls $run/clone_filter/${indiv}_cloneF_1.fq.gz` $run/checkRAD $Enzyme $stacks_dir $LEN " >> $SGE/process_rad_tags_$run.qarray
else
echo "sh $SCRIPT_DIR/process_rad_tags.sh `ls $run/${indiv}_1.fq*` $run/checkRAD $Enzyme $stacks_dir $LEN " >> $SGE/process_rad_tags_$run.qarray
fi
else
echo "sh $SCRIPT_DIR/process_ddrad_tags.sh `ls $run/${indiv}_1.fq*` `ls $run/${indiv}_2.fq*` $run/checkRAD $Enzyme $stacks_dir $LEN " >> $SGE/process_rad_tags_$run.qarray
fi
done

job_id_process_rad_tags=""
for run in $Run
do
id=`qarray -terse -N ${Proj_Name}_process_rad_tags_$run -o $SGE -e $SGE $SGE/process_rad_tags_$run.qarray`
job_id_process_rad_tags=`echo $job_id_process_rad_tags $id`
done

for id in $job_id_process_rad_tags
do
r=0;
while [ $r == 0 ] ; do qstat -j $id >& /dev/null ; r=$? ; sleep 10s ; done 
done
############################ REGENERER LES PAIRES DE SEQUENCE ##############################
# le filtre qualité n'a été lancé que sur la lecture 1. Cela génère 2 fichiers par fastq. Il s'agit maintenant de reconstituer des fichiers pairés.
if [[ $paired = 1 && $double_digest = 0 ]]
then
date
echo -e "\nREORDONNANCEMENT DES PAIRES DE FICHIERS FASTQ"
echo -e "\nREORDONNANCEMENT DES PAIRES DE FICHIERS FASTQ" >&2
cd $RAD_DIR

if [[ $(ls $SGE/ | grep -c recoverMate) -gt 0  ]]
then
rm $SGE/*recoverMate*
fi

cat $INDIV_FILE | while read line
do
run=`echo $line | awk '{print $3}'`
indiv=`echo $line | awk '{print $2}'`
if [[ $dereplication = 1 ]]
then
echo "sh $SCRIPT_DIR/recover_mate.sh `ls $OUT_DIR/$run/checkRAD/${indiv}_cloneF_1*discards` `ls $OUT_DIR/$run/clone_filter/${indiv}_cloneF_2.fq*`" >> $SGE/recoverMateDiscards.qarray
echo "sh $SCRIPT_DIR/recover_mate.sh $OUT_DIR/$run/checkRAD/${indiv}_cloneF_1.fq `ls $OUT_DIR/$run/clone_filter/${indiv}_cloneF_2.fq* `" >> $SGE/recoverMate.qarray
else
echo "sh $SCRIPT_DIR/recover_mate.sh `ls $OUT_DIR/$run/checkRAD/${indiv}_1.fq*discards` `ls $OUT_DIR/$run/${indiv}_2.fq*`" >> $SGE/recoverMateDiscards.qarray
echo "sh $SCRIPT_DIR/recover_mate.sh $OUT_DIR/$run/checkRAD/${indiv}_1.fq `ls $OUT_DIR/$run/${indiv}_2.fq* `" >> $SGE/recoverMate.qarray
fi
done

job_id_recov_dis=`qarray -terse -N ${Proj_Name}_recoverMateDiscards -o $SGE -e $SGE $SGE/recoverMateDiscards.qarray`
job_id_recov_ok=`qarray -terse -N ${Proj_Name}_recoverMate -o $SGE -e $SGE $SGE/recoverMate.qarray`

r=0;
while [ $r == 0 ] ; do qstat -j $job_id_recov_dis >& /dev/null ; r=$? ; sleep 10s ; done
r=0;
while [ $r == 0 ] ; do qstat -j $job_id_recov_ok >& /dev/null ; r=$? ; sleep 10s ; done

else
cat $INDIV_FILE | while read line
do
run=`echo $line | awk '{print $3}'`
indiv=`echo $line | awk '{print $2}'`
echo "gzip $OUT_DIR/$run/checkRAD/${indiv}_*" >> $SGE/gzip_process_radtag.qarray
done

job_id_gzip=`qarray -terse -N ${Proj_Name}_gzip -o $SGE -e $SGE $SGE/gzip_process_radtag.qarray`

r=0;
while [ $r == 0 ] ; do qstat -j $job_id_gzip >& /dev/null ; r=$? ; sleep 10s ; done
fi

############################ CREATION DOSSIER INPUT u/pstacks #############################################
mkdir -p $OUT_DIR/stacks_input

cat $INDIV_FILE | while read line
do
run=`echo $line | awk '{print $3}'`
indiv=`echo $line | awk '{print $2}'`

if [[ $align == 1 ]]
then
if [[ $dereplication == 0 ]]
then
ln -s $OUT_DIR/$run/checkRAD/${indiv}_*.fq.gz $OUT_DIR/stacks_input/
else
for i in `ls $OUT_DIR/$run/checkRAD/${indiv}_cloneF_*.fq.gz`
do
s=`basename $i | sed 's/_cloneF//'`
ln -s $i $OUT_DIR/stacks_input/$s
done 
fi
else
if [[ $double_digest == 1 ]]
then
zcat $OUT_DIR/$run/checkRAD/${indiv}_*.fq.gz > $OUT_DIR/stacks_input/${indiv}.fq
gzip $OUT_DIR/stacks_input/${indiv}.fq
else
if [[ $dereplication == 1 ]]
then
ln -s $OUT_DIR/$run/checkRAD/${indiv}_cloneF_1.fq.gz $OUT_DIR/stacks_input/${indiv}.fq.gz
else
ln -s $OUT_DIR/$run/checkRAD/${indiv}_1.fq.gz $OUT_DIR/stacks_input/${indiv}.fq.gz
fi
fi 
fi
done

############################# STAT #############################
date
echo -e "\nSTATISTICS\n"
echo -e "\nSTATISTICS\n" >&2

echo -e "#run\tpop\tindiv\tdemultiplexed_pairs\tparis_duplicated_removed\tlowqual_read_removed\tbadTAG_read\tfinal_retained_reads" > $STAT_DIR/read_count.txt
for run in $Run
do

	idx=1
	type=`ls $OUT_DIR/$run/*_1.fq*|head -n 1  | awk '{if(match($1,".gz") > 0 ){print "gzfastq"}else{print "fastq"}}'`
	grep $run $INDIV_FILE | while read line
	do
		indiv=`echo $line | cut -d ' ' -f 2`
		pop=`awk -v I=$indiv '{if($1==I){print $2}}' $POP_FILE`
		if [[ "$type" = "gzfastq" ]]
		then
			nb_read=`zcat $OUT_DIR/$run/${indiv}_1* | wc -l |awk '{print $1/4}'`
		else
			nb_read=`wc -l $OUT_DIR/$run/${indiv}_1* | awk '{print $1/4}' `
		fi
		
		if [[ $dereplication == 1 ]]
		then 
			nb_dup=`tail -n 1 $SGE/${Proj_Name}_clone_filter_$run.e*.$idx | awk '{print $12}'`
		else
			nb_dup=0
		fi
		
		nb_lowqual=`tail -n 3 $SGE/${Proj_Name}_process_rad_tags_$run.e*.$idx | grep "low quality" | awk '{print $1}'`
		nb_badRAD=`tail -n 6 $SGE/${Proj_Name}_process_rad_tags_$run.e*.$idx | grep 'ambiguous RAD-Tag drops' | awk '{print $1}'`
		nb_retain=`tail -n 3 $SGE/${Proj_Name}_process_rad_tags_$run.e*.$idx | grep "retained" | awk '{print $1}'`
		
		echo -e "$run\t$pop\t$indiv\t$nb_read\t$nb_dup\t$nb_lowqual\t$nb_badRAD\t$nb_retain" >> $STAT_DIR/read_count.txt
		let idx+=1 
	done
	
	Tdemult=`awk -v R=$run '$1==R' $STAT_DIR/read_count.txt | colstat.sh 4 | awk -v D=$double_digest '{if(D==1){print $4*2}else{print $4}}'`
	Tdup=`awk -v R=$run '$1==R' $STAT_DIR/read_count.txt| colstat.sh 5 | awk -v D=$double_digest '{if(D==1){print $4*2}else{print $4}}'`
	Tlowqual=`awk -v R=$run '$1==R' $STAT_DIR/read_count.txt | colstat.sh 6 | awk '{print $4}'`
	TbadRAD=`awk -v R=$run '$1==R' $STAT_DIR/read_count.txt | colstat.sh 7 | awk '{print $4}'`
	Tretained=`awk -v R=$run '$1==R' $STAT_DIR/read_count.txt | colstat.sh 8 | awk '{print $4}'`
	
	if [[ "$type" = "gzfastq" && -e $OUT_DIR/$run/unmatched_1.fq.gz ]]
	then 
		unmap=`zcat $OUT_DIR/$run/unmatched_1.fq.gz | wc -l | awk '{print $1/4}'`
	else
		if [[ "$type" = "fastq" && -e $OUT_DIR/$run/unmatched_1.fq ]]
		then
			unmap=`wc -l $OUT_DIR/$run/unmatched_1.fq | awk '{print $1/4}'`
		else
			unmap=0
		fi
	fi
	
	if [[ "$type" = "gzfastq" && -e $OUT_DIR/$run/ambiguous_1.fq.gz ]]
	then 
		ambiguous=`zcat $OUT_DIR/$run/ambiguous_1.fq.gz | wc -l | awk '{print $1/4}'`
	else
		if [[ "$type" = "fastq" && -e $OUT_DIR/$run/ambiguous_1.fq ]]
		then
			ambiguous=`wc -l $OUT_DIR/$run/ambiguous_1.fq | awk '{print $1/4}'`
		else
			ambiguous=0
		fi
	fi
	
	if [[ "$type" = "gzfastq" && $(ls $OUT_DIR/$run/ | grep 2rad_1.fq.gz | wc -l ) -gt 0 ]]
	then 
		dbrad=`zcat $OUT_DIR/$run/*2rad_1.fq.gz | wc -l | awk '{print $1/4}'`
	else
		if [[ "$type" = "fastq" && $(ls $OUT_DIR/$run/ | grep 2rad_1.fq.gz | wc -l ) -gt 0 ]]
		then
			dbrad=`wc -l $OUT_DIR/$run/*2rad_1.fq | grep total | awk '{print $1/4}'`
		else
			dbrad=0
		fi
	fi
	
	let Ttot=$Tdemult+$unmap+$dbrad+$ambiguous
	let per_demult=$Tdemult*100/$Ttot
	if [[ $dereplication == 1 ]]
	then
		let per_cloneF=$Tdup*100/$Ttot
	fi
	let per_lowQual=$Tlowqual*100/$Ttot
	let per_badRAD=$TbadRAD*100/$Ttot
	let per_qualFil=$Tretained*100/$Ttot
	let per_unmap=$unmap*100/$Ttot
	let per_ambiguous=$ambiguous*100/$Ttot
	let per_2rad=$dbrad*100/$Ttot
	
	if [[ $double_digest == 1 ]]
	then
		echo ""
		echo "$run :"
		echo "    total reads : $Ttot"
		echo "    wrong barcode : $unmap ($per_unmap %)"
		echo "    ambiguous barcode : $ambiguous ($per_ambiguous %)"
		echo "    2 rad tag : $dbrad ($per_2rad %)"
		echo "    barcoded : $Tdemult ($per_demult %)"
		echo "        low quality read1: $Tlowqual ($per_lowQual %)"
		echo "        bad RADTag: $TbadRAD ($per_badRAD %)"
		echo  "        quality filtered and RAD checked : $Tretained ($per_qualFil %)" 
	else
		echo ""
		echo "$run :"
		echo "    total pairs : $Ttot"
		echo "    wrong barcode : $unmap ($per_unmap %)"
		echo "    ambiguous barcode : $ambiguous ($per_ambiguous %)"
		echo "    2 rad tag : $dbrad ($per_2rad %)"
		echo "    barcoded : $Tdemult ($per_demult %)"
		if [[ $dereplication == 1 ]]
		then
			echo "        clone filtered paired: $Tdup ($per_cloneF %)"
		fi
		echo "        low quality read1: $Tlowqual ($per_lowQual %)"
		echo "        bad RADTag: $TbadRAD ($per_badRAD %)"
		if [[ $dereplication == 1 ]]
		then
			echo "        dereplicated, quality filtered and RAD checked : $Tretained ($per_qualFil %)" 
		else
			echo  "        quality filtered and RAD checked : $Tretained ($per_qualFil %)" 
		fi
	fi
done
echo -e "\nEND OF PREPROCESSING STEP"

date

