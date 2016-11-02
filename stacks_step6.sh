#$ -o [Project_path_dir]/SGE_out/stacks6.out
#$ -e [Project_path_dir]/SGE_out/stacks6.err
#$ -N [Proj_Name]_stacks6
#$ -S /bin/sh

#$ -q unlimitq

#$ -m bea

date
date >&2
Proj_Name=

# VARIABLES DEDIEES
# site de restriction tronqué	==> à adapter pour la double digestion. Impact dans assembly_selection.py et dans la selection des assemblage de bonne qualité
ENZ=""
# SNP filter
QUAL=		# nb indiv * 50 (Qualité moyenne > 50 (sur une echelle de 0 à 100)
AN=			# 50% des individus génotypés = nombre d'individus
# fichier qui contient la référence fasta. Ne rien mettre si vous n'avez pas de référence
GENOME=

# VARIABLES GENERALES
RAD_DIR=[Project_path_dir]
# le chemin complet du dossier qui contient les scripts
SCRIPT_DIR=[Script_path_dir]

# INPUT
# le chemin complet du dossier qui contient tous les fichiers fastq utilisé pour l'analyse stacks: exemple de contenu du dossier: ATL_HARDY_M_1.1.fq ATL_HARDY_M_1.2.fq
SAMPLE_READS_DIR=
# gzfastq or fastq
READS_TYPE="gzfastq"
# le chemin complet du dossier qui contient tous les fichiers de résultats Stacks : example contenu à minima du dossier: ( 4 fichiers par individus) ATL_HARDY_M_1.1.alleles.tsv ATL_HARDY_M_1.1.matches.tsv ATL_HARDY_M_1.1.snps.tsv ATL_HARDY_M_1.1.tags.tsv (et 3 fichiers pour le catalogue) batch_2.catalog.alleles.tsv batch_2.catalog.snps.tsv batch_2.catalog.tags.tsv 
STACKS_RES_DIR=
# le chemin complet du fichier tabulé qui décrit les échantillons : 1 ligne par individus avec colonne 1 un index (1,2,...,N) et colonne 2 un échantillon, exemple de contenu du fichier
#1	ATL_HARDY_M_1
#2	ATL_HARDY_M_2 
INDIV_FILE=
dos2unix $INDIV_FILE
# le chemin complet du fichier qui contient les locus à assembler (1 locus par ligne)
#1
#21
#70
WHITE_LOCUS_LIST=
dos2unix $WHITE_LOCUS_LIST
#################################################### NE PLUS RIEN MODIFIER ENSUITE #############################################################
# OUTPUT
OUT_DIR=$RAD_DIR
SGE=$RAD_DIR/SGE_out
mkdir -p $SGE $OUT_DIR $OUT_DIR/tmp

PAIRS_DIR=$OUT_DIR/locus_pairs
SGE_PAIRS=$SGE/sort_read

SGE_CAP3=$SGE/cap3
ASSEMBLY_DIR=$OUT_DIR/assemblage_cap3

SGE_LOC=$SGE/Genome_localisation
LOC_DIR=$OUT_DIR/Genome_localisation

SGE_SNP=$SGE/GATK_Stacks_compare
SNP_DIR=$OUT_DIR/GATK_Stacks_compare
###################################################################################################
## Reconstitution des fichiers fastq par locus
mkdir -p $PAIRS_DIR $SGE_PAIRS

if [[ -e $SGE_PAIRS/sort_read_pair_aaa.array ]]
then
rm $SGE_PAIRS/sort_read_pair*.array
fi

# split la liste de locus pour une utilisation du cluster politsiquement correcte!
split -l 200 -a 3 $WHITE_LOCUS_LIST $OUT_DIR/tmp/`basename $WHITE_LOCUS_LIST`

cut -f 2 $INDIV_FILE > $OUT_DIR/tmp/white_indiv_list
for locus_list in `ls $OUT_DIR/tmp/$(basename $WHITE_LOCUS_LIST)* `
do
	suffix=`basename $locus_list | awk '{print substr($1,length($1)-2,3)}'`
	echo "perl "$SCRIPT_DIR"/sort_read_pairs_MB.pl -p "$STACKS_RES_DIR" -s "$SAMPLE_READS_DIR" -o "$PAIRS_DIR" -f "$READS_TYPE" -t fastq -W " $OUT_DIR/tmp/white_indiv_list " -w "$locus_list "-r pair" > $SGE_PAIRS/sort_read_pair_$suffix.sh
	
	id_sort_read=`qsub -q unlimitq -terse -l h_vmem=60G -l mem=30G -N ${Proj_Name}_sort_reads_$suffix -o $SGE_PAIRS/ -e $SGE_PAIRS/ $SGE_PAIRS/sort_read_pair_$suffix.sh`
	r=0;
	while [ $r == 0 ] ; do qstat -j $id_sort_read >& /dev/null ; r=$? ; sleep 10s ; done 

	# problème de locus?
	cat $locus_list| while read line
	do
		r1=`wc -l $PAIRS_DIR/${line}_1.fq | awk '{print $1/4}'`
		r2=`wc -l $PAIRS_DIR/${line}_2.fq | awk '{print $1/4}'`
		if [ "$r1" != "$r2" ]
		then
		echo -e "$line\t$r1\t$r2" >> $OUT_DIR/locus_pair.problem
		else
			if [[ "$r1" == "0" ]]
			then
				echo -e "$line\t0" >> $OUT_DIR/locus_pair.problem
			else
				echo -e "$line\t$r1" >> $OUT_DIR/count_read.locus
			fi
		fi
	done

	if [[ -e $OUT_DIR/locus_pair.problem ]]
	then
		awk -v BAD=$OUT_DIR/locus_pair.problem 'BEGIN{while(getline<BAD>0){tab[$1]=1}}{if(tab[$1]!=1){print $1}}' $locus_list >> $OUT_DIR/white_locus_sort_read_ok
	fi

done

if [[ -e $OUT_DIR/white_locus_sort_read_ok ]]
then
WHITE_LOCUS_LIST=$OUT_DIR/white_locus_sort_read_ok
fi

if [[ ! -s $WHITE_LOCUS_LIST ]]
then
echo "sort_read failed for all locus"
fi

##############################################################################################################
########################### 				assemblage cap3 des locus 				       ###################
##############################################################################################################
echo "##############################################################################################################
########################### 				assemblage cap3 des locus 				       ###################"
mkdir -p $SGE_CAP3 $ASSEMBLY_DIR

if [ -e $SGE_CAP3/assembly_aaa.array ]
then 
rm $SGE_CAP3/assembly_*.array
fi

for locus_list in `ls $OUT_DIR/tmp/$(basename $WHITE_LOCUS_LIST)*`
do
	suffix=`basename $locus_list | awk '{print substr($1,length($1)-2,3)}'`
	cat  $locus_list | while read locus
	do
		mkdir -p $ASSEMBLY_DIR/detail/$locus
		echo "sh "$SCRIPT_DIR"/cap3.sh "$ASSEMBLY_DIR"/detail/"$locus" "$PAIRS_DIR"/"${locus}"_1.fq "$PAIRS_DIR"/"${locus}"_2.fq "$locus >> $SGE_CAP3"/assembly_"$suffix".array"
	done
	id_assembly=`qarray -terse -N ${Proj_Name}_assembly_$suffix -o $SGE_CAP3/ -e $SGE_CAP3/ $SGE_CAP3/assembly_$suffix.array`
	r=0;
	while [ $r == 0 ] ; do qstat -j $id_assembly >& /dev/null ; r=$? ; sleep 10s ; done 
done

### STAT ASSEMBLAGE ####

mkdir $ASSEMBLY_DIR/stat
job_id_assembly_stat=`qsub -terse -hold_jid ${Proj_Name}_assembly -N ${Proj_Name}_valid_assembly -o $SGE_CAP3/assembly.stat.out -e $SGE_CAP3/assembly.stat.err $SCRIPT_DIR/assembly_selection.py -e $ENZ -w $WHITE_LOCUS_LIST -a $ASSEMBLY_DIR/detail -o $ASSEMBLY_DIR/stat -l $ASSEMBLY_DIR/stat/assembly_selection.log`

r=0;
while [ $r == 0 ]; do qstat -j $job_id_assembly_stat >& /dev/null; r=$? ;done 

## assemblage 1 contigs > 2 contigs
cut -f 1 $ASSEMBLY_DIR/stat/assembly.stat | grep -v "#" | sort | uniq -c | awk -v DIR=$ASSEMBLY_DIR/stat/ -v OUT=$ASSEMBLY_DIR/stat/assembly_nbcontig.list '{print $1,$2 >> OUT; if ($1==1){n++}else{m++} }END{print n" single contig, "m" >=2 contigs"}'

## sous selection des stat en fonction du nombre de contig
awk -v DIR=$ASSEMBLY_DIR/stat -v COUNT=$ASSEMBLY_DIR/stat/assembly_nbcontig.list 'BEGIN{while(getline<COUNT>0){tab[$2]=$1}}{if($1=="#locus"){print $0 > DIR"/assembly_1contig.stat"; print $0 > DIR"/assembly_sup2contig.stat"; }else{if(tab[$1]==1){print $0 >> DIR"/assembly_1contig.stat";}else{print $0 >> DIR"/assembly_sup2contig.stat";}}}' $ASSEMBLY_DIR/stat/assembly.stat

echo "# stat assembly 1 contig"
## stat assembly 1 contig
echo "# R1 et R2"
## R1 et R2
grep -v "#" $ASSEMBLY_DIR/stat/assembly_1contig.stat | awk '$9>0 && $10>0' | wc -l
echo "# enzyme ok"
## enzyme ok
grep -v "#" $ASSEMBLY_DIR/stat/assembly_1contig.stat | awk '$6>-1' | wc -l
echo "# enzyme abs"
## enzyme abs
grep -v "#" $ASSEMBLY_DIR/stat/assembly_1contig.stat | awk '$6==-1' | cut -f 1

echo "# stat sup 2 contig"
## stat sup 2 contig
echo "# max len + RE"
## max len + RE
grep -v "#" $ASSEMBLY_DIR/stat/assembly_sup2contig.stat | awk '$3==$7 && $6>-1' | wc -l
echo "# max len +RE + R1 + R2"
## max len +RE + R1 + R2
grep -v "#" $ASSEMBLY_DIR/stat/assembly_sup2contig.stat | awk '$3==$7 && $6>-1 && $9>0 && $10>0' | wc -l
echo "# max len + RE ! R1 + R2"
## max len + RE ! R1 + R2
grep -v "#" $ASSEMBLY_DIR/stat/assembly_sup2contig.stat | awk '$3==$7 && $6>-1 && ($9==0 || $10==0) '

echo "# max read + max len "
## max read + max len 
grep -v "#" $ASSEMBLY_DIR/stat/assembly_sup2contig.stat | awk '$3==$7 && $4==$8 '  | wc -l
echo "# max read + max len + RE"
## max read + max len + RE
grep -v "#" $ASSEMBLY_DIR/stat/assembly_sup2contig.stat | awk '$3==$7 && $4==$8 && $6>-1'  | wc -l
echo "# max read + max len + RE + R1 + R2"
## max read + max len + RE + R1 + R2
grep -v "#" $ASSEMBLY_DIR/stat/assembly_sup2contig.stat | awk '$3==$7 && $4==$8 && $6>-1 && $9>0 && $10>0'  | wc -l
echo "# max read | max len + RE"
## max read | max len + RE
grep -v "#" $ASSEMBLY_DIR/stat/assembly_sup2contig.stat | awk '$3==$7 && $4!=$8 && $6 > -1'  | wc -l
grep -v "#" $ASSEMBLY_DIR/stat/assembly_sup2contig.stat | awk '$3!=$7 && $4==$8 && $6 > -1'  | wc -l

echo "# max read | max len + RE"
## max read | max len + RE
grep -v "#" $ASSEMBLY_DIR/stat/assembly_sup2contig.stat | awk '$3==$7 && $4!=$8 && $6 > -1 && $9>0 && $10>0'  | wc -l
grep -v "#" $ASSEMBLY_DIR/stat/assembly_sup2contig.stat | awk '$3!=$7 && $4==$8 && $6 > -1 && $9>0 && $10>0'  | wc -l
#
## assembly ok  : 1 contig +RE +R1 +R2 && sup2contig +max len +max read +RE +R1 +R2 && surp2contig+ max len | max read + RE +R1 +R2
echo "# assembly ok  : 1 contig +RE +R1 +R2 && sup2contig +max len +max read +RE +R1 +R2 && surp2contig+ max len | max read + RE +R1 +R2"
grep "#" $ASSEMBLY_DIR/stat/assembly.stat > $ASSEMBLY_DIR/stat/assembly_ok.stat
head -n 1 $ASSEMBLY_DIR/stat/assembly_1contig.stat > $ASSEMBLY_DIR/stat/assembly_ok.stat
grep -v "#" $ASSEMBLY_DIR/stat/assembly_1contig.stat | awk '$6>-1 && $9>0 && $10>0' >> $ASSEMBLY_DIR/stat/assembly_ok.stat
grep -v "#" $ASSEMBLY_DIR/stat/assembly_sup2contig.stat | awk '$3==$7 && $4==$8 && $6>-1 && $9>0 && $10>0' >> $ASSEMBLY_DIR/stat/assembly_ok.stat
grep -v "#" $ASSEMBLY_DIR/stat/assembly_sup2contig.stat | awk '$3==$7 && $4!=$8 && $6 > -1 && $9>0 && $10>0' >> $ASSEMBLY_DIR/stat/assembly_ok.stat
grep -v "#" $ASSEMBLY_DIR/stat/assembly_sup2contig.stat | awk '$3!=$7 && $4==$8 && $6 > -1 && $9>0 && $10>0' >> $ASSEMBLY_DIR/stat/assembly_ok.stat
grep -v "#" $ASSEMBLY_DIR/stat/assembly_ok.stat | wc -l
# certain locus sont ils détéctés plusieurs fois via différents contig
echo " locus détéctés avec plusieurs \"bon\" assemblage ?"
cut -f 1 $ASSEMBLY_DIR/stat/assembly_ok.stat | sort | uniq -c | awk -v STATOK=$ASSEMBLY_DIR/stat/assembly_ok.stat  '{if($1>1){tab[$2]=1}}END{while(getline<STATOK>0){if($1=="#locus" || tab[$1]==1){print $0}}}'

### MERGE ASSEMBLAGE ### ==> tous les contigs de tous les locus assemblés!! 
mkdir $ASSEMBLY_DIR/merge
for i in ` cut -f 1 $WHITE_LOCUS_LIST  `
do 
cat $ASSEMBLY_DIR/detail/$i/$i.contig.fa >> $ASSEMBLY_DIR/merge/${Proj_Name}_cap3_assembly.fa
done

#~ 
#~ ######################################################################################################################
#~ ########################### 		localisation assemblage sur le génome de ref			##########################
#~ 
echo "######################################################################################################################
########################### 		localisation assemblage sur le génome de ref			##########################"

if [[ -e $GENOME ]]
then

	mkdir -p $LOC_DIR $SGE_LOC

	if [[ -e $SGE_LOC/localisation_aaa.array ]]
	then
	rm $SGE_LOC/localisation*.array
	fi

	if [[ ! -e $GENOME.amb ]]
	then
	ln -s $GENOME $LOC_DIR/.
	GENOME=`basename $GENOME | awk -v DIR=$LOC_DIR '{print DIR"/"$1}'`
	bwa index $GENOME
	fi

	for locus_list in `ls $OUT_DIR/tmp/$(basename $WHITE_LOCUS_LIST)*`
	do
		suffix=`basename $locus_list | awk '{print substr($1,length($1)-2,3)}'`
		cat  $locus_list | while read locus
		do
		CONTIG=$ASSEMBLY_DIR/detail/$locus/$locus.contig.fa
		mkdir -p $LOC_DIR/detail/$locus
		rm $LOC_DIR/detail/$locus/*
		echo sh $SCRIPT_DIR/localisation.sh $GENOME $PAIRS_DIR $CONTIG $LOC_DIR/detail $locus $SCRIPT_DIR $INDIV_FILE >> $SGE_LOC/localisation_$suffix.array
		done
		id_loc=`qarray -terse -N ${Proj_Name}_localisation_$suffix -e $SGE_LOC/ -o $SGE_LOC/ $SGE_LOC/localisation_$suffix.array`
		r=0;
		while [ $r == 0 ] ; do qstat -j $id_loc >& /dev/null ; r=$? ; sleep 10s ; done 
		
	done

	### MERGE BAM ###
	mkdir -p $LOC_DIR/merge
	 
	locus1=` head -n 1 $WHITE_LOCUS_LIST`
	samtools view -H $LOC_DIR/detail/$locus1/${locus1}_pe_RG.bam > $LOC_DIR/merge/${Proj_Name}_genome_read_mapping_RG.sam
	samtools view -H $LOC_DIR/detail/$locus1/${locus1}_contig.bam > $LOC_DIR/merge/${Proj_Name}_genome_contig_mapping.sam

	for i in ` cut -f 1 $WHITE_LOCUS_LIST `
	do
	samtools view $LOC_DIR/detail/$i/${i}_pe_RG.bam >> $LOC_DIR/merge/${Proj_Name}_genome_read_mapping_RG.sam
	samtools view $LOC_DIR/detail/$i/${i}_contig.bam >> $LOC_DIR/merge/${Proj_Name}_genome_contig_mapping.sam
	done

	samtools view -h -S $LOC_DIR/merge/${Proj_Name}_genome_read_mapping_RG.sam -b -o $LOC_DIR/merge/${Proj_Name}_genome_read_mapping_RG.bam
	samtools sort $LOC_DIR/merge/${Proj_Name}_genome_read_mapping_RG.bam $LOC_DIR/merge/${Proj_Name}_genome_read_mapping_RG_sort
	samtools index $LOC_DIR/merge/${Proj_Name}_genome_read_mapping_RG_sort.bam
	rm $LOC_DIR/merge/${Proj_Name}_genome_read_mapping*sam $LOC_DIR/merge/${Proj_Name}_genome_read_mapping_RG.bam

	samtools view -h -S $LOC_DIR/merge/${Proj_Name}_genome_contig_mapping.sam -b -o $LOC_DIR/merge/${Proj_Name}_genome_contig_mapping.bam
	samtools sort $LOC_DIR/merge/${Proj_Name}_genome_contig_mapping.bam $LOC_DIR/merge/${Proj_Name}_genome_contig_mapping_sort
	samtools index $LOC_DIR/merge/${Proj_Name}_genome_contig_mapping_sort.bam
	rm $LOC_DIR/merge/${Proj_Name}_genome_contig_mapping*sam $LOC_DIR/merge/${Proj_Name}_genome_contig_mapping.bam

	### MERGE BED ###

	if [[ -e  $LOC_DIR/merge/${Proj_Name}_cap3_mapping.bed ]]
	then
	rm $LOC_DIR/merge/${Proj_Name}_cap3_mapping.bed
	fi

	cat $WHITE_LOCUS_LIST | while read locus
	do
	cat $LOC_DIR/detail/$locus/$locus.intersect.bed >> $LOC_DIR/merge/${Proj_Name}_cap3_mapping.bed
	done

	# stat BED
	echo -e "locus\tnb_total_localisation\tnb_localisation_readPe\tnb_localisation_contig\tnb_common_localisation\tlocalisations\tnb_localisation_readPe_contigok\tvalid_loc" > $LOC_DIR/merge/${Proj_Name}_cap3_mapping.bed.stat

	for i in `cat $WHITE_LOCUS_LIST`
	do
	exp1=`echo Locus_${i}_`
	nb_loc=`grep -v NULL $LOC_DIR/merge/${Proj_Name}_cap3_mapping.bed | grep -c $exp1  `
	exp2=`echo Locus_${i}_Contig`
	nb_loc_CONTIG=`grep -v NULL $LOC_DIR/merge/${Proj_Name}_cap3_mapping.bed | grep -c $exp2 ` 
	exp3=`echo Locus_${i}_readPe`
	nb_loc_RPE=`grep -v NULL $LOC_DIR/merge/${Proj_Name}_cap3_mapping.bed | grep -c $exp3 ` 
	nb_loc_RPE_AND_CONTIG=`grep -v NULL $LOC_DIR/merge/${Proj_Name}_cap3_mapping.bed | grep $exp2 | grep -c $exp3` 
	loc=`grep -v NULL $LOC_DIR/merge/${Proj_Name}_cap3_mapping.bed | grep $exp2 | grep $exp3 | awk '{loc=loc""$1":"$2"-"$3"("$NF");"}END{if (loc==""){loc= "no_loc;"} print loc}'` 
	nb_loc_RPE_CONTIG_ASSEMBLY_OK=-1
	valid_loc="no_valid_loc;"
	if [[ $(awk -v L=$i '$1==L' $ASSEMBLY_DIR/stat/assembly_ok.stat | wc -l) -gt 0 ]]
	then
	nb_loc_RPE_CONTIG_ASSEMBLY_OK=`grep $exp2 $LOC_DIR/merge/${Proj_Name}_cap3_mapping.bed | grep $exp3 |awk -v L=$i -v STATOK=$ASSEMBLY_DIR/stat/assembly_ok.stat 'BEGIN{while(getline<STATOK>0){if($1==L){tab[$5]=1  }}}{l=split($4,tmp,",") ; for (i=1;i<=l;i++){if(tab[tmp[i]]==1){print $0}}}' | sort -u | wc -l `
	valid_loc=`grep $exp2 $LOC_DIR/merge/${Proj_Name}_cap3_mapping.bed | grep $exp3 |awk -v L=$i -v STATOK=$ASSEMBLY_DIR/stat/assembly_ok.stat 'BEGIN{while(getline<STATOK>0){if($1==L){tab[$5]=1  }}}{l=split($4,tmp,",") ; for (i=1;i<=l;i++){if(tab[tmp[i]]==1){loc=loc""$1":"$2"-"$3"("$NF");"}}}END{if (loc==""){loc= "no_valid_loc;"} print loc}' `
	fi
	echo -e "$i\t$nb_loc\t$nb_loc_RPE\t$nb_loc_CONTIG\t$nb_loc_RPE_AND_CONTIG\t$loc\t$nb_loc_RPE_CONTIG_ASSEMBLY_OK\t$valid_loc" >> $LOC_DIR/merge/${Proj_Name}_cap3_mapping.bed.stat
	done

	grep -v "locus" $LOC_DIR/merge/${Proj_Name}_cap3_mapping.bed.stat | awk 'BEGIN{a=0; r1=0; r1z=0; r1m=0; c=0; cz=0; cm=0; rc=0; rcz=0; rcm=0; r1cok=0; r1cokz=0; common=0}{
	if($2==1){
		a++
	}if ($3==1){
		r1++
	} if($3==0){
		r1z++
	} if($3>1){
		r1m++
	} if($4==1){
		c++
	} if ($4==0){
		cz++
	} if ($4>1){
		cm++
	} if($5==1){
		rc++
	} if ($5==0){
		rcz++
	} if ($5>1){
		rcm++
	} if ($7==0){
		r1cokz++
	} if ($7==1){
		r1cok++
	} if($7>1){
		r1cokm++
	} if($7==-1 && $5==1){
		common++
	}}END{
	print "nombre de locus à 1 seule localisation: "a"\n";
	print "nombre de locus à 0/1/>1 localisation pour les lectures: "r1z" / "r1" / "r1m"\n";
	print "nombre de locus à 0/1/>1 localisation pour les mini-contigs: "cz" / "c" / "cm"\n";
	print "nombre de locus à 0/1/>1 localisation communes entre read et mini-contigs: "rcz" / "rc" / "rcm"\n";
	print "nombre de locus à 0 / 1 / >1 localisation commune entre read et minicontig ok (maxread/maxlen RE R1 R2): "r1cokz" / "r1cok" / "r1cokm"\n" ;
	print "nombre de locus à 1 localisation commune entre read et minicontig avec assemblage de mauvaise qualite: "common"\n"; }' 
else
echo "localisation impossible, car aucun genome de reference specifié"
fi
# fin de test de presence du genome de ref

######################################################################################################################
########################### 		detection des variants GATK + comparaison stacks			######################
echo "######################################################################################################################
########################### 		detection des variants GATK + comparaison stacks			######################"
mkdir -p $SNP_DIR $SGE_SNP

if [[ -e $SGE_SNP/gatk_aaa.array ]]
then
rm $SGE_SNP/*
fi

for locus_list in `ls $OUT_DIR/tmp/$(basename $WHITE_LOCUS_LIST)*`
do
	suffix=`basename $locus_list | awk '{print substr($1,length($1)-2,3)}'`
	cat  $locus_list | while read locus
	do
	if [[ -e $GENOME ]]
	then
		echo sh $SCRIPT_DIR/gatk_snp_calling_ongenome.sh $locus $LOC_DIR/detail/$locus/${locus}_pe_RG.bam $GENOME $QUAL $AN `ls $STACKS_RES_DIR/batch*catalog.snps.tsv*` $SNP_DIR/detail $SCRIPT_DIR >>  $SGE_SNP/gatk_$suffix.array
	else
		REF=$ASSEMBLY_DIR/detail/$locus/$locus.contig.fa
		echo sh $SCRIPT_DIR/gatk_snp_calling_onassembly.sh $locus $PAIRS_DIR/ $REF $INDIV_FILE $QUAL $AN `ls $STACKS_RES_DIR/batch*catalog.snps.tsv*` $SNP_DIR/detail $SCRIPT_DIR >>  $SGE_SNP/gatk_$suffix.array
	fi
	
	mkdir -p $SNP_DIR/detail/$locus
	rm $SNP_DIR/detail/$locus/*
	
	done
	id_gatk=`qarray -terse -l h_vmem=20G -l mem=10G -N ${Proj_Name}_gatk_$suffix -e $SGE_SNP/ -o $SGE_SNP/ $SGE_SNP/gatk_$suffix.array`
	r=0;
	while [ $r == 0 ] ; do qstat -j $id_gatk >& /dev/null ; r=$? ; sleep 10s ; done 
	
done

### MERGE GATK ###
mkdir -p $SNP_DIR/merge
rm $SNP_DIR/merge/*

cat $WHITE_LOCUS_LIST | while read locus
do
if [[ -e $SNP_DIR/detail/$locus/$locus.tsv ]]
then
	if [[ ! -e $SNP_DIR/merge/${Proj_Name}_genome_gatk.tsv ]]
	then
	cat $SNP_DIR/detail/$locus/$locus.tsv > $SNP_DIR/merge/${Proj_Name}_genome_gatk.tsv
	else
	grep -v locus $SNP_DIR/detail/$locus/$locus.tsv >> $SNP_DIR/merge/${Proj_Name}_genome_gatk.tsv
	fi
else
echo -e $locus"\tno_snp_found" >> $SNP_DIR/merge/${Proj_Name}_genome_gatk.tsv
fi

if [[ ! -e $SNP_DIR/merge/${Proj_Name}_cap3_gatk.vcf.gz ]]
then
	if [[ -e $SNP_DIR/detail/$locus/${locus}_filtered.vcf.gz ]]
	then
	cp $SNP_DIR/detail/$locus/${locus}_filtered.vcf.gz $SNP_DIR/merge/${Proj_Name}_genome_gatk.vcf.gz
	fi
else
	if [[ -e $SNP_DIR/detail/$locus/${locus}_filtered.vcf.gz ]]
	then
	vcf-concat $SNP_DIR/merge/${Proj_Name}_cap3_gatk.vcf.gz $SNP_DIR/detail/$locus/${locus}_filtered.vcf.gz > $SNP_DIR/merge/tmp
	mv $SNP_DIR/merge/tmp $SNP_DIR/merge/${Proj_Name}_genome_gatk.vcf
	rm $SNP_DIR/merge/${Proj_Name}_genome_gatk.vcf.gz
	bgzip $SNP_DIR/merge/${Proj_Name}_genome_gatk.vcf
	fi
fi
done

tabix -p vcf $SNP_DIR/merge/${Proj_Name}_genome_gatk.vcf.gz

# Nombre de locus où aucun SNP n'est détecté puisque l'alignement n'est pas de suffisament bonne qualité = absence de vcf
if [[ -e $OUT_DIR/no_vcf ]]
then
rm $OUT_DIR/no_vcf
fi

cat $WHITE_LOCUS_LIST | while read locus
do
if [[ ! -e $SNP_DIR/detail/$locus/${locus}_filtered.vcf.gz ]]
then
echo $locus >> $OUT_DIR/no_vcf
fi
done

if [[ -e $OUT_DIR/no_vcf ]]
then
wc -l $OUT_DIR/no_vcf | awk '{print "Nombre de locus où aucun SNP n est détecté puisque l alignement n est pas de suffisament bonne qualité = absence de vcf : ", $1}'
else
echo "Nombre de locus où aucun SNP n est détecté puisque l alignement n est pas de suffisament bonne qualité = absence de vcf : 0"
fi

# Nombre de locus où aucun SNP n'est détecté avec GATK = vcf vide
if [[ -e $OUT_DIR/empty_vcf ]]
then
rm $OUT_DIR/empty_vcf
fi

cat $WHITE_LOCUS_LIST | while read locus
do
if [[ -e $SNP_DIR/detail/$locus/${locus}_filtered.vcf.gz ]]
then 
zcat $SNP_DIR/detail/$locus/${locus}_filtered.vcf.gz | grep -v "#" | wc -l | awk -v L=$locus '{if ($1==0){print L}}' >> $OUT_DIR/empty_vcf
fi
done

wc -l $OUT_DIR/empty_vcf | awk '{print "Nombre de locus où aucun SNP n est détecté avec GATK = vcf vide: " $1}' 

# Nombre de locus où aucun SNP détectés avec GATK ne passe les filtres de qualité de détection = no PASS
if [[ -e $OUT_DIR/no_PASS ]]
then
rm $OUT_DIR/no_PASS
fi

cat $WHITE_LOCUS_LIST | while read locus
do
if [[ -e $SNP_DIR/detail/$locus/${locus}_filtered.vcf.gz ]]
then 
zcat $SNP_DIR/detail/$locus/${locus}_filtered.vcf.gz | grep -v "#" | awk -v L=$locus 'BEGIN{ok=0; n=0}{n++; if ($7=="PASS"){ok=1}}END{if (n>0 && ok!=1){print L}}' >> $OUT_DIR/no_PASS
fi
done

wc -l $OUT_DIR/no_PASS | awk '{print "Nombre de locus où aucun SNP détectés avec GATK ne passe les filtres de qualité de détection = no PASS: " $1}' 

# Nombre de locus avec SNP de bonne qualité mais sans correspondance avec Stacks
n=`grep -v locus $SNP_DIR/merge/${Proj_Name}_genome_gatk.tsv | awk '$4!="?" && $4!=""' | cut -f 1 | uniq | wc -l`
m=`grep -v locus $SNP_DIR/merge/${Proj_Name}_genome_gatk.tsv | awk '$4=="?"' | cut -f 1 | uniq | wc -l`

echo -e $n" locus avec SNP détecté par stacks.\n"$m" locus sans SNP détecté par Stacks "

# Nombre de locus avec tous les SNP détectés par stacks validé par GATK
if [[ -e $OUT_DIR/ok_stacks_GATK ]]
then
rm $OUT_DIR/ok_stacks_GATK
fi

cat $WHITE_LOCUS_LIST | while read locus
do
gz=`ls $STACKS_RES_DIR/batch*catalog.snps.tsv* | awk '{l=split($1,a,"."); if(a[l]=="gz"){print "gz"}else{print "tsv"}}'`
if [[ $gz = "gz" ]]
then
n=`zcat $STACKS_RES_DIR/batch*catalog.snps.tsv* | awk -v L=$locus '$3==L' | wc -l`
else
n=`cat $STACKS_RES_DIR/batch*catalog.snps.tsv* | awk -v L=$locus '$3==L' | wc -l`
fi
awk -v L=$locus -v N=$n '{if ($1==L && $4!="?" && $2!="no_snp_found"){c++}}END{if(c==N){print L}}' $SNP_DIR/merge/${Proj_Name}_genome_gatk.tsv >> $OUT_DIR/ok_stacks_GATK
done
