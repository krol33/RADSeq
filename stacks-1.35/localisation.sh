#!/bin/sh

GENOME=$1
PAIRS_DIR=$2
CONTIG=$3
LOC_DIR=$4
locus=$5
SCRIPT_DIR=$6
INDIV_FILE=$7

echo "# command line"
echo $0 $*

# coordonnées des lectures
echo -e "\n#coordonnées des lectures"

bwa mem -M $GENOME $PAIRS_DIR/${locus}_1.fq $PAIRS_DIR/${locus}_2.fq > $LOC_DIR/$locus/${locus}_pe_tmp.sam
nb_aln=`samtools view -S $LOC_DIR/$locus/${locus}_pe_tmp.sam | wc -l`

# ajout des read group
sh $SCRIPT_DIR/add_sam_RG.sh $LOC_DIR/$locus/${locus}_pe_tmp.sam $LOC_DIR/$locus/${locus}_pe_tmp_RG.sam $INDIV_FILE

# suppression des multimap
samtools view -Sh -F 0x0100 $LOC_DIR/$locus/${locus}_pe_tmp_RG.sam -bo $LOC_DIR/$locus/${locus}_pe_RG.bam
nb_aln2=`samtools view $LOC_DIR/$locus/${locus}_pe_RG.bam | wc -l`

echo -e "\tremove "`echo $nb_aln - $nb_aln2 | bc -l ` " non primary alignment"
bamToBed -i $LOC_DIR/$locus/${locus}_pe_RG.bam > $LOC_DIR/$locus/${locus}_pe.bed
if [ -s $LOC_DIR/$locus/${locus}_pe.bed ]
then 
	sort -k 1,1 -k 2,2n $LOC_DIR/$locus/${locus}_pe.bed | mergeBed -c 1 -o count -i - | awk -v L=$locus 'BEGIN{OFS="\t"}{print $1,$2,$3,"Locus_"L"_readPe",$4}' > $LOC_DIR/$locus/${locus}_readPe.bed
else
	echo -e "NULL\t0\t1\tLocus_"${locus}"_readPe\t0" > $LOC_DIR/$locus/${locus}_readPe.bed
fi
echo -e "\tnombre de localisations des reads: " `grep -v "NULL" $LOC_DIR/$locus/${locus}_readPe.bed  | wc -l  | awk '{print $1}'`
rm $LOC_DIR/$locus/${locus}_pe.bed

# coordonnées des contigs d'assemblage
echo -e "\n# coordonnées des contigs d'assemblage"
bwa mem -M $GENOME $CONTIG > $LOC_DIR/$locus/${locus}_contig.sam
sam2bed $LOC_DIR/$locus/${locus}_contig.sam $LOC_DIR/$locus/tmp.bed
if [ -s $LOC_DIR/$locus/tmp.bed ]
then
	awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,"1"}' $LOC_DIR/$locus/tmp.bed > $LOC_DIR/$locus/${locus}_contig.bed
else
	echo -e "NULL\t0\t1\tLocus_"${locus}"_Contig\t0" > $LOC_DIR/$locus/${locus}_contig.bed
fi
rm $LOC_DIR/$locus/tmp.bed

echo -e "\tnombre de localisations des "`grep -c '>' $CONTIG | awk '{print $1}' `" contigs: " `grep -v "NULL" $LOC_DIR/$locus/${locus}_contig.bed | wc -l  | awk '{print $1}'`

samtools view -Sh $LOC_DIR/$locus/${locus}_contig.sam -bo $LOC_DIR/$locus/${locus}_contig.bam

# intersection d'alignement entre lecture et contig sur ref
echo -e "\n#intersection d'alignement entre lecture et contig sur ref"
cat $LOC_DIR/$locus/${locus}_readPe.bed $LOC_DIR/$locus/${locus}_contig.bed | sort -k 1,1 -k 2,2n | mergeBed -c 4,5 -o collapse,sum -i - > $LOC_DIR/$locus/$locus.intersect.bed

echo -e "\tnombre de localisations finales de l'intersection des deux alignements: " `grep -v "NULL" $LOC_DIR/$locus/$locus.intersect.bed | wc -l | awk '{print $1}'`
rm $LOC_DIR/$locus/*.sam $LOC_DIR/$locus/${locus}_readPe.bed $LOC_DIR/$locus/${locus}_contig.bed
