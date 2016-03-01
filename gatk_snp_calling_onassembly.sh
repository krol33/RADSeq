#/bin/sh


echo "# command line"
echo $0 $*

LOCUS=$1
PAIRS_DIR=$2
GENOME=$3
INDIV_FILE=$4


PICARD_DIR=/usr/local/bioinfo/src/picard-tools/current
GATK_PATH=/usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar

QUAL=$5
AN=$6

CATALOGUE=$7
SNP_DIR=$8

SCRIPT_DIR=$9


# index bwa
if [[ ! -e $GENOME.amb ]]
then
bwa index $GENOME
fi

# alignement des reads sur le contig
echo -e "\n# alignement des reads sur le contig"
echo -e "\n# alignement des reads sur le contig" >&2
bwa mem -M $GENOME $PAIRS_DIR/${LOCUS}_1.fq $PAIRS_DIR/${LOCUS}_2.fq > $SNP_DIR/$LOCUS/${LOCUS}_pe_tmp.sam
nb_aln=`samtools view -S $SNP_DIR/$LOCUS/${LOCUS}_pe_tmp.sam | wc -l`

# ajout des read group
sh $SCRIPT_DIR/add_sam_RG.sh $SNP_DIR/$LOCUS/${LOCUS}_pe_tmp.sam $SNP_DIR/$LOCUS/${LOCUS}_pe_tmp_RG.sam $INDIV_FILE

#suppression des multimap
samtools view -Sh -F 0x0100 $SNP_DIR/$LOCUS/${LOCUS}_pe_tmp_RG.sam -bo $SNP_DIR/$LOCUS/${LOCUS}_pe_RG.bam
nb_aln2=`samtools view $SNP_DIR/$LOCUS/${LOCUS}_pe_RG.bam | wc -l`
echo -e "\tremove "`echo $nb_aln - $nb_aln2 | bc -l ` " non primary alignment on " $nb_aln " initial alignement"

# sort du bam par locus
echo -e "\n#sort "
echo -e "\n#sort " >&2
samtools sort $SNP_DIR/$LOCUS/${LOCUS}_pe_RG.bam $SNP_DIR/$LOCUS/$LOCUS.RG.sort

# nettoyage
rm $SNP_DIR/$LOCUS/*.sam $SNP_DIR/$LOCUS/${LOCUS}_pe_RG.bam

# filter uniq, mapQ30 alignment
echo -e "\n#filter uniq, mapQ30 alignment and index bam"
echo -e "\n#filter uniq, mapQ30 alignment" >&2
samtools view -F 4 -F 8 -q 30 -h $SNP_DIR/$LOCUS/$LOCUS.RG.sort.bam | awk '{if(substr($1,0,1)=="@"){print $0}else{a=match($0,"XA:");if(a==0){print $0}}}' | samtools view -S - -b -o $SNP_DIR/$LOCUS/$LOCUS.RG.sort.mapQ.bam
samtools index $SNP_DIR/$LOCUS/$LOCUS.RG.sort.mapQ.bam
nb_aln=`samtools view $SNP_DIR/$LOCUS/$LOCUS.RG.sort.bam | wc -l`
nb_aln2=`samtools view $SNP_DIR/$LOCUS/$LOCUS.RG.sort.mapQ.bam | wc -l`
echo -e "\tremove "`echo $nb_aln - $nb_aln2 | bc -l ` " non uniq or low quality alignment on " $nb_aln " initial alignement"

rm $SNP_DIR/$LOCUS/$LOCUS.RG.sort.bam

# fasta and gatk indexes contig ref
echo -e "\n#fasta and gatk indexes contig ref"
echo -e "\n#fasta and gatk indexes contig ref" >&2
if [[ ! -e $GENOME.fai ]]
then
samtools faidx $GENOME
fi
gatk_index=`echo $GENOME | sed 's/\.fa/\.dict/' `
if [[ ! -e $gatk_index ]]
then
java -jar -Xmx4G $PICARD_DIR/CreateSequenceDictionary.jar R=$GENOME O=$gatk_index
fi

# SNP INDEL calling on all read
echo -e "\n#SNP calling on all read"
echo -e "\n#SNP calling on all read" >&2
module load bioinfo/Java7; java -jar -Xmx4G $GATK_PATH  -nt 8 -T UnifiedGenotyper -glm SNP -drf BadMate -R $GENOME -o $SNP_DIR/$LOCUS/$LOCUS.snp.vcf -I $SNP_DIR/$LOCUS/$LOCUS.RG.sort.mapQ.bam
bgzip $SNP_DIR/$LOCUS/$LOCUS.snp.vcf 
module load bioinfo/Java7; java -jar -Xmx4G $GATK_PATH  -nt 8 -T UnifiedGenotyper -glm INDEL -drf BadMate -R $GENOME -o $SNP_DIR/$LOCUS/$LOCUS.indel.vcf -I $SNP_DIR/$LOCUS/$LOCUS.RG.sort.mapQ.bam
bgzip $SNP_DIR/$LOCUS/$LOCUS.indel.vcf 
vcf-concat $SNP_DIR/$LOCUS/$LOCUS.snp.vcf.gz $SNP_DIR/$LOCUS/$LOCUS.indel.vcf.gz | vcf-sort > $SNP_DIR/$LOCUS/$LOCUS.vcf

rm $SNP_DIR/$LOCUS/$LOCUS.snp.vcf* $SNP_DIR/$LOCUS/$LOCUS.indel.vcf*

# filter
echo -e "\n#SNP filter"
echo -e "\n#SNP filter" >&2
module load bioinfo/Java7; java -jar -Xmx4G $GATK_PATH -T VariantFiltration -drf BadMate -R $GENOME -V $SNP_DIR/$LOCUS/$LOCUS.vcf --filterExpression "QUAL < $QUAL || AN < $AN || QD < 2.0 || FS > 60.0 || MQ < 40.0 ||  MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "snp_indel_filter" -o $SNP_DIR/$LOCUS/${LOCUS}_filtered.vcf
bgzip $SNP_DIR/$LOCUS/${LOCUS}_filtered.vcf

rm $SNP_DIR/$LOCUS/$LOCUS.vcf*

# check correspondance avec le catalogue stacks
echo -e "\n#check stacks correspondancies"
echo -e "\n#check stacks correspondancies" >&2
python $SCRIPT_DIR/GATK2tab.py -v $SNP_DIR/$LOCUS/${LOCUS}_filtered.vcf.gz -l $LOCUS -b $SNP_DIR/$LOCUS/$LOCUS.RG.sort.mapQ.bam -c $CATALOGUE
rm $SNP_DIR/$LOCUS/*.idx


