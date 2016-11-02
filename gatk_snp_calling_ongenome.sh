#/bin/sh

echo "# command line"
echo $0 $*


LOCUS=$1
READ_BAM=$2
GENOME=$3


PICARD_DIR=/usr/local/bioinfo/src/picard-tools/current
GATK_PATH=/usr/local/bioinfo/src/GATK/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar

QUAL=$4
AN=$5

CATALOGUE=$6
SNP_DIR=$7

SCRIPT_DIR=$8

rm $SNP_DIR/$LOCUS/*

# sort du bam par locus
echo -e "\n#sort "
echo -e "\n#sort " >&2
samtools sort $READ_BAM $SNP_DIR/$LOCUS/$LOCUS.RG.sort

# filter uniq, mapQ30 alignment
echo -e "\n#filter uniq, mapQ30 alignment and index bam"
echo -e "\n#filter uniq, mapQ30 alignment" >&2
samtools view -F 4 -F 8 -q 30 -h $SNP_DIR/$LOCUS/$LOCUS.RG.sort.bam | awk '{if(substr($1,0,1)=="@"){print $0}else{a=match($0,"XA:");if(a==0){print $0}}}' | samtools view -S - -b -o $SNP_DIR/$LOCUS/$LOCUS.RG.sort.mapQ.bam
rm $SNP_DIR/$LOCUS/$LOCUS.RG.sort.bam

# fasta and gatk indexes contig ref
# récupération des séquence de ref où il y a des lectures
echo -e "\n#extract reference alignment sequence"
echo -e "\n#extract reference alignment sequence" >&2
samtools view $SNP_DIR/$LOCUS/$LOCUS.RG.sort.mapQ.bam | cut -f 3 | sort -u > $SNP_DIR/$LOCUS/ref.list
samtools view -h $SNP_DIR/$LOCUS/$LOCUS.RG.sort.mapQ.bam | awk -v REF=$SNP_DIR/$LOCUS/ref.list 'BEGIN{while(getline<REF>0){tab["SN:"$1]=1}}{if($1=="@SQ"){if(tab[$2]==1){print $0}}else{print $0}}' | samtools view -Sh - -bo $SNP_DIR/$LOCUS/$LOCUS.RG.sort.mapQ_extractRef.bam
samtools index $SNP_DIR/$LOCUS/$LOCUS.RG.sort.mapQ_extractRef.bam
rm $SNP_DIR/$LOCUS/$LOCUS.RG.sort.mapQ.bam

cat $GENOME | fasta_extract.pl $SNP_DIR/$LOCUS/ref.list > $SNP_DIR/$LOCUS/ref.fasta
echo -e "\n#fasta and gatk indexes contig ref"
echo -e "\n#fasta and gatk indexes contig ref" >&2
samtools faidx $SNP_DIR/$LOCUS/ref.fasta
java -jar -Xmx4G $PICARD_DIR/CreateSequenceDictionary.jar R=$SNP_DIR/$LOCUS/ref.fasta O=$SNP_DIR/$LOCUS/ref.dict

# SNP INDEL calling on all read
echo -e "\n#SNP calling on all read"
echo -e "\n#SNP calling on all read" >&2
java -jar -Xmx8G $GATK_PATH  -nt 8 -T UnifiedGenotyper -glm SNP -R $SNP_DIR/$LOCUS/ref.fasta -o $SNP_DIR/$LOCUS/$LOCUS.snp.vcf -I $SNP_DIR/$LOCUS/$LOCUS.RG.sort.mapQ_extractRef.bam
bgzip $SNP_DIR/$LOCUS/$LOCUS.snp.vcf 
java -jar -Xmx8G $GATK_PATH  -nt 8 -T UnifiedGenotyper -glm INDEL -R $SNP_DIR/$LOCUS/ref.fasta -o $SNP_DIR/$LOCUS/$LOCUS.indel.vcf -I $SNP_DIR/$LOCUS/$LOCUS.RG.sort.mapQ_extractRef.bam
bgzip $SNP_DIR/$LOCUS/$LOCUS.indel.vcf 
vcf-concat $SNP_DIR/$LOCUS/$LOCUS.snp.vcf.gz $SNP_DIR/$LOCUS/$LOCUS.indel.vcf.gz | vcf-sort > $SNP_DIR/$LOCUS/$LOCUS.vcf

rm $SNP_DIR/$LOCUS/$LOCUS.snp.vcf* $SNP_DIR/$LOCUS/$LOCUS.indel.vcf*

# filter
echo -e "\n#SNP filter"
echo -e "\n#SNP filter" >&2
java -jar -Xmx4G $GATK_PATH -T VariantFiltration -R $SNP_DIR/$LOCUS/ref.fasta -V $SNP_DIR/$LOCUS/$LOCUS.vcf --filterExpression "QUAL < $QUAL || AN < $AN || QD < 2.0 || FS > 60.0 || MQ < 40.0 ||  MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "snp_indel_filter" -o $SNP_DIR/$LOCUS/${LOCUS}_filtered.vcf
bgzip $SNP_DIR/$LOCUS/${LOCUS}_filtered.vcf

rm $SNP_DIR/$LOCUS/ref* $SNP_DIR/$LOCUS/$LOCUS.vcf*

# check correspondance avec le catalogue stacks
echo -e "\n#check stacks correspondancies"
echo -e "\n#check stacks correspondancies" >&2
python $SCRIPT_DIR/GATK2tab.py -v $SNP_DIR/$LOCUS/${LOCUS}_filtered.vcf.gz -l $LOCUS -b $SNP_DIR/$LOCUS/$LOCUS.RG.sort.mapQ_extractRef.bam -c $CATALOGUE
rm $SNP_DIR/$LOCUS/$LOCUS.RG.sort.mapQ_extractRef.bam* $SNP_DIR/$LOCUS/*.idx


