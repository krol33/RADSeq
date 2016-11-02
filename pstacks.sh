#!/bin/sh
read2=""
ddrad=0

while [ "$1" != "" ]; do
    PARAM=`echo $1 | awk -F= '{print $1}'`
    VALUE=`echo $1 | awk -F= '{print $2}'`
    case $PARAM in
        -h | --help)
            usage
            exit
            ;;
        --read1)
            read1=$VALUE
            ;;
        --rem1)
            rem1=$VALUE
            ;;
        --read2)
            read2=$VALUE
            ;;
        --rem2)
            rem2=$VALUE
            ;;
        --indiv)
            indiv=$VALUE
            ;;
        --id)
            id=$VALUE
            ;;
        --genomeIndex)
            genomeIndex=$VALUE
            ;;            
        --out)
            out=$VALUE
            ;;
        --bin_dir)
            bin_dir=$VALUE
            ;;
        --cov)
            cov=$VALUE
            ;;
        --opt)
            opt=$VALUE
            ;;
        --ddrad)
            ddrad=1
            ;; 
        *)
            echo "ERROR: unknown parameter \"$PARAM\""
            usage
            exit 1
            ;;
    esac
    shift
done

if [[ ! -e $genomeIndex || ! -e $genomeIndex".ann" ]]
then
echo "your genome doesn't exists or is not BWA indexed."
exit 1
fi

# alignement bwa mem : read1 et potentiellement read2
echo -e "BWA alignement :\n Cmd line: \t bwa mem $genomeIndex $read1 $read2 > $out/${indiv}.sam" 
bwa mem $genomeIndex $read1 $read2 > $out/${indiv}.sam
samtools view -Sh $out/${indiv}.sam -bo $out/${indiv}.bam

# alignement potentiel rem1
if [[ -e $rem1 ]]
then
echo " Cmd line: \t bwa mem $genomeIndex $rem1 > $out/${indiv}_rem1.sam" 
bwa mem $genomeIndex $rem1 > $out/${indiv}_rem1.sam
samtools view -Sh $out/${indiv}_rem1.sam -bo $out/${indiv}_rem1.bam
merged="$out/${indiv}_rem1.bam"
fi

# alignement potentiel rem2
if [[ -e $rem2 ]]
then
echo " Cmd line: \t bwa mem $genomeIndex $rem2 > $out/${indiv}_rem2.sam" 
bwa mem $genomeIndex $rem2 > $out/${indiv}_rem2.sam
samtools view -Sh $out/${indiv}_rem2.sam -bo $out/${indiv}_rem2.bam
merged=$merged" $out/${indiv}_rem2.bam"
fi

# merge read1 et? read2 et? rem1 et? rem2
to_filter_bam=$out/${indiv}.bam
if [[ -e $rem1  || -e $rem2 ]]
then
to_filter_bam=$out/${indiv}_merged.bam
echo -e "Merge read1 et? read2 et? rem1 et? rem2\nCmd line : \t samtools merge $to_filter_bam $out/${indiv}.bam $merged "
samtools merge $to_filter_bam $out/${indiv}.bam $merged 
fi

# stat befor filtering
echo "Alignment statistics befor filtering"
samtools flagstat $to_filter_bam

echo "Filter on alignment quality Q20, uniq alignement"
if [[ $ddrad == 1 || ! -e $read2 ]]
then
echo -e " Cmd line : \t samtools view -h $to_filter_bam -q 20 -F 0x0100 -bo $out/${indiv}_uniqMapQ20.bam"
samtools view -h $to_filter_bam -q 20 -F 0x0100 -bo $out/${indiv}_uniqMapQ20.bam
final_bam=$out/${indiv}_uniqMapQ20.bam
else
echo -e " Cmd line : \t samtools view -h $to_filter_bam -q 20 -f 0x40 -F 0x0100 -bo $out/${indiv}_uniqMapQ20_R1.bam"
samtools view -h $to_filter_bam -q 20 -f 0x40 -F 0x0100 -bo $out/${indiv}_uniqMapQ20_R1.bam
final_bam=$out/${indiv}_uniqMapQ20_R1.bam
fi

# stat befor filtering
echo "Alignment statistics after filtering (= input Stacks)"
samtools flagstat $final_bam

# search stack for one individual
echo -e "PSTACKS:\n Cmd line : \t $bin_dir/pstacks -t bam -f $final_bam -i $id -o $out -m $cov $opt"
$bin_dir/pstacks -t bam -f $final_bam -i $id -o $out -m $cov $opt

for i in `ls $out/${indiv}_uniqMapQ20*.tsv.gz`
do
suffix=`basename $i  | sed 's/\./ /g' | awk '{ s=""; for (i=NF-2; i<=NF; i++){s=s"."$i}}END{print s}'`
mv $i $out/$indiv""$suffix
done

for i in `ls $out/${indiv}_*am $out/${indiv}.*am`
do
if [[ $i != $final_bam ]]
then rm $i
fi
done

