#!/bin/sh

dir=$1
read1=$2
read2=$3
id=$4

gz=`echo $read1 | awk '{if(match($1,".gz") > 0 ){print "gzfastq"}else{print "fastq"}}'`

if [[ $gz = "gzfastq" ]]
then
zcat $read1 | awk 'BEGIN{c=1}{if(c==1){print $0".1";c++}else{if(c<4){c++}else{if(c==4){c=1}} print $0 }}' > $dir/$id.fq
zcat $read2 | awk 'BEGIN{c=1}{if(c==1){print $0".2";c++}else{if(c<4){c++}else{if(c==4){c=1}} print $0 }}' >> $dir/$id.fq
else
cat $read1 | awk 'BEGIN{c=1}{if(c==1){print $0".1";c++}else{if(c<4){c++}else{if(c==4){c=1}} print $0 }}' > $dir/$id.fq
cat $read2 | awk 'BEGIN{c=1}{if(c==1){print $0".2";c++}else{if(c<4){c++}else{if(c==4){c=1}} print $0 }}' >> $dir/$id.fq
fi

/home/sigenae/bin/fastq2fastaqual.py  $dir/$id.fq
formcon $dir/$id.fq.fasta 0 600
cap3 $dir/$id.fq.fasta > $dir/$id.cap3.out

sed s%Contig%Locus_${id}_Contig_% $dir/$id.fq.fasta.cap.contigs > $dir/$id.contig.fa 

rm $dir/$id.fq $dir/$id.fq.fasta $dir/$id.*.qual $dir/$id.fq.fasta.con $dir/$id.fq.fasta.cap.contigs $dir/$id*.links $dir/$id*.singlets  $dir/$id*.results  $dir/$id*.info $dir/$id*.out 

if [[ $gz = "fastq" ]]
then
gzip $read1 $read2
fi
