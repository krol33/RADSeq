#!/bin/sh

dir=$1
read1=$2
read2=$3
id=$4


cat $read1 | awk 'BEGIN{c=1}{if(c==1){print $0".1";c++}else{if(c<4){c++}else{if(c==4){c=1}} print $0 }}' > $dir/$id.fq
cat $read2 | awk 'BEGIN{c=1}{if(c==1){print $0".2";c++}else{if(c<4){c++}else{if(c==4){c=1}} print $0 }}' >> $dir/$id.fq
/home/sigenae/bin/fastq2fastaqual.py  $dir/$id.fq
formcon $dir/$id.fq.fasta 0 600
cap3 $dir/$id.fq.fasta > $dir/$id.cap3.out

sed s%Contig%Locus_${id}_Contig_% $dir/$id.fq.fasta.cap.contigs > $dir/$id.contig.fa 

rm $dir/$id.fq $dir/$id.fq.fasta $dir/$id.*.qual $dir/$id.fq.fasta.con $dir/$id.fq.fasta.cap.contigs $dir/$id*.links $dir/$id*.singlets  $dir/$id*.results  $dir/$id*.info $dir/$id*.out 
