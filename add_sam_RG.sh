#!/bin/sh

in_sam=$1
out_sam=$2
indiv=$3

cat $in_sam | awk -v I=$indiv 'BEGIN{while(getline<I>0){tab[$2]=$1} header="True"}{
	if(substr($1,0,1) =="@" && header=="True"){print $0}
	if(substr($1,0,1) !="@" && header=="True"){for (i in tab){print "@RG\tID:"tab[i]"\tSM:"i}; header="False"}
	if(substr($1,0,1) !="@" && header!="True"){read=$1; split(read,l,"|"); print $0"\tRG:Z:"tab[l[2]]}
	}' > $out_sam
	
	



