#!/bin/sh

read1=$1
read2=$2
sample=`basename $read1 | sed -e 's/_1.fq/ /' |awk '{print $1}'`
outFile=`echo $read1 | sed 's/_1.fq/_2.fq/'`

echo $sample, $outFile

R2_type=`echo $read2 | awk '{if(match($1,".gz") > 0 ){print "gzfastq"}else{print "fastq"}}'`

if [[ $R2_type = "gzfastq" ]]
then
zcat $read2 | awk -v NAME=$read1 -v sample=$sample -v out=$outFile '
BEGIN{
	sample="@"sample
	cpt=1
	while (getline < NAME > 0){
		if(cpt==1){
			tab[$1]=1
		}
		if (cpt==4){cpt=0}
		cpt ++
	}
	name=""
	cpt=1
}{
	if (cpt==1){
		if (name!=""){
			print nameTot >> out
			print seq >> out
			print name2 >> out
			print qual >> out
			name=""
			seq=""
			name2=""
			qual=""
		}
		name=$1
		nameTot=$0
	}else{
		if (tab[name]==1){
			if (seq==""){
				seq=$0
			}else{
				if (name2==""){
					name2=$0
				}else{
					if (qual==""){
					qual=$0
					}
				}
			}
		}else{
			name=""
		}
	}
	if (cpt==4){cpt=0}
	cpt ++
}END {
	if (name!=""){
			print nameTot >> out
			print seq >> out
			print name2 >> out
			print qual >> out
			name=""
			seq=""
			name2=""
			qual=""
		}
}
' 
else
cat $read2 | awk -v NAME=$read1 -v sample=$sample -v out=$outFile '
BEGIN{
	sample="@"sample
	cpt=1
	while (getline < NAME > 0){
		if(cpt==1){
			tab[$1]=1
		}
		if (cpt==4){cpt=0}
		cpt ++
	}
	name=""
	cpt=1
}{
	if (cpt==1){
		if (name!=""){
			print nameTot >> out
			print seq >> out
			print name2 >> out
			print qual >> out
			name=""
			seq=""
			name2=""
			qual=""
		}
		name=$1
		nameTot=$0
	}else{
		if (tab[name]==1){
			if (seq==""){
				seq=$0
			}else{
				if (name2==""){
					name2=$0
				}else{
					if (qual==""){
					qual=$0
					}
				}
			}
		}else{
			name=""
		}
	}
	if (cpt==4){cpt=0}
	cpt ++
}END {
	if (name!=""){
			print nameTot >> out
			print seq >> out
			print name2 >> out
			print qual >> out
			name=""
			seq=""
			name2=""
			qual=""
		}
}
'
fi

if [[ $R2_type = "fastq" ]]
then
gzip $read2
fi
gzip $read1 $outFile
