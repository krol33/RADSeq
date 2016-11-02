#!/bin/sh

HAPLOTYPES=$1
ANALYSIS_TYPE=$2
POP_FILE=$3
OUT=$4

cat $HAPLOTYPES | awk -v TYPE=$ANALYSIS_TYPE -v POP=$POP_FILE 'BEGIN{
cpt_start=4; if (TYPE=="genotypes"){cpt_start=5};
while(getline<POP>0){pop[$1]=$2};
OFS="\t";
}{
if ($1=="Catalog"){
for (i=cpt_start; i<=NF ;i++){
index_col[i-1]=$i
homo[$i]=0
hetero[$i]=0
amb[$i]=0
un[$i]=0
}
}else{
for (i=cpt_start-1; i<=NF ;i++){
indiv=index_col[i]
if($i=="-"){
un[indiv]+=1
}else{if($i=="?"){
amb[indiv]+=1
}else{if(match($i,"/")){
hetero[indiv]+=1
}else{
homo[indiv]+=1
}
}
}
}
}
}END{
	print "#pop_num\tindiv\tnum_homo\tnum_hetero\tnum_ambiguous\tnum_unknown"
for (indiv in homo) {
print pop[indiv],indiv,homo[indiv],hetero[indiv],amb[indiv],un[indiv]
}	
}' > $OUT
