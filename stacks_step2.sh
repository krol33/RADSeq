#!/bin/sh

#$ -o [Project_path_dir]/SGE_out/stacks2.out
#$ -e [Project_path_dir]/SGE_out/stacks2.err
#$ -N [Proj_Name]_stacks2
#$ -S /bin/sh

#$ -q unlimitq

#$ -m bea

date
date >&2
Proj_Name=[Proj_Name]

# VARIABLES GENERALES
stacks_dir=/usr/local/bioinfo/src/Stacks/stacks-1.35/bin
libR=[libR_path_dir]
RAD_DIR=[Project_path_dir]
SCRIPT_DIR=[Script_path_dir]

# INPUT
STACKS1_DIR=$RAD_DIR/stacks1
INDIV_FILE=$RAD_DIR/indiv_barcode.txt
dos2unix $INDIV_FILE
POP_FILE=$RAD_DIR/population.map
dos2unix $POP_FILE

# output
CSTACKS_DIR=$RAD_DIR/cstacks
SGE=$RAD_DIR/SGE_out/`basename $CSTACKS_DIR`
STAT_DIR=$RAD_DIR/stat/`basename $CSTACKS_DIR`

mkdir -p $SGE $CSTACKS_DIR $STAT_DIR

# VARIABLE DEDIEE
MISMATCH=1
# voulez vous basé le catalogue sur un autre ? (attention au version) 
KNOWN_CAT=
# travaillez  vous sur des données alignées ou non. oui : align=1, non : align=0
align=1

# INDIVIDUALS SELECTION (1 option au choix, si aucune option tous les individus sont pris en compte)
POP_IN=""
POP_OUT=""
INDIV_IN=""
INDIV_OUT=""

# VARIABLE CLUSTER
MEM="16G"
VMEM="24G"
THREADS=10

###########################################################################################
### catalogue 

### SELECTION DES INDIV 

if [[ -n $POP_IN ]]
then
	if [[ -n $POP_OUT || -n $INDIV_IN || -n $INDIV_OUT ]]
	then
	echo "WRONG1, you must choose one selection option exclusively" >&2
	exit
	fi
else
	if [[ -n $POP_OUT ]]
	then
		if [[ -n $INDIV_IN || -n $INDIV_OUT ]]
		then
		echo "WRONG2, you must choose one selection option exclusively" >&2
		exit
		fi
	else
		if [[ -n $INDIV_IN && -n $INDIV_OUT ]]
		then
		echo "WRONG3, you must choose one selection option exclusively" >&2
		exit
		fi
	fi
fi

args=""

# tous ceux des pop POP_IN
if [[ -n $POP_IN ]]
then
args=`awk -v P="$POP_IN" -v U=$STACKS1_DIR '{if (match(P,$2)){s=s" -s "U"/"$1}}END{print s}' $POP_FILE `
fi

# tous sauf ceux des pop POP_OUT
if [[ -n $POP_IN ]]
then
args=`awk -v P="$POP_OUT" -v U=$STACKS1_DIR '{if (! match(P,$2)){s=s" -s "U"/"$1}}END{print s}' $POP_FILE `
fi

# tous les individus INDIV_IN
if [[ -n $INDIV_IN ]]
then
args=`awk -v P="$INDIV_IN" -v U=$STACKS1_DIR '{if (match(P,$1)){s=s" -s "U"/"$1}}END{print s}' $POP_FILE `
fi

# tous les individus sauf INDIV_OUT
if [[ -n $INDIV_OUT ]]
then
args=`awk -v P="$INDIV_OUT" -v U=$STACKS1_DIR '{if (! match(P,$1)){s=s" -s "U"/"$1}}END{print s}' $POP_FILE `
fi

if [[ ! -n "$args" ]]
then
args=`awk -v U=$STACKS1_DIR '{s=s" -s "U"/"$1}END{print s}' $POP_FILE`
fi

cmd=`echo $stacks_dir/cstacks -b 1 $args -o $CSTACKS_DIR -p $THREADS -n $MISMATCH `

# à partir d'un catalogue 
if [[ -n $KNOWN_CAT ]]
then
cmd=`echo $cmd --catalog $KNOWN_CAT`
fi

# à partir de données alignées
if [[ $align == 1 ]]
then
cmd=`echo $cmd -g`
fi

echo $cmd > $SGE/cstacks.tmp.sh
id_cstacks=`qsub -terse -N ${Proj_Name}_cstacks -l mem=$MEM -l h_vmem=$VMEM -o $SGE/cstacks.out -e $SGE/cstacks.err -pe parallel_smp $THREADS $SGE/cstacks.tmp.sh`
r=0;
while [ $r == 0 ] ; do qstat -j $id_cstacks >& /dev/null ; r=$? ; sleep 10s ;done 

####################################################################################################################
###########################   CSTACK stat : $CSTACKS_DIR ######################################################
echo -e "\n###########################   CSTACK stat : $CSTACKS_DIR ######################################################"
echo -e "\n###########################   CSTACK stat : $CSTACKS_DIR ######################################################" >&2 
echo -e "###########################   CSTACK stat : $CSTACKS_DIR ######################################################"  > $STAT_DIR/summary_statistics.txt

id_gzip_d=`qsub -terse -N ${Proj_Name}_gzip-d_cstacks -o /dev/null -e /dev/null -b y "gzip -d $CSTACKS_DIR/*"`
r=0;
while [ $r == 0 ] ; do qstat -j $id_gzip_d >& /dev/null ; r=$? ; sleep 10s ; done 

# nombre de cluster
echo -e "\nnombre de cluster dans le catalogue et nombre d'individus " >> $STAT_DIR/summary_statistics.txt
wc -l $CSTACKS_DIR/batch_1.catalog.tags.tsv >> $STAT_DIR/summary_statistics.txt

if [[ -e $CSTACKS_DIR/tmp1 ]]
then
rm $CSTACKS_DIR/tmp*
fi

TOT_IND=`grep "Constructing catalog from" $SGE/cstacks.err | awk '{print $4}'`

grep -v "#" $CSTACKS_DIR/batch_1.catalog.tags.tsv | 
awk -v TOTIND=$TOT_IND -v outDist=$STAT_DIR/nb_indiv_per_clust.txt -v out1=$CSTACKS_DIR/tmp1 -v out2=$CSTACKS_DIR/tmp2 -v out50p=$CSTACKS_DIR/tmp50p -v out100p=$CSTACKS_DIR/tmp100p -v POP=$POP_FILE -v IND=$INDIV_FILE 'BEGIN{single=0; deux=0; deux2=0; fifty=0;all=0; while(getline<POP>0){tabP[$1]=$2} while(getline<IND>0){tabI[$1]=$2}; print "#Nb_ind\tnb_cluster" > outDist}
{
a=NF-5
compo=$a
split(compo,tab,",");
for (i in tab){
	split(tab[i],indiv,"_")
	if(tmpI[indiv[1]]!=1){tmpI[indiv[1]]=1; len++}
	pop=tabP[tabI[indiv[1]]]
	if(tmpP[pop]!=1){tmpP[pop]=1; p++}
}
tabCount[len]++;
if(len==1){single++; print $0 >>out1 }; 
if(len>=2){deux++; print $0>>out2; if(p>1){deux2++;}};
if(len>=TOTIND/2){fifty++; print $0>>out50p};
if(len==TOTIND){all++; print $0>>out100p};
len=0; p=0; delete tmpI; delete tmpP;
}END{
for (i=1; i<=TOTIND; i++){
if (tabCount[i] ==""){print i"\t0" >> outDist}
else{print i"\t"tabCount[i] >> outDist}
}
print "\t"single" single individual locus";
print "\t"deux" merge >=2 indiv (dont "deux2" merge >=2 indiv au moins 2 populations différentes)";
print "\t"fifty" merge >=50 pourcent des indiv";
print "\t"all" merge 100 pourcent des  indiv";}'  >> $STAT_DIR/summary_statistics.txt

# nombre de cluster à X individus agrégés
echo -e "\nsee detailed number of cluster in function of number of aggregated individuals available: " $STAT_DIR/nb_indiv_per_clust.txt" and "$STAT_DIR/nb_indiv_per_clust.jpg >> $STAT_DIR/summary_statistics.txt
$libR/plot_dist.R "nb cluster in function of number of aggregated individuals" "number of aggregated individuals" "number of cluster" 1 $STAT_DIR/nb_indiv_per_clust $STAT_DIR/nb_indiv_per_clust.txt

if [[ $align == 0 ]]
then 
# evolution du catalogue en fonction du nombre d'individus agrégés
echo -e "\nplot the distribution of the evolution of the number of catalog's loci in function of the number of sample added into it: "$STAT_DIR/catalog_evol.dist >> $STAT_DIR/summary_statistics.txt

grep -E "tags.tsv.gz|loci in the catalog" $SGE/cstacks.err | awk -v D=$STACKS1_DIR '{sub(D"/","",$0); sub(".tags.tsv.gz","",$NF); if($1=="Parsing" && id==""){id=$NF}else{if(loci=="" && $1 !="Parsing"){loci=$1}else{loci=0; print id"\t"loci; id=$NF; loci=""}}; if(id!="" && loci!=""){print id"\t"loci; id=""; loci=""} }' > $STAT_DIR/catalog_evol.dist
fi

# nombre de SNP dans le catalogue
echo -e "\nnombre de clusters et nombre de SNP" >> $STAT_DIR/summary_statistics.txt
awk '$5=="E"' $CSTACKS_DIR/*.snps.tsv | cut -f 3 | uniq -c > $STAT_DIR/tmp

cat $STAT_DIR/tmp | colstat.sh 1 >> $STAT_DIR/summary_statistics.txt
cat $STAT_DIR/tmp|  distrib_class_percent.sh 1 > $STAT_DIR/catalog.SNP.dist
$libR/plot_dist.R "number of SNP per cluster distribution" "nb SNP/cluster" "proportion" 1 $STAT_DIR/catalog.SNP.dist $STAT_DIR/catalog.SNP.dist
echo -e "\nnombre de clusters provenant d'un seul individus et nombre de SNP" >> $STAT_DIR/summary_statistics.txt
cat $STAT_DIR/tmp| awk -v clustertype=$CSTACKS_DIR/tmp1 'BEGIN{while(getline<clustertype>0){tab[$3]=1}}{if(tab[$2]==1){print $1}}' | colstat.sh >> $STAT_DIR/summary_statistics.txt
echo -e "\nnombre de clusters provenant de plus de 2 individus et nombre de SNP" >> $STAT_DIR/summary_statistics.txt
cat $STAT_DIR/tmp| awk -v clustertype=$CSTACKS_DIR/tmp2 'BEGIN{while(getline<clustertype>0){tab[$3]=1}}{if(tab[$2]==1){print $1}}' | colstat.sh >> $STAT_DIR/summary_statistics.txt
echo -e "\nnombre de clusters provenant de plus 50 % des individus et nombre de SNP" >> $STAT_DIR/summary_statistics.txt
cat $STAT_DIR/tmp|  awk -v clustertype=$CSTACKS_DIR/tmp50p 'BEGIN{while(getline<clustertype>0){tab[$3]=1}}{if(tab[$2]==1){print $1}}' | colstat.sh >> $STAT_DIR/summary_statistics.txt
echo -e "\nnombre de clusters provenant de 100 % des individus et nombre de SNP" >> $STAT_DIR/summary_statistics.txt
cat $STAT_DIR/tmp| awk -v clustertype=$CSTACKS_DIR/tmp100p 'BEGIN{while(getline<clustertype>0){tab[$3]=1}}{if(tab[$2]==1){print $1}}' | colstat.sh >> $STAT_DIR/summary_statistics.txt


# nombre d'haplotype par cluster
echo -e "\nnombre d'happlotypes par cluster avec SNP" >> $STAT_DIR/summary_statistics.txt
awk 'NF==6' $CSTACKS_DIR/batch_1.catalog.alleles.tsv | cut -f 3 | uniq -c | awk '$1==2' | wc -l | awk '{print "\tnombre de cluster avec 2 happlotypes: "$1}' >> $STAT_DIR/summary_statistics.txt
awk 'NF==6' $CSTACKS_DIR/batch_1.catalog.alleles.tsv | cut -f 3 | uniq -c | awk '$1==3' | wc -l| awk '{print "\tnombre de cluster avec 3 happlotypes: "$1}' >> $STAT_DIR/summary_statistics.txt
awk 'NF==6' $CSTACKS_DIR/batch_1.catalog.alleles.tsv | cut -f 3 | uniq -c | awk '$1==4' | wc -l| awk '{print "\tnombre de cluster avec 4 happlotypes: "$1}' >> $STAT_DIR/summary_statistics.txt
awk 'NF==6' $CSTACKS_DIR/batch_1.catalog.alleles.tsv | cut -f 3 | uniq -c | awk '$1>=5' | wc -l| awk '{print "\tnombre de cluster avec >=5 happlotypes: "$1}' >> $STAT_DIR/summary_statistics.txt
awk 'NF==6' $CSTACKS_DIR/batch_1.catalog.alleles.tsv | cut -f 3 | uniq -c | awk '$1>=10' | wc -l| awk '{print "\tnombre de cluster avec >=10 happlotypes: "$1}' >> $STAT_DIR/summary_statistics.txt
awk 'NF==6' $CSTACKS_DIR/batch_1.catalog.alleles.tsv | cut -f 3 | uniq -c | sort -k1,1n | tail -n 1 | awk '{print "\tnombre max d haplotypes: "$1}' >> $STAT_DIR/summary_statistics.txt
awk 'NF==6' $CSTACKS_DIR/batch_1.catalog.alleles.tsv | cut -f 3 | uniq -c | distrib_class_percent.sh  1 > $STAT_DIR/catalog.haplotype.dist
$libR/plot_dist.R "number of haplotype per cluster distribution" "nb haplotype/cluster" "proportion" 1 $STAT_DIR/catalog.haplotype.dist $STAT_DIR/catalog.haplotype.dist


qsub -N ${Proj_Name}_gzip_cstacks -o /dev/null -e /dev/null -b y "gzip $CSTACKS_DIR/*.tsv"
