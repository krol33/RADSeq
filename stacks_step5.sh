#$ -o [Project_path_dir]/SGE_out/stacks5.out
#$ -e [Project_path_dir]/SGE_out/stacks5.err
#$ -N $[Proj_Name]_stacks5
#$ -S /bin/sh

#$ -q unlimitq

#$ -m bea

date
date >&2
Proj_Name=

# VARIABLES GENERALES
# type d'analyse, croisement (ANALYSIS_TYPE="genotypes") ou populations (ANALYSIS_TYPE="populations")
ANALYSIS_TYPE="populations"
RAD_DIR=[Project_path_dir]
SCRIPT_DIR=[Script_path_dir]

# INPUT
INDIV_FILE=$RAD_DIR/indiv_barcode.txt
dos2unix $INDIV_FILE
POP_FILE=$RAD_DIR/${ANALYSIS_TYPE}/pop_extract.map
HAPLOTYPE_DIR=$RAD_DIR/${ANALYSIS_TYPE}
POLYMORPHES_LOCUS=$HAPLOTYPE_DIR/locus_polymorphes.txt

# VARIABLES DEDIEES
# si vous avez une population controle type happloïde doublée, indiquez son numéros ici. Ces individus seront filtrés du tableau d'haplotypage final ainsi que les locus qui sont hétérozygote pour au moins deux de ces individus
POP_CONTROLE="1 2"
# nombre d'individus HD devant être hétéro pour considérer le locus comme duplicat
NB_INDIV=1
# Fichier qui contient les noms des individus que vous ne voulez pas concerver car trop peu génotypés. cf ${ANALYSIS_TYPE}/nb_genotype.txt
BAD_INDIV=$RAD_DIR/${ANALYSIS_TYPE}/bad_indiv

# output
FILTER_DIR=$RAD_DIR/filter_and_stat
SGE=$RAD_DIR/SGE_out/`basename $FILTER_DIR`

mkdir -p $SGE $FILTER_DIR

###########################################################################################
###########################################################################################
############### filtres des locus faux pour les happloïdes doublés
if [[ "$POP_CONTROLE" != "" ]]
then
	# extraction des locus pour lesquels au moins 2 individus control sont hétérozygotes
BAD_LOC=$FILTER_DIR/bad_locus_$NB_INDIV
if [[ -e $FILTER_DIR/control_haplotye.tsv ]] 
then
rm $FILTER_DIR/control_haplotye.tsv
fi

awk -v H=$NB_INDIV -v N="$POP_CONTROLE" -v out_hap=$FILTER_DIR/control_haplotye.tsv -v POP=$POP_FILE  -v TYPE=$ANALYSIS_TYPE 'BEGIN{cpt_start=2; if (TYPE=="genotypes"){cpt_start=3} while(getline<POP>0){c++;if(match(N,$2)){tab[c+cpt_start]=$1}}; 
FS="\t" }
{c=0; hap=""; s=$1 ; for (i=2;i<=cpt_start;i++){s=s"\t"$i} ; for (i=cpt_start+1;i<=NF;i++){ if(tab[i]!=""){ s= s"\t"$i; if(match($i,"/")){c++;hap=hap""tab[i]":"$i"\t"}}} ; print s >> out_hap ; if(c>=H){print $1"\t"c"\t"hap}}' $POLYMORPHES_LOCUS > $BAD_LOC

echo "nombre de locus ayant au moins $NB_INDIV fondateurs hétérozygote ou ambigue"
wc -l $BAD_LOC

# récriture du fichier d'haplotypes polymorphes et du fichier de population
POP_FILE2=`basename $POP_FILE| awk -v D=$FILTER_DIR '{print D"/"$1}'`
if [[ -e $POP_FILE2 ]]
then
rm $POP_FILE2
fi
awk -v BADL=$BAD_LOC -v BADI=$BAD_INDIV -v N="$POP_CONTROLE"  -v POP=$POP_FILE -v OUT=$POP_FILE2 -v TYPE=$ANALYSIS_TYPE 'BEGIN{FS="\t";
cpt_start=2; if (TYPE=="genotypes"){cpt_start=3};  
while(getline<BADI>0){tabBI[$1]=1; }
while(getline<POP>0){c++; if(match(N,$2) || tabBI[$1]==1){tabI[c+cpt_start]=1; }else{print $0 >> OUT}}
while(getline<BADL>0){tabBL[$1]=1}}
{if(tabBL[$1]!=1){s=""; 
g=0; 
for (i=cpt_start+1;i<=NF;i++){
if (tabI[i]!=1){s=s"\t"$i}
else{if($i!="-" && $i!="?"){g++}}
}
if(match($1,"Catalog ID")){
if(TYPE=="genotypes"){print $1"\t"$2"\t"$3""s}
else{print $1"\t"$2""s}
}else{
if(TYPE=="genotypes"){print $1"\t"$2-g"\t"$3""s}
else{print $1"\t"$2-g""s} 
} }}'  $POLYMORPHES_LOCUS > $FILTER_DIR/locus_polymorphes_filter.txt

POLYMORPHES_LOCUS=$FILTER_DIR/locus_polymorphes_filter.txt
POP_FILE=$POP_FILE2
#~ 
#~ 
#~ 
echo "nombre final de loci": `wc -l $POLYMORPHES_LOCUS | awk '{print $1-1}'`
echo "nombre final d'individus ": `wc -l $POP_FILE | awk '{print $1}'`
echo "detail par pop"
cut -f 2 $POP_FILE | sort | uniq -c | awk '{print "\tpop",$2,":",$1}'
fi

# calcul des stats genetique. Génération d'un tableau descriptif des allèles pour toutes les pop et tous les locus+ Generation d'un tableau descriptifs par pop des allèles et des génotypes+ génération d'un tableau descriptif des locus indiquant le nombre de genotypes abs/amb/>2allèles par pop.
prefix=`basename $POLYMORPHES_LOCUS | sed 's/.txt//' `
qsub -N ${Proj_Name}_genetics_stat -o $SGE/genetic_stat.out -e $SGE/genetic_stat.err $SCRIPT_DIR/genetics_stat.py -i $POLYMORPHES_LOCUS -p $POP_FILE -o $FILTER_DIR -n $prefix

