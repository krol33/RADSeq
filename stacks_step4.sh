#$ -o [Project_path_dir]/SGE_out/stacks4.out
#$ -e [Project_path_dir]/SGE_out/stacks4.err
#$ -N [Proj_Name]_stacks4
#$ -S /bin/sh

#$ -q unlimitq

#$ -m bea

date
date >&2
Proj_Name=[Proj_Name]

# VARIABLES DEDIEES
# type d'analyse, croisement (ANALYSIS_TYPE="genotypes") ou populations (ANALYSIS_TYPE="populations")
ANALYSIS_TYPE="populations"

# VARIABLES GENERALES
stacks_dir=/usr/local/bioinfo/src/Stacks/stacks-1.44/bin
RAD_DIR=[Project_path_dir]
SCRIPT_DIR=[Script_path_dir]

# INPUT
STACKS1_DIR=$RAD_DIR/stacks1
CSTACKS_DIR=$RAD_DIR/cstacks
SSTACKS_DIR=$RAD_DIR/sstacks
INDIV_FILE=$RAD_DIR/indiv_barcode.txt
dos2unix $INDIV_FILE
POP_FILE=$RAD_DIR/population.map
dos2unix $POP_FILE

# si besoins de générer le tableau d'haplotypage uniquement sur une selection de population en particulier indiquez entre "" les numéros de population séparer par un espace.
EXTRACT_POP=""
# genotypes/populations options  : 
# pour populations : 
# 	autre que -P -M -b -k -t --fstats (notamment les options de format de sortie  --genepop --fasta --structure --phase --beagle --plink --phylip --vcf ..., de  boostrap, de fusion et phasage, de filtre (attention à la version de stacks utilisé)
# pour genotypes : 
# 	autre que -b 1 -P $input -t GEN -s -o joinmap
OPT="--genepop --fasta "

# SGE option
MEM="16G"
VMEM="32G"
THREADS=10

# output
OUT_DIR=$RAD_DIR/${ANALYSIS_TYPE}
SGE=$RAD_DIR/SGE_out/`basename $OUT_DIR`
STAT_DIR=$RAD_DIR/stat/`basename $OUT_DIR`

mkdir -p $SGE $OUT_DIR $STAT_DIR

###########################################################################################

if [[ -n $EXTRACT_POP ]]
then
	if [ -e $OUT_DIR/pop_extract.map ]
	then 
	rm $OUT_DIR/pop_extract.map
	fi
	for p in $EXTRACT_POP
	do
		awk -v P=$p '$2==P' $POP_FILE >> $OUT_DIR/pop_extract.map
	done
	POP_FILE=$OUT_DIR/pop_extract.map
	cat $POP_FILE | while read line
	do
		indiv=`echo $line | awk '{print $1}'`
		ln -s $STACKS1_DIR/${indiv}* $OUT_DIR/
		ln -s $SSTACKS_DIR/${indiv}* $OUT_DIR/
	done
	ln -s $CSTACKS_DIR/* $OUT_DIR/
else
	ln -s $STACKS1_DIR/* $OUT_DIR/
	ln -s $CSTACKS_DIR/* $OUT_DIR/
	ln -s $SSTACKS_DIR/* $OUT_DIR/
fi

echo qsub -N ${Proj_Name}_$ANALYSIS_TYPE -o $SGE/$ANALYSIS_TYPE.out -e $SGE/$ANALYSIS_TYPE.err -l h_vmem=$VMEM -l mem=$MEM -pe parallel_smp $THREADS $SCRIPT_DIR/$ANALYSIS_TYPE.sh $OUT_DIR $stacks_dir $POP_FILE $THREADS \"$OPT\"
id=`qsub -b y -terse -N ${Proj_Name}_$ANALYSIS_TYPE -o $SGE/$ANALYSIS_TYPE.out -e $SGE/$ANALYSIS_TYPE.err -l h_vmem=$VMEM -l mem=$MEM -pe parallel_smp $THREADS $SCRIPT_DIR/$ANALYSIS_TYPE.sh $OUT_DIR $stacks_dir $POP_FILE $THREADS \"$OPT\"| awk -F "." '{print $1}'`


r=0;
while [ $r == 0 ] ; do qstat -j $id >& /dev/null ; r=$? ; sleep 10s ; done 

#############################################################################################
###########################   SSTACK stat : $OUT_DIR_DIR ######################################################
echo -e "\n###########################   SSTACK stat : $OUT_DIR ######################################################"
echo -e "\n###########################   SSTACK stat : $OUT_DIR ######################################################" >&2 

python $SCRIPT_DIR/stacks_summary.py --stacks-prog $ANALYSIS_TYPE --res-dir $OUT_DIR --pop-map $POP_FILE --summary $STAT_DIR/${ANALYSIS_TYPE}_summary.html

mkdir $OUT_DIR/input
mv $OUT_DIR/*tags.tsv* $OUT_DIR/input
mv $OUT_DIR/*snps.tsv* $OUT_DIR/input
mv $OUT_DIR/*alleles.tsv* $OUT_DIR/input
mv $OUT_DIR/*matches.tsv* $OUT_DIR/input
mv $OUT_DIR/*models.tsv* $OUT_DIR/input

# recherche des génotypes ambigues, remplace "-" par "?"
if [[ $ANALYSIS_TYPE == "populations" ]]
then
python $SCRIPT_DIR/retreat_haplotypes.tsv.py -i $OUT_DIR/batch_1.haplotypes.tsv -s $OUT_DIR/input -o $OUT_DIR/batch_1.haplotypes_amb.tsv
else
python $SCRIPT_DIR/retreat_haplotypes.tsv.py -i $OUT_DIR/batch_1.haplotypes_1.tsv -s $OUT_DIR/input -o $OUT_DIR/batch_1.haplotypes_amb.tsv
fi

# réordonnancement du tableaux en fonction de l'ordre fourni dans le fichier map
head -n 1 $OUT_DIR/batch_1.haplotypes_amb.tsv | awk -v OUT=$OUT_DIR/batch_1.haplotypes_amb_ord.tsv -v IN=$OUT_DIR/batch_1.haplotypes_amb.tsv -v POP=$POP_FILE '{l=split($0,tab,"\t");
for(i=1;i<=l;i++){
idx[tab[i]]=i
}
first=1}END{
while(getline<POP>0){
if (first==1){
system("cut -f 1,2,"idx[$1]" "IN" > "OUT)
first=0
}else{
system( "cut -f "idx[$1]" "IN" | paste "OUT" - > "OUT"tmp; mv "OUT"tmp "OUT)
}
}}'
rm $OUT_DIR/batch_1.haplotypes_amb.tsv

# séparation des locus monomorphes des locus polymorphes
grep -E 'Catalog|consensus' $OUT_DIR/batch_1.haplotypes_amb_ord.tsv > $OUT_DIR/locus_monomorphes.txt
	# pour les locus polymorphes ajout des couvertures total et allelique à chaque genotype individuel
grep -v consensus $OUT_DIR/batch_1.haplotypes_amb_ord.tsv > $OUT_DIR/tmp_poly
datamash transpose < $OUT_DIR/tmp_poly > $OUT_DIR/tmp_poly_transposed
rm $OUT_DIR/tmp_poly
python $SCRIPT_DIR/add_cov_to_tab.py -i $OUT_DIR/tmp_poly_transposed -d $OUT_DIR/input -o $OUT_DIR/tmp_poly_cov_transposed
if [[ $? == 0 ]]
then
        rm $OUT_DIR/tmp_poly_transposed
        datamash transpose < $OUT_DIR/tmp_poly_cov_transposed > $OUT_DIR/locus_polymorphes_cov.txt
        rm $OUT_DIR/tmp_poly_cov_transposed
        POLY=$OUT_DIR/locus_polymorphes_cov.txt
else
        rm $OUT_DIR/tmp_poly_transposed $OUT_DIR/tmp_poly_cov_transposed
        grep -v consensus $OUT_DIR/batch_1.haplotypes_amb_ord.tsv > $OUT_DIR/locus_polymorphes.txt
        POLY=$OUT_DIR/locus_polymorphes.txt
fi

# génération du fichier fasta des locus polymorphes (séquence consensus du catalog)
zcat $OUT_DIR/input/batch_1.catalog.tags.tsv.gz | awk -v L=$POLY 'BEGIN{while(getline<L>0){tab[$1]=1}}{if(tab[$3]==1){print ">Locus_"$3; print $9 }}'> $OUT_DIR/locus_polymorphes.fa

#		statistique excell resumé 
# allel_summarized : #Locus	nb_SNP	nb_Allel	code_Allel	nb_indiv_genotyped
# allel_detailed :allel_summarized+	pop1_nb_indiv_genotyped	pop1nb_allele_diff	pop1_nb_a1	pop1_nb_aXX	...	pop1_nb_a5	pop2_nb_indiv_genotyped	pop2nb_allele_diff	pop2_nb_a1	pop2_nb_aX	...	pop2_nb_a5
# genotype_detailed :allel_summarized+	pop1_nb_indiv_genotyped	pop1_nb_genotype_diff	pop1_nb_g1/1	pop1_nb_g1/2	...	pop1_nb_g2/2	pop1_nb_g2/3	...	pop1_nb_g5/5
python $SCRIPT_DIR/genetics_stat.py -d -g -i $POLY -p $POP_FILE -o $OUT_DIR -n all_locus_polymorphes

