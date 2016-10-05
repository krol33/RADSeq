#!/bin/sh

input=$1
bin_dir=$2
map=$3
script_dir=$4
libR=$5
STAT_DIR=$6
OPT=$7

echo $bin_dir/populations -P $input -M $map -b 1 -k -t 10 --fstats $OPT
$bin_dir/populations -P $input -M $map -b 1 -k -t 10 --fstats $OPT

echo "Stacks populations finished" >&2

mkdir $input/input
mv $input/*tags.tsv* $input/input
mv $input/*snps.tsv* $input/input
mv $input/*alleles.tsv* $input/input
mv $input/*matches.tsv* $input/input
NB_IND=`ls $input/input | grep -c matches.tsv`

# recherche des génotypes ambigues, remplace "-" par "?"
python $script_dir/retreat_haplotypes.tsv.py -i $input/batch_1.haplotypes.tsv -s $input/input -o $input/batch_1.haplotypes_amb.tsv

# réordonnancement du tableaux en fonction de l'ordre fourni dans le fichier map
head -n 1 $input/batch_1.haplotypes_amb.tsv | awk -v OUT=$input/batch_1.haplotypes_amb_ord.tsv -v IN=$input/batch_1.haplotypes_amb.tsv -v POP=$map '{l=split($0,tab,"\t");
for(i=1;i<=l;i++){
	idx[tab[i]]=i
	first=1
}}END{
	while(getline<POP>0){
		if (first==1){
			system("cut -f 1,2,"idx[$1]" "IN" > "OUT)
			first=0
		}else{
			system( "cut -f "idx[$1]" "IN" | paste "OUT" - > "OUT"tmp; mv "OUT"tmp "OUT)
		}
	}
}'
rm $input/batch_1.haplotypes_amb.tsv

# calcul du taux de Hz/Ho zygotie par individus (#pop_num	indiv	num_homo	num_hetero	num_ambiguous	num_unknown)
sh $script_dir/zygosity_rate.sh $input/batch_1.haplotypes_amb_ord.tsv populations $map $STAT_DIR/zygosity_rate.txt
echo "taux d'homozygotie" > $STAT_DIR/summary_statistics.txt
grep -v '#' $STAT_DIR/zygosity_rate.txt  | awk '{sum=$3+$4; printf("%s\t%.2f\n", $2, $3*100/sum)}' | colstat.sh 2 >> $STAT_DIR/summary_statistics.txt
echo "taux d'hheterozygotie" >> $STAT_DIR/summary_statistics.txt
grep -v '#' $STAT_DIR/zygosity_rate.txt  | awk '{sum=$3+$4; printf("%s\t%.2f\n", $2, $4*100/sum)}' | colstat.sh 2 >> $STAT_DIR/summary_statistics.txt

# séparation des locus monomorphes des locus polymorphes
grep -E 'Catalog|consensus' $input/batch_1.haplotypes_amb_ord.tsv > $input/locus_monomorphes.txt
	# pour les locus polymorphes ajout des couvertures total et allelique à chaque genotype individuel
grep -v consensus $input/batch_1.haplotypes_amb_ord.tsv > $input/tmp_poly
datamash transpose < $input/tmp_poly > $input/tmp_poly_transposed
rm $input/tmp_poly
python $script_dir/add_cov_to_tab.py -i $input/tmp_poly_transposed -d $input/input -o $input/tmp_poly_cov_transposed
rm $input/tmp_poly_transposed 
datamash transpose < $input/tmp_poly_cov_transposed > $input/locus_polymorphes_cov.txt
rm $input/tmp_poly_cov_transposed

# génération du fichier fasta des locus polymorphes (séquence consensus du catalog)
zcat $input/input/batch_1.catalog.tags.tsv.gz | awk -v L=$input/locus_polymorphes_cov.txt 'BEGIN{while(getline<L>0){tab[$1]=1}}{if(tab[$3]==1){print ">Locus_"$3; print $9 }}'> $input/locus_polymorphes.fa

# statistiques sur le nombre de locus genotypé, genotypé pour au moins 50% des individus
wc -l $input/batch_1.haplotypes_amb_ord.tsv | awk '{print "nombre de locus genotypé: "$1-1}'   >> $STAT_DIR/summary_statistics.txt
grep -v Catalog $input/batch_1.haplotypes_amb_ord.tsv | awk -v TOT=$NB_IND '{if($2>TOT/2){n++}}END{print n" locus genotypés pour au moins 50% des individus"}'  >> $STAT_DIR/summary_statistics.txt
# 		détail sur les locus polymorphes
grep -v Catalog $input/locus_polymorphes_cov.txt | awk -v TOT=$NB_IND '{t++; if($2>TOT/2){n++}}END{print "nombre de locus polymorphes: "t; print n" locus polymorphes genotypés pour au moins 50% des individus"}'  >> $STAT_DIR/summary_statistics.txt
#		statistique excell resumé 
# allel_summarized : #Locus	nb_SNP	nb_Allel	code_Allel	nb_indiv_genotyped
# allel_detailed :allel_summarized+	pop1_nb_indiv_genotyped	pop1nb_allele_diff	pop1_nb_a1	pop1_nb_aXX	...	pop1_nb_a5	pop2_nb_indiv_genotyped	pop2nb_allele_diff	pop2_nb_a1	pop2_nb_aX	...	pop2_nb_a5
# genotype_detailed :allel_summarized+	pop1_nb_indiv_genotyped	pop1_nb_genotype_diff	pop1_nb_g1/1	pop1_nb_g1/2	...	pop1_nb_g2/2	pop1_nb_g2/3	...	pop1_nb_g5/5
python $script_dir/genetics_stat.py -d -g -i $input/locus_polymorphes_cov.txt -p $map -o $input -n all_locus_polymorphes
awk '{if($2==1 && $3==2){c++}}END{print c,"locus polymorphe à 1 SNP 2 alleles" }' $input/all_locus_polymorphes.alleles_detailed_stat.tsv  >> $STAT_DIR/summary_statistics.txt

wc -l $input/locus_monomorphes.txt | awk '{print "nombre de locus monomorphes : "$1-1}'   >> $STAT_DIR/summary_statistics.txt

awk -v TOT=$NB_IND '{if($1!="Catalog"){tab[$2]++}}END{for (i=1; i<= TOT; i++){if(tab[i]!=""){print i"\t"tab[i]}else{print i"\t0" }}}' $input/batch_1.haplotypes_amb_ord.tsv > $STAT_DIR/nb_indiv_genotyped_per_locus.txt
awk -v TOT=$NB_IND '{if($1!="Catalog"){tab[$2]++}}END{for (i=1; i<= TOT; i++){if(tab[i]!=""){print i"\t"tab[i]}else{print i"\t0" }}}' $input/locus_polymorphes_cov.txt > $STAT_DIR/nb_indiv_genotyped_per_locus_polymorphes.txt
echo -e "\nsee detailed number of cluster in function of number of genotyped individuals available: " $STAT_DIR/nb_indiv_genotyped_per_locus.txt", "$STAT_DIR/nb_indiv_genotyped_per_locus_polymorphes.txt " and "$STAT_DIR/nb_indiv_genotyped_per_locus.jpg  >> $STAT_DIR/summary_statistics.txt
$libR/plot_dist.R "nb cluster in function of number of genotyped individuals" "number of genotyped individuals" "number of cluster" 1 $STAT_DIR/nb_indiv_genotyped_per_locus $STAT_DIR/nb_indiv_genotyped_per_locus.txt $STAT_DIR/nb_indiv_genotyped_per_locus_polymorphes.txt
