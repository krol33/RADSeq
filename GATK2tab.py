#!/usr/local/bioinfo/src/python/current/bin/python
# -*- coding: utf-8 -*- 

############################# IMPORT ###################################

import os, sys 
import logging
from itertools import groupby
import argparse
import gzip
import glob
import pysam

############################# FUNCTION ###################################
def is_gzip( file ):
    """
    @return: [bool] True if the file is gziped.
    @param file : [str] Path to processed file.
    """
    is_gzip = None
    FH_input = gzip.open( file )
    try:
        FH_input.readline()
        is_gzip = True
    except:
        is_gzip = False
    finally:
        FH_input.close()
    return is_gzip

def revcom(s):
     basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
     letters = list(s)
     letters = [basecomplement[base] for base in letters[::-1]]
     return ''.join(letters)
     
def parse_bam(bam):
	bamfile = pysam.AlignmentFile(bam, "rb")
	aln={}
	for read in bamfile.fetch() :
		if read.is_read2 :
			continue
		strand = "R" if read.is_reverse else "F"
		# premiere position +1 car pysam enleve 1 pour avoir une indexation qui commence à 0. Or dans le VCF la position est en index qui commence à 1.
		if strand == "R":
			start = read.get_reference_positions(full_length=True)[-1]+1 if not read.get_reference_positions(full_length=True)[-1] == None else None
		else:
			start = read.get_reference_positions(full_length=True)[0]+1 if not read.get_reference_positions(full_length=True)[0] == None else None
		ref = bamfile.references[read.reference_id]
		if not start == None :
			if not ref in aln:
				aln[ref]=[start]
			elif not start in aln[ref]:
				aln[ref].append(start)
	# aln={'scaffold_654': [341066]}
	return aln

def record_stacks_snp (stacks_file, locus, stacks_list) :
			
	if not is_gzip(stacks_file):
		handle = open( stacks_file )
	else:
		handle = gzip.open( stacks_file )
	
	for line in handle.readlines():
		l=line.strip().split()
		if l[2] == locus :
			# stacks 1.19 : 0	1	1	23	0	A	T
			# stacks 1.35 : 0	2	1	53	E	0	T	G	-	- ==> jusqu'à 4 alleles
			allele=l[5:] if len(l)==7 else [n for n in l[6:] if n!="-"]
			allele.sort()
			stacks_list.append(str(int(l[3]))+"_"+"_".join(allele))
	
	# result
	# ['81_A_G']

def parse_GATK(vcf, locus, aln_dic, stacks_list):
	# locus
	# 100011
	# aln_dic
	# {'scaffold_654': [341066]}
	# stacks_list
	# ['81_A_G']

	if not is_gzip(vcf):
		handleIN = open( vcf )
	else:
		handleIN = gzip.open( vcf )

	#~ snp_dic={"position":"","equivalent_stacks":"", "ref":"", "alt":"" }
	snp_dic={}
	correspondance_dic={}
	snp_count=0
	for line in handleIN.readlines():
		if not line.startswith("##"):
			snp_list=line.strip().split()
			# enregistrement de l'ordre des individus
			if line.startswith("#CHROM") :	
				i=9
				for indiv in snp_list[9:] :
					correspondance_dic[i]=indiv
					i+=1
			
			snp_filter=snp_list[6]
			if snp_filter =="." or snp_filter =="PASS":
				snp_count+=1
				seq=snp_list[0]
				pos=int(snp_list[1])
				if not seq in snp_dic :
					snp_dic[seq]={pos:{"equivalent_stacks":"", "ref":"", "alt":"", "haplotype" :{} }}
				else:
					snp_dic[seq][pos]={"equivalent_stacks":"", "ref":"", "alt":"", "haplotype" :{}}
					
				# enregistrement des genotypes
				allele=snp_list[4].split(",")
				allele.append(snp_list[3])
				allele.sort()
				#recherche du snp qui correspond à la référence stacks
				equivalent_stacks="?"
				if seq in aln_dic:
					for aln in aln_dic[seq] : 
						convert_start = pos-aln
						a=allele
						if convert_start < 0:
							a=[revcom(al) for al in a]
							a.sort()
						#~ print str(abs(convert_start))+"_"+"_".join(a) ,"in", stacks_list
						if str(abs(convert_start))+"_"+"_".join(a) in stacks_list:
							equivalent_stacks=str(abs(convert_start))

				snp_dic[seq][pos]["equivalent_stacks"]=equivalent_stacks
				
				#enregistrement des haplotypes	
				snp_dic[seq][pos]["ref"]=snp_list[3]
				snp_dic[seq][pos]["alt"]=snp_list[4]
				
				i=9
				for genotype in snp_list[9:] :
					indiv=correspondance_dic[i]
					gen_list=genotype.split(":")
					# genotype inconnu ou avec genotype quality <=9					# genotype quality = différence de Phred score entre le meilleur (toujours 0) et le second meilleur (échelle de 0 à 99 puisqu'au dessus de 99 on est sûre du genotype donc pas la peine d'aller plus loin.
					if gen_list[0] == "./." or int(gen_list[3]) <= 9:									
						snp_dic[seq][pos]["haplotype"][indiv]="."
					
					# genotype hetero avec genotype quality > 9
					elif gen_list[0] == "0/1" or gen_list[0]=="1/0" :								
						snp_dic[seq][pos]["haplotype"][indiv]="H"
					
					# genotype avec genotype quality > 9
					else :																			
						snp_dic[seq][pos]["haplotype"][indiv]=genotype[0]				
					i+=1
	handleIN.close()
	
	handleOUT=open(os.path.join(os.path.dirname(vcf),locus+".tsv") ,"w")
	string="locus\tseq_ref\tpositions\tequivalent_stacks\tallel_ref\tallel_alt\thaplotype_count\t"+"\t".join(sorted(correspondance_dic.values()))+"\n"
	if len(snp_dic) == 0 :
		string=string+locus+"\tno_snp_found"
	else :
		for seq in snp_dic :
			for pos in snp_dic[seq]:
				string=string+locus+"\t"+seq+"\t"+str(pos)+"\t"+snp_dic[seq][pos]["equivalent_stacks"]+"\t"+snp_dic[seq][pos]["ref"]+"\t"+snp_dic[seq][pos]["alt"]+"\t"
				hap_count=""
				for hap,group in groupby(sorted(snp_dic[seq][pos]["haplotype"].values())):
					if hap_count=="" :
						hap_count=hap+":"+str(len(list(group)))
					else:
						hap_count+=";"+hap+":"+str(len(list(group)))
				string=string+hap_count
				for indiv in sorted(correspondance_dic.values()):
					string=string+"\t\""+snp_dic[seq][pos]["haplotype"][indiv]+"\""
				string=string+"\n"
	handleOUT.write(string)
	handleOUT.close()
	

############################# MAIN ###################################
parser = argparse.ArgumentParser(description="")
group_input = parser.add_argument_group('Inputs')
group_input.add_argument('-v', '--vcf', default=None, required=True, help='GATK snp calling on read 1 (vcf[.gz])')
group_input.add_argument('-l', '--locus', default=None, required=True, help='Locus id')
group_input.add_argument('-c', '--catalog-snp', default=None, required=True, help='Stacks catalog SNP file (gz or not)')
group_input.add_argument('-b', '--bam', default=None, required=True, help='Read 1 alignment (bam)')

args = parser.parse_args()

dic_aln = parse_bam(args.bam)
print "dic_aln",dic_aln
snp_stacks=[]
record_stacks_snp(args.catalog_snp, args.locus, snp_stacks)
print "snp_stacks", snp_stacks

parse_GATK(args.vcf,args.locus,dic_aln,snp_stacks )
