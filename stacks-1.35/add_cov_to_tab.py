#!/usr/bin/env python2.7
# -*-coding:Latin-1 -*

## IMPORT
import os
import glob
import gzip
import argparse
from datetime import datetime

# FUNCTION
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
    
def stacks_matches(indiv,matches_dir,locus_list):
    
    indiv_matches = glob.glob(os.path.join(matches_dir,indiv+".matches.tsv*"))
    if len(indiv_matches) != 1 :
        print "Error finding matches file for "+indiv+"\n";
        exit(1)
    
    if not is_gzip(indiv_matches[0]):
        handle_matches = open( indiv_matches[0] )
    else:
        handle_matches = gzip.open( indiv_matches[0] )
    
    dic_cov={}
    for line in handle_matches:
        #line = SQL_ID    batch_ID    locus_ID    sample_ID    sampl_locus_ID    genotype    coverage
        l=line.strip().split()
        if l[2] in locus_list :
            if not l[2] in dic_cov :
                dic_cov[l[2]]={l[5]:int(l[6])}
            else :
                dic_cov[l[2]][l[5]]=int(l[6])
            
    handle_matches.close()
    return dic_cov

# check de la cohérence entre les sorties matches et les sorties genotype
# list_alleles_matches : ['CGA','TGT']
# genotype : 'CA/TT'
# supprime la position commune à tous les matches par individus si elles sont exactes.
# faire la correspondance entre sortie matches et genotype
# si problème de correspondance alors mettre genotype à ?
def check_matches_genotype(list_alleles_matches, genotype):
	
	is_coherent = True
	alleles=genotype.split("/")

	for a in alleles:
		if not a in list_alleles_matches :
			is_coherent=False
			break
	#~ print "\tis_coherent" , is_coherent
	#~ if not is_coherent :
		#~ tmp_matches=dict()
		#~ rm_pos=list()
		#~ # recherche des positions communes entre genotype d'un individus
		#~ for pos in range(0,len(list_alleles_matches[0])):
			#~ a=list()
			#~ for allel in list_alleles_matches :
				#~ if allel[pos] not in a:
					#~ a.append(allel[pos])
			#~ if len(a) == 1:
				#~ rm_pos.append(pos)
		#print "\t\t",rm_pos
		#~ # creation d'un dictionnaire de correspondance entre matches et alleles sans les positions non variante dans 1 individus
		#~ for allel in list_alleles_matches:
			#~ a=""
			#~ for pos in range(len(allel)):
				#~ if not pos in rm_pos:
					#~ a=a+allel[pos]
			#~ tmp_matches[a]=allel
		#print "\t\t",tmp_matches
		#~ # correction des genotypes pour prendre en compte toutes les positions. Si une incohérence, alors genotype = "?"
		#~ for i in range(len(alleles)):
			#~ if alleles[i] in tmp_matches :
				#~ alleles[i]=tmp_matches[alleles[i]]
			#~ else:
				#~ genotype = "?"
				#~ break
				#~ 
		#~ if genotype != "?":
			#~ genotype = "/".join(alleles)
		#print "\t\t",genotype
	#~ 
	#~ return genotype
	return is_coherent

## MAIN
parser = argparse.ArgumentParser(description="add coverage to haplotype file")

# Inputs
group_input = parser.add_argument_group('Inputs')
group_input.add_argument('-i', '--input', default=None, required=True, help="transposed haplotype file")
group_input.add_argument('-d', '--matches-dir', default=None, required=True, help='Path to the sstacks matches files directory')

group_output = parser.add_argument_group('Outputs')
group_output.add_argument('-o', '--out', required=True, help='Path to output haplotype + allel coverage file')
args = parser.parse_args()

handle_in=open(args.input)
locus_list=list()
handle_out=open(args.out,"w")
 
for line in handle_in:
    if locus_list == [] :
        locus_list=line.strip().split("\t")[1:]
        handle_out.write(line)
        handle_out.write(handle_in.next())
    else:
        start=datetime.today()
        
        indiv=line.split("\t")[0]
        print indiv+" in treatment..."
        allelel_list=line.strip().split("\t")[1:]
        
        indiv_cov_list = stacks_matches(indiv,args.matches_dir,locus_list)
        
        for idx,genotype in enumerate(allelel_list):
			locus = locus_list[idx]
			#~ print locus
			#~ print "\t1",genotype
			if locus in indiv_cov_list : 
				#~ print "\t",indiv_cov_list[locus]
				#~ genotype = check_matches_genotype(indiv_cov_list[locus].keys(), genotype)
				is_coherent = check_matches_genotype(indiv_cov_list[locus].keys(), genotype)
			#~ print "\t2",genotype
			if genotype != "-" and genotype !="?" and is_coherent:
				alleles = genotype.split("/")      
				genotype_cov = "/".join(alleles)+":"+str(sum(indiv_cov_list[locus].values()))+":"+"/".join([str(indiv_cov_list[locus][a]) for a in alleles])
				#~ print "\t\t",genotype_cov
				allelel_list[idx]=genotype_cov
			else:
				allelel_list[idx]=genotype
        handle_out.write(indiv+"\t"+"\t".join(allelel_list)+"\n")
        elapsed_time = datetime.today() - start
        print "\t## ended in ",elapsed_time
                
                
handle_in.close()
handle_out.close()                





