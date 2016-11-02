#!/usr/bin/env python2.7
# -*- coding: utf-8 -*- 

############################# IMPORT ###################################

from Bio import SeqIO
from Bio.Sequencing import Ace
import os, sys 
import logging
import re
import argparse
import glob
############################# FUNCTION ###################################
def error(message):
	print message
	usage()
	sys.exit(2)

def revcom(s):
     basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
     letters = list(s)
     letters = [basecomplement[base] for base in letters[::-1]]
     return ''.join(letters)
     

def fasta_enz_search(locus_fastafile,locus,enzyme,dic_out) :
	dic_out[locus]={"nb_contig":0, "tot_read":"?", "assembly":{}}
	
	handle = open(locus_fastafile, "rU")
	for record in SeqIO.parse(handle, "fasta") :
		dic_out[locus]["nb_contig"]+=1
		start=record.seq.find(enzyme,0 , 20)
		if start==-1:
			start=record.seq.find(revcom(enzyme),len(record.seq)-20)
			if start >-1 :
				start+=len(enzyme)-1
		
		dic_out[locus]["assembly"][record.id]={"length":len(record.seq),"nb_read":"?","enz":start,"read1":"?","read2":"?"}
		
	handle.close()

	maxlen=0
	tab=[]
	for contig in sorted(dic_out[locus]["assembly"]) :
		tab.append(contig+"\t"+`dic_out[locus]["assembly"][contig]["enz"]`+"\t"+`dic_out[locus]["assembly"][contig]["length"]`+"\t"+dic_out[locus]["assembly"][contig]["nb_read"]+"\t"+dic_out[locus]["assembly"][contig]["read1"]+"\t"+dic_out[locus]["assembly"][contig]["read2"]+"\n")
		
		if int(dic_out[locus]["assembly"][contig]["length"]) > maxlen:
			maxlen=int(dic_out[locus]["assembly"][contig]["length"])
			
	header=locus+"\t"+`dic_out[locus]["tot_read"]`+"\t"+str(maxlen)+"\t?\t"
	dic_out[locus]["bilan"]=header+header.join(tab)		

def ace_enz_search(locus_acefile,locus,enzyme,dic_out) :
	
	handle=open(locus_acefile)
	contigs = Ace.parse(handle)
	dic_out[locus]={"nb_contig":0, "tot_read":0, "assembly":{}}
	for contig in contigs :
		read1=0
		read2=0
		dic_out[locus]["nb_contig"]+=1
		dic_out[locus]["tot_read"]+=contig.nreads
		# contig length sans gap non reporté dans le fasta
		length=len(contig.sequence)-contig.sequence.count("*")
		# enzyme start
		start=contig.sequence.replace("*","").find(enzyme,0,20)
		if start==-1:
			start=contig.sequence.replace("*","").find(revcom(enzyme),length-20)
			if start >-1 :
				start+=len(enzyme)-1
		# nb read1 et read 2
		for r in contig.af:
			if r.name.endswith(".1"):
				read1+=1
			elif r.name.endswith(".2"):
				read2+=1
				
		dic_out[locus]["assembly"]["Locus_"+locus+"_"+contig.name.replace("Contig","Contig_")]={"length":length,"nb_read":str(contig.nreads),"enz":start,"read1":str(read1),"read2":str(read2)}
	handle.close()
	if len(dic_out[locus]["assembly"].keys())==0:
		logging.warning("No assembly: "+locus+"less than 6 reads per contig")
	else:
		maxlen=0
		maxread=0
		tab=[]
		for contig in sorted(dic_out[locus]["assembly"]) :
			tab.append(contig+"\t"+`dic_out[locus]["assembly"][contig]["enz"]`+"\t"+`dic_out[locus]["assembly"][contig]["length"]`+"\t"+dic_out[locus]["assembly"][contig]["nb_read"]+"\t"+dic_out[locus]["assembly"][contig]["read1"]+"\t"+dic_out[locus]["assembly"][contig]["read2"]+"\n")
			
			if int(dic_out[locus]["assembly"][contig]["length"]) > maxlen:
				maxlen=int(dic_out[locus]["assembly"][contig]["length"])
			if int(dic_out[locus]["assembly"][contig]["nb_read"]) > maxread:
				maxread=int(dic_out[locus]["assembly"][contig]["nb_read"])
				
		header=locus+"\t"+`dic_out[locus]["tot_read"]`+"\t"+str(maxlen)+"\t"+str(maxread)+"\t"
		dic_out[locus]["bilan"]=header+header.join(tab)
						

		
############################# MAIN ###################################
parser = argparse.ArgumentParser(description="compute cap3/velvet assembly statistics")
parser.add_argument('-e', '--enzyme', default=None, required=True, help="troncated enzyme sequence")
#~ parser.add_argument('-f', '--format', default="CAP3", choices=['CAP3','Velvet'], help='CAP3 or Velvet')
parser.add_argument('-l', '--log', default="assembly_selection.log", help='Path to the assembly selection log file')

# Inputs
group_input = parser.add_argument_group('Inputs')
group_input.add_argument('-w', '--white-list', required=True, help='Path to the locus file list')
group_input.add_argument('-a', '--assembly-dir', required=True, help='Path to the assembly directory (containing one directory per loci')
group_output = parser.add_argument_group('Outputs')
group_output.add_argument('-o', '--output-dir', required=True, help='Path to output directory')
args = parser.parse_args()
		

logging.basicConfig(filename=args.log,level=logging.DEBUG)

# enzyme start in assembly
logging.info("#### Researche enzyme : "+args.enzyme+", start site in assembly ###\n")
logging.info("...\n")
enz_start_dic={}

for locus in open(args.white_list) : 
	if os.path.exists(args.assembly_dir+"/"+locus.strip()+"/"+locus.strip()+".fq.fasta.cap.ace"):
			ace_enz_search(args.assembly_dir+"/"+locus.strip()+"/"+locus.strip()+".fq.fasta.cap.ace",locus.strip(),args.enzyme,enz_start_dic)
	else:
		fasta_list=[n for n in glob.glob(args.assembly_dir+"/"+locus.strip()+"/"+locus.strip()+"*.fa") if os.path.getsize(n) > 0 ]
		for fasta in fasta_list :
			fasta_enz_search(fasta,locus.strip(),args.enzyme,enz_start_dic)
FO=open(args.output_dir+"/assembly.stat","w")
FO.write("#locus\ttot_input_read\tcontig_max_length\tcontig_max_read\tcontig\tstart_enz\tlength\tcontig_read\tnb_read1\tnb_read2\n")
for locus in enz_start_dic:
	FO.write(enz_start_dic[locus]["bilan"])
	
#~ if(args.format=="VELVET"):
	#~ fasta_enz_search(args.white_list,args.assembly_dir,args.enzyme,enz_start_dic)
	#~ print enz_start_dic
#~ elif(args.format=="CAP3"):
	#~ for locus in open(args.white_list) : 
		#~ ace_enz_search(args.assembly_dir+"/"+locus.strip()+"/"+locus.strip()+".fq.fasta.cap.ace",locus.strip(),args.enzyme,enz_start_dic)
	#~ logging.info("#### WRITE args.enzyme : "+args.enzyme+", START POS : "+args.output_dir+"/assembly.stat ###\n")
	#~ logging.info("...\n")
	#~ FO=open(args.output_dir+"/assembly.stat","w")
	#~ FO.write("#locus\ttot_input_read\tcontig_max_length\tcontig_max_read\tcontig\tstart_enz\tlength\tcontig_read\tnb_read1\tnb_read2\n")
	#~ for locus in enz_start_dic:
		#~ FO.write(enz_start_dic[locus]["bilan"])
