#!/usr/bin/python
# -*- coding: utf-8 -*- 
########################################################################
## IMPORT
import sys, getopt, os,glob
from os import listdir
from os.path import isfile, join, islink
import gc
import gzip

########################################################################
## FUNCTIONS
def usage():
    print 'retreat_haplotypes.tsv.py -i <input-haplotype> -o <output-haplotype> -s <sstacks_results_dir>'
    sys.exit(2)
def research_ambiguous_haplotypes(sstacks_dir,dic_out):
    
    for f in glob.glob(sstacks_dir+'/*.matches.tsv*') :
        name=os.path.basename(f).split(".matches")[0]
        dic_out[name]={}
        tmpd={}
        
        if f.endswith(".gz"):
            handle = gzip.open( f, "r" )
        else :
            handle=open(f,"r")
            
        for l in handle.readlines():
            tmplist=l.split()
            locus=tmplist[2]
            if not locus in tmpd:
                tmpd[locus]=1
            else:
                tmpd[locus]+=1
                if not locus in dic_out[name] :
                    dic_out[name][locus]=2
                else :
                    dic_out[name][locus]+=1
        
        tmpd.clear()
        gc.collect()    
                

def rewrite_happlotypes(infile,dic_amb,outfile):
    FO=open(outfile,"w")
    
    FI=open(infile,"r")
    header=FI.readline().strip()
    nameOut=header.split("\t")[1:]
    c=0
    
    string_out=header+"\n"
    i=1
    
    for l in FI.readlines():
        tmplist=l.strip().split()
        locus=tmplist[0]
        s=locus
        for i in xrange(0,len(nameOut)):
            sample=nameOut[i]
            if not sample in dic_amb:                     # colonne de comptage d'individus génotypé (et de fréquence pour le programme genotype)
                s+="\t"+tmplist[i+1] 
            elif tmplist[i+1]!="-":                        # si genotype connu on le conserve
                s+="\t"+tmplist[i+1] 
            else : 
                if locus in dic_amb[sample] :
                    s+="\t?"
                else :
                    s+="\t-"
        string_out+=s+"\n"
        i+=1
        if i==2000:
            c+=1
            FO.write(string_out)
            string_out=""
            i=0
            print "\trewrite haplotypes for  "+`c*2000`+" locus "
    
    if not string_out == "":
        print "\trewrite haplotypes for  "+`c*2000+i`+" locus"
        FO.write(string_out)
        string_out=""
        i=0

# lecture du fichier d'haplotypage de Stacks et enregistrement dans un dictionnaire
# sortie = dict[locus][sample]=haplotypes dict[locus][Cnt]=nb_indiv_haplotypé
def read_haplotypes(infile,dicOut):
    nameOut=[]
    FI=open(infile,"r")
    header=FI.readline().strip()
    nameOut=header.split("\t")[1:]
    for locus in FI.readlines():
        tmplist=locus.split()
        dicOut[int(tmplist[0])]=dict(zip(nameOut,tmplist[1:]))

    FI.close()
    #~ print nameOut
    #~ print dicOut.keys()
    #~ print dicOut["10258"]
    return nameOut

# lectures des fichiers matches.tsv (sorties sstacks) et enregistrement du nombre de cluster d'un individus correspondant à un locus du catalog
# sortie = dict[sample][locus]=nb cluster
def read_sstacks(indir,dicOut) :
    #~ print [ f for f in listdir(indir) if isfile(join(indir,f))]
    matchesfiles = [ f for f in listdir(indir) if (isfile(join(indir,f)) or islink(join(indir,f))) and f.endswith(".matches.tsv") ]
    #~ print matchesfiles
    for f in matchesfiles:
        FI=open(os.path.realpath(indir+"/"+f),"r")
        sample=f.replace(".matches.tsv","")
        d={}
        for l in FI.readlines():
            tmplist=l.split()
            locus=int(tmplist[2])
            if not locus in d:
                d[locus]=1
            else:
                d[locus]+=1
        dicOut[sample]=d
        FI.close()
    #~ print dicOut.keys()
    #~ print dicOut["EFFICACE_11_1"]["10258"]

#~ def research_ambiguities_haplotypes(haplotypes_dict, header, matches_dict, outfile) :
    #~ FO=open(outfile,"w")
#~ 
    #~ FO.write("Catalog ID\t"+"\t".join(header)+"\n")
    #~ string=""
    #~ i=0
    #~ amb=0
    #~ c=0
    #~ for locus in sorted(haplotypes_dict.keys()):
        #~ #print c*2000+i
        #~ s=`locus`
        #~ for sample in header :
            #~ if not sample in matches_dict.keys():                # colonne de comptage d'individus génotypé (et de fréquence pour le programme genotype)
                #~ s+="\t"+haplotypes_dict[locus][sample]
            #~ else:
                #~ if not haplotypes_dict[locus][sample]=="-" :    # si genotype connu on le conserve
                    #~ s+="\t"+haplotypes_dict[locus][sample]
                #~ else:                                           # genotype inconnu pour l'individus
                    #~ if locus in matches_dict[sample].keys():        #inconnu car ambigue
                        #~ s+="\t?"
                        #~ amb+=1
                    #~ else :                                          # reellement inconnu
                        #~ s+="\t-"
        #~ string+=s+"\n"
        #~ i+=1
        #~ if i==2000:
            #~ c+=1
            #~ FO.write(string)
            #~ string=""
            #~ i=0
            #~ print "\tresearch ambiguous haplotypes for  "+`c*2000`+" locus on "+`len(haplotypes_dict.keys())`
    #~ 
    #~ if not string == "":
        #~ print "\tresearch ambiguous haplotypes for  "+`c*2000+i`+" locus on "+`len(haplotypes_dict.keys())`
        #~ FO.write(string)
        #~ string=""
        #~ i=0
    #~ print "Found "+`amb`+" ambiguous haplotypes\n"
    #~ FO.close()
    
########################################################################
## MAIN
SSTACKS_DIR = ''
HAPLOTYPES_IN= ''
HAPLOTYPES_OUT= ''

try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:s:o:",["input=","sstacks-dir=","output="])
except getopt.GetoptError:
    usage()
for opt, arg in opts:
    if opt == '-h':
        usage()
    elif opt in ("-i", "--input"):
        HAPLOTYPES_IN = arg
    elif opt in ("-o", "--output"):
        HAPLOTYPES_OUT = arg
    elif opt in ("-s", "--sstacks-dir"):
        SSTACKS_DIR = arg
if HAPLOTYPES_IN =='' or not os.path.exists(HAPLOTYPES_IN):
    print "infile :"+HAPLOTYPES_IN+" ,does not exist"
    usage()
if SSTACKS_DIR =='' or not os.path.exists(SSTACKS_DIR):
    print "sstacks results directory :"+SSTACKS_DIR+" ,does not exist"
    usage()

if HAPLOTYPES_OUT=='':
    tmplist=HAPLOTYPES_IN.split('.')
    HAPLOTYPES_OUT='.'.join(tmplist[0:-1])+"_retreat."+tmplist[-1]
    print "output haplotypes files will be : "+HAPLOTYPES_OUT

print "record amb matches"
amb_matches={}
research_ambiguous_haplotypes(SSTACKS_DIR,amb_matches)

print "rewrite haplotypes"
rewrite_happlotypes(HAPLOTYPES_IN,amb_matches,HAPLOTYPES_OUT)


#~ haplotypes_dict={}
#~ header=read_haplotypes(HAPLOTYPES_IN,haplotypes_dict)
#~ print "Found "+`len(haplotypes_dict.keys())`+" loci\n"
 #~ 
#~ matches_dict={}
#~ read_sstacks(SSTACKS_DIR,matches_dict)
#~ print "Found "+`len(matches_dict.keys())`+" samples\n"
#~ research_ambiguities_haplotypes(haplotypes_dict, header, matches_dict, HAPLOTYPES_OUT)
