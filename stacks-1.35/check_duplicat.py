#!/usr/bin/python

import sys, getopt, os
########################################################################
## FUNCTIONS
def usage():
    print 'check_duplicat.py -i <input-haplotype> -p <poplation-map> -n "name of individuals duplicated separate by space" -s <min_incoherence_accepted> -a <analysis type genotypes/populations> -o <output_dir>'
    sys.exit(2)

def check_dup_genotype(names, genotypes_in, pop_file, start, min_inc, out):
    dic_dup_name=dict.fromkeys(names.split())

    for i,line in enumerate(open(pop_file,"r")):
        for dup in dic_dup_name:
            if dic_dup_name[dup]==None:
                    dic_dup_name[dup]=[]
            if dup in line:
                dic_dup_name[dup].append(i+start)

    one_dup=0
    for dup in dic_dup_name:
        if len(dic_dup_name[dup]) == 1:
            one_dup+=1
    if one_dup == len(dic_dup_name):
        print str(one_dup)+" on "+str(len(dic_dup_name))+" individual are present only one time"
        usage()

    dic_gen_dup=dict.fromkeys(names.split())
    FH_rm=open(out,"w")
    FH_rm.write("#locus\tnb_incoherence\tgenotype des duplicats\n")
    for line in open(genotypes_in,"r") :
        if not line.startswith("Catalog"):
            loc=line.split()[0]
            b=0
            dic_gen_dup=dict.fromkeys(names.split())
            for dup in dic_dup_name:
                if dic_gen_dup[dup]==None:
                    dic_gen_dup[dup]=[]
                for col in dic_dup_name[dup]:
                    if not line.split()[col] in dic_gen_dup[dup] :
                        dic_gen_dup[dup].append(line.split()[col])
                gen_known=[n for n in dic_gen_dup[dup] if n != "-" ]
                #~ print gen_known
                if len(gen_known) > 1 :
                    b+=1
            if b > min_inc :
                s=loc+"\tincoherence:"+`b`+"\t"
                for dup in dic_gen_dup:
                    s+=dup+":"+";".join(dic_gen_dup[dup])+"\t"
                FH_rm.write(s+"\n")
    FH_rm.close()
            
########################################################################
## MAIN
POP_FILE=""
ANALYSIS_TYPE=""
POLYMORPHES_LOCUS=""
DUP_NAME=""
SEUIL=1
OUT_DIR=""

try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:p:a:n:s:o:",["input=","population-map=","analysis-type=","dup-name=","min-incoherence=","out-dir="])
except getopt.GetoptError:
    usage()
for opt, arg in opts:
    if opt == '-h':
        usage()
    elif opt in ("-i", "--input"):
        POLYMORPHES_LOCUS = arg
    elif opt in ("-p", "--population-map"):
        POP_FILE = arg
    elif opt in ("-a", "--analysis-type"):
        ANALYSIS_TYPE = arg
    elif opt in ("-n", "--dup-name"):
        DUP_NAME = arg        
    elif opt in ("-s", "--min-incoherence"):
        SEUIL = int(arg)
    elif opt in ("-o", "--out-dir"):
        OUT_DIR = arg
    
if POLYMORPHES_LOCUS =='' or not os.path.exists(POLYMORPHES_LOCUS):
    print "haplotype file :"+POLYMORPHES_LOCUS+" ,does not exist"
    usage()
if POP_FILE =='' or not os.path.exists(POP_FILE):
    print "population map file :"+POP_FILE+" ,does not exist"
    usage()
if not os.path.exists(OUT_DIR):
    print "output directory :"+OUT_DIR+" ,does not exist"
    usage()
    
start=3
if ANALYSIS_TYPE not in ["populations","genotypes"]:
    print "Analysis type must be either \"populations\" or \"genotypes\" "
    usage()
elif ANALYSIS_TYPE=="populations":
    start=2

if DUP_NAME == "":
    print "You must at least give one name of duplicated individual"
    usage()

if OUT_DIR != "":
    FH_rm = OUT_DIR+"/dup_check_locus_rm.txt"
else :
    FH_rm = os.path.dirname(POLYMORPHES_LOCUS)+"/dup_check_locus_rm.txt" if os.path.dirname(POLYMORPHES_LOCUS) != "" else "./dup_check_locus_rm.txt"
    
check_dup_genotype(DUP_NAME, POLYMORPHES_LOCUS, POP_FILE, start, SEUIL, FH_rm)
