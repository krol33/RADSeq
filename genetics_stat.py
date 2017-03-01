#!/usr/bin/env python2.7
# -*-coding:Latin-1 -*
########################################################################
## IMPORT
import sys, getopt, os
from os import listdir
from os.path import isfile, join, islink
import argparse
import time 
from collections import Counter
########################################################################
## FUNCTIONS
##############
def prevent_shell_injections( argparse_namespace ):
    """
    @summary: Raises an exception if one parameter contains a backquote or a semi-colon.
    @param argparse_namespace: [Namespase] The result of parser.parse_args().
    """
    for param_name in argparse_namespace.__dict__.keys():
        param_val = getattr(argparse_namespace, param_name)
        if issubclass(param_val.__class__, list):
            new_param_val = list()
            for val in param_val:
                if ';' in val.encode('utf8') or '`' in val.encode('utf8') or '|' in val.encode('utf8'):
                    raise Exception( "';' and '`' are unauthorized characters." ) 
        elif param_val is not None and issubclass(param_val.__class__, str):
            if ';' in param_val.encode('utf8') or '`' in param_val.encode('utf8') or '|' in param_val.encode('utf8'):
                raise Exception( "';' and '`' are unauthorized characters." )
##############
def readPop (popFile, dic) :
    FI=open(popFile,"r")
    for l in FI:
        tmpList=l.strip().split()
        dic["ind"][tmpList[0]]=tmpList[1]
        if tmpList[1] in dic["pop"].keys():
            dic["pop"][tmpList[1]].append(tmpList[0])
        else:
            dic["pop"][tmpList[1]]=[tmpList[0]]
    
    FI.close()
    return len(dic["pop"].keys())
##############
def readGenotype(infile,dicPop, multiallelic, ambiguous, dic):
    """
    @summary: lectures des genotypes et stockage pour ceux qui sont connues 
    @param infile [str] : path to haplotype.tsv file
    @param dicPop [dict] : dictionnary with key = individuals, value = pop index
    @param min_count [int] : on allel is kept if its count is >= to min_count
    @param multiallelic [bool] : genotype is still kept if multiallelic ?
    @param ambiguous [bool] : take ambiguous genotype (?) into account
    @param dicGenOut [dict] : output dictionnary with key = locus, value = {"ind": {ind:genotype},"code_chiffre"={chiffre:allel},"code_allel"={allel:chiffre} }   
    """
    maxAllel=0
    hapFile=open(infile,"r")
    header=hapFile.readline().strip()
    if "Seg Dist" in header:
        nameInd=header.split("\t")[3:]
    else:
        nameInd=header.split("\t")[2:]
    for l in hapFile.readlines():
        tmpDic={} 
        tmpList=l.strip().split("\t")
        locus=int(tmpList[0])
        genotypes_list=[g.split(":")[0] for g in tmpList[2:]] if not "Seg Dist" in header else [g.split(":")[0] for g in tmpList[3:]]
        dic[locus]={"ind":dict(zip(nameInd,genotypes_list))}    # enregistrement des genotypes de chaque individus pour chaque locus : dic[locus]["ind"][indiv]=genotype

        #~ print dic[locus]
        locusAllelList=Counter()
        #~ enregistrement des alleles pour calculer le nombre max d'alleles/genotypes
        for ind in nameInd:    # par des locus
            pop=dicPop[ind]
            
            indAllelList=[a for a in dic[locus]["ind"][ind].split('/')]    # lecture du genotype et split des alleles. On ne va analyser que les genotypes connus et au max biallèlique
            #~ if len(indAllelList) <=2 and dic[locus]["ind"][ind]!="-" and dic[locus]["ind"][ind]!="?" :
            if dic[locus]["ind"][ind] =="-":
                dic[locus]["ind"].pop(ind)
            elif ambiguous == False and dic[locus]["ind"][ind] =="?" :
                dic[locus]["ind"].pop(ind)  
            elif multiallelic == False and len(indAllelList) > 2:
                dic[locus]["ind"].pop(ind)
            else :
                if len(indAllelList)==1 : 
                    locusAllelList.update(indAllelList*2)
                else:
                    locusAllelList.update(indAllelList)

        #~ print locusAllelList
        if maxAllel < len(locusAllelList):
                maxAllel = len(locusAllelList)
        dic[locus]["allele_count"]=locusAllelList
        dic[locus]["code_chiffre"]=dict(zip(list(xrange(1,len(locusAllelList)+1)),[ a[0] for a in locusAllelList.most_common() ]))
        dic[locus]["code_allel"]=dict(zip([ a[0] for a in locusAllelList.most_common() ],list(xrange(1,len(locusAllelList)+1))))
    
    hapFile.close()
    return maxAllel+1

##############
# dicPop={"ind":{},"pop":{}}
# dicGenotype: 
# dicGenotype[locus]: "ind": dictionnaire de génotype par individus
# dicGenotype[locus]: "code_chiffre": dictionnaire des code alleles clé=chiffre, valeur = allele
# dicGenotype[locus]: "code_allel": dictionnaire des code alleles clé=allel, valeur = chiffre
def write_alleles_stat(out,detailed_out, dicGenotype,dicPop,maxAllel):
    
    FO=open(out,"w")

    # header des fichiers
    stringAllel="Locus\tnb_SNP\tnb_Allel\tcode_Allel\tnb_indiv_genotyped\tallele_count\t"
    
    
    if detailed_out :
        tmpallelHeader=[]
        for i in xrange(1,maxAllel):
            tmpallelHeader.append("nb_a"+`i`)
    
        for idx,pop in enumerate(sorted(dicPop["pop"].keys())):
            stringAllel+="pop"+pop+"_nb_indiv_genotyped\tpop"+pop+"nb_allele_diff\t"+"\t".join(["pop"+pop+"_"+a for a in tmpallelHeader ])+"\t"
    
    stringAllel+="\n"
    
    # parcours des locus et remplissage des tableaux de comptage
    nb_locus=0
    for locus in sorted(dicGenotype.keys()):
        nb_locus+=1
        nbAllel=len(dicGenotype[locus]["code_chiffre"].keys())
        nbSNP=0
        if nbAllel ==1 and dicGenotype[locus]["code_chiffre"][1]=="?":
            continue
        if nbAllel > 0:
            nbSNP=len(dicGenotype[locus]["code_chiffre"][1]) if dicGenotype[locus]["code_chiffre"][1] != "?" else len(dicGenotype[locus]["code_chiffre"][2])
        codeAllel=" ".join(`code`+":"+dicGenotype[locus]["code_chiffre"][code] for code in sorted(dicGenotype[locus]["code_chiffre"]) )
        countAllel=" ".join(`idx+1`+":"+`allele[1]` for idx,allele in enumerate(dicGenotype[locus]["allele_count"].most_common()) )
        
        nb_ind=0
        # initialisation des dictionnaires de comptage d'allele
        dicPopAllel={}
        for p in dicPop["pop"].keys():
            dicPopAllel[p]=dict(zip(list(xrange(1,maxAllel)), [0]*maxAllel))
            dicPopAllel[p]["nb_ind"]=0
        # parcours des génotypes et incrémentation des tableaux
        for ind in dicGenotype[locus]["ind"]:
            p=dicPop["ind"][ind]
            gen=dicGenotype[locus]["ind"][ind]
            nb_ind+=1
            dicPopAllel[p]["nb_ind"]+=1
            if detailed_out :
                alleles=[dicGenotype[locus]["code_allel"][a] for a in gen.split("/")]
                alleles.sort()
                if len(alleles)==1:
                    dicPopAllel[p][alleles[0]]+=2
                else:
                    dicPopAllel[p][alleles[0]]+=1
                    dicPopAllel[p][alleles[1]]+=1

        # ecriture des fichiers
        string=`locus`+"\t"+`nbSNP`+"\t"+`nbAllel`+"\t"+codeAllel+"\t"+`nb_ind`+"\t"+countAllel
        stringAllel+=string
        if detailed_out :
            for idx,pop in enumerate(sorted(dicPop["pop"].keys())):
                nb_all_pop=len([a for a in dicPopAllel[pop] if dicPopAllel[pop][a] > 0 and not a =="nb_ind"])
                allel_count=""
                for a in sorted(dicPopAllel[pop].keys()):
                    if not a == "nb_ind":
                        allel_count += "\t"+`dicPopAllel[pop][a]`

                # fichier allele
                stringAllel+="\t"+`dicPopAllel[pop]["nb_ind"]`+"\t"+`nb_all_pop`+allel_count
        stringAllel+="\n"
        
        if nb_locus==2000:
            FO.write(stringAllel)
            nb_locus=0
            stringAllel=""

    if nb_locus !=0:
        FO.write(stringAllel)
        FO.close()

def write_genotypes_stat(prefixout,dicGenotype,dicPop,maxAllel):
    
    # header des fichiers
    summary="Locus\tnb_SNP\tnb_Allel\tcode_Allel\tnb_indiv_genotyped\tgenotype_count\t"
    stringGenotypeList=[summary]*len(dicPop["pop"].keys())
    
    tmpgenotypHeader=[]
    for i in xrange(1,maxAllel):
        for j in xrange(i,maxAllel):
            tmpgenotypHeader.append("nb_g"+`i`+"/"+`j`)
    
    for idx,pop in enumerate(sorted(dicPop["pop"].keys())):
        stringGenotypeList[idx]+="pop"+pop+"_nb_indiv_genotyped\tpop"+pop+"_nb_genotype_diff\t"+"\t".join(["pop"+pop+"_"+g for g in tmpgenotypHeader ])+"\n"
    
    # parcours des locus et remplissage des tableaux de comptage
    nb_locus=0
    for locus in sorted(dicGenotype.keys()):
        nb_locus+=1
        nbAllel=len(dicGenotype[locus]["code_chiffre"].keys())
        nbSNP=0
        if nbAllel==1 and dicGenotype[locus]["code_chiffre"][1]=="?":
            continue
        if nbAllel > 0:
            nbSNP=len(dicGenotype[locus]["code_chiffre"][1]) if dicGenotype[locus]["code_chiffre"][1] != "?" else len(dicGenotype[locus]["code_chiffre"][2])
        codeAllel=" ".join(`code`+":"+dicGenotype[locus]["code_chiffre"][code] for code in sorted(dicGenotype[locus]["code_chiffre"]) )
        nb_ind=0
        # initialisation des dictionnaires de comptage d'allele et de genotype
        dicPopGenotype={}
        for p in dicPop["pop"].keys():
            dicPopGenotype[p]={}
            for i in xrange(1,maxAllel):
                for j in xrange(i,maxAllel):
                    dicPopGenotype[p][`i`+"/"+`j`]=0
        
        genotype_count=Counter()
        # parcours des génotypes et incrémentation des tableaux
        for ind in dicGenotype[locus]["ind"]:
            p=dicPop["ind"][ind]
            gen=dicGenotype[locus]["ind"][ind]
            nb_ind+=1
            alleles=[dicGenotype[locus]["code_allel"][a] for a in gen.split("/")]
            alleles.sort()
            if len(alleles)==1:
                genCode=`alleles[0]`+"/"+`alleles[0]`
                dicPopGenotype[p][genCode]+=1
            else:
                genCode=`alleles[0]`+"/"+`alleles[1]`
                dicPopGenotype[p][genCode]+=1
            genotype_count.update([genCode])
        
        genCount=" ".join([ g[0]+":"+`g[1]` for g in genotype_count.most_common()])
        # ecriture des fichiers
        string=`locus`+"\t"+`nbSNP`+"\t"+`nbAllel`+"\t"+codeAllel+"\t"+`nb_ind`+"\t"+genCount
        for i in xrange(0,len(stringGenotypeList)):
            stringGenotypeList[i]+=string

        for idx,pop in enumerate(sorted(dicPop["pop"].keys())):
            nb_ind_pop=sum(dicPopGenotype[pop].values())
            nb_gen_pop=len([g for g in dicPopGenotype[pop] if dicPopGenotype[pop][g] > 0])
            genotype_count=""
            for a1 in xrange(1,maxAllel):
                for a2 in xrange(a1,maxAllel):
                    genotype_count+="\t"+`dicPopGenotype[pop][`a1`+"/"+`a2`]`
            
            # fichier genotype
            stringGenotypeList[idx]+="\t"+`nb_ind_pop`+"\t"+`nb_gen_pop`+genotype_count+"\n"

        
        if nb_locus==2000:
            for idx,pop in enumerate(sorted(dicPop["pop"].keys())):
                FO=open(prefixout+".pop"+pop+"_genotypes_statistics.tsv","a")
                FO.write(stringGenotypeList[idx])
                FO.close()
            nb_locus=0
            stringGenotypeList=[""]*len(dicPop["pop"].keys())
    if nb_locus !=0:
        for idx,pop in enumerate(sorted(dicPop["pop"].keys())):    
            FO=open(prefixout+".pop"+pop+"_genotypes_statistics.tsv","a")
            FO.write(stringGenotypeList[idx])
            FO.close()
    
########################################################################
## MAIN
parser = argparse.ArgumentParser(description="generate resume or detailed allele stastics and genotypes statistics")
#params
parser.add_argument('-m', '--multiallelic', default=False, action="store_true", help="take also account multiallelic genotypes (more thant 2 alleles), default = False")
parser.add_argument('-a', '--ambiguous', default=False, action="store_true", help="take also account of ambiguous genotype (?) as known allel ")
parser.add_argument('-s', '--summarized', default=False, action="store_true", help="summarized alleles statistics")
parser.add_argument('-d', '--detailed', default=False, action="store_true", help="detailed alleles statistics")
parser.add_argument('-g', '--genotype-stat', default=False, action="store_true", help="detailed genetics statistics")
parser.add_argument('-n', '--name', type=str, default=str(time.time()), help="prefix")
parser.add_argument('--debug', default=False, action='store_true', help="Keep temporary files to debug program.")
# Inputs
group_input = parser.add_argument_group('Inputs')
group_input.add_argument('-p', '--population-map', required=True, help='population map.')
group_input.add_argument('-i', '--haplotype-tsv', required=True, help='Stacks haplotype tsv file')
group_input.add_argument('-o', '--output-directory', default='', help='Output directory')
args = parser.parse_args()
prevent_shell_injections(args)

if args.output_directory != '' and not os.path.exists(args.output_directory):
    raise Exception("output directory :"+args.output_directory+" ,does not exist")
elif args.output_directory == '' :
    output_directory = os.path.dirname(os.path.abspath(args.haplotype_tsv))
else:
    output_directory = args.output_directory
    
if args.haplotype_tsv=='' or not os.path.exists(args.haplotype_tsv):
    print "haplotype file :"+args.haplotype_tsv+" ,does not exist\n"
    usage()
if args.population_map =='' or not os.path.exists(args.population_map):
    print "population map file :"+args.population_map+" ,does not exist\n"
    usage()

if not args.detailed and not args.summarized and not args.genotype_stat:
    print "you must choose at least one output type betweeb -r -d -g options\n"
    usage()
####
# read population map
print "read population map: " + args.population_map
dicPop={"ind":{},"pop":{}}
nb_pop=readPop( args.population_map ,dicPop)

###
# dicGenotype: 
# dicGenotype[locus]: "ind": dictionnaire de génotype par individus
# dicGenotype[locus]: "code_chiffre": dictionnaire des code alleles clé=chiffre, valeur = allele
# dicGenotype[locus]: "code_allel": dictionnaire des code alleles clé=allele, valeur = chiffre
print "read haplotype file : "+args.haplotype_tsv
dicGenotype={}
maxAllel = readGenotype(args.haplotype_tsv,dicPop["ind"],args.multiallelic, args.ambiguous, dicGenotype)

###
print "outputs statistcis files will be :"
outAllele_summarized=""
outAllele_detailed=""
if args.summarized:
    outAllele_summarized = os.path.join(output_directory,args.name+".alleles_summarized_stat.tsv")
    print "\t",outAllele_summarized
if args.detailed:
    outAllele_detailed =  os.path.join(output_directory,args.name+".alleles_detailed_stat.tsv") if not args.detailed else os.path.join(args.output_directory,args.name+".alleles_detailed_stat.tsv")
    print "\t",outAllele_detailed

if args.genotype_stat :
    for pop in sorted(dicPop["pop"].keys()):
        outGen=os.path.join(output_directory,args.name+".pop"+pop+"_genotypes_statistics.tsv")
        print "\t"+outGen
        if os.path.exists(outGen):
            os.remove(outGen)

print maxAllel," maximum allel observed in one population"
if args.summarized:
    write_alleles_stat(outAllele_summarized,False, dicGenotype,dicPop,maxAllel)
if args.detailed:
    write_alleles_stat(outAllele_detailed,True, dicGenotype,dicPop,maxAllel)
if args.genotype_stat : 
    write_genotypes_stat(os.path.join(output_directory,args.name), dicGenotype,dicPop,maxAllel)
