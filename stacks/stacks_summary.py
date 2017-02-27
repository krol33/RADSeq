# -*- coding: utf-8 -*-
#
# Copyright (C) 2016 INRA
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Maria Bernard - Sigenae Team'
__copyright__ = 'Copyright (C) 2016 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'maria.bernard@inra.fr'

import sys
import re
import os
import argparse
import glob
import json
import gzip
import numpy
import operator


def parse_log(log):
    """
    @summary : parse Stacks log file, and store cstacks informations.
    @param log : [path] Path to Stacks denovo map log file
    @output cstacks_dict : [dict] dictionnary with one key (to be update in parse_cstacks) of the evolution of catalog: cstacks_dict["evol_cat"]=[int, int, int]
    @output first_stacks_title : name of the first stacks program launched, either pstacks or ustacks.
    """
    print "\t... "+ log +" ..."
    cstacks_dict={"evol_cat":[[0,0]]}
    FH_log = open(log)

    end_cstacks=False
    first_stacks_title = ""

    for line in FH_log.readlines():
        if "ustacks" in line and first_stacks_title =="":
            first_stacks_title = "ustacks"
        elif  "pstacks" in line and first_stacks_title =="":
            first_stacks_title = "pstacks"

        # CSTACKS # store the evolution of the number of loci in the catalog
        if not end_cstacks :
            if "loci were newly added to the catalog." in line:
                num_cat_ind = cstacks_dict["evol_cat"][-1][0] +1
                num_cat_loc = cstacks_dict["evol_cat"][-1][1] + int(line.split()[0])
                cstacks_dict["evol_cat"].append([num_cat_ind,num_cat_loc])
        # END CSTACKS
        if "sstacks" in line:
            break

    #Cstacks not launched
    if len(cstacks_dict["evol_cat"]) == 1:
        cstacks_dict=dict()
    return cstacks_dict , first_stacks_title

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

def parse_ustacks(res_dir,ustacks_dict):
    """
    @summary : Parse usatcks output files and store information to summerize
    @param res_dir : [path] Path to ustacks output directory
    @output ustacks_dict : [dict] dictionnary of list per sample :
    ustacks_dict={
                 "titles" : ["Nb used reads", "Nb clusters", "Nb validated clusters","Nb reads clustered","Nb secondary reads clustered", "Mean cluster coverage", "Nb polymorphic clusters","Nb reads on polymorphic clusteres","Mean polymorphic cluster coverage"],
                 "samples" : {
                     "PopA_02": [60, 3, 3, 60, 5, 20.0, 1, 12, 12.0],
                     "PopA_01": [66, 3, 3, 66, 3, 22.0, 1, 18, 18.0],
                     "PopB_03": [92, 3, 3, 92, 7, 30.67, 1, 38, 38.0]
                    }
                }
    """
    files=[ f for f in glob.glob(os.path.join(res_dir,'*tags.tsv*')) if not "catalog" in f]
    if len(files) == 0:
        print "\t... Ustacks not launched"
        return

    ustacks_dict["titles"]=["Nb used reads", "Nb clusters", "Nb validated clusters","Nb reads clustered","Nb secondary reads clustered", "Average cluster coverage", "Nb polymorphic clusters","Nb reads on polymorphic clusters","Average polymorphic cluster coverage"]
    ustacks_dict["samples"]=dict()

    for file in files:
        print "\t... "+file+" ..."
        sample=os.path.basename(file).replace(".tags.tsv"," ").split()[0]
        ustacks_dict["samples"][sample]=[0, 0, 0, 0, 0, 0, 0, 0, 0]

        FH_in = None
        if not is_gzip(file):
            FH_in = open( file )
        else:
            FH_in = gzip.open( file )

        # white=False
        # poly=False

        for line in FH_in.readlines():
            if line.startswith("#"):
                continue
            # count clusters
            elif "consensus" in line:
                white=False     # not blacklisted clusters
                poly=False      # polymorphic clusters
                ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Nb clusters")]+=1
                if line.split()[8]!="1":
                    white=True
                    ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Nb validated clusters")]+=1
            # count number of whitelisted polymorphic clusters
            elif "model" in line and white and "E" in line.split()[-1]:
                poly=True
                ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Nb polymorphic clusters")]+=1
            # coverage
            elif not "model" in line and not "consensus" in line :
                ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Nb used reads")] += 1
                # coverage on whitelisted clusters
                if white :
                    ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Nb reads clustered")]+=1
                    # coverage on polymorphic clusters
                    if poly:
                        ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Nb reads on polymorphic clusters")]+=1
                    # secondary reads coverage on whitelisted clusters
                    if "secondary" in line :
                        ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Nb secondary reads clustered")]+=1

        # compute mean coverage on whitelisted clusters
        if ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Nb validated clusters")] > 0 :
            ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Average cluster coverage")] = round(float(ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Nb reads clustered")])/ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Nb validated clusters")],2)
        else :
            ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Average cluster coverage")] = 0.0
        # compute mean coverage on whitelisted polymorphic clusters
        if ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Nb polymorphic clusters")] > 0 :
            ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Average polymorphic cluster coverage")] = round(float(ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Nb reads on polymorphic clusters")])/ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Nb polymorphic clusters")],2)
        else:
            ustacks_dict["samples"][sample][ustacks_dict["titles"].index("Average polymorphic cluster coverage")] = 0.0

def parse_cstacks(res_dir, cstacks_dict):
    """
    @summary : Parse csatcks output files and store information to summerize
    @param res_dir : [path] Path to cstacks output directory
    @param cstacks_dict : [dict] dictionnary initiated in parse_log cstacks_dict["evol_cat"]=[[nb_indiv_added, nb_clust]]
    @output cstacks_dict : [dict] updated dictionnary with
                {
                 "min_nb_indiv_per_cluster": [[1, 4], [2, 3], [3, 3]],
                 "evol_cat": [[0, 0], [2, 3], [3, 3]]
                 "summary": {
                    "titles" : ["total","single loci", "loci merging at least 2 individuals", "loci merging at least 50% individuals", "loci merging 100% individuals"];
                    "details" : [
                        ["loci" ,102632, 10161, 92471, 69459, 21367],
                        ["polymorphic loci", 32222, 369, 31853, 27577, 8061],
                        ["SNPs" , 47451, 543, 46908, 40551, 11946]
                    ]
                 },
                 "nb_haplotype_distribution": [["1", 0],["2", 11029], ["3", 1376], ["4", 288], [">=5", 132], [">=10", 2]],
                 "nb_snp_distribution": [["1", 9338], ["2", 2416], ["3", 677], ["4", 0], [">=5", 90], [">=10", 1]]
                }
    """
    catalog_files=glob.glob(os.path.join(res_dir,"*catalog*"))
    catalog_files.sort(reverse=True)
    nb_agregated = []
    catalog_dict= {}
    # Cstacks not launched
    if len(catalog_files) == 0:
        print "\t... Cstacks not launched"
        return

    for file in catalog_files:
        print "\t... "+file+" ..."

        if "tags" in file :
            # open tags file
            FH_in = None
            if not is_gzip(file):
                FH_in = open( file )
            else:
                FH_in = gzip.open( file )
            # parse tags file and count number of samples agregated per cluster
            for line in FH_in:
                if line.startswith("#"):
                    continue
                samples=list()
                locus_id = line.split()[2]
                # parse sample agregated
                for s in line.strip().split()[7].split(','):
                    # add new sample identifier
                    if not s.split('_')[0] in samples:
                        samples.append(s.split("_")[0])
                catalog_dict[locus_id] = {"nb_ind" : len(samples)}
            FH_in.close()

        elif "snps" in file:
            # open snps file
            FH_in = None
            if not is_gzip(file):
                FH_in = open( file )
            else:
                FH_in = gzip.open( file )
            # parse snps file and count number of snps per cluster
            for line in FH_in:
                if line.startswith("#"):
                    continue
                locus_id = line.split()[2]
                if line.split()[4]=="E":
                    if "nb_snp" in catalog_dict[locus_id]:
                        catalog_dict[locus_id]["nb_snp"] += 1
                    else:
                        catalog_dict[locus_id]["nb_snp"] = 1
            FH_in.close()

        elif "alleles" in file:
            # open alleles file
            FH_in = None
            if not is_gzip(file):
                FH_in = open( file )
            else:
                FH_in = gzip.open( file )
            # parse alleles file and count number of haplotype per cluster
            for line in FH_in:
                if line.startswith("#"):
                    continue
                locus_id = line.split()[2]
                if len(line.split()) == 6:
                    if "nb_hap" in catalog_dict[locus_id]:
                        catalog_dict[locus_id]["nb_hap"] += 1
                    else:
                        catalog_dict[locus_id]["nb_hap"] = 1
            FH_in.close()

    # parse cstacks dictionnary
    # cluster distribution in function of nb indiv agregated
    cstacks_dict["min_nb_indiv_per_cluster"]=list()

    # extract information of nb_ind per cluster
    nb_indiv_clust = [ catalog_dict[locus_id]["nb_ind"] for locus_id in catalog_dict]
    # count occurences for each number of indiv ==> [[nb_ind,nb_clust]] reverse sorted by number of individuals
    nb_indiv_clust_count = [ [c,nb_indiv_clust.count(c)] for c in set(nb_indiv_clust)]
    nb_indiv_clust_count.sort(reverse=True)
    c=0
    for n,cl in nb_indiv_clust_count:
        cstacks_dict["min_nb_indiv_per_cluster"].insert(0,[n,c + cl] )
        c += cl

    # cstacks summary
    cstacks_dict["summary"] = dict()
    cstacks_dict["summary"]["titles"] = ["total","single loci", "loci merging at least 2 individuals", "loci merging at least 50% individuals", "loci merging 100% individuals"]
    cstacks_dict["summary"]["details"] = [["loci" , len(catalog_dict),0,0,0,0], ["polymorphic loci",0,0,0,0,0], ["SNPs",0,0,0,0,0]]
    # cstacks haplotype and snp distributions
    cstacks_dict["nb_haplotype_distribution"]=[["1",0],["2",0],["3",0],["4",0],[">=5",0],[">=10",0]]
    cstacks_dict["nb_snp_distribution"]=[["1",0],["2",0],["3",0],["4",0],[">=5",0],[">=10",0]]
    for locus_id in catalog_dict:
        # selection of polymorphic cluster
        if "nb_snp" in catalog_dict[locus_id]:
            # count total number of polymorphic loci and snp
            cstacks_dict["summary"]["details"][1][cstacks_dict["summary"]["titles"].index("total")+1] += 1
            cstacks_dict["summary"]["details"][2][cstacks_dict["summary"]["titles"].index("total")+1] += catalog_dict[locus_id]["nb_snp"]
            # update snp distribution
            if catalog_dict[locus_id]["nb_snp"]==1:
                cstacks_dict["nb_snp_distribution"][0][1]+=1
            elif catalog_dict[locus_id]["nb_snp"]==2:
                cstacks_dict["nb_snp_distribution"][1][1]+=1
            elif catalog_dict[locus_id]["nb_snp"]==3:
                cstacks_dict["nb_snp_distribution"][2][1]+=1
            elif catalog_dict[locus_id]["nb_snp"]==2:
                cstacks_dict["nb_snp_distribution"][3][1]+=1
            elif catalog_dict[locus_id]["nb_snp"]>=5:
                cstacks_dict["nb_snp_distribution"][4][1]+=1
                if catalog_dict[locus_id]["nb_snp"]>=10:
                    cstacks_dict["nb_snp_distribution"][5][1]+=1
            # update haplotype distribution
            if catalog_dict[locus_id]["nb_hap"]==2:
                cstacks_dict["nb_haplotype_distribution"][1][1]+=1
            elif catalog_dict[locus_id]["nb_hap"]==3:
                cstacks_dict["nb_haplotype_distribution"][2][1]+=1
            elif catalog_dict[locus_id]["nb_hap"]==4:
                cstacks_dict["nb_haplotype_distribution"][3][1]+=1
            elif catalog_dict[locus_id]["nb_hap"]>=5:
                cstacks_dict["nb_haplotype_distribution"][4][1]+=1
                if catalog_dict[locus_id]["nb_hap"]>=10:
                    cstacks_dict["nb_haplotype_distribution"][5][1]+=1
        # loci singleton
        if catalog_dict[locus_id]["nb_ind"]==1:
            # count number of total singleton loci
            cstacks_dict["summary"]["details"][0][cstacks_dict["summary"]["titles"].index("single loci")+1] += 1
            # count number of polymorphic loci and snp on singleton loci
            if "nb_snp" in catalog_dict[locus_id]:
                cstacks_dict["summary"]["details"][1][cstacks_dict["summary"]["titles"].index("single loci")+1] += 1
                cstacks_dict["summary"]["details"][2][cstacks_dict["summary"]["titles"].index("single loci")+1] += catalog_dict[locus_id]["nb_snp"]
        # loci merging at least 2 samples
        if catalog_dict[locus_id]["nb_ind"]>=2:
            # count number of total loci merging at least 2 individuals
            cstacks_dict["summary"]["details"][0][cstacks_dict["summary"]["titles"].index("loci merging at least 2 individuals")+1] += 1
            # count number of polymorphic loci and snp on loci merging at least 2 individuals
            if "nb_snp" in catalog_dict[locus_id]:
                cstacks_dict["summary"]["details"][1][cstacks_dict["summary"]["titles"].index("loci merging at least 2 individuals")+1] += 1
                cstacks_dict["summary"]["details"][2][cstacks_dict["summary"]["titles"].index("loci merging at least 2 individuals")+1] += catalog_dict[locus_id]["nb_snp"]
        # loci merging at least 50% of samples
        if catalog_dict[locus_id]["nb_ind"]>= float(cstacks_dict["evol_cat"][-1][0])*50/100:
            # count number of total loci merging at least 50% individuals
            cstacks_dict["summary"]["details"][0][cstacks_dict["summary"]["titles"].index("loci merging at least 50% individuals")+1] += 1
            # count number of polymorphic loci and snp on loci merging at least 50 individuals
            if "nb_snp" in catalog_dict[locus_id]:
                cstacks_dict["summary"]["details"][1][cstacks_dict["summary"]["titles"].index("loci merging at least 50% individuals")+1] += 1
                cstacks_dict["summary"]["details"][2][cstacks_dict["summary"]["titles"].index("loci merging at least 50% individuals")+1] += catalog_dict[locus_id]["nb_snp"]
        # loci merging 100% of samples
        if catalog_dict[locus_id]["nb_ind"] == cstacks_dict["evol_cat"][-1][0]:
            # count number of total loci merging 100% individuals
            cstacks_dict["summary"]["details"][0][cstacks_dict["summary"]["titles"].index("loci merging 100% individuals")+1] += 1
            # count number of polymorphic loci and snp on loci merging 100% individuals
            if "nb_snp" in catalog_dict[locus_id]:
                cstacks_dict["summary"]["details"][1][cstacks_dict["summary"]["titles"].index("loci merging 100% individuals")+1] += 1
                cstacks_dict["summary"]["details"][2][cstacks_dict["summary"]["titles"].index("loci merging 100% individuals")+1] += catalog_dict[locus_id]["nb_snp"]

def parse_sstacks(res_dir, pop_map,sstacks_stat):
    """
    @summary : Parse ssatcks output files and store information to summerize
    @param res_dir : [path] Path to sstacks output directory
    @param pop_map : [path] Path to the populaiton map file (if available)
    @output sstacks_dict : [dict] dictionnary describing number of catalog loci retrieved in sample, and corresponding number of haplotype. (Pop id added if available)
                {
                    "titles": ["Samples", "population id", "matching loci", "matching haplotype"],
                    "samples": {
                        "EFFICACE_113": ["2", 50822, 54830],
                        "EFFICACE_209": ["2", 48143, 52354],
                        "EFFICACE_43": ["1", 48481, 52541]
                    }
                }
    """
    matches_files=glob.glob(os.path.join(res_dir,"*matches*"))
    if len(matches_files) == 0 :
        print "\t... Sstacks not launched"
        return
    # parse pop_map file
    pop_dic={}
    if pop_map != None:
        Fh_pop= open(pop_map)
        for line in Fh_pop:
            pop_dic[line.split()[0]] = line.strip().split()[1]
        Fh_pop.close()

    sstacks_stat["titles"]=["Samples","matching loci","matching haplotype"]
    if len(pop_dic) > 0 :
        sstacks_stat["titles"].insert(1,"population id")

    sstacks_stat["samples"]=dict()
    # How many matches per sample
    for file in matches_files:
        # open matches file
        print "\t... "+file+" ..."
        sample=os.path.basename(file).replace(".matches.tsv"," ").split()[0]
        dic_sstacks = dict()
        FH_in = None
        if not is_gzip(file):
            FH_in = open( file )
        else:
            FH_in = gzip.open( file )
        # parse matches file
        for line in FH_in:
            if line.startswith("#"):
                continue
            # catalog loci id
            loci_id = line.split()[2]
            if loci_id in dic_sstacks:
                # update number of haplotype for on catalog loci
                dic_sstacks[loci_id] += 1
            else:
                # add a new catalog loci
                dic_sstacks[loci_id] = 1
        FH_in.close()
        # summarise for each sample
        sstacks_stat["samples"][sample]=[len(dic_sstacks), sum(dic_sstacks.values())]
        # update with pop id
        if len(pop_dic)>0:
            sstacks_stat["samples"][sample].insert(0,pop_dic[sample])

def parse_haplotype(res_dir, pop_map,haplotype_stat):
    """
    @summary : Parse population/genotype output file and store information to summerize
    @param res_dir : [path] Path to population/genotype output directory
    @param pop_map : [path] Path to the populaiton map file (if available)
    @output return Stacks program launched, either Populations or Genotypes
    @output haplotype_stat : [dict] dictionnary describing number of catalog loci genotyped, zygotie rate per sample (Pop id added if available), and number of (polymorphic) loci kept in function of minimum number of individuals genotyped
                {
                    "zygosity_detail": {
                        "titles": ["Sample", "population id", "missing genotype", "heterozygote genotype", "homozygote genotype", "heterozygosity rate (%)"], "samples": {
                            "EFFICACE_113": ["2", 24172, 3710, 45804, 7.49],
                            "EFFICACE_209": ["2", 26887, 3752, 43047, 8.02],
                            "EFFICACE_43": ["1", 26615, 3594, 43477, 7.64]
                        }
                    },
                    "summary": {
                        "titles": ["Nb loci genotyped", "Nb loci genotyped by at least 50% of individuals", "Nb polymorphic loci genotyped", "Nb polymorphic loci genotyped by at least 50% of individuals", "Nb polymorphic loci with 1 SNP and 2 alleles"],
                        "values": [73686, 47261, 8729, 7333, 7946]
                    },
                    "min_nb_indiv_per_cluster": {
                        "all loci": [["1", 73686], ["2", 47261], ["3", 22437]],
                        "polymorphic loci": [["1", 8729], ["2", 7333], ["3", 4317]]
                    }
                }
    """
    haplotype_file=glob.glob(os.path.join(res_dir,"batch_1.haplotypes*.tsv*"))
    if len(haplotype_file) == 0:
        print "\t... Populations/Genotypes not launched"
        return
    else:
        print "\t... "+haplotype_file[0]+" ..."

    ## if populations Stacks program launched first sample is in the 3rd column, if genotypes Stacks program in the 4th
    if os.path.basename(haplotype_file[0]).startswith("batch_1.haplotypes_"):
        start = 3
    else:
        start = 2

    # cluster distribution in function of nb indiv agregated
    haplotype_stat["min_nb_indiv_per_cluster"]=dict()
    haplotype_stat["min_nb_indiv_per_cluster"]["all loci"]=list()
    haplotype_stat["min_nb_indiv_per_cluster"]["polymorphic loci"]=list()

    # extract information of nb_ind per cluster
    nb_indiv_clust = list()
    nb_indiv_poly_clust = list()
    cnt = 0
    FH_in = open(haplotype_file[0])
    # parse haplotype file
    for line in FH_in:
        if "Catalog" in line:
            samples = line.strip().split("\t")[start:]
            # initialize zygosity table and summary
            haplotype_stat["zygosity_detail"]=dict()
            haplotype_stat["zygosity_detail"]["titles"]= ["Sample","missing genotype", "heterozygote genotype", "homozygote genotype", "heterozygosity rate (%)"]
            haplotype_stat["zygosity_detail"]["samples"]={sample : [0, 0, 0, 0.0] for sample in samples}

            haplotype_stat["summary"] = dict()
            haplotype_stat["summary"]["titles"]= ["Nb loci genotyped", "Nb loci genotyped by at least 50% of individuals","Nb polymorphic loci genotyped", "Nb polymorphic loci genotyped by at least 50% of individuals","Nb polymorphic loci with 1 SNP and 2 alleles"]
            haplotype_stat["summary"]["values"]= [0, 0, 0, 0, 0]

            # retrieved the column index for the number of individual genotyped
            cnt = line.split("\t").index("Cnt")
        else :
            # number of individuals genotyped
            nb_indiv_clust.append(int(line.split()[cnt]))

            # update summary with number of (polymorphic) cluster genotyped (for at least 50% of individuals)
            haplotype_stat["summary"]["values"][0] += 1
            if int(line.split()[cnt]) > 50.0*len(samples)/100:
                haplotype_stat["summary"]["values"][1] += 1
            if "/" in line:
                nb_indiv_poly_clust.append(int(line.split()[cnt]))
                haplotype_stat["summary"]["values"][2] += 1
                if int(line.split()[cnt]) > 50.0*len(samples)/100:
                    haplotype_stat["summary"]["values"][3] += 1

            # parse haplotype
            haplotypes = line.strip().split("\t")[start:]
            alleles = list()
            for idx,hap in enumerate(haplotypes):
                sample = samples[idx]
                # count unknown haplotype per sample
                if hap == "-":
                    haplotype_stat["zygosity_detail"]["samples"][sample][0] += 1
                # count heterozygote haplotype per sample
                elif "/" in hap:
                    for a in hap.split("/"):
                        if not a in alleles:
                            alleles.append(a)
                    haplotype_stat["zygosity_detail"]["samples"][sample][1] += 1
                # count homozygote haplotype per sample
                else:
                    haplotype_stat["zygosity_detail"]["samples"][sample][2] += 1
                    if not hap == "consensus" and not hap in alleles:
                        alleles.append(hap)
            # count number of polymorphic cluster with 1 SNP and 2 alleles
            if len(alleles) == 2 and len(alleles[0]) == 1 :
                haplotype_stat["summary"]["values"][4] += 1
    FH_in.close()

    # add heterozygot rate per sample in zigosity details
    for sample in haplotype_stat["zygosity_detail"]["samples"]:
        Hz_rate = round(haplotype_stat["zygosity_detail"]["samples"][sample][1] * 100.0 / sum(haplotype_stat["zygosity_detail"]["samples"][sample][1:3]),2)
        haplotype_stat["zygosity_detail"]["samples"][sample][3] = Hz_rate

    # Add pop ID information:
    if pop_map != None:
        FH_pop= open(pop_map)
        haplotype_stat["zygosity_detail"]["titles"].insert(1,"population id")
        for line in FH_pop:
            if line.split()[0] in haplotype_stat["zygosity_detail"]["samples"]:
                haplotype_stat["zygosity_detail"]["samples"][line.split()[0]].insert(0,line.strip().split()[1])
        FH_pop.close()

    # compute cumulative distribution of the number of cluster kept with a minimum number of individuals
    # count occurences for each number of indiv ==> [[nb_ind,nb_clust]] reverse sorted by number of individuals
    nb_indiv_clust_count = [ [c,nb_indiv_clust.count(c)] for c in set(nb_indiv_clust)]
    nb_indiv_clust_count.sort(reverse=True)
    c=0
    for n,cl in nb_indiv_clust_count:
        haplotype_stat["min_nb_indiv_per_cluster"]["all loci"].insert(0,[n,c + cl] )
        c += cl

    nb_indiv_poly_clust_count = [ [c,nb_indiv_poly_clust.count(c)] for c in set(nb_indiv_poly_clust)]
    nb_indiv_poly_clust_count.sort(reverse=True)
    c=0
    for n,cl in nb_indiv_poly_clust_count:
        haplotype_stat["min_nb_indiv_per_cluster"]["polymorphic loci"].insert(0,[n,c + cl] )
        c += cl

    # return Stacks program launched for tab title in HTML output file
    if start == 3:
        return "Genotypes"
    else:
        return "Populations"


def summarize_results( prog, out_dir, log, pop_map, summary_html ):
    """
    @summary : parse log and stacks results files and generate summerize results HTML file, based on template (stacks_summary_tpl.html)
    @param prog : [str] Stacks prog launched
    @param out_dir : [Path] Stacks output directory
    @param log : [path] path to Stacks log file (only for cstacks, ref_map.pl and denovo_map.pl)
    @param pop_map : [path] path to Stacks input pop_map file (only for populations and potentially for ref_map.pl and denovo_map.pl)
    @param summary_html: [path] Path to HTML output file
    """
    Fh_tpl = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "stacks_summary_tpl.html"))
    FH_summary_tpl = open( summary_html,"w" )

    # one of the individuals Stacks step
    first_step =""
    # either ustacks or pstacks (if launched)
    first_stacks_title=""

    ### PARSE LOG
    if not log == None :
        print "parsing log"
        cstacks_stat,first_stacks_title = parse_log(log)
    else:
        cstacks_stat = dict()

    ### PARSE USTACKS/PSTACKS RES
    u_pstacks_stat = dict()
    if prog in ["ustacks", "pstacks", "ref_map.pl", "denovo_map.pl"] :
        print "\nparsing ustacks results"
        parse_ustacks(out_dir,u_pstacks_stat)
        if u_pstacks_stat != dict():
            first_step = "u-pstacks"
        if first_stacks_title == "":
            first_stacks_title = prog

    ### PARSE CSTACKS RES
    if prog in ["cstacks", "ref_map.pl", "denovo_map.pl"] :
        print "\nparsing cstacks results"
        parse_cstacks(out_dir, cstacks_stat)
        if first_step =="" and cstacks_stat != dict():
            first_step = "cstacks"

    ### PARSE SSTACKS RES
    sstacks_stat = dict()
    if prog in ["sstacks", "ref_map.pl", "denovo_map.pl"] :
        print "\nparsing sstacks results"
        parse_sstacks(out_dir, pop_map, sstacks_stat)
        if first_step == "" and sstacks_stat != dict():
            first_step = "sstacks"

    ### PARSE POPULATION/GENOTYPE RES
    haplotype_stat = dict()
    if prog in ["populations", "genotypes", "ref_map.pl", "denovo_map.pl"] :
        print "\nparsing populations/genotypes results"
        haplotype_type = parse_haplotype(out_dir, pop_map, haplotype_stat)
        if first_step == "" and haplotype_stat != dict():
            first_step = "haplotype"
    if first_step == "":
        raise Exception("Something wrong append, there is no Stacks ouput file to parse\n")

    ### WRITE SUMMARY.HTML
    for line in Fh_tpl.readlines():
        # USTACKS / PSTACKS
        # add USTACKS / PSTACKS data and uncomment var u_pstacks
        if "### STACKS_DETAIL ###" in line and u_pstacks_stat != dict():
            line = line.replace("//","")
            line = line.replace( "### STACKS_DETAIL ###", json.dumps(u_pstacks_stat) )
        # udpate first stacks step tab title and activate tab
        elif "### FIRST_STEP_TITLE ###" in line and u_pstacks_stat != dict():
            line = line.replace("<!--","").replace("-->","")
            if first_step == "u-pstacks" :
                line = line.replace("<li","<li class=\"active\"")
            line = line.replace("### FIRST_STEP_TITLE ###", first_stacks_title[0].upper()+first_stacks_title[1:])
        elif "<div id=\"u-pstacks-stat\"" in line and first_step == "u-pstacks" :
            line = line.replace("tab-pane disabled" ,"tab-pane active")
        # CSTACKS
        # add CSTACKS data and uncomment var cstacks
        elif "### CSTACKS_DETAIL ###" in line and cstacks_stat != dict():
            line = line.replace("//","")
            line = line.replace( "### CSTACKS_DETAIL ###", json.dumps(cstacks_stat) )
        # activate cstacks tab
        elif "Cstacks Stat</a></li>" in line and cstacks_stat != dict():
            line = line.replace("<!--","").replace("-->","")
            if first_step == "cstacks" :
                line = line.replace("<li","<li class=\"active\"")
        elif "<div id=\"cstacks-stat\"" in line and first_step == "cstacks" :
            line = line.replace("tab-pane disabled" ,"tab-pane active")
        # SSTACKS
        # add SSTACKS data and uncomment var sstacks
        elif "### SSTACKS_DETAIL ###" in line and sstacks_stat != dict():
            line = line.replace("//","")
            line = line.replace( "### SSTACKS_DETAIL ###", json.dumps(sstacks_stat) )
        # activate sstacks tab
        elif "Sstacks Stat</a></li>" in line and sstacks_stat != dict():
            line = line.replace("<!--","").replace("-->","")
            if first_step == "sstacks" :
                line = line.replace("<li","<li class=\"active\"")
        elif "<div id=\"sstacks-stat\"" in line and first_step == "sstacks" :
            line = line.replace("tab-pane disabled" ,"tab-pane active")
        # GENOTYPES / POPULATIONS
        # add GENOTYPES / POPULATIONS data and uncomment var haplotype
        elif "### HAPLOTYPE_TAB_TITLE ###" in line and haplotype_stat != dict():
            line = line.replace("<!--","").replace("-->","")
            line = line.replace( "### HAPLOTYPE_TAB_TITLE ###", haplotype_type )
            if first_step == "haplotype" :
                line = line.replace("<li","<li class=\"active\"")
        # udpate last stacks step tab title and activate tab
        elif "### HAPLOTYPE_DETAIL ###" in line and haplotype_stat != dict():
            line = line.replace("//","")
            line = line.replace( "### HAPLOTYPE_DETAIL ###", json.dumps(haplotype_stat) )
        elif "<div id=\"haplotype-stat\"" in line and first_step == "haplotype" :
            line = line.replace("tab-pane disabled" ,"tab-pane active")

        FH_summary_tpl.write(line)

    FH_summary_tpl.close()
    Fh_tpl.close()
