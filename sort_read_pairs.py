#!/usr/bin/env python2.7
## -*-coding:Latin-1 -*
# Copyright (C) 2014 INRA
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
from test.test_support import args_from_interpreter_flags

__author__ = 'Maria Bernard - SIGENAE'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'maria.bernard@jouy.inra.fr'

import threading
import multiprocessing
import glob
import time
import argparse
import re
from sort_read_pairs_utils import *

##################################################################################################################################################
#
# CLASS
#
##################################################################################################################################################



##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

def store_read_files(sample_dir, sample_file_list, read_type, read_format, read1_files, read2_files, sample_names, log_file) :
    """
    @summary: Create a dictionnary of file from sample_dir in read_files.
    @param sample_dir: [str] input directory path.
    @param sample_file_list: [str] if not None, only select sample in this file.
    @param read_type : [str] either single or pair
    @param read_format : [str] either fastq or gzfastq
    @param read1_files: [list] output read1 file list order by sample : [sample1_1.fq, sample2_1.fq ... ]
    @param read2_files: [list] output read2 file list order by sample : [sample1_2.fq, sample2_2.fq ... ]
    @param sample_names: [list] output sample name list
    @param log_file : [str] path to log file
    """
    tmp_dict=dict()
    sample_select = list()
    if not sample_file_list is None:
         FH = open(sample_file_list,"r")
         sample_select = [n.strip() for n in FH.readlines()]
    
    list_dir=list()
    if read_format == "gzfastq":
        list_dir = glob.glob(os.path.join(sample_dir,"*.f*q.gz"))
    else:
        list_dir = glob.glob(os.path.join(sample_dir,"*.f*q"))
    
    for f in list_dir:
        match = re.search("(?P<sample_name>\S+)(.|_)R*(?P<num_read>(1|2))\.f.*q\.gz",os.path.basename(f)) if read_format == "gzfastq" else re.search("(?P<sample_name>\S+)(.|_)R*(?P<num_read>(1|2))\.f.*q",os.path.basename(f))
        if match : 
            sample_name = match.group("sample_name")
            num_read = match.group("num_read")
            if sample_name in sample_select or len(sample_select)==0:
                if not sample_name in sample_names:
                    sample_names.append(sample_name)
                    
                if sample_name in tmp_dict:
                    tmp_dict[sample_name][num_read] = f
                else:
                    tmp_dict[sample_name]={num_read : f}
        else :
            match = re.search("(?P<sample_name>\S+)\.f.*q\.gz",os.path.basename(f)) if read_format == "gzfastq" else re.search("(?P<sample_name>\S+)\.f.*q",os.path.basename(f))
            if match:
                sample_name = match.group("sample_name")
                if sample_name in sample_select or len(sample_select)==0:
                    if not sample_name in sample_names:
                        sample_names.append(sample_name)
                        
                    if sample_name in tmp_dict:
                        tmp_dict[sample_name]["1"] = f
                    else:
                        tmp_dict[sample_name]={"1" : f}
                    
    # check sample input files
    log = Logger(log_file)
    err=False            
    for sample in sample_names : 
        if not "1" in tmp_dict[sample]:
            log.write("sample "+sample+" missing read 1 file\n")
            err=True
        else:
            read1_files.append(tmp_dict[sample]['1'])
        if read_type == "pair" :
            if not "2" in tmp_dict[sample]:
                log.write("sample "+sample+" missing read 2 file\n")
                err=True
            else:
                read2_files.append(tmp_dict[sample]['2'])           
    log.close()    
    if err : 
        raise ValueError("some samples are missing read files\nSee "+log_file+" for details\n")
            
def store_stacks_files(stacks_dir, sample_names, tags_files,matches_files, log_file):
    """
    @summary: Create a dictionnary of file from stacks_dir in stacks_files.
    @param : stacks_dir : [str] input directory path.
    @param sample_names: [list] sample name list
    @param tags_files: [list] output list order by sample : [sample1.tags.tsv, sample2.tags.tsv, ... ]
    @param matches_files: [list] output list order by sample : [sample1.matches.tsv, sample2.matches.tsv, ... ]
    @param log_file : [str] path to log file
    """
    tmp_dict = dict()
    list_dir=list()
    list_dir = glob.glob(os.path.join(stacks_dir,"*.tsv*"))
    
    for f in list_dir:
        match = re.search("(?P<sample_name>\S+)\.(?P<stacks_file>(tags|matches))\.tsv.*",os.path.basename(f))
        if match:
            sample_name = re.sub("\.1$","",match.group("sample_name"))
            stacks_file = match.group("stacks_file")
            if sample_name in sample_names:
                if sample_name in tmp_dict:
                    tmp_dict[sample_name][stacks_file]=f
                else:
                    tmp_dict[sample_name]={stacks_file : f }
    
    # check sample stacks input file
    log = Logger(log_file)
    err=False             
    for sample in sample_names:
        if not sample in tmp_dict:
            log.write("sample "+sample+" missing stacks files\n")
            err=True
        else:
            if not "tags" in tmp_dict[sample]:
                log.write("sample "+sample+" missing tags stacks files\n")
                err=True
            else:
                tags_files.append(tmp_dict[sample]["tags"])        
            if not "matches" in tmp_dict[sample]:
                log.write("sample "+sample+" missing matches stacks files\n")
                err=True
            else:
                matches_files.append(tmp_dict[sample]["matches"])
    log.close()    
    if err : 
        raise Exception("some samples are missing stacks files\nSee log_file for details\n")
       
def store_locus_id(locus_file_list, locus_id, locus_list):
    """
    @summary: store de locus id to treat
    @param locus_file_list or locus_id: [str] path the locus file list or one locus id
    @parma locus_list: [list] output locus list
    """
    if not locus_file_list is None:
        FH = open(locus_file_list)
        for line in FH.readlines():
            locus_list.append(line.strip())
        FH.close()
    elif not locus_id is None:
        locus_list.append(locus_id)

def parse_multiple_samples(sample_names, read1_files, read2_files, tags_files, matches_files, log_files, locus_list, tmp_files, args) :
    """
    @summary: parse each sample to divide its reads into split locus file.
    @param: sample_names : [list] list of sample to treat
    @param read1_files: [list] list of read1 path files
    @param read2_files: [list] list of read2 path files
    @param tags_files: [list] list of tags path files
    @param matches_files: [list] list of matches path files
    @param log_files: [list] list of log path files
    @param locus_list : [list] of locus id to treat
    @param locus_read1_files / locus_read2_files : [dict] dictionnary of read1/read2 files for each locus: key = locus, value = list of output path files
    @param tmp_files : [TmpFiles] to manage temporary files
    @param args : [Namespace] Global parameters.
    """
    for idx in range(len(sample_names)):
        if args.read_in_type == "pair" : 
            parse_sample(sample_names[idx], read1_files[idx], read2_files[idx], tags_files[idx], matches_files[idx], log_files[idx], locus_list, tmp_files, args)
        else:
            parse_sample(sample_names[idx], read1_files[idx], None, tags_files[idx], matches_files[idx], log_files[idx], locus_list, tmp_files, args)
 
def check_multimap(matches_dict, sample_name, log_file):
    """
    @summary: check that only one sample locus correspond to one catalog (remove if multimap). In contrario, you can have multiple sample locus that match the same catalog locus
    @param matches_dict: [dict] matches dictionnary : keys = sample locus, values = list of catalog locus
    @param sample_name: [str] sample name
    @param log_file : [str] path to log file   
    """
    rm_locus_list= list()
    tmp_dict=dict()
    
    for sample_locus in matches_dict:
        if matches_dict[sample_locus] in tmp_dict :
            tmp_dict[matches_dict[sample_locus]].append(sample_locus)
        else :  
            tmp_dict[matches_dict[sample_locus]] = [sample_locus]

    c=0
    for cat_locus in tmp_dict:

        if len(tmp_dict[cat_locus]) > 1:
            for sample_locus in tmp_dict[cat_locus]:
                matches_dict.pop(sample_locus)
                c+=1
    if c>0:
        Logger.static_write(log_file, "\t\tWARNING sample "+sample_name+" has "+str(c)+" of its locus that map same catalog locus\n")

def parse_readFile(sample_name, read_file, num_read, reads_dict , tmp_files, args):
    """
    @summary: parse each sample to divide its reads into split locus file.
    @param: sample_name : [str] sample name
    @param read_file: [str] read1 path file
    @param num_read ; [str] either 1 or 2 for read 1 or read 2
    @param reads_dict: [dict] dictionnary of reads cat locus : keys = read_id, values = [cat_locus sample_locus] 
    @param tmp_files : [TmpFiles] to manage temporary files
    @param args : [Namespace] Global parameters.
    """
#     print read_file
#     first=True
    FH_seq = SequenceFileReader.factory( read_file )
    cat_locus_dict=dict()
    for record in FH_seq:
#         if first:
#             print record.id, reads_dict.keys()[0]
#             first=False
        if len(record.id.split("_")) == 5 and re.search("_(1|2)$",record.id) :
                record.id = record.id[:-2]
        if record.id in reads_dict:
            cat_locus = reads_dict[record.id][0]
            sample_locus = reads_dict[record.id][1]
            id = cat_locus+"|"+sample_name+"|"+sample_locus+"|"+record.id
            record.id = id
            if cat_locus in cat_locus_dict : 
                cat_locus_dict[cat_locus].append(record)
            else:
                cat_locus_dict[cat_locus]=[record]
            #
            #    Ecriture des record si 200 lectures déjà stockées pour le locus en cours
            #
            if len(cat_locus_dict[cat_locus]) > 200:
                if not os.path.join(args.outdir,tmp_files.prefix+"_"+cat_locus+"_"+sample_name+"_"+num_read+"."+args.read_out_format) in tmp_files.files:
                    f = tmp_files.add(args.outdir,cat_locus+"_"+sample_name+"_"+num_read+"."+args.read_out_format)
                if args.read_out_format=="fasta":
                    FH_out_seq = FastaIO(f,"a")
                elif args.read_out_format=="fastq":
                    FH_out_seq = FastqIO(f,"a")
                for r in cat_locus_dict[cat_locus]:
                    FH_out_seq.write(r)
                FH_out_seq.close()
                cat_locus_dict[cat_locus] = []
    #
    #    Ecriture des record restant encore dans cat_locus_dict
    #
    for cat_locus in cat_locus_dict:
        if not os.path.join(args.outdir,tmp_files.prefix+"_"+cat_locus+"_"+sample_name+"_"+num_read+"."+args.read_out_format) in tmp_files.files:
            f=tmp_files.add(cat_locus+"_"+sample_name+"_"+num_read+"."+args.read_out_format)
        if args.read_out_format=="fasta":
            FH_out_seq = FastaIO(f,"a")
        elif args.read_out_format=="fastq":
            FH_out_seq = FastqIO(f,"a")
        for r in cat_locus_dict[cat_locus]:
            FH_out_seq.write(r)
        FH_out_seq.close()
        cat_locus_dict[cat_locus] = []
        
def parse_sample(sample_name, read1_file, read2_file, tags_file, matches_file, log_file, locus_list, tmp_files, args) :
    """
    @summary: parse each sample to divide its reads into split locus file.
    @param: sample_name : [str] sample name
    @param read1_file: [str] read1 path file
    @param read2_file: [str] read2 path file
    @param tags_file: [str] tags path file
    @param matches_file: [str] matches path file
    @param log_file : [str] path to log output file
    @param locus_list : [list] of locus id to treat
    @param tmp_files : [TmpFiles] to manage temporary files
    @param args : [Namespace] Global parameters.
    """
    # dictionnary of matches locus : keys = sample locus, values = list of catalog locus: len= 1 after check_multimap
    Logger.static_write(log_file, "## Parsing sample : "+sample_name+"\n")
    start_time = time.strftime("%d %b %Y %H:%M:%S", time.localtime())
    Logger.static_write(log_file,"\tparse matches result...")
    matches_dict = dict()
    
    if is_gzip(matches_file):
        matches_handle = gzip.open( matches_file )
    else:
        matches_handle = open(matches_file)
        
    for line in matches_handle.readlines():
        cat_locus = line.strip().split()[2]
        sample_locus = line.strip().split()[4]
        if len(locus_list)==0 or cat_locus in locus_list:
            matches_dict[sample_locus] = cat_locus
    Logger.static_write(log_file,str(len(matches_dict))+" locus correspondancies stored\n")        
    check_multimap ( matches_dict,sample_name, log_file)
    if len(matches_dict) == 0 :
        Logger.static_write(log_file,"\tNo need to go forward for this sample\n")
    
    else:
        # dictionnary of reads cat locus : keys = read_id, values = [cat_locus sample_locus]
        reads_dict = dict()
        Logger.static_write(log_file,"\tparse tags result...")
        if is_gzip(tags_file):
            tags_handle = gzip.open( tags_file )
        else:
            tags_handle = open(tags_file)
            
        for line in tags_handle.readlines():
            sample_locus=line.strip().split()[2]
            if sample_locus in matches_dict:
                cat_locus = matches_dict[sample_locus]
                if re.search("(primary|secondary)",line):
                    read_id = line.strip().split("\t")[8].split()[0]
                    if len(read_id.split("_"))==5 and read_id.endswith("_1"):
                        read_id = read_id[:-2]
                    reads_dict[read_id] = [cat_locus,sample_locus]
        if len(reads_dict) == 0 and len(matches_dict)!=0:
            raise Exception ("Problem with sample "+sample_name+" "+str(len(matches_dict))+" locus correspondancies but no reads found\n")
        
        Logger.static_write(log_file,str(len(reads_dict))+" reads stored\n")
        
        # parse read1 file and split reads into locus files
        # store record list into cat_locus_dict, when len record > 200 > write into file
        Logger.static_write(log_file,"\tparsing sample read 1 file ... \n")
        parse_readFile(sample_name, read1_file, "1", reads_dict , tmp_files, args)
    
        # parse read2 file if not None and split reads into locus files
        # store record list into cat_locus_dict, when len record > 200 > write into file
        if not read2_file is None:
            Logger.static_write(log_file,"\tparsing sample read 2 file ... \n")
            parse_readFile(sample_name, read2_file, "2", reads_dict , tmp_files, args)

    Logger.static_write(log_file,"\tstart time : "+start_time+"; end time : "+ time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + "\n\n")
            
def get_fastq_nb_seq( fastq_file ):
    """
    @summary: Returns the number of sequences in fastq_file.
    @param fastq_file: [str] Path to the fastq file processed.
    @return: [int] The number of sequences.
    """
    FH_input = None
    if not is_gzip(fastq_file):
        FH_input = open( fastq_file )
    else:
        FH_input = gzip.open( fastq_file )
    nb_line = 0
    for line in FH_input:
        nb_line += 1
    FH_input.close()
    nb_seq = nb_line/4
    return nb_seq

def get_fasta_nb_seq( fasta_file ):
    """
    @summary: Returns the number of sequences in fasta_file.
    @param fasta_file: [str] Path to the fasta file processed.
    @return: [int] The number of sequences.
    """
    FH_input = None
    if not is_gzip(fasta_file):
        FH_input = open( fasta_file )
    else:
        FH_input = gzip.open( fasta_file )
    nb_seq = 0
    for line in FH_input:
        if line.startswith(">"):
            nb_seq += 1
    FH_input.close()
    return nb_seq
      

def concatenate_readFiles(locus, sample_names, num_read, prefix, args) :
    """
    @summary: concatenante individuals locus_sample read files into one locus file
    @param locus: [str] locus to treat
    @param sample_names: [list] list of sample names (to treat all locus in the same order)
    @param num_read : [str] to select read1 or read2 file
    @param args [NameSpace]: program's argument, notably outdir or read out format    
    """
    count = 0
    filenames = list()
    for sample in sample_names:
        list_dir = glob.glob(os.path.join(args.outdir,prefix+"_"+locus+"_"+sample+"_"+num_read+"."+args.read_out_format))
        if len(list_dir) ==1:
            f=list_dir[0]
            nb_seq = get_fasta_nb_seq(f) if args.read_out_format == "fasta" else get_fastq_nb_seq(f)
            filenames.append(f)
            count += nb_seq
    
    with open(os.path.join(args.outdir,locus+"_"+num_read+"."+args.read_out_format), 'w') as FH_out_seq:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    FH_out_seq.write(line)
    return count

def summariseResults(sample_names, locus_list, prefix, args) : 
    """
    @summary: concatenate locus sample read files into one, and calculate number of reads per locus
    @param sample_names : [list] list of sample name to treat each locus in the same order
    @param locus_list: [list] list of locus to treat. if empty treat all locus
    @param args : [NameSpace] program's argument, notably outdir or read out format    
    """
    
    if len(locus_list) == 0:
        list_dir=list()
        list_dir = glob.glob(os.path.join(args.outdir,prefix+"_*_*_1."+args.read_out_format))
        for f in list_dir:
            locus = f.split("_")[0]
            if not locus in locus_list:
                locus_list.append(locus)

    Logger.static_write(args.log_file,"##Final concatenation \n\n")
    start_time = time.strftime("%d %b %Y %H:%M:%S", time.localtime())
    for locus in locus_list :
        count1 = concatenate_readFiles(locus, sample_names, "1" , prefix, args)
        if args.read_in_type == "pair" :
            count2=concatenate_readFiles(locus, sample_names, "2" , prefix, args)
            if not count1 == count2:
                raise Exception ("Not the same number of reads between reads 1 and reads 2 in locus :"+locus+"\n")
            Logger.static_write(args.log_file,"\tlocus "+locus+" has "+str(count1)+" pair of reads\n")
        else:
            Logger.static_write(args.log_file,"\tlocus "+locus+" has "+str(count1)+" reads\n")
            
    Logger.static_write(args.log_file,"start time : "+start_time+";  end time : "+ time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + "\n")   

def log_append_files( log_file, appended_files ):
    """
    @summary: Append content of several log files in one log file.
    @param log_file: [str] The log file where contents of others are appended.
    @param appended_files: [list] List of log files to append.
    """
    FH_log = Logger( log_file )
    FH_log.write( "\n" )
    for current_file in appended_files:
        FH_input = open(current_file)
        for line in FH_input:
            FH_log.write( line )
        FH_input.close()
        FH_log.write( "\n" )
    FH_log.close()  
##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Sort reads by locus instead of samples")
    parser.add_argument('-r', '--read-in-type', type=str, choices=['single', 'pair'], default="pair",  help="Read input type, either single or pair (default)")
    parser.add_argument('-t', '--read-in-format', type=str, choices=["fastq", "gzfastq"], default="fastq",  help="Read input format, either fastq (default) or gzfastq ")
    parser.add_argument('-f', '--read-out-format', type=str, choices=['fasta', 'fastq'], default="fastq",  help="Read output format, either fasta or fastq (default)")
    parser.add_argument('-n', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used." )
    parser.add_argument('--debug', default=False, action='store_true', help="Keep temporary files to debug program.")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('--stacks-dir', required=True, help='path to Stacks output results')
    group_input.add_argument('--sample-dir', required=True, help='Path to sample reads')
    group_input.add_argument('-s','--sample-list', type=str, required=False, help='Path to sample file list')
    locus_input = group_input.add_mutually_exclusive_group()
    locus_input.add_argument('-l','--locus-list', type=str, required=False, help='Path to locus file list (exclusive with -i/--locus_id)')
    locus_input.add_argument('-i','--locus-id', type=str, required=False, help='one locus ID (exclusive with -l/--locus_list)')
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('--outdir', required=True, help='Path to output directory')
    group_output.add_argument('--log-file', default=sys.stdout, help='The list of commands executed.')
    args = parser.parse_args()
    prevent_shell_injections(args)    

    Logger.static_write(args.log_file, "## Application\nSoftware: " + os.path.basename(sys.argv[0]) + " (version: " + str(__version__) + ")\nCommand: " + " ".join(sys.argv) + "\n\n")

	# Check parameters
    if not os.path.exists(args.stacks_dir):	raise argparse.ArgumentTypeError( "'--stacks-dir option value doe not exist" )
    if not os.path.exists(args.sample_dir):	raise argparse.ArgumentTypeError( "'--sample-dir option value doe not exist" )
    if not args.sample_list is None and not os.path.exists(args.sample_list):    raise argparse.ArgumentTypeError( "'--sample-list option value doe not exist" )
    if not args.locus_list is None and not os.path.exists(args.locus_list):    raise argparse.ArgumentTypeError( "'--locus-list option value doe not exist" )
    if not os.path.exists(args.outdir):    raise argparse.ArgumentTypeError( "'--outdir option value doe not exist" )
    
    tmp_files = TmpFiles( args.outdir )
    
    # Process
    try:
        read1_files=list()
        read2_files=list()
        sample_names=list()
        tags_files=list()
        matches_files=list()
        locus_list=list()
        log_files = list()
        
        ##�store input file by sample
        # reads
        Logger.static_write(args.log_file, "Check sample input reads ... \n")
        start_time = time.strftime("%d %b %Y %H:%M:%S", time.localtime())
        store_read_files(args.sample_dir, args.sample_list, args.read_in_type, args.read_in_format, read1_files, read2_files, sample_names, args.log_file )
        if args.read_in_type == "pair":
            Logger.static_write(args.log_file, "\t"+str(len(sample_names))+" pair of sample files will be treated\n")
        else:
            Logger.static_write(args.log_file, "\t"+str(len(sample_names))+" single sample files will be treated\n")
        if args.debug :
                Logger.static_write(args.log_file,"\tdetails : "+", ".join(sample_names)+"\n")
        Logger.static_write(args.log_file,"\tstart time : "+start_time+";  end time : "+ time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + "\n")
        
        if len(sample_names) == 0 :
            raise Exception("Error: No sample read files found or precised\n")
        
        #stacks results
        Logger.static_write(args.log_file, "\nCheck stacks input results ... \n")
        start_time = time.strftime("%d %b %Y %H:%M:%S", time.localtime())
        store_stacks_files(args.stacks_dir, sample_names, tags_files, matches_files, args.log_file)
        Logger.static_write(args.log_file,"\tstart time : "+start_time+";  end time : "+ time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + "\n")
        
        ##�store input locus list (if precised)
        Logger.static_write(args.log_file, "\nSelect locus to treat ... \n")        
        if not args.locus_list is None or not args.locus_id is None:
            start_time = time.strftime("%d %b %Y %H:%M:%S", time.localtime())
            store_locus_id(args.locus_list, args.locus_id, locus_list)
            Logger.static_write(args.log_file, "\t"+str(len(locus_list))+" locus will be treated ... \n")
            if args.debug :
                Logger.static_write(args.log_file,"\tdetails : "+", ".join(locus_list)+"\n")
            Logger.static_write(args.log_file,"\tstart time : "+start_time+";  end time : "+ time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + "\n")
        else:
            Logger.static_write(args.log_file, "\tall locus will be treated ... \n")
        
        ## store temporary log files
        log_files = [tmp_files.add(current_sample + '_log.txt') for current_sample in sample_names]
    
        # Process
        nb_processses_used = min( len(read1_files), args.nb_cpus )
        if nb_processses_used == 1:
            parse_multiple_samples(sample_names, read1_files, read2_files, tags_files, matches_files, log_files, locus_list, tmp_files, args)
        else:
            processes = [{'process':None,'sample_names': [], 'read1_files':[], 'read2_files':[], 'tags_files':[],'matches_files':[], 'log_files':[]} for idx in range(nb_processses_used)]
            # Set processes
            for idx in range(len(sample_names)):
                process_idx = idx % nb_processses_used
                processes[process_idx]['sample_names'].append(sample_names[idx])
                processes[process_idx]['read1_files'].append(read1_files[idx])
                if args.read_in_type == "pair":
                    processes[process_idx]['read2_files'].append(read2_files[idx])
                processes[process_idx]['tags_files'].append(tags_files[idx])
                processes[process_idx]['matches_files'].append(matches_files[idx])
                processes[process_idx]['log_files'].append(log_files[idx])
            # Launch processes
            for current_process in processes:
                if idx == 0: # First process is threaded with parent job
                    current_process['process'] = threading.Thread( target=parse_multiple_samples, 
                                                                   args=(current_process['sample_names'], current_process['read1_files'], 
                                                                         current_process['read2_files'], current_process['tags_files'], 
                                                                         current_process['matches_files'], current_process['log_files'], 
                                                                         locus_list, tmp_files, args) )
                else: # Others processes are processed on different CPU
                    current_process['process'] = multiprocessing.Process( target=parse_multiple_samples, 
                                                                   args=(current_process['sample_names'], current_process['read1_files'], 
                                                                         current_process['read2_files'], current_process['tags_files'], 
                                                                         current_process['matches_files'], current_process['log_files'], 
                                                                         locus_list, tmp_files, args) )
                current_process['process'].start()
            # Wait processes end
            for current_process in processes:
                current_process['process'].join()
            # Check processes status
            for current_process in processes:
                if issubclass(current_process['process'].__class__, multiprocessing.Process) and current_process['process'].exitcode != 0:
                    raise Exception( "Error in sub-process execution." )

        # Write summary
        log_append_files( args.log_file, log_files )
        summariseResults( sample_names, locus_list, tmp_files.prefix, args) 
        
        
        
    # Remove temporary files
    finally:
        if not args.debug:
            tmp_files.deleteAll()	
            
            list_dir = glob.glob(os.path.join(tmp_files.tmp_dir,tmp_files.prefix+"*."+args.read_out_format))
            if len(list_dir) > 0:
                for f in list_dir:
                    os.remove(f)