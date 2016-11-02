#!/usr/bin/env python2.7
# -*-coding:Latin-1 -*
########################################################################
## IMPORT
import argparse
from random import randint

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

def make_indiv_list(indiv_list, nb_indiv):
	full_list=indiv_list
	random_list=list()
	
	for i in range(0,nb_indiv):
		r = randint(0,len(full_list))
		random_list.append(full_list[r])
		full_list.pop(r)
		
	return random_list 
########################################################################
## MAIN
parser = argparse.ArgumentParser(description="generate resume or detailed allele stastics and genotypes statistics")
#params
parser.add_argument('-s', '--snp', default=False, help="snp diversity")
parser.add_argument('-a', '--allel', default=False , help="allel diversity")
parser.add_argument('--debug', default=False, action='store_true', help="Keep temporary files to debug program.")
# Inputs
group_input = parser.add_argument_group('Inputs')
group_input.add_argument('-p', '--population-map', required=True, help='population map.')
group_input.add_argument('-i', '--haplotype-tsv', required=True, help='Stacks haplotype tsv file')
group_input.add_argument('-o', '--output-directory', default='', help='Output directory')
args = parser.parse_args()
prevent_shell_injections(args)