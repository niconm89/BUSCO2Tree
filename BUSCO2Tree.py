#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:13:22 2022

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
"""
#%% Imports
#import os
import argparse
from time import time

#%% Class definition
#%% Functions definition
def validate_steps(STEPS):
	for n in STEPS:
		if n > 4 or n < 1:
			raise ValueError("You must enter step numbers between 1 and 4...")
	if not checkConsecutive(STEPS):
		raise ValueError("You must enter consecutive step numbers (e.g. '1 2 3 4', '1 2', '1 2 3', '2 3 4', '3 4')")
#end
def checkConsecutive(l):
    return sorted(l) == list(range(min(l), max(l)+1))
#end
#%% Menu -> is executed when the script is called independently
def usage():
	parser = argparse.ArgumentParser(
		description='''01_find_singlecopy.py recover ID sequences from a BUSCO output directory.''',
		epilog="""End of the help""")

	parser.add_argument('-s', '--steps', nargs='+', type=int, required=True, help='Number of the steps to run. Only consecutive steps are possible to select. Options are: 1: find common single copy BUSCOs; 2: align common single copy BUSCOs; 3: create phylogenetic matrix from BUSCOs; 4: Do phylogenetic Analysis.')
	parser.add_argument('-b', '--buscodir', type=str, required=False, help='Path to the BUSCO output directory. Individual genome directory must be place within this path. Directory names will be take as genomes names.')
	parser.add_argument('-o', '--outdir', type=str, required=False, help='Path to save output files. If it does not exists, it will be created.')
	parser.add_argument('-d', '--odb', type=int, required=False, default=10, help='OrthoDB version. This parameter is only use to complete the path of the single copy sequences directory.')
	parser.add_argument('-l', '--lineage', type=str, required=False, default='lepidoptera', help='OrthoDB lineage. This parameter is only use to complete the path of the single copy sequences directory.')
	parser.add_argument('-c', '--config', metavar='<config file>', type=str, required=False, help='Config file for alignment setting. You can find a template config file in the docs directory. If no file is provided, the alignments will be done using default parameters.')
	parser.add_argument('-a', '--aligndir', type=str, required=False, help='Path to the directory containing the alignments in fasta format.')
	parser.add_argument('-f', '--outformat', type=str, required=False, choices=['nexus','phylip'], default='phylip', help='Format for the phylogenetic matrix.')
	parser.add_argument('-P', '--prefix', type=str, required=False, default="iqtree", help='Prefix to name the output dataset and results.')
	#parser.add_argument('-b', '--bootstrap', type=int, required=False, default="1000", help='Number of bootstrap replicate to run.')
	parser.add_argument('-t', '--threads', type=int, required=False, default="16", help='Number of threads to use in IQ-Tree.')
	
	return parser.parse_args()
#end
#%% Main program
if __name__ == '__main__':
	start = time()
	args = usage()	
	#print(args)
	#first we validate that steps introduced are ok...
	validate_steps(args.steps)
	#now we see which steps we have to do...
	
	
	print(f'Time taken to run: {time() - start} seconds.')
#%% End




