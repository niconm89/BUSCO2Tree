#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:13:22 2022

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
"""
#%% Imports
#import os
import argparse
import os
from time import time
import find_singlecopy_BUSCOs as step1
#from scripts import align_BUSCOs

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
def BUSCO2Tree(STEPS, BUSCODIR, OUTDIR, ODB, LINEAGE):
    step1_dir = ""
    step2_dir = ""
    step3_dir = ""
    step4_dir = ""
    for step in STEPS:
        if step == 1: #Finding single-copy BUSCOs...
            print("Finding single-copy BUSCOs...")
            #01_single-copy 02_alignments 02_Matrix 03_Tree
            try:
                if not os.path.isdir(OUTDIR):
                    os.mkdir(OUTDIR)
                step1_dir = os.path.join(OUTDIR, "01_single-copy")
                os.mkdir(step1_dir) #creating output dir OUTDIR/01_single-copy
            except:
                    raise RuntimeError("Output directory can not be created. Check your path!")
            #first we find common single-copy BUSCOs among genomes
            genomes_names, common_busco_ids = step1.find_singlecopy(BUSCODIR, step1_dir, ODB, LINEAGE)
            #now we create multifasta files with common single-copy BUSCOs
            step1.create_busco_fasta(BUSCODIR, step1_dir, ODB, LINEAGE, genomes_names, common_busco_ids, common_busco_dir="common_busco_sequences")
        #print("Finding single-copy BUSCOs...")
        if step == 2: #Finding single-copy BUSCOs...
            print("Aligning common single-copy BUSCOs...")
            step2_dir 
#end


#%% Menu -> is executed when the script is called independently
def usage():
	parser = argparse.ArgumentParser(
		description='''BUSCO2Tree.py helps the processing of BUSCO outputs to create a phylogenetic tree.''',
		epilog="""End of the help""")

	parser.add_argument('-s', '--steps', nargs='+', type=int, required=True, help='Number of the steps to run. Only consecutive steps are possible to select. Options are: 1: find common single copy BUSCOs; 2: align common single copy BUSCOs; 3: create phylogenetic matrix by concatenating BUSCO alignmets; 4: Generate the phylogenetic tree.')
	parser.add_argument('-i', '--inputdir', type=str, required=False, help='Path to the directory where individual BUSCO outputs are located (you can select the BUSCO output when running in batch mode). Directories names will be taken as the genomes names.')
	parser.add_argument('-o', '--outdir', type=str, required=False, help='Path to the directory where the output of each step will be saved. If it does not exists, it will be created.')
	parser.add_argument('-d', '--odb', type=int, required=False, default=10, help='OrthoDB version. This parameter is only use to complete the path of the single copy sequences directory.')
	parser.add_argument('-l', '--lineage', type=str, required=False, default='lepidoptera', help='OrthoDB lineage. This parameter is only use to complete the path of the single copy sequences directory.')
	parser.add_argument('-c', '--config', metavar='<config file>', type=str, required=False, help='Config file for alignment setting. You can find a template config file in the docs directory. If no file is provided, the alignments will be done using default parameters.')
	parser.add_argument('-a', '--aligndir', type=str, required=False, help='Path to the directory containing the alignments in fasta format.')
	parser.add_argument('-f', '--outformat', type=str, required=False, choices=['nexus','phylip'], default='phylip', help='Format for the phylogenetic matrix.')
	parser.add_argument('-P', '--prefix', type=str, required=False, default="iqtree", help='Prefix to name the output dataset and results.')
	parser.add_argument('-b', '--bootstrap', type=int, required=False, default="1000", help='Number of bootstrap replicate to run.')
	parser.add_argument('-t', '--threads', type=int, required=False, default="16", help='Number of threads to use.')
	
	return parser.parse_args()
#end
#%% Main program
if __name__ == '__main__':
    start = time() #time 0
    args = usage()	
    #print(args)
    #first we validate that steps introduced are ok...
    print("Validating steps selected...")
    print(args.steps)
    validate_steps(args.steps)
    #now we see which steps we have to do...
    print("\t\tSteps " + str(args.steps) + " selected.")
    #validating parameters
    print("Validating parameters...")
    #DO SOMETHING TO VALIDATE PARAMETERS
    print("Running BUSCO2Tree...")
    BUSCO2Tree(args.steps, args.inputdir, args.outdir, args.odb, args.lineage)
    print("BUSCO2Tree has finished.")
    print(f'Time taken to run: {time() - start} seconds.')
#%% End




