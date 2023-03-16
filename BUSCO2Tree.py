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
#import find_singlecopy_BUSCOs as step1
#import align_BUSCOs as step2
#import create_matrix as step3
#import phylogenetic_analysis as step4
from scripts import find_singlecopy_BUSCOs as step1
from scripts import align_BUSCOs as step2
from scripts import create_matrix as step3
from scripts import phylogenetic_analysis as step4

#%% Functions definition
def validate_steps(STEPS):
    for n in STEPS:
        if n > 4 or n < 1:
            raise ValueError("You must enter step numbers between 1 and 4...")
    if not checkConsecutive(STEPS):
        raise ValueError("You must enter consecutive step numbers (e.g. '1 2 3 4', '1 2', '1 2 3', '2 3 4', '3 4')")
#end
def validate_params(ARGUMENTS):
	mandatory_args = {}
	if 1 in ARGUMENTS.steps:
		if not BUSCODIR:
			raise RuntimeError("The buscodir parameter was not introduced.")
	if 2 in ARGUMENTS.steps:
		if not arguments.fastadir:
			raise RuntimeError("The fastadir parameter was not introduced.")
	if 3 in ARGUMENTS.steps:
		if not ARGUMENTS.aligndir:
			raise RuntimeError("The aligndir parameter was not introduced.")
	if 4 in ARGUMENTS.steps:
		if not ARGUMENTS.matrix and not ARGUMENTS.partitions:
			raise RuntimeError("You must introduce the matrix and partitions files when selecting step 4.")
			
#end
def checkConsecutive(l):
    return sorted(l) == list(range(min(l), max(l)+1))
#end
def BUSCO2Tree(STEPS, BUSCODIR, OUTDIR, ODB, LINEAGE, FASTADIR, CONFIG, COMMAND, TRIM, TRIMPARAMS, ALIGNDIR, FORMAT, MATRIX, PARTITIONS, SEQTYPE, PREFIX, BOOTSTRAP, THREADS):
	#BUSCO2Tree(args.steps, args.inputdir, args.outdir, args.odb, args.lineage, 
	#args.aligndir, args.config, args.command, args.trim, args.trimparams, 
	#args.aligndir, args.outformat)
	#args.matrix, args.partitions, args.outdir, args.seqtype, args.prefix, args.bootstrap, args.threads)
	try:
		if not os.path.isdir(OUTDIR):
			os.mkdir(OUTDIR)
	except:
		raise RuntimeError("Output directory can not be created. Check your path!")
	step1_dir = ""
	step2_dir = ""
	step3_dir = ""
	step4_dir = ""
	for step in STEPS:
		if step == 1: #Finding single-copy BUSCOs...
			print("Finding single-copy BUSCOs...")
			#01_single-copy 02_alignments 02_Matrix 03_Tree
			step1_dir = os.path.join(OUTDIR, "01_single-copy")
			os.mkdir(step1_dir) #creating output dir OUTDIR/01_single-copy
            #first we find common single-copy BUSCOs among genomes
			genomes_names, common_busco_ids = step1.find_singlecopy(BUSCODIR, step1_dir, ODB, LINEAGE)
			#now we create multifasta files with common single-copy BUSCOs
			step1.create_busco_fasta(BUSCODIR, step1_dir, ODB, LINEAGE, genomes_names, common_busco_ids, common_busco_dir="common_busco_sequences")
		#print("")
		if step == 2: #Aligning BUSCOs multifasta files...
			print("Aligning common single-copy BUSCOs...")
			step2_dir = os.path.join(OUTDIR, "02_alignments")
			os.mkdir(step2_dir) #creating output dir OUTDIR/02_alignments
			
			if COMMAND:
				print("MAFFT will be excecuted using a command introduced by the user.")
				step2.align_command_mafft(FASTADIR, step2_dir, COMMAND)
			else:
				print("MAFFT will be excecuted using a configuration file.")
				step2.align_config_mafft(FASTADIR, step2_dir, CONFIG)
			print("Alignments completed...")
			if TRIM: #trimming alingments
				print("Trimming alignments...")
				step2.trim_alns(step2_dir, TRIMPARAMS)
		if step == 3:
			print("Generating the phylogenetic matrix.")
			step3_dir = os.path.join(OUTDIR, "03_matrix")
			os.mkdir(step3_dir) #creating output dir OUTDIR/03_matrix
			step3.cat_alignments(ALIGNDIR, step3_dir, FORMAT)
		if step == 4:
			print("Creating the phylogenetic tree with IQTree.")
			step4_dir = os.path.join(OUTDIR, "04_phylogenetic_tree")
			os.mkdir(step4_dir) #creating output dir OUTDIR/04_phylogenetic_tree
			step4.model_partitions(args.matrix, args.partitions, step4_dir, args.seqtype, args.prefix, args.bootstrap, args.threads)
			#step4.run_IQTree_manual(args.matrix, args.partitions, args.outdir, args.command)
#end


#%% Menu -> is executed when the script is called independently
def usage():
	parser = argparse.ArgumentParser(
		description='''BUSCO2Tree.py helps the processing of BUSCO outputs to create a phylogenetic tree.''',
		epilog="""End of the help""")
	#general arguments
	parser.add_argument('-s', '--steps', nargs='+', type=int, required=True, help='Number of the steps to run. Only consecutive steps are possible to select. Options are: 1: find common single copy BUSCOs; 2: align common single copy BUSCOs; 3: create phylogenetic matrix by concatenating BUSCO alignmets; 4: Generate the phylogenetic tree.')
	parser.add_argument('-o', '--outdir', type=str, required=True, help='Path to the directory where the output of each step will be saved. If it does not exists, it will be created.')
	parser.add_argument('-t', '--threads', type=int, required=False, default=4, help='Number of threads to use.')
	#step 1
	parser.add_argument('-b', '--buscodir', type=str, required=False, help='Path to the directory where individual BUSCO outputs are located (you can select the BUSCO output when running in batch mode). Directories names will be taken as the genomes names.')
	parser.add_argument('-d', '--odb', type=int, required=False, default=10, help='OrthoDB version. This parameter is only use to complete the path of the single copy sequences directory.')
	parser.add_argument('-l', '--lineage', type=str, required=False, default='lepidoptera', help='OrthoDB lineage. This parameter is only use to complete the path of the single copy sequences directory.')
	#step 2
	parser.add_argument('-f', '--fastadir', type=str, required=False, help='Path to the directory containing the alignments in fasta format.')
	parser.add_argument('-cnf', '--config', metavar='<config file>', type=str, required=False, help='Config file for alignment setting. You can find a template config file in the docs directory. If no file is provided, the alignments will be done using default parameters.')
	parser.add_argument('-cmd', '--command', metavar='<command>', type=str, required=False, help='Command (between quote marks) to use in the MAFFT alingment. The names of input and out files must be avoided, e.g. "mafft --thread 8 --unalignlevel 0.1 --leavegappyregion --ep 0.12 --globalpair --maxiterate 1000".')
	parser.add_argument('--trim', action='store_true', required=False, help='Trim alignments using TrimAl in automated mode, which must be available in the path.')
	parser.add_argument('--trimparams', required=False, help='TrimAl parameters to use. They must be passed in quotes marks and avoiding input and output files, e.g. "-gt 0.3 -nogaps -phylip".')
	#step 3
	parser.add_argument('-a', '--aligndir', type=str, required=False, help='Path to the directory containing the alignments in fasta format (step 3).')
	parser.add_argument('-fmt', '--format', type=str, required=False, choices=['nexus','phylip'], default='phylip', help='Format for the phylogenetic matrix generated by concatenating BUSCOs alignment files (step 3).')
	#step 4
	parser.add_argument('-m', '--matrix', type=str, required=True, help='Path to the phylogenetic matrix file.')
	parser.add_argument('-p', '--partitions', type=str, required=True, help='Path to the partitions file in nexus format.')
	parser.add_argument('-P', '--prefix', type=str, required=False, default="iqtree", help='Prefix to name the output dataset and results.')
	parser.add_argument('-B', '--bootstrap', type=int, required=False, default="1000", help='Number of bootstrap replicate to run.')
	parser.add_argument('-st', '--seqtype', type=str, required=False, choices=["DNA","AA"], default="AA", help='Type of sequence.')
	
	return parser.parse_args()
#end
#%% Main program
if __name__ == '__main__':
    start = time() #time 0
    args = usage()	
    #print(args)
    #first we validate that steps introduced are ok...
    print("Validating steps selected...")
    validate_steps(args.steps)
	print("Validating parameters...")
    validate_params(args)
    #now we see which steps we have to do...
    print("\t\tSteps " + str(args.steps) + " selected.")
    #validating parameters
    print("Validating parameters...")
    #DO SOMETHING TO VALIDATE PARAMETERS
    print("Running BUSCO2Tree...")
    BUSCO2Tree(args.steps, args.buscodir, args.outdir, args.odb, args.lineage, args.fastadir, args.config, args.command, args.trim, args.trimparams, args.aligndir, args.format, args.matrix, args.partitions, args.seqtype, args.prefix, args.bootstrap, args.threads)
    print("BUSCO2Tree has finished.")
    print(f'Time taken to run: {time() - start} seconds.')
#%% End

