#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:13:22 2022

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
"""
#%% Imports
import argparse
import os
from time import time
from scripts import find_singlecopy_BUSCOs as step1
from scripts import align_BUSCOs as step2
from scripts import create_matrix as step3
from scripts import phylogenetic_analysis as step4

#%% Function definition
def validate_steps(STEPS):
    for n in STEPS:
        if n > 4 or n < 1:
            raise ValueError("The available steps range from 1 to 4!!!")
    if not checkConsecutive(STEPS):
        raise ValueError("Only consecutive step numbers are allowed, e.g. '1 2 3 4', '1 2', '1 2 3', '2 3 4', '3 4'.")    
#end validate_steps
def validate_params(ARGUMENTS):
	'It receives the object of arguments and validates mandatory arguments to run the pipeline.'
	mandatory_args = {}
	if 1 in ARGUMENTS.steps:
		if not ARGUMENTS.buscodir:
			raise RuntimeError("The directory containing BUSCO outputs was not defined.")
	if 2 in ARGUMENTS.steps and 1 not in ARGUMENTS.steps:
		if not ARGUMENTS.fastadir:
			raise RuntimeError("The directory contining common BUSCO groups in fasta format was not defined.")
	if 3 in ARGUMENTS.steps and 2 not in ARGUMENTS.steps:
		if not ARGUMENTS.aligndir:
			raise RuntimeError("The directy containing alignments was not defined.")
	if 4 in ARGUMENTS.steps and 3 not in ARGUMENTS.steps:
		if not ARGUMENTS.matrix:
			raise RuntimeError("The matrix file was not defined.")
		if not ARGUMENTS.partitions:
			raise RuntimeError("The partitions file was not defined.")
#end validate_params
def checkConsecutive(l):
    'It receives a list of steps selected to run and return TRUE or FALSE whether the steps are consecutive or not. On the last case, the program will exit.'
    return sorted(l) == list(range(min(l), max(l)+1))
#end
def BUSCO2Tree(STEPS, BUSCODIR, OUTDIR, ODB, LINEAGE, FASTADIR, CONFIG, COMMAND, TRIM, TRIMPARAMS, ALIGNDIR, FORMAT, MATRIX, PARTITIONS, SEQTYPE, PREFIX, BOOTSTRAP, THREADS):
	'''
	Main function that receives several arguments and run this pipeline.
	BUSCO2Tree(args.steps, args.inputdir, args.outdir, args.odb, args.lineage, 
	args.aligndir, args.config, args.command, args.trim, args.trimparams, 
	args.aligndir, args.outformat)
	Additional arguments: args.matrix, args.partitions, args.outdir, args.seqtype, args.prefix, args.bootstrap, args.threads
	'''
	try:
		if not os.path.isdir(OUTDIR):
			os.mkdir(OUTDIR)
	except:
		raise RuntimeError("Output directory can not be created. Check your path!")
	step1_dir, step2_dir, step3_dir, step4_dir = ["","","",""]
	for step in STEPS:
		#01_single-copy 02_alignments 02_Matrix 03_Tree
		if step == 1: #Finding single-copy BUSCOs...
			print("1. Looking for single-copy BUSCO groups in the lineage %s obtained from the ODB v%d database..." % (args.lineage,args.odb))
			step1_dir = os.path.join(OUTDIR, "01_single-copy")
			os.mkdir(step1_dir) #creating output dir OUTDIR/01_single-copy
            #first we find common single-copy BUSCOs among genomes
			genomes_names, common_busco_ids = step1.find_singlecopy(BUSCODIR, step1_dir, ODB, LINEAGE)
			#now we create multifasta files with common single-copy BUSCOs
			step1.create_busco_fasta(BUSCODIR, step1_dir, ODB, LINEAGE, genomes_names, common_busco_ids, common_busco_dir="common_busco_sequences")
		if step == 2: #Aligning BUSCOs multifasta files...
			print("2. Aligning common single-copy BUSCOs...")
			step2_dir = os.path.join(OUTDIR, "02_alignments")
			os.mkdir(step2_dir) #creating output dir OUTDIR/02_alignments
			if 1 in STEPS: #se hace en cadena
				FASTADIR = os.path.join(step1_dir, "common_busco_sequences")
			if COMMAND:
				print("MAFFT will be excecuted using a user-defined command.")
				step2.align_command_mafft(FASTADIR, step2_dir, COMMAND)
			else:
				print("MAFFT will be excecuted using a configuration file.")
				step2.align_config_mafft(FASTADIR, step2_dir, CONFIG)
			print("Alignments completed...")
			if TRIM: #trimming alingments
				print("Trimming alignments to remove poorly aligned regions...")
				step2.trim_alns(step2_dir, TRIMPARAMS)
				print("Poorly aligned regiones have been removed.")
			else:
				if TRIMPARAMS: #trimming alingments
					print("Trimming alignments to remove poorly aligned regions...")
					step2.trim_alns(step2_dir, TRIMPARAMS)
					print("Poorly aligned regiones have been removed.")
		if step == 3: #Generating the phylogenetic matrix
			print("3. Generating the phylogenetic matrix.")
			step3_dir = os.path.join(OUTDIR, "03_matrix")
			os.mkdir(step3_dir) #creating output dir OUTDIR/03_matrix
			if 2 in STEPS and TRIM:
				ALIGNDIR =  os.path.join(step2_dir, "trimAl")
			elif 2 in STEPS and TRIMPARAMS:
				ALIGNDIR =  os.path.join(step2_dir, "trimAl")
			step3.cat_alignments(ALIGNDIR, step3_dir, FORMAT)
		if step == 4: #Building the phylogenetic tree
			print("4. Creating the phylogenetic tree with IQTree.")
			step4_dir = os.path.join(OUTDIR, "04_phylogenetic_tree")
			os.mkdir(step4_dir) #creating output dir OUTDIR/04_phylogenetic_tree
			if 3 in STEPS:
				files_in_step3_dir = os.listdir(step3_dir) #['phylomatrix.phylip', 'busco_coords.partitions.tsv']
				for file in files_in_step3_dir:
					if '.phy' in file or '.nex' in file:
						MATRIX = os.path.join(step3_dir, file)
					else:
						PARTITIONS = os.path.join(step3_dir, "busco_coords.partitions.tsv")
			step4.model_partitions(MATRIX, PARTITIONS, step4_dir, SEQTYPE, PREFIX, BOOTSTRAP, THREADS)
#end
#%% Menu -> is executed when this script is called as main program
def usage():
	parser = argparse.ArgumentParser(
		description='''BUSCO2Tree.py helps to process the BUSCO outputs over several genomes/proteomes/transcriptomes to create a species tree.''',
		epilog="""End of the help""")
	#general arguments
	general_arguments = parser.add_argument_group('General arguments')
	general_arguments.add_argument('-s', '--steps', nargs='+', type=int, required=True, help='Steps to run. Only single or consecutive steps are allowed. Available options are: 1) find common single copy BUSCOs; 2) align common single copy BUSCOs; 3) create phylogenetic matrix by concatenating BUSCO alignmets; 4) Generate the phylogenetic tree.')
	general_arguments.add_argument('-o', '--outdir', type=str, required=False, default='B2T_output', help='Path to the directory where the output of each step will be saved. If it does not exists, it will be created. Default: B2T_output')
	general_arguments.add_argument('-t', '--threads', type=int, required=False, default=4, help='Number of threads to use. Default: 4')
	#step 1
	step1_arguments = parser.add_argument_group('Step1: Find common single-copy BUSCO groups among BUSCO outputs')
	step1_arguments.add_argument('-b', '--buscodir', type=str, required=False, help='Path to the directory where individual BUSCO outputs are located (you can select the BUSCO output when running in batch mode). Directories names will be taken as the genomes names.')
	step1_arguments.add_argument('-d', '--odb', type=int, required=False, default=10, help='OrthoDB version. This parameter is only use to complete the path of the single copy sequences directory. Default: 10')
	step1_arguments.add_argument('-l', '--lineage', type=str, required=False, default='eukaryota', help='OrthoDB lineage. This parameter is only use to complete the path of the single copy sequences directory. Default: eukaryota')
	#step 2
	step2_arguments = parser.add_argument_group('Step2: Align common BUSCO groups')
	step2_arguments.add_argument('-f', '--fastadir', type=str, required=False, help='Path to the directory containing the BUSCO groups in fasta format.')
	step2_arguments.add_argument('-cnf', '--config', metavar='<config file>', type=str, required=False, help='Config file for alignment setting. Users can find a config file template in the docs directory. If no file is provided, the alignments will be done using default parameters.')
	step2_arguments.add_argument('-cmd', '--command', metavar='<command>', type=str, required=False, help='MAFFT parameters to apply to each alingment. The parmeters must be defined in a command line style, between quote marks, and the avoiding the names of in/output files, e.g. "--unalignlevel 0.1 --leavegappyregion --ep 0.12 --globalpair --maxiterate 1000".')
	step2_arguments.add_argument('--trim', action='store_true', required=False, help='Trim alignments using trimAl in automated mode, which must be available in the path.')
	step2_arguments.add_argument('--trimparams', required=False, help='trimAl parameters to apply when removing poorly aligned regions. The parameters must be defined in quotes marks and avoiding in/output files, e.g. "-gt 0.3".')
	#step 3
	step3_arguments = parser.add_argument_group('Step3: Create phylogenetic matrix and partition files')
	step3_arguments.add_argument('-a', '--aligndir', type=str, required=False, help='Path to the directory containing the alignments in fasta format (step 3).')
	step3_arguments.add_argument('-fmt', '--format', type=str, required=False, choices=['nexus','phylip'], default='phylip', help='Format for the phylogenetic matrix generated by concatenating BUSCOs alignment files (step 3). Default: phylip')
	#step 4
	step4_arguments = parser.add_argument_group('Step4: Genetare phylogenetic tree')
	step4_arguments.add_argument('-m', '--matrix', type=str, required=False, help='Path to the phylogenetic matrix file.')
	step4_arguments.add_argument('-p', '--partitions', type=str, required=False, help='Path to the partitions file in nexus format.')
	step4_arguments.add_argument('-P', '--prefix', type=str, required=False, default="iqtree", help='Prefix to name the output dataset and results. Default: iqtree')
	step4_arguments.add_argument('-B', '--bootstrap', type=int, required=False, default="1000", help='Number of bootstrap replicate to run. Default: 1000')
	step4_arguments.add_argument('-st', '--seqtype', type=str, required=False, choices=["DNA","AA"], default="AA", help='Type of sequence. Default: AA')
	
	return parser.parse_args()
#end
#%% Main program
if __name__ == '__main__':
    start = time() #time 0
    args = usage()
    print("0.1. Validating the selected steps...") #first we validate that steps introduced are ok
    validate_steps(args.steps)
    print("0.2. Validating parameters...") #now we validate parameters
    validate_params(args)
    print("\t\tSteps " + str(args.steps) + " selected.") #now we see which steps were requested
    print("0.3. Running BUSCO2Tree pipeline...")
    BUSCO2Tree(args.steps, args.buscodir, args.outdir, args.odb, args.lineage, args.fastadir, args.config, args.command, args.trim, args.trimparams, args.aligndir, args.format, args.matrix, args.partitions, args.seqtype, args.prefix, args.bootstrap, args.threads)
    print("BUSCO2Tree has finished.")
    print(f'Time taken to run: {time() - start} seconds.')
#%% End