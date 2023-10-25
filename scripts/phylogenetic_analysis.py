#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 00:37:16 2022

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
"""
#%% Imports
import argparse
import os
import subprocess
from time import time	

# Importing necessary libraries. argparse for command line arguments, os for operating system related tasks,
# subprocess for running commands in a subshell, time for timing operations.

#%% Functionts definition
def model_partitions(MATRIXFILE, PARTITIONFILE, OUTDIR, SEQTYPE, PREFIX, BOOTSTRAP, THREADS):
	#iqtree -s supermatrix.aln.faa.phy -p partitions-scheme.txt -m MFP --seqtype AA --prefix Drosophila -B 1000 --mem 50 -T 8
	# This function runs the IQ-TREE software with the given parameters.
	# It takes as input the matrix file, partition file, output directory, sequence type, prefix for the output files,
	# number of bootstrap replicates and number of threads to use.
	cwd = os.getcwd()
	try:
		if not os.path.isdir(OUTDIR):
			os.mkdir(OUTDIR)
		os.chdir(OUTDIR)
	except:
		raise ValueError("Output directory %s can not be created." % OUTDIR)
	# The function first checks if the output directory exists, if not it creates it.
	# If the directory cannot be created, it raises an error.
	IQTree_MFP = ""
	matrixfile = MATRIXFILE
	if not os.path.isabs(MATRIXFILE):
		matrixfile = os.path.join(cwd, MATRIXFILE)
	partitionfile = PARTITIONFILE
	if not os.path.isabs(PARTITIONFILE):
		partitionfile = os.path.join(cwd, PARTITIONFILE)
	# The function then checks if the matrix and partition files are absolute paths, if not it makes them absolute.
	if BOOTSTRAP < 1000:
		raise ValueError("Number of bootstrap replicates must be an integer >=1000.")
	# It checks if the number of bootstrap replicates is at least 1000, if not it raises an error.
	IQTree_MFP = "iqtree -s " + matrixfile + " -p " + partitionfile + " -m MFP --seqtype " + SEQTYPE + " --prefix " + PREFIX + " -B " + str(BOOTSTRAP) + " -T " + str(THREADS)
	#run_iqtree = subprocess.call([IQTree_MFP], shell=True, stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
	subprocess.call([IQTree_MFP], shell=True, stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
	# The function then constructs the IQ-TREE command and runs it in a subshell.
#end
def gene_trees(MATRIXFILE, PREFIX, THREADS):
	IQTree_GENE = "iqtree -s " + MATRIXFILE + " -S " + PREFIX + ".best_scheme.nex" + " --prefix loci" + " -T " + str(THREADS)
	#run_iqtree = subprocess.call([IQTree_MFP], shell=True, stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
	subprocess.call([IQTree_GENE], shell=True, stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
	# The function then constructs the IQ-TREE command and runs it in a subshell.
#end
def concordance_factors(MATRIXFILE, PREFIX, THREADS):
	cwd = os.getcwd()
	#os.chdir(OUTDIR)
	# The function first checks if the output directory exists, if not it creates it.
	# If the directory cannot be created, it raises an error.
	IQTree_gCF = "iqtree -t " + PREFIX + ".treefile" + " --gcf " + "loci.treefile" + " --prefix " + PREFIX + ".treefile.gCF" + " -T " + str(THREADS)
	#iqtree2 -t $prefix.treefile --gcf loci.treefile --prefix $prefix.treefile.gCF
	IQTree_sCF = "iqtree -te " + PREFIX + ".treefile.gCF.cf.tree" + " -s " + MATRIXFILE + " -p " + PREFIX + ".best_scheme.nex" + " -blfix -scf 100 " + " --prefix " + PREFIX + ".treefile.sCF" + " -T " + str(THREADS)
	#iqtree2 -te $prefix.treefile.gCF.cf.tree -s matrix.phylip -p $prefix.best_scheme.nex -blfix --scf 100 --prefix $prefix.treefile.sCF -T $cpu
	subprocess.call([IQTree_gCF], shell=True, stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
	subprocess.call([IQTree_sCF], shell=True, stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
	# The function then constructs the IQ-TREE command and runs it in a subshell.
#end

#%% Menu -> Usage
def usage():
	parser = argparse.ArgumentParser(
		description='''phylogenetic_analysis.py creates a species tree using a phylogenetic matrix and a partition file (optional).''',
		epilog="""End of the help""")
	parser.add_argument('-m', '--matrix', type=str, required=True, help='Path to the matrix file.')
	parser.add_argument('-p', '--partitions', type=str, required=True, help='Path to the partitions file in nexus format.')
	parser.add_argument('-s', '--seqtype', type=str, required=False, choices=["DNA","AA"], default="AA", help='Sequence type. Default: AA')
	parser.add_argument('-o', '--outdir', type=str, required=False, default='04_phylogenetic_tree', help='Path to save results. If the output directory does not exists, it will be created. Default: 04_phylogenetic_tree')
	parser.add_argument('-P', '--prefix', type=str, required=False, default="iqtree", help='Prefix for the output dataset and results. Default: iqtree')
	parser.add_argument('-b', '--bootstrap', type=int, required=False, default="1000", help='Number of bootstrap replicates. Default: 1000')
	parser.add_argument('-t', '--threads', type=int, required=False, default="8", help='Number of threads to run in IQ-Tree. Default: 8')
	parser.add_argument('-gt', '--genetrees', action='store_true', required=False, help='Generate loci (gene) trees using iqtree')
	parser.add_argument('-cf', '--concordance', action='store_true', required=False, help='Calculate concordance factors using iqtree. The gene trees estimation (-gt or --genetress) must also be set.')
	return parser.parse_args()
	
	# This function sets up the command line arguments for the script using argparse.
	# It returns the parsed arguments.

#%% Main program
if __name__ == '__main__':
	args = usage()
	start = time() #time 0
	print("Starting the phylogenetic analysis with IQ-Tree...")
	cwd = os.getcwd()
	model_partitions(args.matrix, args.partitions, args.outdir, args.seqtype, args.prefix, args.bootstrap, args.threads)
	if args.genetrees:
		if not os.path.isabs(args.matrix):
			matrixfile = os.path.join(cwd, args.matrix)
		else:
			matrixfile = args.matrix
		if args.concordance:
			print("Estimating gene trees using IQ-Tree...")
			gene_trees(matrixfile, args.prefix, args.threads)
			print("Calculating concordance factors using IQ-Tree...")
			concordance_factors(matrixfile, args.prefix, args.threads)
		else:
			print("Estimating gene trees using IQ-Tree...")
			gene_trees(args.matrix, args.prefix, args.threads)
	print("IQ-Tree has finished.")
	print(f'Time taken to run: {time() - start} seconds.')
	
	# This is the main part of the script. It first parses the command line arguments, then starts the timer,
	# runs the model_partitions function with the parsed arguments, then prints the time taken to run.

#end
