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
#from Bio import SeqIO
from time import time	

#%% Functionts definition
def model_partitions(MATRIXFILE, PARTITIONFILE, OUTDIR, SEQTYPE, PREFIX, BOOTSTRAP, THREADS):
	cwd = os.getcwd()
	try:
		if not os.path.isdir(OUTDIR):
			os.mkdir(OUTDIR)
		os.chdir(OUTDIR)
		print(os.getcwd())
	except:
		raise ValueError("Output directory %s can not be created." % OUTDIR)
	#iqtree -s supermatrix.aln.faa.phy -p partitions-scheme.txt -m MFP #basic example
    #--seqtype AA --prefix Drosophila -B 1000 --mem 50 -T 8
	IQTree_MFP = ""
	if os.path.isabs(MATRIXFILE):
		print("path absoluto")
		if os.path.isabs(PARTITIONFILE):
			IQTree_MFP = "iqtree -s " + MATRIXFILE + " -p " + PARTITIONFILE + " -m MFP --seqtype " + SEQTYPE + " --prefix " + PREFIX + " -B " + str(BOOTSTRAP) + " -T " + str(THREADS)
		else:
			partitionfile = os.path.join(cwd, PARTITIONFILE)
			IQTree_MFP = "iqtree -s " + MATRIXFILE + " -p " + partitionfile + " -m MFP --seqtype " + SEQTYPE + " --prefix " + PREFIX + " -B " + str(BOOTSTRAP) + " -T " + str(THREADS)
	else:
		print("path matrix relativo")
		matrixfile = os.path.join(cwd, MATRIXFILE)
		if os.path.isabs(PARTITIONFILE):
			IQTree_MFP = "iqtree -s " + matrixfile + " -p " + PARTITIONFILE + " -m MFP --seqtype " + SEQTYPE + " --prefix " + PREFIX + " -B " + str(BOOTSTRAP) + " -T " + str(THREADS)
		else:
			print("path partitions relativo")
			partitionfile = os.path.join(cwd, PARTITIONFILE)
			IQTree_MFP = "iqtree -s " + matrixfile + " -p " + partitionfile + " -m MFP --seqtype " + SEQTYPE + " --prefix " + PREFIX + " -B " + str(BOOTSTRAP) + " -T " + str(THREADS)
	run_iqtree = subprocess.call([IQTree_MFP], shell=True, stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
#end
'''
def run_IQTree_manual(MATRIXFILE, PARTITIONFILE, OUTDIR, COMMAND):
	cwd = os.getcwd()
	try:
		if not os.path.isdir(OUTDIR):
			os.mkdir(OUTDIR)
			os.chdir(OUTDIR)
	except:
		raise ValueError("Output directory %s can not be created." % OUTDIR)
	IQTree_command = ""
	if os.path.isabs(MATRIXFILE):
		if os.path.isabs(PARTITIONFILE):
			IQTree_command = "iqtree -s " + MATRIXFILE + " -p " + PARTITIONFILE + COMMAND
		else:
			partitionfile = os.path.join(cwd, PARTITIONFILE)
			IQTree_command = "iqtree -s " + MATRIXFILE + " -p " + partitionfile + COMMAND
	else:
		matrixfile = os.path.join(cwd, MATRIXFILE)
		if os.path.isabs(PARTITIONFILE):
			IQTree_command = "iqtree -s " + matrixfile + " -p " + PARTITIONFILE + COMMAND
		else:
			partitionfile = os.path.join(cwd, PARTITIONFILE)
			IQTree_command = "iqtree -s " + matrixfile + " -p " + partitionfile + COMMAND
	run_iqtree = subprocess.call([IQTree_command], shell=True, stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
#end
'''
#%% Menu -> Usage
def usage():
	parser = argparse.ArgumentParser(
		description='''03_create_matrix.py creates a phylogenetic matrix with coordinates by concatenating alignment files in fasta format.''',
		epilog="""End of the help""")

	parser.add_argument('-m', '--matrix', type=str, required=True, help='Path to the phylogenetic matrix file.')
	parser.add_argument('-p', '--partitions', type=str, required=True, help='Path to the partitions file in nexus format.')
	parser.add_argument('-s', '--seqtype', type=str, required=False, choices=["DNA","AA"], default="AA", help='Sequence type.')
	parser.add_argument('-o', '--outdir', type=str, required=True, help='Full/relative path where the phylogenetic analysis results will be saved. If it does not exists, it will be created.')
	parser.add_argument('-P', '--prefix', type=str, required=False, default="iqtree", help='Prefix for the output dataset and results.')
	parser.add_argument('-b', '--bootstrap', type=int, required=False, default="1000", help='Number of bootstrap replicates.')
	parser.add_argument('-t', '--threads', type=int, required=False, default="8", help='Number of threads to use in IQ-Tree.')
	parser.add_argument('--manual', type=int, action='store_true', required=False, help='a')

	return parser.parse_args()

#%% Main program
if __name__ == '__main__':
	args = usage()
	start = time() #time 0
	if args.manual:
		print("Starting the phylogenetic analysis with IQ-Tree in manual mode...")
		run_IQTree_manual(args.matrix, args.partitions, args.outdir, args.command)
	else:
		print("Starting the phylogenetic analysis with IQ-Tree in automated mode") #using arguments: %s\t%s\t%s\t%n\t%n")
		model_partitions(args.matrix, args.partitions, args.outdir, args.seqtype, args.prefix, args.bootstrap, args.threads)
		print("Finishing IQ-Tree...")
		print(f'Time taken to run: {time() - start} seconds.')
#end