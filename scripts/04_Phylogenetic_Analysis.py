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
    try:
        if not os.path.isdir(OUTDIR):
            os.mkdir(OUTDIR)
            os.chdir(OUTDIR)
    except:
        raise ValueError("Output directory %s can not be created." % OUTDIR)
		
	#iqtree -s supermatrix.aln.faa.phy -p partitions-scheme.txt -m MFP #basic example
    #--seqtype AA --prefix Drosophila -B 1000 --mem 50 -T 8
    IQTree_MFP = "iqtree -s " + MATRIXFILE + " -p " + PARTITIONFILE + " -m MFP --seqtype " + SEQTYPE + " --prefix " + PREFIX + " -B " + str(BOOTSTRAP) + " -T " + str(THREADS)
    subprocess.call([IQTree_MFP], shell=True)	
#end
def run_IQTree_manual(MATRIXFILE, PARTITIONFILE, OUTDIR, COMMAND):
	cwd = os.getcwd()
	try:
		if not os.path.isdir(OUTDIR):
			os.mkdir(OUTDIR)
			os.chdir(OUTDIR)
	except:
		raise ValueError("Output directory %s can not be created." % OUTDIR)
	IQTree_command = "iqtree -s " + cwd + "/" + MATRIXFILE + " -p " + cwd + "/" + PARTITIONFILE + COMMAND
	run_iqtree = subprocess.call([IQTree_command], shell=True, stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
#end
#%% Menu -> Usage
def usage():
	parser = argparse.ArgumentParser(
		description='''03_create_matrix.py creates a phylogenetic matrix with coordinates by concatenating alignment files in fasta format.''',
		epilog="""End of the help""")

	parser.add_argument('-m', '--matrix', type=str, required=True, help='Full path to the phylogenetic matrix file.')
	parser.add_argument('-p', '--partitions', type=str, required=True, help='Full path to the partitions file in nexus format.')
	parser.add_argument('-s', '--seqtype', type=str, required=False, choices=["DNA","AA"], default="AA", help='Sequence type.')
	parser.add_argument('-o', '--outdir', type=str, required=True, help='Full/relative path where the phylogenetic analysis results will be saved. If it does not exists, it will be created.')
	parser.add_argument('-P', '--prefix', type=str, required=False, default="iqtree", help='Prefix for the output dataset and results.')
	parser.add_argument('-b', '--bootstrap', type=int, required=False, default="1000", help='Number of bootstrap replicates.')
	parser.add_argument('-t', '--threads', type=int, required=False, default="16", help='Number of threads to use in IQ-Tree.')
	#parser.add_argument('--manual', action='store_true', required=False, help='IQTree in manual mode.')
	#parser.add_argument('-c', '--command', required=False, help='IQTree parameters to use. They must be passed in quotes marks and avoiding input and output files/directories, e.g. "-m MFP --seqtype AA --prefix Drosophila -B 1000 --mem 50 -T 16".')

	return parser.parse_args()

#%% Main program
if __name__ == '__main__':
	args = usage()
	start = time() #time 0
	#if args.manual:
	#	print("Step1: Starting the phylogenetic analysis with IQ-Tree in manual mode")
	#	run_IQTree_manual(args.matrix, args.partitions, args.outdir, args.command)
	#else:
	print("Step1: Starting the phylogenetic analysis with IQ-Tree in automated mode") #using arguments: %s\t%s\t%s\t%n\t%n")
	model_partitions(args.matrix, args.partitions, args.outdir, args.seqtype, args.prefix, args.bootstrap, args.threads)
	print("Finishing IQ-Tree...")
	print(f'Time taken to run: {time() - start} seconds.')
#end