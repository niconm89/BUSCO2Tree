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
#%% Functionts definition
def model_partitions(MATRIXFILE, PARTITIONFILE, OUTDIR, SEQTYPE, PREFIX, BOOTSTRAP, THREADS):
	#iqtree -s supermatrix.aln.faa.phy -p partitions-scheme.txt -m MFP --seqtype AA --prefix Drosophila -B 1000 --mem 50 -T 8
	cwd = os.getcwd()
	try:
		if not os.path.isdir(OUTDIR):
			os.mkdir(OUTDIR)
		os.chdir(OUTDIR)
		print(os.getcwd())
	except:
		raise ValueError("Output directory %s can not be created." % OUTDIR)
	IQTree_MFP = ""
	matrixfile = MATRIXFILE
	if not os.path.isabs(MATRIXFILE):
		matrixfile = os.path.join(cwd, MATRIXFILE)
	partitionfile = PARTITIONFILE
	if not os.path.isabs(PARTITIONFILE):
		partitionfile = os.path.join(cwd, PARTITIONFILE)
	IQTree_MFP = "iqtree -s " + matrixfile + " -p " + partitionfile + " -m MFP --seqtype " + SEQTYPE + " --prefix " + PREFIX + " -B " + str(BOOTSTRAP) + " -T " + str(THREADS)
	'''
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
	'''
	run_iqtree = subprocess.call([IQTree_MFP], shell=True, stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
#end
#%% Menu -> Usage
def usage():
	parser = argparse.ArgumentParser(
		description='''phylogenetic_analysis.py creates a species tree using a phylogenetic matrix and a partition file (optional).''',
		epilog="""End of the help""")
	parser.add_argument('-m', '--matrix', type=str, required=True, help='Path to the matrix file.')
	parser.add_argument('-p', '--partitions', type=str, required=True, help='Path to the partitions file in nexus format.')
	parser.add_argument('-s', '--seqtype', type=str, required=False, choices=["DNA","AA"], default="AA", help='Sequence type. Default: AA')
	parser.add_argument('-o', '--outdir', type=str, required=True, help='Path to save results. If the output directory does not exists, it will be created.')
	parser.add_argument('-P', '--prefix', type=str, required=False, default="iqtree", help='Prefix for the output dataset and results. Default: iqtree')
	parser.add_argument('-b', '--bootstrap', type=int, required=False, default="1000", help='Number of bootstrap replicates. Default: 1000')
	parser.add_argument('-t', '--threads', type=int, required=False, default="8", help='Number of threads to run in IQ-Tree. Default: 8')
	return parser.parse_args()
#%% Main program
if __name__ == '__main__':
	args = usage()
	start = time() #time 0
	print("Starting the phylogenetic analysis with IQ-Tree...")
	model_partitions(args.matrix, args.partitions, args.outdir, args.seqtype, args.prefix, args.bootstrap, args.threads)
	print("IQ-Tree has finished.")
	print(f'Time taken to run: {time() - start} seconds.')
#end