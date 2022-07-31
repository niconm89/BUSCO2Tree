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
def model_partitions(MATRIXFILE, PARTITIONFILE, OUTDIR, PREFIX, BOOTSTRAP, THREADS):
	#iqtree -s supermatrix.aln.faa.phy -p partitions-scheme.txt -m MFP 
	#--seqtype AA --prefix Drosophila -B 1000 --mem 50 -T 8
	try:
		if not os.path.isdir(OUTDIR):
			os.mkdir(OUTDIR)
	except:
		raise ValueError("Output directory %s can not be created." % OUTDIR)
		
	print("Starting the phylogenetic analysis with IQ-Tree")#using arguments: %s\t%s\t%s\t%n\t%n")
	IQTree_MFP = "iqtree -s " + MATRIXFILE + " -p " + PARTITIONFILE + " -m MFP --seqtype AA --prefix " + PREFIX + " -B " + str(BOOTSTRAP) + " -T " + str(THREADS)
	print(IQTree_MFP)
	#process = subprocess.Popen(IQTree_MFP.split(), stdout=subprocess.PIPE)
	#output, error = process.communicate()
	subprocess.call([IQTree_MFP], shell=True)	
	print("Finishing IQ-Tree")
#end

#%% Menu -> Usage
def usage():
	parser = argparse.ArgumentParser(
		description='''03_create_matrix.py creates a phylogenetic matrix with coordinates by concatenating alignment files in fasta format.''',
		epilog="""End of the help""")

	parser.add_argument('-m', '--matrix', type=str, required=True, help='Path to the phylogenetic matrix file to run.')
	parser.add_argument('-p', '--partitions', type=str, required=True, help='Partitions file in nexus format.')
	parser.add_argument('-o', '--outdir', type=str, required=True, help='Path to save the phylogenetic analysis results. If it does not exists, it will be created.')
	parser.add_argument('-P', '--prefix', type=str, required=False, default="iqtree", help='Prefix to name the output dataset and results.')
	parser.add_argument('-b', '--bootstrap', type=int, required=False, default="1000", help='Number of bootstrap replicate to run.')
	parser.add_argument('-t', '--threads', type=int, required=False, default="16", help='Number of threads to use in IQ-Tree.')
	parser.add_argument('-s', '--seqtype', type=str, required=False, choices=["NUC","AA"], default="AA", help='Sequence type')
	#seqtype not implemented

	return parser.parse_args()

#%% Main program
if __name__ == '__main__':
	args = usage()
	start = time()
	model_partitions(args.matrix, args.partitions, args.outdir, args.prefix, args.bootstrap, args.threads)
	print("Phylogenetic matrix created...")
	print(f'Time taken to run: {time() - start} seconds.')
#
'''
iqtree -s ../05_Matrix_Partitions/phylomatrix.phylip -p ../05_Matrix_Partitions/busco_coords.partitions.tsv -m MFP \
--seqtype AA --prefix Pyraloidea -B 1000 -T 8
'''