#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 10:08:55 2022

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
"""
#%% Imports
import argparse
import os
#from Bio import SeqIO
from Bio import AlignIO
from time import time	

#%% Functions definition
def cat_alignments(ALIGNDIR, OUTDIR, FORMAT):
	try:
		aln_files = os.listdir(ALIGNDIR)
	except:
		raise ValueError("Directory %s containing alignments can not be found." % ALIGNDIR)
	try:
		if not os.path.isdir(OUTDIR):
			os.mkdir(OUTDIR)
	except:
		raise ValueError("Output directory %s can not be created." % ALIGNDIR)
	#start reading alignments and concatenating into a dictionart
	matrix = {} #phylogenetic matrix. Each taxon me will be a key and the value will be the sequence made by concatenating the corresponding sequence of all alignments.
	coordinates = {}
	matrixlen = 0
	for n, file in enumerate(aln_files):
		infile = os.path.join(ALIGNDIR, file)
		alignment = AlignIO.read(infile, "fasta")
		for i, taxon in enumerate(alignment):
			if i == 0:
				coordinates[taxon.name.split('|')[1]] = [matrixlen+1, matrixlen+len(str(taxon.seq))]
				matrixlen += len(str(taxon.seq))
			if n == 0:
				matrix[taxon.name.split('|')[0]] = str(taxon.seq)
			else:
				matrix[taxon.name.split('|')[0]] += str(taxon.seq)
		
	#once the matrix is ready, we need yo save it in the outformat selected.
	do_matrix(matrix, OUTDIR, FORMAT)
	#save coordinates
	coords_out = os.path.join(OUTDIR, "busco_coords.partitions.tsv")
	with open(coords_out, 'wt') as COORDS:
		COORDS.write("#nexus\nbegin sets;\n")
		for buscoid, x_y in coordinates.items():
			COORDS.write("\tcharset " + buscoid + '\t=\t' + str(x_y[0]) + '-' + str(x_y[1]) + ';\n')
		COORDS.write("end;\n")
#end
def do_matrix(matrix, OUTDIR, out_format):
	#header first
	ntaxa = len(matrix)
	name_seqs = list(matrix.keys())
	nchars = len(matrix[name_seqs[0]])
	out_matrix_path = ''
	#set output format for the matrix
	if out_format == 'nexus':
		out_matrix_path = os.path.join(OUTDIR, "phylomatrix.nexus")
		header = "#nexus\nbegin data;\n" + ' '*8 + "DIMENSIONS  NTAX=" + str(ntaxa) + " NCHAR=" + str(nchars) + ";\n" + ' '*8 + "FORMAT DATATYPE = DNA GAP = - MISSING = ?;\n" + ' '*8 + "MATRIX\n"
		tail = ";\nend;\n"
	elif out_format == 'phylip':
		out_matrix_path = os.path.join(OUTDIR, "phylomatrix.phylip")
		header = ' '*8 + str(ntaxa) + ' ' + str(nchars) + '\n'
	else:
		raise ValueError("The output matrix format is not recognized. Possible options are nexus and phylip." % out_format)
	#start writing the matrix
	with open(out_matrix_path,'wt') as OUT:
		OUT.write(header)
		for taxon, sequence in matrix.items():
			taxon_line = taxon + ' ' + str(sequence) + '\n'
			OUT.write(taxon_line)
		if out_format == 'nexus':
			OUT.write(tail)
#end
#%% Menu -> Usage
def usage():
	parser = argparse.ArgumentParser(
		description='''03_create_matrix.py creates a phylogenetic matrix with coordinates by concatenating alignment files in fasta format.''',
		epilog="""End of the help""")

	parser.add_argument('-a', '--aligndir', type=str, required=True, help='Path to the directory containing the alignments in fasta format.')
	parser.add_argument('-o', '--outdir', type=str, required=True, help='Path to save the phylogenetic matrix and the coordinates. If it does not exists, it will be created.')
	parser.add_argument('-f', '--outformat', type=str, required=False, choices=['nexus','phylip'], default='nexus', help='Format for the phylogenetic matrix.')

	return parser.parse_args()

#%% Main program
if __name__ == '__main__':
	args = usage()
	print("Starting the cancatenation of alignments files...")
	start = time()
	cat_alignments(args.aligndir, args.outdir, args.outformat)
	print("Phylogenetic matrix created...")
	print(f'Time taken to run: {time() - start} seconds.')
#