#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 10:08:55 2022

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
"""
#%% Imports
import argparse
import os
from Bio import AlignIO
from time import time	
#%% Functions definition
def cat_alignments(ALIGNDIR, OUTDIR, FORMAT):
	'It receives the paths of a directory containing alignment files, an directory to save results and the matrix format. It returns the phylogenetic matrix and the coordinates of each locus (BUSCO group).'
	try:
		aln_files = os.listdir(ALIGNDIR)
	except:
		raise ValueError("The directory %s containing alignments can not be found." % ALIGNDIR)
	try:
		if not os.path.isdir(OUTDIR):
			os.mkdir(OUTDIR)
	except:
		raise ValueError("The output directory %s can not be created." % OUTDIR)
	matrix = {} #phylogenetic matrix. Taxa are keys and sequences are the values.
    #Sequences will made by concatenating the alignments for each taxon.
	coordinates = {} #it saves the star and end positions of each locus
	matrixlen = 0
	for n, file in enumerate(aln_files): #reading alignments
		infile = os.path.join(ALIGNDIR, file)
		alignment = AlignIO.read(infile, "fasta")
		for i, taxon in enumerate(alignment):
			if i == 0:
				coordinates[taxon.name.split('|')[1]] = [matrixlen+1, matrixlen+len(str(taxon.seq))]
				matrixlen += len(str(taxon.seq))
			if n == 0: #defining sequence in the first alignment
				matrix[taxon.name.split('|')[0]] = str(taxon.seq)
			else: #concatenating sequences after the first alignment
				matrix[taxon.name.split('|')[0]] += str(taxon.seq)
	do_matrix(matrix, OUTDIR, FORMAT) #create matrix
	do_partitions(OUTDIR, coordinates) #create partition
#end
def do_partitions(OUTDIR, COORDINATES):
	coords_out = os.path.join(OUTDIR, "busco_coords.partitions.tsv") #saving coordinates of each locus
	with open(coords_out, 'wt') as COORDS:
		COORDS.write("#nexus\nbegin sets;\n")
		for buscoid, x_y in COORDINATES.items():
			COORDS.write("\tcharset " + buscoid + '\t=\t' + str(x_y[0]) + '-' + str(x_y[1]) + ';\n')
		COORDS.write("end;\n")	
#end
def do_matrix(matrix, OUTDIR, out_format):
	'Receive a matrix in a dictionary, the output directory and matrix format to generate the phylogenetic matrix file.'
	ntaxa = len(matrix) #number of taxa
	name_seqs = list(matrix.keys()) #name of taxa
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
		raise ValueError("The output format %s is not recognized. Available options are nexus or phylip." % out_format)
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
		description='''create_matrix.py creates a phylogenetic matrix by concatenating alignments and a loci coordinates file.''',
		epilog="""End of the help""")
	parser.add_argument('-a', '--aligndir', type=str, required=True, help='Path to the directory containing the alignment files in fasta format.')
	parser.add_argument('-o', '--outdir', type=str, required=True, help='Path to save the phylogenetic matrix and the loci coordinate file. If it does not exists, it will be created.')
	parser.add_argument('-f', '--outformat', type=str, required=False, choices=['nexus','phylip'], default='phylip', help='Format for the phylogenetic matrix.')
	return parser.parse_args()
#end
#%% Main program
if __name__ == '__main__':
	args = usage() #reading arguments
	start = time() #time 0
	print("Starting the concatenation of alignments...")
	cat_alignments(args.aligndir, args.outdir, args.outformat)
	print("Phylogenetic matrix have been created.")
	print(f'Time taken to run: {time() - start} seconds.')
#