#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 10:08:55 2022

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)

This script is used to create a phylogenetic matrix by concatenating alignments and a loci coordinates file.
"""
#%% Imports
import argparse
import os
from Bio import AlignIO
from time import time	

#%% Functions definition

# Function to concatenate alignments and create a phylogenetic matrix and coordinates of each locus
def cat_alignments(ALIGNDIR, OUTDIR, FORMAT):
	"""
	This function receives the paths of a directory containing alignment files, an directory to save results and the matrix format.
	It returns the phylogenetic matrix and the coordinates of each locus (BUSCO group).
	"""
	# Check if the directory containing alignments exists
	try:
		aln_files = os.listdir(ALIGNDIR)
	except:
		raise ValueError("The directory %s containing alignments can not be found." % ALIGNDIR)
	# Check if the output directory exists, if not, create it
	try:
		if not os.path.isdir(OUTDIR):
			os.mkdir(OUTDIR)
	except:
		raise ValueError("The output directory %s can not be created." % OUTDIR)
	# Initialize the phylogenetic matrix and coordinates dictionary
	matrix = {} # Phylogenetic matrix. Taxa are keys and sequences are the values.
	coordinates = {} # Dictionary to save the start and end positions of each locus
	matrixlen = 0 # Initialize the length of the matrix
	# Loop through each alignment file
	for n, file in enumerate(aln_files):
		infile = os.path.join(ALIGNDIR, file)
		# Skip if the file is a directory
		if os.path.isdir(infile):
			continue
		# Read the alignment file
		alignment = AlignIO.read(infile, "fasta")
		# Loop through each taxon in the alignment
		for i, taxon in enumerate(alignment):
			# For the first taxon, save the start and end positions of the locus
			if i == 0:
				coordinates[taxon.name.split('|')[1]] = [matrixlen+1, matrixlen+len(str(taxon.seq))]
				matrixlen += len(str(taxon.seq))
			# For the first alignment, define the sequence
			if n == 0:
				matrix[taxon.name.split('|')[0]] = str(taxon.seq)
			# For subsequent alignments, concatenate the sequences
			else:
				matrix[taxon.name.split('|')[0]] += str(taxon.seq)
	# Create the matrix and partition
	do_matrix(matrix, OUTDIR, FORMAT)
	do_partitions(OUTDIR, coordinates)

# Function to create the partition
def do_partitions(OUTDIR, COORDINATES):
	"""
	This function receives the output directory and coordinates of each locus.
	It creates a file with the coordinates of each locus.
	"""
	# Define the output file path
	coords_out = os.path.join(OUTDIR, "busco_coords.partitions.tsv")
	# Write the coordinates to the file
	with open(coords_out, 'wt') as COORDS:
		COORDS.write("#nexus\nbegin sets;\n")
		for buscoid, x_y in COORDINATES.items():
			COORDS.write("\tcharset " + buscoid + '\t=\t' + str(x_y[0]) + '-' + str(x_y[1]) + ';\n')
		COORDS.write("end;\n")	

# Function to create the phylogenetic matrix
def do_matrix(matrix, OUTDIR, out_format):
	"""
	This function receives a matrix in a dictionary, the output directory and matrix format.
	It generates the phylogenetic matrix file.
	"""
	# Get the number of taxa and the length of the sequences
	ntaxa = len(matrix)
	name_seqs = list(matrix.keys())
	nchars = len(matrix[name_seqs[0]])
	out_matrix_path = ''
	# Set the output format for the matrix
	if out_format == 'nexus':
		out_matrix_path = os.path.join(OUTDIR, "matrix.nexus")
		header = "#nexus\nbegin data;\n" + ' '*8 + "DIMENSIONS  NTAX=" + str(ntaxa) + " NCHAR=" + str(nchars) + ";\n" + ' '*8 + "FORMAT DATATYPE = DNA GAP = - MISSING = ?;\n" + ' '*8 + "MATRIX\n"
		tail = ";\nend;\n"
	elif out_format == 'phylip':
		out_matrix_path = os.path.join(OUTDIR, "matrix.phylip")
		header = ' '*8 + str(ntaxa) + ' ' + str(nchars) + '\n'
	else:
		raise ValueError("The output format %s is not recognized. Available options are nexus or phylip." % out_format)
	# Start writing the matrix
	with open(out_matrix_path,'wt') as OUT:
		OUT.write(header)
		for taxon, sequence in matrix.items():
			taxon_line = taxon + ' ' + str(sequence) + '\n'
			OUT.write(taxon_line)
		if out_format == 'nexus':
			OUT.write(tail)

# Function to parse command line arguments
def usage():
	"""
	This function parses command line arguments.
	It returns the parsed arguments.
	"""
	parser = argparse.ArgumentParser(
		description='''create_matrix.py creates a phylogenetic matrix by concatenating alignments and a loci coordinates file.''',
		epilog="""End of the help""")
	parser.add_argument('-a', '--aligndir', type=str, required=True, help='Path to the directory containing the alignment files in fasta format.')
	parser.add_argument('-o', '--outdir', type=str, required=False, default='03_matrix', help='Path to save the phylogenetic matrix and the loci coordinate file. If it does not exists, it will be created. Default: 03_matrix')
	parser.add_argument('-f', '--outformat', type=str, required=False, choices=['nexus','phylip'], default='phylip', help='Format for the phylogenetic matrix.')
	return parser.parse_args()

# Main program
if __name__ == '__main__':
	# Parse command line arguments
	args = usage()
	# Start the timer
	start = time()
	print("Starting the concatenation of alignments...")
	# Call the function to concatenate alignments and create the matrix
	cat_alignments(args.aligndir, args.outdir, args.outformat)
	print("Phylogenetic matrix have been created.")
	# Print the time taken to run the script
	print(f'Time taken to run: {time() - start} seconds.')
