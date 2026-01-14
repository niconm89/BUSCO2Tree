#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 10:08:55 2022

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)

This script is used to create a phylogenetic matrix by concatenating alignments and a loci coordinates file.
"""
#%% Imports
import argparse
from pathlib import Path
from Bio import AlignIO
from time import time	

#%% Functions definition

# Function to concatenate alignments and create a phylogenetic matrix and coordinates of each locus
def cat_alignments(ALIGNDIR, OUTDIR, FORMAT, DATATYPE="PROTEIN"):
	"""
	This function receives the paths of a directory containing alignment files, an directory to save results and the matrix format.
	It returns the phylogenetic matrix and the coordinates of each locus (BUSCO group).
	"""
	# Check if the directory containing alignments exists
	try:
		aln_files = sorted([f for f in Path(ALIGNDIR).iterdir() if f.is_file()])

	except:
		raise ValueError("The directory %s containing alignments can not be found." % ALIGNDIR)
	# Check if the output directory exists, if not, create it
	try:
		if not Path(OUTDIR).is_dir():
			Path(OUTDIR).mkdir(parents=True, exist_ok=True)
	except:
		raise ValueError("The output directory %s can not be created." % OUTDIR)
	# Initialize the phylogenetic matrix and coordinates dictionary
	matrix = {} # Phylogenetic matrix. Taxa are keys and sequences are the values.
	coordinates = {} # Dictionary to save the start and end positions of each locus
	matrixlen = 0 # Initialize the length of the matrix
	expected_taxa = None

	# Loop through each alignment file
	for n, file in enumerate(aln_files):
		infile = file
		# Skip if the file is a directory
		if Path(infile).is_dir():
			continue
		# Read the alignment file
		alignment = AlignIO.read(infile, "fasta")
		# Validación: todos los alignments deben contener exactamente los mismos taxa
		present_taxa = set()
		for taxon in alignment:
			parts = taxon.id.split('|')
			if len(parts) < 2:
				raise ValueError(
					f"Invalid FASTA header in {infile}. Expected 'genome|busco_id' but got: '{taxon.id}'"
				)
			present_taxa.add(parts[0])

		if n == 0:
			expected_taxa = present_taxa
		else:
			if present_taxa != expected_taxa:
				missing = sorted(expected_taxa - present_taxa)
				extra = sorted(present_taxa - expected_taxa)
				raise RuntimeError(
					f"Taxa mismatch in {infile}. Missing: {missing if missing else 'none'}; "
					f"Extra: {extra if extra else 'none'}"
				)
		# Loop through each taxon in the alignment
		for i, taxon in enumerate(alignment):
			parts = taxon.id.split('|')
			if len(parts) < 2:
				raise ValueError(
					f"Invalid FASTA header in {infile}. Expected 'genome|busco_id' but got: '{taxon.id}'"
				)

			genome_id = parts[0]
			busco_id = parts[1]
			seq = str(taxon.seq)

			# For the first taxon, save the start and end positions of the locus
			if i == 0:
				coordinates[busco_id] = [matrixlen + 1, matrixlen + len(seq)]
				matrixlen += len(seq)
			# For the first alignment, define the sequence
			if n == 0:
				matrix[genome_id] = seq
			# For subsequent alignments, concatenate the sequences
			else:
				if genome_id not in matrix:
					raise RuntimeError(
						f"Taxon '{genome_id}' is missing from the first alignment but present in {infile}. "
						"All alignments must contain the same taxa."
					)
				matrix[genome_id] += seq
	# Create the matrix and partition
	do_matrix(matrix, OUTDIR, FORMAT, DATATYPE)
	do_partitions(OUTDIR, coordinates)

# Function to create the partition
def do_partitions(OUTDIR, COORDINATES):
	"""
	This function receives the output directory and coordinates of each locus.
	It creates a file with the coordinates of each locus.
	"""
	# Define the output file path
	coords_out = Path(OUTDIR) / "partitions.nex"
	# Write the coordinates to the file
	with open(coords_out, 'wt') as COORDS:
		COORDS.write("#nexus\nbegin sets;\n")
		for buscoid, x_y in COORDINATES.items():
			COORDS.write("\tcharset " + buscoid + '\t=\t' + str(x_y[0]) + '-' + str(x_y[1]) + ';\n')
		COORDS.write("end;\n")	

# Function to create the phylogenetic matrix
def do_matrix(matrix, OUTDIR, out_format, datatype="PROTEIN"):
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
		out_matrix_path = Path(OUTDIR) / "matrix.nex"
		header = ("#nexus\nbegin data;\n"
		+ " " * 8
		+ f"DIMENSIONS  NTAX={ntaxa} NCHAR={nchars};\n"
		+ " " * 8
		+ f"FORMAT DATATYPE = {datatype} GAP = - MISSING = ?;\n"
		+ " " * 8
		+ "MATRIX\n"
		)
		tail = ";\nend;\n"
	elif out_format == 'phylip':
		out_matrix_path = Path(OUTDIR) / "matrix.phylip"
		header = f"{ntaxa} {nchars}\n"
	else:
		raise ValueError("The output format %s is not recognized. Available options are nexus or phylip." % out_format)
	# Start writing the matrix
	with open(out_matrix_path,'wt') as OUT:
		OUT.write(header)
		for taxon, sequence in sorted(matrix.items()):
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
	parser.add_argument("--datatype", choices=["DNA", "PROTEIN"], default="PROTEIN",
                    help="NEXUS DATATYPE. Use PROTEIN for BUSCO amino acid alignments.")
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
	cat_alignments(args.aligndir, args.outdir, args.outformat, args.datatype)
	print("Phylogenetic matrix have been created.")
	# Print the time taken to run the script
	print(f'Time taken to run: {time() - start} seconds.')
