#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 11:33:46 2022

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
"""
#%% Imports
import os
import argparse
#import check_paths

#%% Function definition
def check_paths(BUSCODIR, OUTDIR):
	if not os.path.isdir(BUSCODIR):
		print("Error: BUSCO output directory does not exist. Check your paths! " + BUSCODIR)
		exit(1)

	if '/' in OUTDIR:
		elements_path = OUTDIR.split('/')
		parentdir = '/'.join(elements_path[:-1])
		if not os.path.isdir(parentdir):
			print("Error: Output directory can not be found or created in the path you have passed. Parent directory does not exist! " + parentdir)
			exit(1)
#end
def read_fasta(infasta):
	if not check_paths.filepath(infasta):
		print("Error: Input fasta file does not exist. Check your paths! " + infasta)
		exit(1)

	seqs_dict = {}
	with open(infasta, "rt") as filein:
		fasta = filein.readlines()
		for line in fasta:
			line = line.rstrip() #let's discard the newline at the end (if any)
			#now we need to distinguish header from sequence
			if line[0] == '>': #identify the id of each sequence
				#words = line.split() #in case you want to avoid spaces in seq name
				#name = words[0][1:] #keep only the first word of the seq ID
				name = line[1:] #avoid ">"
				seqs_dict[name] = ''
			else:
				seqs_dict[name] = seqs_dict[name] + line

	seqs_list = list(seqs_dict.items()) #seqs_list is a list of tuples (id,seq)
	return seqs_list
#end
def create_BUSCO_FASTA(BUSCODIR, OUTDIR, ODB, LINEAGE, GENOMES, COMMONIDS):
	check_paths(BUSCODIR, OUTDIR) #check if BUSCO directory exists and output directory can be created
	outputfolder = os.path.join(OUTDIR, "busco_alignments")
	if not os.path.isdir(outputfolder):
		os.mkdir(outputfolder)
	for busco in COMMONIDS:
		busco_fasta_out = os.path.join(outputfolder, busco + ".faa")
		filenames = []
		with open(busco_fasta_out, 'wt') as FASTA:
			for n,name in enumerate(GENOMES):
				GENOMEDIR = os.path.join(BUSCODIR, name, "run_" + LINEAGE + "_" + ODB, "busco_sequences/single_copy_busco_sequences")
				busco_file = os.path.join(GENOMEDIR, busco + '.faa')
				filenames.append(busco_file)
			for fname in filenames:
				with open(fname) as infile:
					FASTA.write(infile.read()) #if fasta files are big, then they must be read line by line
#end		
#%% Menu -> is executed when the script is called independently
def main():
	parser = argparse.ArgumentParser(
		description='''find_singlecopy.py recover ID sequences from a BUSCO output directory.''',
		epilog="""End of the help""")

	parser.add_argument('-b', '--busco', type=str, required=True, help='Path to the file containing the list of common BUSCO single copy among genomes.')
	parser.add_argument('-o', '--outdir', type=str, required=True, help='Path to save output files. If it does not exists, it will be created.')
	parser.add_argument('-m', '--mode', type=str, required=False, choices=['genome', 'protein', 'transcript'], default='genome', help='BUSCO run mode. The genome mode is only currently available and set as default.')
	parser.add_argument('-d', '--odb', type=str, required=False, default='odb10', help='OrthoDB version. This parameter is only use to complete the path of the single copy sequences directory.')
	parser.add_argument('-l', '--lineage', type=str, required=True, help='OrthoDB lineage. This parameter is only use to complete the path of the single copy sequences directory.')

	args=parser.parse_args()
	#find_singlecopy(args.buscodir, args.outdir, args.odb, args.lineage)
	create_BUSCO_matrix(args.busco, args.outdir, args.odb, args.lineage)

#%% Main program
if __name__ == '__main__':
    main()

#%% End