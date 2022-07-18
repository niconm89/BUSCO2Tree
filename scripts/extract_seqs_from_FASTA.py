#!/home/nmoreyra/Soft/miniconda3/bin/python

"""
Created on Tue Mar 17 09:11:56 2021

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
"""

#%% Imports
import argparse
import os
from Bio import SeqIO
import random #comment this line if you are not using the random extraction mode and don't want to install the library
import gzip

#%% Function definitions
def check_paths(inputfile,outputfile):
	if not os.path.isfile(inputfile):
		print("Error: Input file does not exist! " + inputfile)
		exit(1)
	
	if '/' in outputfile:
		elements_path = outputfile.split('/')
		outputfolder = '/'.join(elements_path[:-1])
		if not os.path.isdir(outputfolder):
			print("Error: Output folder does not exist! " + outputfolder)
			exit(1)
#end
def check_mode(arguments):
	if arguments.mode == 'byIDs':
		if arguments.IDsfile:
			if not os.path.isfile(arguments.IDsfile):
				print("Error: IDs file does not exist! " + arguments.IDsfile)
				exit(1)
		else:
			print("Extraction mode was set to byIDs but no IDs file was input. Aborting...")
			exit(1)
	else:
		if not arguments.nseqs:
			print("Unless you are using byIDs extraction mode, you must define a the number of sequences you want to recover (-n, --nseqs). Aborting...")
			exit(1)
#end
def read_IDs(IDsfile):
	'receive a file with one ID per line and return a list'
	IDs_list = []
	with open(IDsfile,'rt') as inputIDs:
		lines = inputIDs.readlines()
		for line in lines:
			IDs_list.append(line.rstrip())
	return IDs_list
#
def check_gzipped(inputfile):
	'Return true if file is gzipped. Otherwise return false'
	# The first two bytes of a gzip file are: 1f8b
	gzipped = False
	GZIP_MAGIC_NUMBER = b'\x1f\x8b'
	with open(inputfile,'rb') as INPUT:
		twobytes = INPUT.read(2)
		if twobytes == GZIP_MAGIC_NUMBER:
			gzipped = True
	return gzipped
#end
def read_fasta(inputfile):
	sequences = {}
	if check_gzipped(inputfile):
		with gzip.open(inputfile, "rt") as FASTA:
			count = 0
			for record in SeqIO.parse(FASTA, "fasta"):
				count +=1
				seqname = str(record.id)
				seqprot = str(record.seq)
				if seqname in sequences:
					print("There have repeates sequence IDs in the inputfile. Aborting...")
					exit(0)
				sequences[seqname] = seqprot
	else:
		with open(inputfile, "rt") as FASTA:
			for record in SeqIO.parse(FASTA, "fasta"):
				seqname = str(record.id)
				seqprot = str(record.seq)
				if seqname in sequences:
					print("There have repeates sequence IDs in the inputfile. Aborting...")
					exit(0)
				sequences[seqname] = seqprot
	return sequences
#end
def extract_nseqs(inputfile, outputfile, nseqs, mode, IDsfile):
	'main function that manage the differente types of extraction.'
	sequences = read_fasta(inputfile)
	subseqs = {}
	
	if mode == 'head' or mode == 'tail':
		subseqs = extract_head_tail(sequences, nseqs, mode)
	else:
		if mode == 'random':
			subseqs = extract_random(sequences, nseqs)
		else:
			IDs_list = read_IDs(IDsfile)
			subseqs = extract_byIDs(sequences, IDs_list)
	save_subseqs(subseqs, outputfile)
#end
def save_subseqs(subseqs, outputfile):
	'get a dictionary with ids (key) and sequences (values) and save them in a outputfile.'
	with open(outputfile, 'wt') as OUT:
		for k,v in subseqs.items():
			OUT.write(">" + str(k) + "\n" + str(v) + "\n")
#end
def extract_head_tail(sequences, nseqs, mode):
	subsequences = {}
	if mode == 'head':
		for ID in list(sequences)[0:nseqs]:
			subsequences[ID] = sequences[ID]
	else:
		pos = len(sequences) - nseqs
		for ID in list(sequences)[pos:]:
			subsequences[ID] = sequences[ID]
	return subsequences
#end
def extract_random(sequences, nseqs):
	subsequences = {}
	lista_IDs = list(sequences.keys())
	for n in range(0, nseqs):
		IDseq = random.choice(lista_IDs)
		lista_IDs.remove(IDseq)
		subsequences[IDseq] = sequences[IDseq]
	return subsequences
#end
def extract_byIDs(sequences, IDs_list):
	subsequences = {}
	for ID in IDs_list:
		if ID in sequences:
			subsequences[ID] = sequences[ID]
		else:
			print("Sequence ID " + str(ID) + " not found. ID avoided...")
	return subsequences
#end

#%%
def main():
	parser = argparse.ArgumentParser(
		description='''extract_seqs_from_FASTA.py recover sequences from a FASTA file.''',
		epilog="""End of the help""")

	parser.add_argument('-f', '--infile', type=str, required=True, help='Full path to the FASTA file to convert. It can be gzipped.')
	parser.add_argument('-o', '--outfile', type=str, required=True, help='Full path where the new FASTA file will be placed. Output file will not be compressed.')
	parser.add_argument('-m', '--mode', type=str, required=True, choices=['head', 'tail', 'random', 'byIDs'], help='Extraction mode. head: first -n sequences; tail: last -n sequences; random: randomly selection of -n sequences; ID_file: file containing a sequence ID per line.')
	parser.add_argument('-n', '--nseqs', type=int, required=False, help='Numbers of sequences to recover from the input.')
	parser.add_argument('-d', '--IDsfile', type=str, required=False, help='Full path to the file with IDs to extract. Sequence descriptions must be excluded.')

	args=parser.parse_args()

	check_paths(args.infile, args.outfile) #check if input file and output folder exist
	check_mode(args)

	extract_nseqs(args.infile, args.outfile, args.nseqs, args.mode, args.IDsfile)

#%% Main program
if __name__ == '__main__':
    main()

