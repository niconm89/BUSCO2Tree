#!/usr/bin/env python3
"""
Created on Tue Nov 07 12:11:56 2020

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
"""
#%% Imports
import argparse

#%%Function definitions
def output_name(fullpath):
	'Parse file fullpath and extract just the name of the file.'
	last_part = fullpath
	if '/' in fullpath: #true if the fullpath to input file was introduced
		parts = fullpath.split('/') #split fullpath
		last_part = parts[-1] #get the last one (filename)
	inname = last_part.split('.')[0] #keep name only, exclude file extention
	return inname
	
def split_fasta(seqs_list,nparts,inname,outpath):
	'''
	This function splits list of seqs in nparts and output nparts files using the prefix.
	'''

	num_seqs = len(seqs_list) #number of sequences in the list
	if num_seqs >= nparts:
		seqs_per_part = int(num_seqs/nparts) #number of sequences to be inserted in each part
	else:
		print("Input file must have at least the same number of sequences than amount of parts introduced.")

	ini = 0 #start index
	for part in range(1,nparts+1):
		end = ini + seqs_per_part #slice to extract
		if part < nparts:
			subseqs = seqs_list[ini:end]
		else:
			subseqs = seqs_list[ini:] #in the last part keep the rest of the sequences (in case rest of divition was not zero)
		with open(outpath + '/' + inname + "_part" + str(part) + ".fasta",'wt') as OUT:
			for i in range(0,len(subseqs)):
					OUT.write(">"+subseqs[i][0]+"\n"+subseqs[i][1]+"\n")
		ini += seqs_per_part

def generate_seqsdict(infasta):
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

def main():
	parser = argparse.ArgumentParser(
		description='''split_fasta.py a fasta file in several fasta files.''',
		epilog="""End of the help""")

	parser.add_argument('-f', '--file', type=str, required=True, help='path to the file to split')
	parser.add_argument('-n', '--nparts', type=int, required=False, default=4, help='Numbers of parts the input will be divided. By defaults, n=4.')
	parser.add_argument('-o', '--output', type=str, required=True, help='path where outputs will be located')

	args=parser.parse_args()

	seqs_list = generate_seqsdict(args.file)
	inname = output_name(args.file)
	split_fasta(seqs_list,args.nparts,inname,args.output)

#%% Main program
if __name__ == '__main__':
    main()


