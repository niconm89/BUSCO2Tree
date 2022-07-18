#!/usr/bin/env python3
"""
Created on Tue Nov 07 12:11:56 2020

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
"""
#%% Imports
import argparse
import os

#%% Function definitions
def convert_aln2fasta(infile, outfile):
	with open(infile,'rt') as inputfile:
		dict_records = {}
		ALN = inputfile.readlines()
		for i, line in enumerate(ALN):
			if line[0] == '' or line[0] == ' ' or line[0] == '\n':
				pass
			else:
				if i > 0:
					line = line.rstrip()
					columns = " ".join(line.split())
					elements = columns.split(' ')
					if elements[0] not in dict_records and len(elements) > 1:
						dict_records[elements[0]] = elements[1]
					else:
						dict_records[elements[0]] += elements[1]
	with open(outfile, 'wt') as OUT:
		for record in sorted(dict_records.keys()):
			OUT.write(">"+record+"\n"+str(dict_records[record])+"\n")

def main():
	parser = argparse.ArgumentParser(
		description='''clustalw2fasta.py conversts a clustalw alignment file (.aln) into a fasta file.''',
		epilog="""End of the help""")

	parser.add_argument('-f', '--infile', type=str, required=True, help='Full path to the clustalw alignment file to convert')
	parser.add_argument('-o', '--outfile', type=str, required=True, help='Full path where the fasta file will be placed')

	args=parser.parse_args()

	convert_aln2fasta(args.infile,args.outfile)

#%% Main program
if __name__ == '__main__':
    main()
