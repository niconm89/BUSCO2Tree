#!/usr/bin/python

import argparse
import sys
import os
from collections import Counter
#from Bio import SeqIO
from Bio import AlignIO

def count_nucs(dictionary_of_records,position):
	list_nucs = []
	for ids,seq in dictionary_of_records.items():
		list_nucs.append(seq[position])
	counter = Counter(list_nucs)
	mas_comun = counter.most_common()
	return mas_comun[0][0]

##Menu and help
def getOptions(args=sys.argv[1:]):
        parser = argparse.ArgumentParser(description="Script command help.")
        parser.add_argument("-d", "--directory", help="Path to the folder whith the genes in fasta format.")
        parser.add_argument("-o", "--output", help="Your destination output folder.")
        parser.add_argument("-v", "--verbose",dest='verbose',action='store_true', help="Verbose mode.")
        options = parser.parse_args(args)
        return options

options = getOptions(sys.argv[1:])

lista_IUPAC = ['r','y','s','w','k','m','b','d','h','v','n']

#dict_IUPAC
'''
A	Adenine
C	Cytosine
G	Guanine
T (or U)	Thymine (or Uracil)
R	A or G
Y	C or T
S	G or C
W	A or T
K	G or T
M	A or C
B	C or G or T
D	A or G or T
H	A or C or T
V	A or C or G
N	any base
. or -	gap
'''

for codon in os.listdir(options.directory):
	codon = codon.rstrip()
	if codon.endswith('codon'):
		gene_file = os.path.join(options.directory, codon)
		print(gene_file)
		if os.stat(gene_file).st_size > 0:
			alignment = AlignIO.read(open(gene_file), 'fasta')

			aln = {}
			for record in alignment:
				aln[record.id] = record.seq


			for idseq,seq in aln.items():
				new_seq = ''
				for i,letter in enumerate(seq):
					if letter in lista_IUPAC:
						mas_comun = count_nucs(aln,i)
						new_seq = new_seq + str(mas_comun)
					else:
						new_seq = new_seq + str(seq[i])
				aln[idseq] = new_seq

			codonout = os.path.join(options.output,codon)
			with open(codonout, 'wt') as OUT:
				for idseq,seq in aln.items():
					OUT.write('>'+str(idseq)+'\n'+str(seq)+'\n')
