#!/usr/bin/python

import argparse
import sys
import os
from Bio import SeqIO

##Menu and help
def getOptions(args=sys.argv[1:]):
	parser = argparse.ArgumentParser(description="Script command help.")
	parser.add_argument("-d", "--directory", help="Path to the folder whith the genes in fasta format.")
	parser.add_argument("-o", "--output", help="Your destination output folder.")
	#parser.add_argument("-n", "--number", type=int, help="A number.")
	parser.add_argument("-v", "--verbose",dest='verbose',action='store_true', help="Verbose mode.")
	options = parser.parse_args(args)
	return options

options = getOptions(sys.argv[1:])
if options.verbose:
	print("Verbose mode on")
else:
	print("Verbose mode off")

print(options.directory)
print(options.output)
#print(options.number)
##End

#Beginning of the scripts
species = ['Dald','Dari','Dato','Dbrb','Dbuz','Dhyd','DkoeA','DkoeB','Dmel','Dmoj','Dnav','Drep','Dvir']

final_matrix_file = os.path.join(str(options.output),"supermatrix.phylip")
phylip = open(final_matrix_file,"w")

coordinates = {}

seqlen = 0

count = 0
for name in species:
	out = open(str(options.output)+"/"+name+"_BUSCO_groups.fasta","w")
	out.write("\n>"+name+"\t")
	suma = 1
	sequence = ''
	for gen in os.listdir(options.directory):
		gen_file = os.path.join(options.directory,gen.rstrip())
		BUSCO = SeqIO.parse(open(gen_file),'fasta')
		coordinates[gen.rstrip()] = []
		for fasta in BUSCO:
			gen_sequence = str(fasta.seq)
			gen_id = str(fasta.id)
			if gen_id == name:
				seqlen = len(gen_sequence)
				inicio_codon = suma
				suma += seqlen
				fin_codon = suma -1
				sequence += gen_sequence
				
				coordinates.update({ gen.rstrip() : [inicio_codon,fin_codon] })
				
				
	out.write(sequence)
	if count == 0:
		phylip.write(str(len(species))+" "+str(len(sequence))+"\n")
		count = count + 1
	phylip.write(name+"\n"+sequence+"\n")
	out.close()
	
	
	
with open(os.path.join(str(options.output),'partitions.tsv'), 'wt') as coords:
	#coords.write("#nexus\nbegin sets;\n")
	for k,v in coordinates.items():
		coords.write(k + "\t" + "=" + "\t" + str(v[0]) + "-" + str(v[1]) + ";\n")

