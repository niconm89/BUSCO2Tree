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

phylip = open(str(options.output)+"/supermatrix.phylip","w")

count = 0
for name in species:
	out = open(str(options.output)+"/"+name+"_BUSCO_groups.fasta","w")
	out.write("\n>"+name+"\t")
	print(name)
	suma = 0
	#BUSCOs_groups = open("/home/nmoreyra/Documents/Software/busco_usecases/phylogenomics/workdir/divergence_times/lista_BUSCOs.txt","r")
	sequence = ''
	for gen in os.listdir(options.directory):
	        #print(gen)
		gen_file = options.directory+"/"+gen.rstrip()
		#print(gen_file)
		BUSCO = SeqIO.parse(open(gen_file),'fasta')
		for fasta in BUSCO:
			gen_sequence = str(fasta.seq)
			gen_id = str(fasta.id)
			if gen_id == name:
				suma += len(gen_sequence)
				sequence += gen_sequence
	out.write(sequence)
	if count == 0:
		phylip.write(str(len(species))+" "+str(len(sequence))+"\n")
		count = count + 1
	phylip.write(name+"\n"+sequence+"\n")
	print(suma)
	out.close()
