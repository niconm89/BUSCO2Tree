#!/home/nmoreyra/Software/miniconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 9 11:03:02 2022

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
"""
#%% Imports
import os
import argparse
from time import time

#%% Functions definition
def check_paths(BUSCODIR, OUTDIR):
	if not os.path.isdir(BUSCODIR):
		print("Error: BUSCO output directory %s does not exist. Check your paths! " % BUSCODIR)
		exit(1)

	if '/' in OUTDIR:
		elements_path = OUTDIR.split('/')
		parentdir = '/'.join(elements_path[:-1])
		if not os.path.isdir(parentdir):
			print("Error: Output directory %s can not be found or created in the path you have passed. Parent directory does not exist! " % parentdir)
			exit(1)
#end
def find_singlecopy(BUSCODIR, OUTDIR, ODB, LINEAGE):
	gene_dict = {} #each gene with number of appearences among genomes
	genomes_dirs = os.listdir(BUSCODIR) #read each genome directory
	#ngenomes = len(genomes_dirs) #this does not work if other files are located in the directory
	ngenomes = 0
	genomes_names = []
	common_busco_ids = []
	for genome in genomes_dirs:
		if genome != "batch_summary.txt" and genome != "logs" and genome != "list_common_busco.txt":
			genomes_names.append(genome)
			ngenomes += 1 #number of genomes will be accumulating when read a genome directory
	        #we set the directory with single copy BUSCO groups found in genome
			singlecopypath = os.path.join(BUSCODIR, genome, "run_" + LINEAGE + "_" + "odb" + str(ODB), "busco_sequences/single_copy_busco_sequences")
	        #we read all .faa files
			singlecopyfiles = os.listdir(singlecopypath) #read files in single copy busco sequences directory
	
			for busco in singlecopyfiles: #busco will get the file name of each busco group
				if busco.endswith(".faa"): #check that all files have the same format
					#now we fill the gene dictionary
					if busco not in gene_dict: #if busco was not previously found, we add it to the dictionary
						gene_dict[busco] = 1
					else: # if busco is in the dictionary, we only add one occurrence
						gene_dict[busco] += 1

	#now, after counting all busco groups in each genome, we need to extract only those that are present in all genomes.
	if not os.path.isdir(OUTDIR): #check if output directory exists
		os.mkdir(OUTDIR) #if not, it is created
	outfile = os.path.join(OUTDIR,"list_common_busco.txt") #name of outputfile
	with open(outfile,'wt') as BUSCO:	#writing results in one single file
		for busco,n in gene_dict.items():
			if n == ngenomes: #check if busco is equal to the number of genomes (common)
				BUSCO.write(busco[:-4] + "\n")
				common_busco_ids.append(busco[:-4])
		
	#finally we save the name of the genomes
	namesfile = os.path.join(OUTDIR,"genomes_names.txt") #name of outputfile
	with open(namesfile,'wt') as NAMES:
		for name in genomes_names:
			NAMES.write(name + "\n")

	return genomes_names, common_busco_ids
#end
def create_busco_fasta(BUSCODIR, OUTDIR, ODB, LINEAGE, genomes_names, common_busco_ids, common_busco_dir="common_busco_sequences"):
	out_seq_dir = os.path.join(OUTDIR, common_busco_dir)
	if not os.path.isdir(out_seq_dir):#create directory with fasta files
		os.mkdir(os.path.join(OUTDIR, common_busco_dir))
	for busco_id in common_busco_ids:
		busco_multifasta = os.path.join(out_seq_dir, busco_id + '.faa')
		with open(busco_multifasta, 'wt') as MULTIFASTA:
			for genome in genomes_names:
				busco_fasta_file = os.path.join(BUSCODIR, genome, "run_" + LINEAGE + "_" + "odb" + str(ODB), "busco_sequences/single_copy_busco_sequences", busco_id + '.faa')
				#CHECK FILE EXISTS
				with open(busco_fasta_file, 'rt') as FASTA:
					for line in FASTA:
						if line[0] == '>':
							MULTIFASTA.write(">" + genome + "|" + busco_id + '\n')
						else:
							MULTIFASTA.write(line)
					#MULTIFASTA.write(FASTA.read()) #esto tiene el error de que no se pueden distinguir las especies y las secuencias tendrán todas el mismo nombre
#end
#%% Menu -> is executed when the script is called independently
def usage():
	parser = argparse.ArgumentParser(
		description='''01_find_singlecopy.py recover ID sequences from a BUSCO output directory.''',
		epilog="""End of the help""")

	parser.add_argument('-b', '--buscodir', type=str, required=True, help='Path to the BUSCO output directory. Individual genome directory must be place within this path. Directory names will be take as genomes names.')
	parser.add_argument('-o', '--outdir', type=str, required=True, help='Path to save output files. If it does not exists, it will be created.')
	parser.add_argument('-d', '--odb', type=int, required=False, default=10, help='OrthoDB version. This parameter is only use to complete the path of the single copy sequences directory.')
	parser.add_argument('-l', '--lineage', type=str, required=False, default='lepidoptera', help='OrthoDB lineage. This parameter is only use to complete the path of the single copy sequences directory.')

	return parser.parse_args()
#end
#%% Main program
if __name__ == '__main__':
	args = usage()
	print("Checking paths...")
	start = time()
	check_paths(args.buscodir, args.outdir) #check if BUSCO directory exists and output directory can be created
	
	#find common single copy BUSCO groups
	print("Finding single copy BUSCOs")
	genomes_names, common_busco_ids = find_singlecopy(args.buscodir, args.outdir, args.odb, args.lineage)
	
	#create a directory to put busco sequences in multifasta files
	print("Creating single copy multifasta files")
	#common_busco_output_dir = "common_busco_sequences"
	create_busco_fasta(args.buscodir, args.outdir, args.odb, args.lineage, genomes_names, common_busco_ids)#, common_busco_output_dir)
	print("Common busco search finished...")
	print(f'Time taken to run: {time() - start} seconds.')
#%% End







