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
    'Check paths of BUSCO nad output directories.'
    if not os.path.isdir(BUSCODIR):
        raise ValueError("Error: BUSCO output directory %s does not exist. Check your paths! " % BUSCODIR)

    if '/' in OUTDIR:
        elements_path = OUTDIR.split('/')
        parentdir = '/'.join(elements_path[:-1])
        if not os.path.isdir(parentdir):
            raise ValueError("Output directory %s can not be found or created in the path you have passed. Parent directory does not exist! " % parentdir)
#end
def find_singlecopy(BUSCODIR, OUTDIR, ODB, LINEAGE):
	'Find common BUSCO groups among several genomes (BUSCO output directories) and creates a multi FASTA file per BUSCO group.'
	gene_dict = {} #each gene with number of appearences among genomes
	genomes_dirs = os.listdir(BUSCODIR) #read each genome directory
	#ngenomes = len(genomes_dirs) #this does not work if other files are located in the directory
	ngenomes = 0
	genomes_names = []
	common_busco_ids = []
	for genome in genomes_dirs:
        #in case your BUSCO output directory contain other files/dirs, you must remove them or modify this if sentence.
		if genome != "batch_summary.txt" and genome != "logs" and genome != "list_common_busco.txt":
			genomes_names.append(genome)
			ngenomes += 1 #number of genomes will be accumulating when read a genome directory
	        #setting the directory with all single copy BUSCO groups found in the genome
			singlecopypath = os.path.join(BUSCODIR, genome, "run_" + LINEAGE + "_" + "odb" + str(ODB), "busco_sequences/single_copy_busco_sequences")
	        #reading all *.faa files
			singlecopyfiles = os.listdir(singlecopypath)
	
			for busco in singlecopyfiles: #the busco var will get the file name of each busco group
				if busco.endswith(".faa"): #check that all files have the same format
					#now we fill the gene dictionary
					if busco not in gene_dict: #if busco was not previously found, we add it to the dictionary
						gene_dict[busco] = 1
					else: # if busco is in the dictionary, we only add one occurrence
						gene_dict[busco] += 1

	#now, after counting all busco groups in each genome, we need to extract only those that are present in all genomes (ngenomes ocurrences).
	if not os.path.isdir(OUTDIR): #check if output directory exists
		os.mkdir(OUTDIR) #if not, it is created
	outfile = os.path.join(OUTDIR,"list_common_busco.txt") #name of outputfile
	with open(outfile,'wt') as BUSCO:	#writing results in one single file
		for busco,n in gene_dict.items():
			if n == ngenomes: #check if busco was found in all genomes (it is common)
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
	'Creates a multi FASTA files por each common BUSCO group.'
	out_seq_dir = os.path.join(OUTDIR, common_busco_dir) #output directory to save the fasta files
	if not os.path.isdir(out_seq_dir): #create directory if it does not exist
		os.mkdir(os.path.join(OUTDIR, common_busco_dir))
	for busco_id in common_busco_ids:
		busco_multifasta = os.path.join(out_seq_dir, busco_id + '.faa') #BUSCO groupo multi FASTA file
		with open(busco_multifasta, 'wt') as MULTIFASTA:
			for genome in genomes_names: #reading the busco protein of each genome to save it within the multifasta file
				busco_fasta_file = os.path.join(BUSCODIR, genome, "run_" + LINEAGE + "_" + "odb" + str(ODB), "busco_sequences/single_copy_busco_sequences", busco_id + '.faa')
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
		description='''01_find_singlecopy.py recover common ID sequences from several BUSCO output directories and creates multi FASTA files of each common BUSCO group.''',
		epilog="""End of the help""")

	parser.add_argument('-b', '--buscodir', type=str, required=True, help='Path to the BUSCO output directory. The genome directories must be placed within this path. The names of the directories will be taken as the names of the genomes.')
	parser.add_argument('-o', '--outdir', type=str, required=True, help='Directory path where output files will be saved. If it does not exists, it will be created.')
	parser.add_argument('-d', '--odb', type=int, required=False, default=10, help='OrthoDB version. This parameter is only use to complete the path of the single copy sequences directory.')
	parser.add_argument('-l', '--lineage', type=str, required=False, default='lepidoptera', help='OrthoDB lineage. This parameter is only use to complete the path of the single copy sequences directory.')

	return parser.parse_args()
#end
#%% Main program
if __name__ == '__main__':
	args = usage()
	start = time() #time 0
	print("Step 0: Checking paths...") #check if busco and output directories exist
	check_paths(args.buscodir, args.outdir) #check if BUSCO directory exists and output directory can be created
	
	#find common single copy BUSCO groups
	print("Step 1: Finding single copy BUSCOs...")
	genomes_names, common_busco_ids = find_singlecopy(args.buscodir, args.outdir, args.odb, args.lineage)
	
	#create a directory to place the BUSCO sequence of each genome into a multifasta file
	print("Step 2: Creating single copy multifasta files...")
	#common_busco_output_dir = "common_busco_sequences"
	create_busco_fasta(args.buscodir, args.outdir, args.odb, args.lineage, genomes_names, common_busco_ids)
	print("The search for common BUSCO groups is complete...")
	print(f'Time taken to run: {time() - start} seconds.') #time out
#%% End