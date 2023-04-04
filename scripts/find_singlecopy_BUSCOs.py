#!/usr/bin/env python3
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
def check_paths(BUSCODIR, OUTDIR,out_seq_dir=False):
    'Check paths of BUSCO and output directories. Creates output directory if it does not exist.'
    if not os.path.isdir(BUSCODIR):
        raise ValueError("Error: BUSCO output directory %s does not exist." % BUSCODIR)
    if not os.path.isdir(OUTDIR):#check if output directory exists
        try:
            os.mkdir(OUTDIR)#if not, it is created
        except:
            raise ValueError("Error: Output directory does not exist and can not be created in the path %s." % OUTDIR)
    if out_seq_dir:
        if not os.path.isdir(out_seq_dir): #create directory if it does not exist
            try:
                os.mkdir(out_seq_dir)
            except:
                raise ValueError("Error: Output directory to save sequences does not exist and can not be created at %s." % out_seq_dir)
#end
def find_singlecopy(BUSCODIR, OUTDIR, ODB, LINEAGE):
	'Find common BUSCO groups among several genomes (BUSCO output directories) and creates a multi FASTA file per BUSCO group.'
	check_paths(BUSCODIR, OUTDIR) #check if BUSCO directory exists and output directory can be created
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
			singlecopyfiles = os.listdir(singlecopypath)#reading all *.faa files
			for busco in singlecopyfiles: #the busco var get the file name of each busco group
				if busco.endswith(".faa"): #check that all files have the same format
					if busco not in gene_dict: #fill the gene dictionary, if busco is not in the dictionary, we add it
						gene_dict[busco] = 1
					else: # if busco is in the dictionary, add one occurrence
						gene_dict[busco] += 1
	outfile = os.path.join(OUTDIR,"list_common_busco.txt") #name of outputfile
	with open(outfile,'wt') as BUSCO: #writing results in one single file
		for busco,n in gene_dict.items():
			if n == ngenomes: #check if busco was found in all genomes (common)
				BUSCO.write(busco[:-4] + "\n")
				common_busco_ids.append(busco[:-4])
	namesfile = os.path.join(OUTDIR,"genomes_names.txt") #name of outputfile
	with open(namesfile,'wt') as NAMES:#we save the name of the genomes
		genomes_names = sorted(genomes_names)
		for name in genomes_names:
			NAMES.write(name + "\n") 
	return genomes_names, common_busco_ids
#end
def create_busco_fasta(BUSCODIR, OUTDIR, ODB, LINEAGE, genomes_names, common_busco_ids, common_busco_dir="common_busco_sequences"):
	'Creates a multi FASTA file por each common BUSCO group.'
	out_seq_dir = os.path.join(OUTDIR, common_busco_dir) #output directory to save the fasta files
	check_paths(BUSCODIR, OUTDIR, out_seq_dir) #check if BUSCO and output directories exist (using this function make sense when it is imported)
	for busco_id in common_busco_ids:
		busco_multifasta = os.path.join(out_seq_dir, busco_id + '.faa') #BUSCO groupo multi FASTA file
		with open(busco_multifasta, 'wt') as MULTIFASTA:#read the busco protein of each genome to save it within the multifasta file
			for genome in genomes_names:
				busco_fasta_file = os.path.join(BUSCODIR, genome, "run_" + LINEAGE + "_" + "odb" + str(ODB), "busco_sequences/single_copy_busco_sequences", busco_id + '.faa')
				with open(busco_fasta_file, 'rt') as FASTA:
					#MULTIFASTA.write(FASTA.read()) #esto tiene el error de no poder distinguir las especies ya que las secuencias tendrán el mismo nombre
					for line in FASTA:
						if line[0] == '>':
							MULTIFASTA.write(">" + genome + "|" + busco_id + '\n')
						else:
							MULTIFASTA.write(line)
#end
#%% Menu -> is executed when the script is called independently
def usage():
	parser = argparse.ArgumentParser(
		description='''find_singlecopy.py recover common ID sequences from several BUSCO output directories and creates a multi FASTA file of each common BUSCO group.''',
		epilog="""End of the help""")
	parser.add_argument('-b', '--buscodir', type=str, required=True, help='Path to the BUSCO output directory. Whithin this directory must be placed the directories of each genome (BUSCO output). The names of the directories will be used as genome names.')
	parser.add_argument('-o', '--outdir', type=str, required=False, default='01_single-copy' , help='Path to the directory to save output files. If it does not exists, it will be created. Default: 01_single-copy')
	parser.add_argument('-d', '--odb', type=int, required=False, default=10, help='OrthoDB version. This parameter is only use to complete the path of the single copy sequences directory. Please note that the BUSCO2tree pipeline was developed to process the results of BUSCO v5, so some path and name of directories within BUSCO output may change in past and future versions.')
	parser.add_argument('-l', '--lineage', type=str, required=False, default='eukaryota', help='OrthoDB lineage. This parameter is only use to complete the path of the single copy sequences directory.Please note that the BUSCO2tree pipeline was developed to process the results of BUSCO v5, so some path and name of directories within BUSCO output may change in past and future versions.')
	return parser.parse_args()
#end
#%% Main program
if __name__ == '__main__':
	args = usage()
	start = time() #time 0
	print("Step 1.1: Looking for single-copy BUSCO groups in the lineage %s obtained from the ODB v%d database..." % (args.lineage,args.odb))
	genomes_names, common_busco_ids = find_singlecopy(args.buscodir, args.outdir, args.odb, args.lineage)
	print("Step 1.2: Creating single-copy BUSCO multifasta files...")
	create_busco_fasta(args.buscodir, args.outdir, args.odb, args.lineage, genomes_names, common_busco_ids)
	print("The search for common BUSCO groups has been completed.")
	print(f'Time taken to run: {time() - start} seconds.') #time out
#%% End