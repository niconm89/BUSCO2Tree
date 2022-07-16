#!/home/nmoreyra/Software/miniconda3/envs/orthomcl/bin/python

#%% Imports
from os import path as path_check
import os
import shutil
import sys
from Bio import SeqIO

#%% Function definition
def convert_trans2prot(convertion):
	print("Creating convertion prot to trans dictionary.\n")
	converted = {}
	with open(convertion,'rt') as convtab:
		table = convtab.readlines()
		for line in table:
			line = line.rstrip()
			columns = line.split(" ")
			if len(columns) == 2:
				converted[columns[1]] = columns[0]
	return converted
#end

def dict_trans(transcripts_path,biopython=False):
	'transcript path is the location where all species transcripts are.'
	print("Creating transcriptomes dictionary.\n")
	seqs = {}
	with open(transcripts_path,"rt") as transfile:
		if biopython == False:
			transcriptome = transfile.readlines()
			for line in transcriptome:
				line = line.rstrip() #let's discard the newline at the end (if any)
				if len(line) > 0 and line[0] == '>': #identify the id of each sequence
					words = line.split()
					name = words[0][1:]
					seqs[name] = ''
				else:
					seqs[name] = seqs[name] + line
		else:
			TRANS = SeqIO.parse(transfile,'fasta')
			for i,tran in enumerate(TRANS):
				seqs[str(tran.id)] = str(tran.seq)
	return seqs
# end


#%% Main program

workdir = '/home/nmoreyra/Documents/Genomics/03_orthologs_search/00_Phylogenomics_BUSCO-groups/03_divergence-times_1pergene/00_transcripts'
busco_folder = '/home/nmoreyra/Documents/Genomics/03_orthologs_search/00_Phylogenomics_BUSCO-groups/01_extracted_1pergene'
transcripts_folder = '/home/nmoreyra/Documents/Genomics/03_orthologs_search/01_orthogroups/transcripts_data/01_modified-with-P450'
transcripts_path = '/home/nmoreyra/Documents/Genomics/03_orthologs_search/01_orthogroups/transcripts_data/01_modified-with-P450/transcripts.all.fasta'
convertion = '/home/nmoreyra/Documents/Genomics/03_orthologs_search/01_orthogroups/transcripts_data/01_modified-with-P450/convertion-table.trans2prot'
busco_seqs = '/home/nmoreyra/Documents/Genomics/03_orthologs_search/00_Phylogenomics_BUSCO-groups/busco4_odb10/results_odb10_1prot'
species_list = ["Dald","Dari","Dato","Dbrb","Dbuz","Dhyd","DkoeA","DkoeB","Dmel","Dmoj","Dnav","Drep","Dvir"]

#generate dictionary with convertions like FBpp2FBtr
convertion_table = convert_trans2prot(convertion)

#generate dictionary with transcripts of species
transcriptome_dict = dict_trans(transcripts_path, True)

for busco in os.listdir(busco_folder):
	if busco[-5:] == '7.faa':
		print(busco)
		busco_transcript = os.path.join(workdir,busco)
		with open(busco_transcript,'wt') as OUT:
			for species in species_list:
				folder = 'run_' + species + '_diptera_odb10/busco_sequences/single_copy_busco_sequences'
				busco_species = os.path.join(busco_seqs,folder,busco)
				with open(busco_species, 'rt') as IN:
					busco_seq = IN.readlines()
					for line in busco_seq:
						line = line.rstrip()
						if line[0] == '>':
							description = line[1:]
							prot_name = description.split(' ')[0]
							trans_name = convertion_table[prot_name]
							OUT.write('>'+trans_name+'\n'+transcriptome_dict[trans_name]+'\n')


