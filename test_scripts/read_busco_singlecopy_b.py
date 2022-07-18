#!/home/nmoreyra/Software/miniconda3/bin/python

import os
import csv
import sys #argv

ODB = "odb10"
LINEAGE = "lepidoptera"
#BUSCOPATH = "/home/nmoreyra/Cactoblastis/02_Assembly/03_BUSCO/02_Phylogenomics/02_BUSCO"
BUSCOPATH = sys.argv[1] #user must pass the busco directory
OUTDIR = sys.argv[2] #user must pass the output directory

gene_dict = {} #acá se van a ir poniendo las apariciones de todos los genes en cada genoma
genomes_dirs = os.listdir(BUSCOPATH)
#ngenomes = len(genomes_dirs)
ngenomes = 0
genomes_names = []
for genome in genomes_dirs:
	if genome != "batch_summary.txt" and genome != "logs" and genome != "common_busco.txt":
		genomes_names.append(genome)
		ngenomes += 1
		#~/Cactoblastis/02_Assembly/03_BUSCO/02_Phylogenomics/02_BUSCO/Cbuc.l3.pd.p.fa/run_lepidoptera_odb10/full_table.tsv
		fulltablepath = os.path.join(BUSCOPATH, genome, "run_" + LINEAGE + "_" + ODB, "full_table.tsv")
		#'''
		with open(fulltablepath, 'rt') as tsvfile:
			tsvreader = csv.reader(tsvfile, delimiter="\t")
			for line in tsvreader:
				if line[0][0] != '#':
					buscoid,status = line[0:2]
					if status == "Complete": #test time for if line[1] == "complete": and work with line[0] for buscoid
						if buscoid not in gene_dict:
							gene_dict[buscoid] = 1
						else:
							gene_dict[buscoid] += 1
		'''
		with open(fulltablepath, 'rt') as tsvfile:
			tsvreader = tsvfile.readlines()
			for line in tsvreader:
				if line[0] != '#':
					buscoid,status = line.split("\t")[0:2]
					if status == "complete":
						if buscoid not in gene_dict:
							gene_dict[buscoid] = 1
						else:
							gene_dict[buscoid] += 1
		'''					
#now, after reading all busco in each genome, we need to extract only those busco genes that are present in all genomes.
if not os.path.isdir(OUTDIR):
	os.mkdir(OUTDIR)
outfile = os.path.join(OUTDIR,"common_busco.txt")
common_busco = []
for busco,n in gene_dict.items():
	if n == ngenomes:
		common_busco.append(busco)

with open(outfile,'wt') as OUT:
	for busco in common_busco:
		OUT.write(busco + "\n")

#finally we save the name of the genomes
namesfile = os.path.join(OUTDIR,"genomes_names.txt") #name of outputfile
with open(namesfile,'wt') as NAMES:
	for name in genomes_names:
		NAMES.write(name + "\n")

