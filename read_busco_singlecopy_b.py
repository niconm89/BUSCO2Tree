#!/home/nmoreyra/Software/miniconda3/bin/python

import os
import csv

ODB = "odb10"
LINEAGE = "lepidoptera"
BUSCOPATH = "/home/nmoreyra/Cactoblastis/02_Assembly/03_BUSCO/02_Phylogenomics/02_BUSCO"
WD = BUSCOPATH

gene_dict = {} #acá se van a ir poniendo las apariciones de todos los genes en cada genoma
genomes_dirs = os.listdir(BUSCOPATH)
#ngenomes = len(genomes_dirs)
ngenomes = 0
for genome in genomes_dirs:
	if genome != "batch_summary.txt" and genome != "logs" and genome != "common_busco.txt":
		ngenomes += 1
		#~/Cactoblastis/02_Assembly/03_BUSCO/02_Phylogenomics/02_BUSCO/Cbuc.l3.pd.p.fa/run_lepidoptera_odb10/full_table.tsv
		fulltablepath = os.path.join(BUSCOPATH, genome, "run_" + LINEAGE + "_" + ODB, "full_table.tsv")
		
		with open(fulltablepath, 'rt') as tsvfile:
			tsvreader = csv.reader(tsvfile, delimiter="\t")
			for line in tsvreader:
				buscoid,status = line[0:2]
				if status == "complete": #test time for if line[1] == "complete": and work with line[0] for buscoid
					if buscoid not in gene_dict:
						gene_dict[buscoid] = 1
					else:
						gene_dict[buscoid] += 1
		'''
		with open(fulltablepath, 'rt') as tsvfile:
			tsvreader = tsvfile.readlines()
			for line in tsvreader:
				if line != '#':
					buscoid,status = line.split("\t")[0:2]
					if status == "complete":
						if buscoid not in gene_dict:
							gene_dict[buscoid] = 1
						else:
							gene_dict[buscoid] += 1
		'''					
#now, after reading all busco in each genome, we need to extract only those busco genes that are present in all genomes.
outfile = os.path.join(WD,"common_busco.txt")
common_busco = []
for busco,n in gene_dict.items():
	if n == ngenomes:
		common_busco.append(busco)

with open(outfile,'wt') as OUT:
	for busco in common_busco:
		OUT.write(busco + "\n")

		

