#!/home/nmoreyra/Software/miniconda3/bin/python

import os


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
		singlecopypath = os.path.join(BUSCOPATH, genome, "run_" + LINEAGE + "_" + ODB, "busco_sequences/single_copy_busco_sequences")
		singlecopyfiles = os.listdir(singlecopypath) #read files in single copy busco sequences directory

		for busco in singlecopyfiles:
			if busco.endswith(".faa"): #just in case some weird file appears in the directory
				#now we fill the gene dictionary
				if busco not in gene_dict:
					gene_dict[busco] = 1
				else:
					gene_dict[busco] += 1

#now, after reading all busco in each genome, we need to extract only those busco genes that are present in all genomes.
outfile = os.path.join(WD,"common_busco.txt")
common_busco = []
for busco,n in gene_dict.items():
	if n == ngenomes:
		common_busco.append(busco)

with open(outfile,'wt') as OUT:
	for busco in common_busco:
		OUT.write(busco + "\n")
		
		
		
		
		
		
		
		

