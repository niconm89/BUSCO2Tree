#!/usr/bin/env python3
# coding: utf-8
"""
.. module::
   :synopsis: Extract the best BUSCOs sequence for which a single copy complete sequence is found in all species
.. moduleauthor:: Mathieu Seppey <mathieu.seppey@unige.ch>
.. versionadded:: 1.0

From BUSCO run folders, this code extracts complete single copy BUSCO genes

You will need to do something or edit the script where you find !!! >

"""

import os
import fetch_best_sequence
import re

# Config here
# !!! > you need to untar unzip the files in the run folder before running this script, if you used the -z option of BUSCO
# - single_copy_busco_sequences
# - augustus_output/extracted_proteins

wd = '/workdir/'  # !!! > place all BUSCO run folders here
output = '/workdir/extracted/'  # !!! > output of this script will be here. Create the folder before running the script.
genesets = '/workdir/proteins' # !!! > complete protein set for each species, in fasta, should be placed here

if wd[-1] != '/':
    wd += '/'

if output[-1] != '/':
    output += '/'

if genesets[-1] != '/':
    genesets += '/'

# !!! > Write the following mapping manually. 
# key = what is after run_ in your BUSCO result folders. 
# 1st entry, the name you want to see in the tree
# If protein, write the name of the fasta containing the proteome as 3rd entry.
# If transcriptome, just write tran, no 3rd entry, the sequence exists in the busco run folder
# run_for_phylogeny_TEST => 'for_phylogeny_TEST': ['TEST', 'prot', 'protein_file_for_TEST.faa']

#mapping_name = {
#'': ['', 'prot',''],
#'': ['', 'prot',''],
#'': ['', 'prot',''],
#'': ['', 'tran'],
#'': ['', 'tran'],
#'': ['', 'tran']
#}

mapping_name = {
'for_phylogeny_TEST': ['TEST', 'prot', 'TEST.faa'],
'for_phylogeny_TEST_tran': ['TEST_tran', 'tran'],
}

# First, clean the output files.
for file in os.listdir(output):
    if '.faa' in file:
        os.remove('%s%s' % (output, file))

species = {}

# Then init a dict for each species and each kind of prediction (C,D,F,M)
for content in os.listdir(wd):
    if content.startswith('run_') and content.split('run_')[-1] in list(mapping_name.keys()):
        species.update({'_'.join(content.split('_')[1:]): {'Complete': [], 'Duplicated': [],
                                                           'Fragmented': [], 'Missing': []}})

# Determine which BUSCO to keep
complete = set([])
for a_species in species:
    try:
        full_table = open('%srun_%s/full_table_%s.tsv' % (wd, a_species, a_species), 'r')
    except FileNotFoundError as e:
        continue
    complete_tmp = []
    for full_table_line in full_table:
        try:
            if 'Complete' in full_table_line.split()[1]: #  !!! >  if you don't want the best duplicated
            #if 'Complete' in full_table_line.split()[1] or 'Duplicated' in full_table_line.split()[1]: # if you want to include the best scoring duplicated. In our experience, it is less reliable, but we have not really benchmarked that.
                complete_tmp.append(full_table_line.split()[0])
        except IndexError:
            pass
    if complete:
        complete = complete.intersection(set(complete_tmp))
    else:
        complete = set(complete_tmp)

# Now using this list, retrieve all sequences, one file for each BUSCO. When duplicate, keep the best score
for busco in complete:
    output_file = open('%s%s.faa' % (output, busco), 'a')
    for species in mapping_name:
        if mapping_name[species][1] == 'prot':
            seq = fetch_best_sequence.fetch(busco, 'prot', '%srun_%s' % (wd, species),
                                            '%s' % species,
                                            '%s%s' % (genesets, mapping_name[species][2]))
        else:
            seq = fetch_best_sequence.fetch(busco, 'tran', '%srun_%s' % (wd, species),
                                            '%s' % species)

        search_string = r'^>.*\n'
        seq = re.sub(search_string, '>%s\n' % mapping_name[species][0], seq)
        output_file.write('%s\n' % seq)

