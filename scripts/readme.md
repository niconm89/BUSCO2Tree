This folder contains scripts to extract BUSCO genes from several BUSCO run folders (! having used the same BUSCO dataset) and produce a species phylogeny.

- `run.sh`: the main script to proceed to the phylogenomics analysis. The final output is final.nwk
- `extract_buscos_phylo.py`: the main script to identify and extract BUSCO sequences
- `fetch_best_sequence.py`: an additional script required by the previous one
- `superalignment.py`: it fuses all single gene alignement together

You will need to adapt several scripts for your own use, each files are commented.
For feedbacks and help, please use the main BUSCO gitlab project: https://gitlab.com/ezlab/busco
