# BUSCO2Tree User Guide

---
**Author**
========

[Nicolás Nahuel Moreyra](https://github.com/niconm89/Curriculum_Vitae)

* Instituto de Ecología, Genética y Evolución de Buenos Aires (IEGEBA), CONICET, Argentina
* Departamento de Ecología, Genética y Evolución (EGE), Facultad de Ciencias Exactas y Naturales (FCEN), Universidad de Buenos Aires (UBA), Argentina
+ E-mail: niconm89@gmail.com |  nmoreyra@ege.fcen.uba.ar
+ Twitter: [@nicomoreyra](https://twitter.com/NicoMoreyra)

---
**About BUSCO2Tree**
========

BUSCO2Tree is a pipeline designed to automate downstream phylogenomic analyses for generating a species tree after running BUSCO on multiple genomes/proteomes/transcriptomes. This guide provides an example data set obtained from a recent study [1] where four lepidopteran genomes were assessed using BUSCO v5.5.0 to obtain completeness scores and dipteran single-copy orthologs. These common BUSCOs among all genomes serve as the input for this pipeline.

---
**Table of Contents**
========

-   [Author](#author)
-   [What is BUSCO2Tree?](#about)
-   [Installation Guide](#installation)
    -   [Software Dependencies](#dependencies)
        -   [Required Tools](#mandatory-tools-and-python-modules)
        -   [Optional Tools](#optional-tools)
        -   [Installing Dependencies with Anaconda](#installing-dependencies-with-anaconda)
    -   [Components of BUSCO2Tree](#busco2tree-components)
-   [Example Data](#example-data)
-   [How to Run BUSCO2Tree](#running-BUSCO2tree)
    -   [Command Line Options](#description-of-command-line-options)
    -   [Steps in the BUSCO2Tree Pipeline](#different-BUSCO2tree-pipeline-steps)
-   [Understanding the Output of BUSCO2Tree](#output-of-BUSCO2tree)
-   [License Information](#license-of-BUSCO2tree)


---
**Installation Guide**
========

```
git clone https://github.com/niconm89/BUSCO2Tree.git
cd BUSCO2Tree
./BUSCO2Tree.py -h
```

# Dependencies

Supported software versions
---------------------------

At the time of release, this pipeline was tested with:

-   [BUSCO v5.4.6](https://busco.ezlab.org/)

-   [Python v3.11.0](https://www.python.org/downloads/)

-   [Mafft v7.520](https://mafft.cbrc.jp/alignment/software/)

-   [trimAl v1.4.1](http://trimal.cgenomics.org/downloads)

-   [IQTree v2.2.0.3](http://www.iqtree.org/)

-------
Bioinformatics software dependencies
------------------------------------

BUSCO2Tree calls upon various bioinformatics software tools and uses python modules that are not part of this pipeline. Some tools are mandatory while others are optional. Please install all tools that are required for running BUSCO2Tree in the mode of your choice.

## Mandatory tools and Python modules
Running BUSCO2Tree requires a Linux-system (only tested on this platform) with `bash` and `Python3` (and some modules), but depending on user requirements it may be usefull to install `TrimAl` and `IQTree` as well (see [Running BUSCO2tree](#running-BUSCO2tree)).

### Python3

On Ubuntu, Python3 is usually installed by default, and `python3` will be in your `$PATH` variable, and BUSCO2Tree will automatically locate it. In case your system is not configured in this way, you may include the directory containing the python3 binary to the `$PATH` or just call the pipeline's main script BUSCO2Tree.py using a command like:

`/path/to/python3 /path/to/BUSCO2Tree.py -h`

Note: You may want to simplify

### Python modules

BUSCO2Tree requires the following Python modules to be installed (except to Biopython, the remiaining modules are commonly installed by default with Python):

-   argparse
-   os
-   time
-   subprocess
-   Biopython

### Biopython
On Ubuntu, install Python3 package manager with:

    `sudo apt-get install python3-pip`

Then, install Biopython with:

    `sudo pip3 install biopython`

In case you do not have root permitions, you can follow the steps shown in [Installing dependencies with Anaconda](#installing-dependencies-with-anaconda).

### trimAl
trimAl automatically removes spurious sequences or poorly aligned regions from a multiple sequence alignment.
Download trimAl from [http://trimal.cgenomics.org/downloads](http://trimal.cgenomics.org/downloads) or install it using [Anaconda](#installing-dependencies-with-anaconda).

### IQTree
Download IQTree from [http://www.iqtree.org/#download](http://www.iqtree.org/#download) or install it using [Anaconda](#installing-dependencies-with-anaconda).

##### **Important: modification of the $PATH variable**

Adding directories containing binaries and scripts to your `$PATH` variable enables your system to automatically locate these tools. It is a requirement for running BUSCO2Tree to do this, because BUSCO2Tree will try to find them from the location of the `$PATH` variable. To do this:

a) For your current `bash` session (you will have to repeat these steps after you close your terminal or finish the current bash sessions), type:

```
PATH="/home/nmoreyra/Soft/miniconda3/bin:/home/nmoreyra/Soft/bin:$PATH"
export PATH
```

b) For all your BASH sessions, add the above lines to a linux startup script (e.g.`~/.bashrc`). 

## Optional tools

### Mafft
It is not necessary to install mafft since it is executed using `MafftCommandline` module available in Biopython. However, in case you need to make run Mafft by your own, you can download it from [here](https://mafft.cbrc.jp/alignment/software/) or install it using [Anaconda](#installing-dependencies-with-anaconda).

-------
Installing dependencies with Anaconda
------------------------------------
On Ubuntu, if you do not have root permissions on the Linux machine, you can install the modules with by setting up an **Anaconda** (<https://www.anaconda.com/distribution/>) environment as follows:

```
conda create -n BUSCO2Tree -c anaconda python #it will create an environment with python including os, argparse, time and subprocess modules
conda activate BUSCO2Tree
conda install -c conda-forge biopython
conda install -c bioconda mafft
conda install -c bioconda trimal 
conda install -c bioconda iqtree
```

Alternatively, user may install dependencies from each software's website.

-------
# BUSCO2Tree components

BUSCO2Tree is a collection of five Python scripts. The main script that will be called in order to run the pipeline is [BUSCO2Tree.py](/BUSCO2Tree.py), and additional four Python scripts that employ different stages of the analysis are:

-   Step 1: [find_singlecopy_BUSCOs.py](/scripts/find_singlecopy_BUSCOs.py)

-   Step 2: [align_BUSCOs.py](/scripts/align_BUSCOs.py)

-   Step 3: [create_matrix.py](/scripts/create_matrix.py)

-   Step 4: [phylogenetic_analysis.py](/scripts/phylogenetic_analysis.py)

It is recommended to give executable permissions to all python scripts that are part of BUSCO2Tree as this facilitates the use of the entire pipeline or each script individually. This should already be the case if you download BUSCO2Tree from GitHub.

It is important that the `x` in `-rwxr-xr-x` is present for each script when you run `ls -l BUSCO2Tree.py scripts/*py` . If that is not the case, you may run the following command in order to change file attributes:

    chmod a+x BUSCO2Tree.py scripts/*py

You may find it helpful to add the directories in which `BUSCO2Tree.py` and supplementary python scripts reside to your `$PATH` environment variable. For a single bash session, enter:

```
PATH=/your_path_to_galba/:$PATH
export PATH
```

To make this `$PATH` modification available to all bash sessions, add the above lines to a startup script (e.g.`~/.bashrc`).

---
**Example data**
========
The example data provide here was obtained from a study that was part of my PhD, in which I and my collegues search for genomic innovations that could be associated with cactophily and host plant specialization. The study is [Phylogenomics provides insights into the evolution of cactophily and host plant shifts in *Drosophila*](https://doi.org/10.1016/j.ympev.2022.107653). I employed the genomes of the thirteen species included in this study and run [BUSCO v5.4.6](https://busco.ezlab.org/) on each one in search of [BUSCO groups for the eukaryota lineage obtained from ODB v10](https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2020-09-10.tar.gz).

To test BUSCO2Tree, please decompress the [BUSCO results zip file](/example/BUSCO_results.zip) placed in the [example](/example) directory.

---
**Running BUSCO2Tree**
========
## Description of command line options
~~~~
BUSCO2Tree.py --help

usage: BUSCO2Tree.py [-h] -s STEPS [STEPS ...] [-o OUTDIR] [-t THREADS] [-b BUSCODIR] [-d ODB] [-l LINEAGE] [-f FASTADIR] [-cnf <config file>] [-cmd <command>] [--trim] [--trimparams TRIMPARAMS] [-a ALIGNDIR] [-fmt {nexus,phylip}] [-m MATRIX] [-p PARTITIONS] [-P PREFIX] [-B BOOTSTRAP] [-st {DNA,AA}]

options:
  -h, --help            show this help message and exit

General arguments:
  -s STEPS [STEPS ...], --steps STEPS [STEPS ...]
                        Steps to run. Only single or consecutive steps are allowed. Available options are: 1) find common single copy BUSCOs; 2) align
                        common single copy BUSCOs; 3) create phylogenetic matrix by concatenating BUSCO alignmets; 4) Generate the phylogenetic tree.
  -o OUTDIR, --outdir OUTDIR
                        Path to the directory where the output of each step will be saved. If it does not exists, it will be created. Default:
                        B2T_output
  -t THREADS, --threads THREADS
                        Number of threads to use. Default: 4

Step1: Find common single-copy BUSCO groups among BUSCO outputs:
  -b BUSCODIR, --buscodir BUSCODIR
                        Path to the directory where individual BUSCO outputs are located (you can select the BUSCO output when running in batch mode).
                        Directories names will be taken as the genomes names.
  -d ODB, --odb ODB     OrthoDB version. This parameter is only use to complete the path of the single copy sequences directory. Default: 10
  -l LINEAGE, --lineage LINEAGE
                        OrthoDB lineage. This parameter is only use to complete the path of the single copy sequences directory. Default: eukaryota

Step2: Align common BUSCO groups:
  -f FASTADIR, --fastadir FASTADIR
                        Path to the directory containing the BUSCO groups in fasta format.
  -cnf <config file>, --config <config file>
                        Config file for alignment setting. Users can find a config file template in the docs directory. If no file is provided, the
                        alignments will be done using default parameters.
  -cmd <command>, --command <command>
                        MAFFT parameters to apply to each alingment. The parmeters must be defined in a command line style, between quote marks, and
                        the avoiding the names of in/output files, e.g. "--unalignlevel 0.1 --leavegappyregion --ep 0.12 --globalpair --maxiterate
                        1000".
  --trim                Trim alignments using trimAl in automated mode, which must be available in the path.
  --trimparams TRIMPARAMS
                        trimAl parameters to apply when removing poorly aligned regions. The parameters must be defined in quotes marks and avoiding
                        in/output files, e.g. "-gt 0.3".

Step3: Create phylogenetic matrix and partition files:
  -a ALIGNDIR, --aligndir ALIGNDIR
                        Path to the directory containing the alignments in fasta format (step 3).
  -fmt {nexus,phylip}, --format {nexus,phylip}
                        Format for the phylogenetic matrix generated by concatenating BUSCOs alignment files (step 3). Default: phylip

Step4: Genetare phylogenetic tree:
  -m MATRIX, --matrix MATRIX
                        Path to the phylogenetic matrix file.
  -p PARTITIONS, --partitions PARTITIONS
                        Path to the partitions file in nexus format.
  -P PREFIX, --prefix PREFIX
                        Prefix to name the output dataset and results. Default: iqtree
  -B BOOTSTRAP, --bootstrap BOOTSTRAP
                        Number of bootstrap replicate to run. Default: 1000
  -st {DNA,AA}, --seqtype {DNA,AA}
                        Type of sequence. Default: AA

End of the help
~~~~

## Different BUSCO2tree pipeline steps
The commands shown below were tested using example data.
### Single-steps commands
#### Step 1. Find common BUSCO groups among set of genomes
`BUSCO2Tree.py --steps 1 --buscodir Drosophila_proteins --lineage eukaryota --outdir B2T --threads 8`

#### Step 2. Align common BUSCO groups
`BUSCO2Tree.py --steps 2 --fastadir B2T/01_single-copy/common_busco_sequences --outdir B2T --threads 8`

**Note**: it is possible to activate alignment trimming by setting `--trim` parameter. In addition, the parameter `--trimparams` can be set when it is necessary to use a timAl user-defined command.

#### Step 3. Generate phylogenetic matrix and partitions files
`BUSCO2Tree.py --steps 3 --aligndir B2T/02_alignments --outdir B2T --threads 8`

**Note**: if step 2 was run using the --trim/trimparams parameters, trimmed alignments will be placed in BT2/02_alignments/trimAl.

#### Step 4. Build phylogenetic tree
`BUSCO2Tree.py --steps 4 --matrix B2T/03_matrix/matrix.phylip --partitions B2T/03_matrix/busco_coords.partitions.tsv --outdir B2T --threads 8`

---

### Multiple-steps commands
#### Running the complete pipeline.
`BUSCO2Tree.py -s 1 2 3 4 -b Drosophila_proteins --trim -o B2T -t 8`

Note: this command takes ~23 minutes.

#### Step 1 & 2.
`BUSCO2Tree.py -s 1 2 -b Drosophila_proteins --config /path/to/BUSCO2Tree/example/mafft_config.txt --trim --trimparams "-gt 0.3" -o B2T -t 8`

#### Step 1, 2 & 3.
`BUSCO2Tree.py -s 1 2 3 -b Drosophila_proteins -cnf mafft_config.txt --trim --trimparams "-gt 0.3" --format nexus -o B2T -t 8`

#### Step 2 & 3.
`BUSCO2Tree.py -s 2 3 -f B2T/01_single-copy/common_busco_sequences --trim --trimparams "-gt 0.3" -o B2T -t 8`

#### Step 2, 3 & 4.
`BUSCO2Tree.py -s 2 3 4 -f B2T/01_single-copy/common_busco_sequences --trim --trimparams "-gt 0.3" -B 500 -P TEST234 -o B2T -t 8`

#### Step 3 & 4.
`BUSCO2Tree.py -s 3 4 -a B2T/02_alignments --trim -B 1500 -P TEST34 -o B2T -t 8`

---
**Output of BUSCO2tree**
========
Running the complete pipeline will generate an output directory containing four main folder, one per each step:
~~~~
B2T/
├── 01_single-copy
│   ├── common_busco_sequences
│   |   ├── #common BUSCO groups fasta files
│   ├── genomes_names.txt
│   └── list_common_busco.txt
├── 02_alignments
│   ├── #common BUSCO groups alignments
│   ├── trimAl
│   │   ├── #trimmed alignments
├── 03_matrix
│   ├── busco_coords.partitions.tsv
│   └── matrix.phylip
└── 04_phylogenetic_tree
    ├── iqtree.log
    └── iqtree.model.gz
~~~~

- **B2T/01_single-copy/common_busco_sequences**: it contains common BUSCO groups among species in fasta files.
- **B2T/01_single-copy/genomes_names.txt**: file with names of the species placed in the directory set for --buscodir.
- **B2T/01_single-copy/list_common_busco.txt**: file with IDs of common BUSCO group.
- **B2T/02_alignments**: it contains alignments of common BUSCO groups in fasta format.
- **B2T/02_alignments/trimAl**: it contains trimmed alignments in fasta format.
- **B2T/03_matrix/busco_coords.partitions.tsv**: partition file in nexus format. Each partition corresponds to a BUSCO group.
- **B2T/03_matrix/matrix.phylip**: phylogenetic matrix file in phylip format (nexus is also allowed).
- **B2T/04_phylogenetic_tree**: it contains the IQtree outputs.

When running BUSCO2Tree step by step (or using combinations of them), it is recommended to place the results (--outdir) in the same main directory where subdirectories containing files (fasta, aligment, and matrix/partitions files) are. For instance, following the tree directory shown above, if we want to run step 3:
~~~~
BUSCO2Tree.py --steps 3 --aligndir B2T/02_alignments/trimAl --outdir B2T --threads 8
~~~~
Check that --aligndir in set with a directry placed within the B2T directory and --outdir is thus set with B2T. As a result, a directory named "03_matrix" will be created within B2T.

---
License
========
[GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007](/LICENSE).
