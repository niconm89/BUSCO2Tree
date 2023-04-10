# BUSCO2Tree User Guide

**This GIT README is being created!!!!**

---
**Author**
========

[Nicolás Nahuel Moreyra](https://github.com/niconm89/Curriculum_Vitae)

* Instituto de Ecología, Genética y Evolución de Buenos Aires (IEGEBA), CONICET, Argentina
* Departamento de Ecología, Genética y Evolución (EGE), Facultad de Ciencias Exactas y Naturales (FCEN), Universidad de Buenos Aires (UBA), Argentina
+ E-mail: niconm89@gmail.com |  nmoreyra@ege.fcen.uba.ar
+ Twitter: [@nicomoreyra](https://twitter.com/NicoMoreyra)

---
**About**
========

This pipeline was developed to automatize the downstream phylogenomic analyses to obtain a species tree after running BUSCO over several genomes/proteomes/transcriptomes.
We provide one data example to test funcionallity. Data was obtained from a study we have recently published [1]. Briefly, thirteen Drosophila genomes were assessed using BUSCO v5.4.6 (REF) to obtain completeness scores and dipteran single-copy orthologs. Common BUSCOs among all genomes were the input of this pipeline.

---
Contents
========

-   [Author](#author)
-   [About BUSCO2tree?](#about)
-   [Installation](#installation)
    -   [Dependencies](#dependencies)
        -   [Mandatory tools](#mandatory-tools-and-python-modules)
        -   [Optional tools](#optional-tools)
    -   [BUSCO2tree](#busco2tree-components)
-   [Running BUSCO2tree](#running-BUSCO2tree)
    -   [BUSCO2tree pipeline modes](#different-BUSCO2tree-pipeline-modes)
    -   [Description of command line options](#description-of-command-line-options)
-   [Output of BUSCO2tree](#output-of-BUSCO2tree)
-   [Example data](#example-data)
    -   [Data description](#data-description)
    -   [Testing BUSCO2tree with genomes](#testing-braker-with-genomes)
    -   [Testing BUSCO2tree with proteins](#testing-braker-with-proteins)
    -   [Testing BUSCO2tree with transcripts](#testing-braker-with-transcripts)


---

Installation
========

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
BUSCO2Tree components
------------------------------------

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

Example data
========
The example data provide here was obtained from a study that was part of my PhD, in which I and my collegues search for genomic innovations that could be associated with cactophily and host plant specialization. The study is [Phylogenomics provides insights into the evolution of cactophily and host plant shifts in *Drosophila*](https://doi.org/10.1016/j.ympev.2022.107653). I employed the genomes of the thirteen species included in this study and run [BUSCO v5.4.6](https://busco.ezlab.org/) on each one in search of [BUSCO groups for the eukaryota lineage obtained from ODB v10](https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2020-09-10.tar.gz).

To test BUSCO2Tree, please decompress the [BUSCO results zip file](/example/BUSCO_results.zip) placed in the [example](/example) directory.

---
Running BUSCO2Tree
========

## Tests with example data
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
`BUSCO2Tree.py --steps 4 -m B2T/03_matrix/matrix.phy -p B2T/03_matrix/busco_coords.partitions.nexus --outdir B2T --threads 8`

---

### Multiple-steps commands
#### Step 1 & 2.
`BUSCO2Tree.py -s 1 2 -b Drosophila_proteins --config /path/to/BUSCO2Tree/example/mafft_config.txt --trim --trimparams "-gt 0.3" -o B2T -t 8`

#### Step 1, 2 & 3.
`BUSCO2Tree.py -s 1 2 3 -b Drosophila_proteins -cnf mafft_config.txt --trim --trimparams "-gt 0.3" --format nexus -o B2T -t 8`

#### Step 2 & 3.
`BUSCO2Tree.py -s 2 3 -b B2T/01_single-copy -cnf mafft_config.txt --trim --trimparams "-gt 0.3" -o B2T -t 8`

#### Step 2, 3 & 4.
`BUSCO2Tree.py -s 2 3 4 -b B2T/01_single-copy --trim --trimparams "-gt 0.3" -B 500 -P TEST234 -o B2T -t 8`

#### Step 3 & 4.
`BUSCO2Tree.py -s 3 4 -a B2T/02_alignments --trim -B 500 -P TEST34 -o B2T -t 8`

#### Running the complete pipeline.
`BUSCO2Tree.py -s 1 2 3 4 -b Drosophila_proteins --trim -o B2T -t 8`


