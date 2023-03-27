# BUSCO2Tree User Guide

**This GIT README is being created!!!!**

**Author**

[Nicolás Nahuel Moreyra](https://github.com/niconm89/Curriculum_Vitae)

* Instituto de Ecología, Genética y Evolución de Buenos Aires (IEGEBA), CONICET, Argentina
* Departamento de Ecología, Genética y Evolución (EGE), Facultad de Ciencias Exactas y Naturales (FCEN), Universidad de Buenos Aires (UBA), Argentina
+ E-mail: niconm89@gmail.com |  nmoreyra@ege.fcen.uba.ar
+ Twitter: [@nicomoreyra](https://twitter.com/NicoMoreyra)

---
**About BUSCO2Tree**

This pipeline was developed to automatize the downstream phylogenomic analyses to obtain a species tree after running BUSCO over several genomes/proteomes/transcriptomes.
We provide one data example to test funcionallity. Data was obtained from a study we have recently published [1]. Briefly, thirteen Drosophila genomes were assessed using BUSCO v5.4.6 (REF) to obtain completeness scores and dipteran single-copy orthologs. Common BUSCOs among all genomes were the input of this pipeline.

---
Contents
========

-   [Authors](#Author)
-   [What is BUSCO2tree?](#About-BUSCO2tree)
-   [Installation](#installation)
    -   [Dependencies](#dependencies)
            -   [Mandatory tools](#mandatory-tools)
            -   [Optional tools](#optional-tools)
    -   [BUSCO2tree](#BUSCO2tree)
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

## Dependencies
The long way around is to manually install all dependencies of GALBA.

Supported software versions
---------------------------

At the time of release, this pipeline was tested with:

-   BUSCO v5.4.6

-   Python v3.11.0

-   Mafft v7.520

-   TrimAl v1.4.1

-   IQTree v2.2.0.3


-------
Bioinformatics software dependencies
------------------------------------

BUSCO2Tree calls upon various bioinformatics software tools and uses python modules that are not part of this pipeline. Some tools are mandatory while others are optional. Please install all tools that are required for running BUSCO2Tree in the mode of your choice.

### Mandatory tools and Python modules
Running BUSCO2Tree requires a Linux-system with `bash`, Python, several python modules, TrimAl and IQTree.

#### Python3

On Ubuntu, Python3 is usually installed by default, `python3` will be in your `$PATH` variable, by default, and GALBA will automatically locate it. However, you have the option to specify the `python3` binary location in two other ways:

1.  Export an environment variable `$PYTHON3_PATH`, e.g. in your `~/.bashrc` file:

        export PYTHON3_PATH=/path/to/python3/

2.  Specify the command line option `--PYTHON3_PATH=/path/to/python3/` to `galba.pl`.

#### Python modules

BUSCO2Tree requires the following Python modules to be installed:

-   argparse
-   os
-   time
-   Biopython
-   subprocess

#### Biopython
If Biopython is installed, GALBA can generate FASTA-files with coding sequences and protein sequences predicted by AUGUSTUS and generate track data hubs for visualization of a GALBA run with MakeHub <sup name="a8">[R8](#f8)</sup>.
These are optional steps. The first can be disabled with the command-line flag `--skipGetAnnoFromFasta`, the second can be activated by using the command-line options `--makehub --email=your@mail.de`, Biopython is not required if neither of these optional steps shall be performed.

On Ubuntu, install Python3 package manager with:

    `sudo apt-get install python3-pip`

Then, install Biopython with:

    `sudo pip3 install biopython`

#### TrimAl

#### IQTree

Download AUGUSTUS from its master branch at <https://github.com/Gaius-Augustus/Augustus>. Unpack AUGUSTUS and install AUGUSTUS according to AUGUSTUS `README.TXT`. ***Do not use outdated AUGUSTUS versions from other sources, e.g. Debian package or the Bioconda package! GALBA highly depends in particular on an up-to-date Augustus/scripts directory, and other sources are often lagging behind.***

You should compile AUGUSTUS on your own system in order to avoid problems with versions of libraries used by AUGUSTUS. Compilation instructions are provided in the AUGUSTUS `README.TXT` file (`Augustus/README.txt`).

In order to make the variable available to all Bash sessions, add the above line to a startup script, e.g., `~/.bashrc`.

##### Important:

##### Modification of $PATH

Adding directories of binaries and scripts to your `$PATH` variable enables your system to locate these tools,
automatically. It is a requirement for running BUSCO2Tree to do this, because BUSCO2Tree will try to find them from the location of the `$PATH` variable we recommend to add them to your `$PATH` variable. 

a) For your current bash session, type:

```
PATH="/home/nmoreyra/Soft/miniconda3/bin:/home/nmoreyra/Soft/bin:$PATH"
export PATH
```

b) For all your BASH sessions, add the above lines to a startup script (e.g.`~/.bashrc`).

### Optional tools

#### Mafft

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

Subsequently, install other software "as usual" in your conda environment.

-------
BUSCO2Tree components
------------------------------------

BUSCO2Tree is a collection of Python scripts. The main script that will be called in order to run the pipeline is `BUSCO2Tree.py`. Additional Python components are:

-   [find_singlecopy_BUSCOs.py]()

-   [align_BUSCOs.py]()

-   [create_matrix.py]()

-   [phylogenetic_analysis.py]()

All scripts (files ending with `*.py`) that are part of this pipeline must be executable in order to run BUSCO2Tree. This should already be the case if you download BUSCO2Tree from GitHub.

It is important that the `x` in `-rwxr-xr-x` is present for each script. If that is not the case, run

    `chmod a+x *.pl *.py`

in order to change file attributes.

You may find it helpful to add the directory in which GALBA perl scripts reside to your `$PATH` environment variable. For a single bash session, enter:

```
    PATH=/your_path_to_galba/:$PATH
    export PATH
```

To make this `$PATH` modification available to all bash sessions, add the above lines to a startup script (e.g.`~/.bashrc`).

---

Example data
========
[BUSCO groups for the eukaryota lineage obtained from ODB v10](https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2020-09-10.tar.gz)

