# BUSCO2Tree User Guide

**Author of BUSCO2Tree:**

Nicolás Nahuel Moreyra
* Instituto de Ecología, Genética y Evolución de Buenos Aires (IEGEBA), CONICET, Argentina
* Departamento de Ecología, Genética y Evolución (EGE), Facultad de Ciencias Exactas y Naturales (FCEN), Universidad de Buenos Aires (UBA), Argentina
+ E-mail: niconm89@gmail.com |  nmoreyra@ege.fcen.uba.ar
+ Twitter: [@nicomoreyra](https://twitter.com/NicoMoreyra)

---
**About BUSCO2Tree**

This pipeline was developed to automatize the downstream phylogenomic analyses to obtain a species tree after running BUSCO. 
We provide data example to test funcionallity. Data was obtained from a study we have recently published [1]. Briefly, thirteen Drosophila genomes were assessed using BUSCO v5.4.2 to obtain completeness scores and Eukaryote single-copy orthologs (BUSCO groups). Common BUSCOs among all genomes were the input of this pipeline.

---

Contents
========

-   [Authors](#author-of-BUSCO2tree)
-   [What is BUSCO2tree?](#what-is-BUSCO2tree)
-   [Overview of modes for running BUSCO2tree](#overview-of-modes-for-running-BUSCO2tree)
-   [Installation](#installation)
    -   [Supported software versions](#supported-software-versions)
    -   [BUSCO2tree](#BUSCO2tree)
        -   [Perl pipeline dependencies](#perl-pipeline-dependencies)
        -   [BRAKER components](#braker-components)
        -   [Bioinformatics software dependencies](#bioinformatics-software-dependencies)
            -   [Mandatory tools](#mandatory-tools)
            -   [Optional tools](#optional-tools)
-   [Running BUSCO2tree](#running-BUSCO2tree)
    -   [BUSCO2tree pipeline modes](#different-BUSCO2tree-pipeline-modes)
        -   [BUSCO2tree with RNA-Seq data](#BUSCO2tree-with-rna-seq-data)
        -   [BUSCO2tree with proteins of any evolutionary distance](#BUSCO2tree-with-proteins-of-any-evolutionary-distance)
        -   [BUSCO2tree with proteins of short evolutionary distance](#BUSCO2tree-with-proteins-of-short-evolutionary-distance)
    -   [Description of selected BUSCO2tree command line options](#description-of-selected-BUSCO2tree-command-line-options)
        -   [--ab_initio](#--ab_initio)
        -   [--augustus_args=--some\_arg=bla](#--augustus_args--some_argbla)
        -   [--cores=INT](#--coresint)
        -   [--fungus](#--fungus)
        -   [--softmasking](#--softmasking)
        -   [--useexisting](#--useexisting)
        -   [--crf](#--crf)
        -   [--lambda=int](#--lambdaint)
        -   [--UTR=on](#--utron)
        -   [--addUTR=on](#--addutr)
        -   [--stranded=+,-,.,...](#--stranded-)
	    -   [--makehub --email=your@mail.de](#--makehub---emailyourmailde)
-   [Output of BUSCO2tree](#output-of-BUSCO2tree)
-   [Example data](#example-data)
    -   [Data description](#data-description)
    -   [Testing BUSCO2tree with genomes](#testing-braker-with-genomes)
    -   [Testing BUSCO2tree with proteins](#testing-braker-with-proteins)
