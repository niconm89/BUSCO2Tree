# BUSCO2Tree User Guide

**This GIT README is in being created!!!!**

**Author of BUSCO2Tree**
[Nicolás Nahuel Moreyra](https://github.com/niconm89/Curriculum_Vitae)

* Instituto de Ecología, Genética y Evolución de Buenos Aires (IEGEBA), CONICET, Argentina
* Departamento de Ecología, Genética y Evolución (EGE), Facultad de Ciencias Exactas y Naturales (FCEN), Universidad de Buenos Aires (UBA), Argentina
+ E-mail: niconm89@gmail.com |  nmoreyra@ege.fcen.uba.ar
+ Twitter: [@nicomoreyra](https://twitter.com/NicoMoreyra)

---
**About BUSCO2Tree**

This pipeline was developed to automatize the downstream phylogenomic analyses to obtain a species tree after running BUSCO over several genomes/proteomes/transcriptomes.
We provide one data example to test funcionallity. Data was obtained from a study we have recently published [1]. Briefly, thirteen Drosophila genomes were assessed using BUSCO v5.4.2 to obtain completeness scores and dipteran single-copy orthologs. Common BUSCOs among all genomes were the input of this pipeline.

---

Contents
========

-   [Authors](#author-of-BUSCO2tree)
-   [What is BUSCO2tree?](#what-is-BUSCO2tree)
-   [Overview of modes for running BUSCO2tree](#overview-of-modes-for-running-BUSCO2tree)
-   [Installation](#installation)
    -   [Supported software versions](#supported-software-versions)
    -   [BUSCO2tree](#BUSCO2tree)
        -   [Dependencies](#dependencies)
            -   [Mandatory tools](#mandatory-tools)
            -   [Optional tools](#optional-tools)
-   [Running BUSCO2tree](#running-BUSCO2tree)
    -   [BUSCO2tree pipeline modes](#different-BUSCO2tree-pipeline-modes)
    -   [Description of selected BUSCO2tree command line options](#description-of-selected-BUSCO2tree-command-line-options)
-   [Output of BUSCO2tree](#output-of-BUSCO2tree)
-   [Example data](#example-data)
    -   [Data description](#data-description)
    -   [Testing BUSCO2tree with genomes](#testing-braker-with-genomes)
    -   [Testing BUSCO2tree with proteins](#testing-braker-with-proteins)
