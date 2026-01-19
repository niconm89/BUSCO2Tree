# BUSCO2Tree

## Author

**Nicolás Nahuel Moreyra, PhD**  
Instituto de Ecología, Genética y Evolución de Buenos Aires (IEGEBA), CONICET  
Departamento de Ecología, Genética y Evolución (EGE),  
Facultad de Ciencias Exactas y Naturales, Universidad de Buenos Aires (UBA)

Contact  
- Email: niconm89@gmail.com  
- GitHub: https://github.com/niconm89  
- LinkedIn: https://www.linkedin.com/in/nicolasmoreyra/

---

## Overview

**BUSCO2Tree** is a modular phylogenomics pipeline designed to automate downstream analyses starting from BUSCO results obtained from multiple genomes, proteomes, or transcriptomes.

From an evolutionary perspective, the pipeline focuses on the use of **shared single-copy orthologs** as a robust genomic substrate for phylogenetic inference. By restricting analyses to conserved genes that are present in single copy across all taxa, BUSCO2Tree minimizes the confounding effects of gene duplication, loss, and paralogy, which are common sources of noise in genome-scale phylogenetic analyses.

The pipeline identifies **shared single-copy BUSCO orthologs**, aligns them, concatenates the resulting alignments into a partitioned phylogenetic matrix, and infers a species tree. This strategy enables the reconstruction of evolutionary relationships across species while retaining locus-level information that can later be used to explore gene tree heterogeneity, model fit, and patterns of evolutionary constraint.

Each step can be executed independently or as part of an end-to-end workflow, allowing users to adapt the pipeline to different comparative genomics and evolutionary biology questions.

---

## Pipeline steps

BUSCO2Tree consists of one orchestrator script and four main steps:

1. **Find shared single-copy BUSCOs**  
   Extracts BUSCO groups present as single copy across all datasets.

2. **Align BUSCO orthologs**  
   Performs multiple sequence alignments using MAFFT, optionally trimming alignments with trimAl.

3. **Build phylogenetic matrix**  
   Concatenates alignments into a supermatrix and generates partition files.

4. **Infer phylogenetic tree**  
   Runs IQ-TREE to estimate the species tree, substitution models, bootstrap support, and concordance factors.

Each step can be executed independently or in combination.

---

## Repository structure

```
BUSCO2Tree/
├── BUSCO2Tree.py
├── steps/
│   ├── find_singlecopy.py
│   ├── align_BUSCOs.py
│   ├── create_matrix.py
│   └── phylogenetic_analysis.py
|── docker/
│   └── Dockerfile/
├── example/
│   ├── Drosophila_proteins_tiny/
│   └── config/
│       └── config_mafft.yaml
├── nextflow/
│   ├── main.nf
│   ├── nextflow.config
│   └── README_nextlow.md
|── scripts/
│   └── run_docker.sh
├── tests/
│   └── smoke_test.sh
├── environment.yml
├── CITATION.cff
└── README.md
```

---

## Installation

BUSCO2Tree can be executed in three supported ways:

1. Using a local installation of the required dependencies (standalone Python execution).
2. Using Docker for fully reproducible, containerized runs ([Docker execution](#docker-execution)).
3. Using Nextflow on top of the Docker image for scalable and fully reproducible workflows ([Nextflow execution](#nextflow-execution)).

For most users, Docker or Nextflow execution is recommended.
See the corresponding sections below for details.

BUSCO2Tree can be executed either using a local installation of the dependencies
or via Docker or Nextflow for fully reproducible runs ).

### Clone repository

```bash
git clone https://github.com/niconm89/BUSCO2Tree.git
cd BUSCO2Tree
```

### Create Conda environment

```bash
conda env create -f environment.yml
conda activate busco2tree
```

### Verify installation

```bash
python BUSCO2Tree.py -h
```

---

## Dependencies

BUSCO2Tree has been tested with:

- Python ≥ 3.9  
- BUSCO v5 (results only; BUSCO itself is not run by this pipeline)  
- MAFFT  
- trimAl  
- IQ-TREE2  

All dependencies are listed in `environment.yml` and must be available in the system `$PATH`.

---

## Configuration files

### MAFFT configuration (YAML)

BUSCO2Tree supports optional YAML configuration files for alignment settings.

Example: `example/config/config_mafft.yaml`

```yaml
align:
  mafft_path: mafft
  align_method: auto
  threads: 4
```

Notes:
- All parameters are optional.
- If not provided, default values are used.
- When running via the orchestrator, the `--threads` CLI argument overrides YAML settings.

---

## Example data

The example data provided in this repository were obtained from a study carried out during my PhD, aimed at identifying genomic innovations associated with cactophily and host plant specialization in *Drosophila*.

The study is:

*Phylogenomics provides insights into the evolution of cactophily and host plant shifts in* **Drosophila**  
https://doi.org/10.1016/j.ympev.2022.107653

For this work, I employed the genomes of the thirteen *Drosophila* species included in the study and ran **BUSCO v5.4.6** on each dataset to identify conserved orthologs. BUSCO analyses were performed using the **eukaryota lineage dataset from OrthoDB v10 (eukaryota_odb10)**.

The directory `example/Drosophila_proteins_tiny/` contains a reduced, curated subset of BUSCO outputs derived from these analyses and is intended exclusively for testing and demonstration purposes. This dataset is used by the reproducible smoke test provided with the pipeline.

---

## Reproducible smoke test

A full end-to-end reproducible test is provided using the example BUSCO dataset.

### Run smoke test

```bash
conda activate busco2tree
./tests/smoke_test.sh
```

This test:
- runs Steps 1–3 unconditionally
- runs Step 4 if `iqtree2` is available
- validates expected outputs
- writes results to `.smoke_out/` (ignored by git)

---

## Running BUSCO2Tree

BUSCO2Tree expects BUSCO v5 results as input and runs downstream steps including
BUSCO filtering, multiple sequence alignment, matrix construction, and phylogenetic
inference.

### Step 1: Find shared single-copy BUSCOs

```bash
BUSCO2Tree.py   -s 1   -b <BUSCO_results_directory>   -l eukaryota   -d 10   -o B2T_output
```

### Step 2: Align BUSCO orthologs

```bash
BUSCO2Tree.py   -s 2   -f B2T_output/01_single-copy/common_busco_sequences   -c example/config/config_mafft.yaml   -o B2T_output
```

Optional trimming:

```bash
--trim
```

or custom trimAl parameters:

```bash
--trimparams "-gt 0.3"
```

### Step 3: Build phylogenetic matrix

```bash
BUSCO2Tree.py   -s 3   -a B2T_output/02_alignments/00_raw_aln   -o B2T_output   --format phylip
```

### Step 4: Infer phylogenetic tree

```bash
BUSCO2Tree.py   -s 4   -m B2T_output/03_matrix/matrix.phylip   -p B2T_output/03_matrix/partitions.nex   -o B2T_output   -B 1000
```

### Full pipeline (Steps 1–4)

```bash
BUSCO2Tree.py   -s 1 2 3 4   -b <BUSCO_results_directory>   -c example/config/config_mafft.yaml   --trim   -o B2T_output   -t 8
```

---

### Docker execution

BUSCO2Tree can be executed using Docker, which provides a fully reproducible
environment including all required dependencies (MAFFT, trimAl, IQ-TREE).

#### Build the Docker image

From the root of the repository:

```bash
docker build -t busco2tree:latest -f docker/Dockerfile .
```

#### Run BUSCO2Tree using Docker

A helper script is provided to simplify execution:

```bash
chmod +x scripts/run_docker.sh
```

Example run:

```bash
./scripts/run_docker.sh \
  -b /path/to/BUSCO_results_directory \
  -t 8
```

Notes:

- The BUSCO input directory is mounted read-only inside the container.
- All outputs are written to the current working directory.
- The Docker image uses IQ-TREE v3 for phylogenetic inference.

---

### Nextflow execution

BUSCO2Tree includes an initial Nextflow workflow (DSL2)  available under `nextflow/` that orchestrates steps 1-4 using the BUSCO2Tree Docker image.
See: `nextflow/README_nextflow.md`.

Example Run:

```bash
nextflow run nextflow/main.nf -profile docker \
  --buscodir /path/to/BUSCO_results_dir \
  --outdir B2T_output_nf \
  --steps 1,2,3,4
```

Optional flags:
```bash
--config <config_mafft.yaml>: MAFFT YAML config for Step2
--trim / --trimparams "<trimal params>": enable trimming
--format {phylip,nexus}: matrix format (Step3)
--seqtype {AA,DNA}: sequence type
--prefix, --bootstrap, --threads: IQ-TREE options
```

Notes:

- The Docker image must not define an ENTRYPOINT (Nextflow needs to execute arbitrary commands).

- Step3 writes 03_matrix/metadata.yaml, which Step4 reads to locate the matrix and partitions.
---

## Output structure

```
B2T_output/
├── 01_single-copy
│   ├── common_busco_sequences/
│   ├── genomes_names.txt
│   └── list_common_busco.txt
├── 02_alignments
│   ├── 00_raw_aln/
│   └── 01_trim_aln/
├── 03_matrix
│   ├── matrix.phylip
│   └── partitions.nex
└── 04_phylogenetic_tree
    ├── *.treefile
    ├── *.iqtree
    └── *.log
```

---

## License

This project is distributed under the **GNU General Public License v3.0**.  
See the `LICENSE` file for details.
