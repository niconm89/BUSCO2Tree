# BUSCO2Tree
**Reproducible phylogenomics from BUSCO outputs.** From orthology inference to partitioned species trees, at scale.

---

## Overview
BUSCO2Tree is a modular phylogenomics pipeline that converts **BUSCO v5 results** into a **partitioned supermatrix** and a **species tree**, with optional gene trees and concordance metrics.  
It is designed for **reproducible, scalable, and production-ready comparative genomics** workflows.

The pipeline focuses on **shared single-copy orthologs**, providing a high-confidence signal for phylogenetic inference while minimizing noise from paralogy, duplication, and gene loss.

---

## Why BUSCO2Tree
BUSCO2Tree was designed with a focus on reproducibility, scalability, and methodological clarity across heterogeneous computing environments:

- **Reproducibility first**  
  Identical inputs yield identical outputs via containerized execution.

- **Production-oriented architecture**  
  Clear step boundaries, deterministic outputs, and explicit I/O contracts.

- **Scalable by design**  
  Runs locally, in Docker, or orchestrated with Nextflow on servers or HPC.

- **Scientifically grounded**  
  Built around best practices in orthology-based phylogenomics.

---

## Typical use cases
BUSCO2Tree supports workflows such as:

- Orthology-driven phylogenomic inference
- Assembly and annotation quality control
- Comparative genomics feature discovery
- Partition-aware model selection
- Gene-tree discordance and concordance analyses
- Evolutionary benchmarking of new genomes or proteomes

---

## Pipeline summary
Converts BUSCO v5 results into a partitioned supermatrix and a species tree based on shared single-copy orthologs.

**Input**  
A directory containing multiple BUSCO v5 result folders, one per dataset (genome, proteome, or transcriptome). Directory names are treated as sample identifiers.

**Steps**
1. Identification of shared single-copy BUSCOs  
2. Multiple sequence alignment (MAFFT) with optional trimming (trimAl)  
3. Concatenation into a supermatrix with partitions  
4. Species tree inference using IQ-TREE2  
   Optional gene trees and concordance factors

**Supported sequence types**  
- Amino acids (AA)  
- Nucleotides (DNA)

---

## Technology stack
Core tools:
- Python (pipeline logic)
- MAFFT (multiple sequence alignment)
- trimAl (alignment trimming, optional)
- IQ-TREE2 (phylogenetic inference)

Execution modes:
- Local execution (Conda)
- Docker (fully reproducible runtime)
- Nextflow DSL2 (workflow orchestration)

---

## Repository structure
```
.
├── BUSCO2Tree.py
├── steps/
│   ├── find_singlecopy.py
│   ├── align_BUSCOs.py
│   ├── create_matrix.py
│   └── phylogenetic_analysis.py
├── docker/
│   └── Dockerfile
├── nextflow/
│   ├── main.nf
│   ├── nextflow.config
│   └── README_nextflow.md
├── example/
├── tests/
├── environment.yml
├── CITATION.cff
└── README.md
```

---

## Installation

### Cloning the repo:
```bash
git clone https://github.com/niconm89/BUSCO2Tree.git
cd BUSCO2Tree
```

### Option 1: Local installation (Conda)
```bash
conda env create -f environment.yml
conda activate busco2tree

python BUSCO2Tree.py -h
```

### Option 2: Docker (recommended)
```bash
docker build -t busco2tree:latest -f docker/Dockerfile .
```

Run via helper script:
```bash
chmod +x scripts/run_docker.sh

./scripts/run_docker.sh \
  -b /path/to/BUSCO_results_directory \
  -t 8
```

Notes:
- BUSCO inputs are typically mounted read-only
- Outputs are written to the working directory
- All dependencies are bundled in the image

### Option 3: Nextflow (scalable execution)
```bash
nextflow run nextflow/main.nf -profile docker \
  --buscodir /path/to/BUSCO_results_directory \
  --outdir B2T_output_nf \
  --steps 1,2,3,4
```

Important:
The Docker image intentionally avoids a fixed ENTRYPOINT to allow Nextflow to execute arbitrary commands inside the container.

---

## Usage (Python CLI)

### Step 1: Find shared single-copy BUSCOs
```bash
BUSCO2Tree.py -s 1 \
  -b <BUSCO_results_directory> \
  -l eukaryota \
  -d 10 \
  -o B2T_output
```

### Step 2: Align orthologs
```bash
BUSCO2Tree.py -s 2 \
  -f B2T_output/01_single-copy/common_busco_sequences \
  -c example/config/config_mafft.yaml \
  -o B2T_output
```

Optional trimming:
```bash
--trim --trimparams "-gt 0.3"
```

### Step 3: Build supermatrix and partitions
```bash
BUSCO2Tree.py -s 3 \
  -a B2T_output/02_alignments/00_raw_aln \
  --format phylip \
  -o B2T_output
```

### Step 4: Species tree inference
```bash
BUSCO2Tree.py -s 4 \
  -m B2T_output/03_matrix/matrix.phylip \
  -p B2T_output/03_matrix/partitions.nex \
  -B 1000 \
  -o B2T_output
```

### End-to-end execution
```bash
BUSCO2Tree.py -s 1 2 3 4 \
  -b <BUSCO_results_directory> \
  -c example/config/config_mafft.yaml \
  --trim \
  -o B2T_output \
  -t 8
```

---

## Outputs
```
B2T_output/
├── 01_single-copy/
│   ├── common_busco_sequences/
│   ├── genomes_names.txt
│   └── list_common_busco.txt
├── 02_alignments/
│   ├── 00_raw_aln/
│   └── 01_trim_aln/
├── 03_matrix/
│   ├── matrix.phylip
│   └── partitions.nex
└── 04_phylogenetic_tree/
    ├── *.treefile
    ├── *.iqtree
    └── *.log
```

---

## Design principles
- Deterministic execution and explicit outputs
- Clear separation between orchestration and scientific logic
- Container-first reproducibility
- Transparent, inspectable intermediate results
- Minimal assumptions about infrastructure

---

## Example data
The repository includes reduced BUSCO example datasets derived from a published phylogenomics study, intended for testing and demonstration purposes.

---

## Roadmap
- Versioned Docker image releases
- Nextflow profiles for common HPC schedulers
- CI-based smoke tests (Docker + Nextflow)
- Input/output schema validation
- Runtime and memory benchmarking

---

## Citation
If you use BUSCO2Tree in academic work, please cite this repository and the underlying tools (BUSCO, MAFFT, trimAl, IQ-TREE2).  
A `CITATION.cff` file is included.

---

## License
MIT License.

---

## Author
**Nicolás Nahuel Moreyra, PhD**  
Bioinformatician | Comparative Genomics | Omics Data Integration  

GitHub: https://github.com/niconm89  
Email: niconm89@gmail.com  
LinkedIn: https://www.linkedin.com/in/nicolasmoreyra/
