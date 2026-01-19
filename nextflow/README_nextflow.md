# BUSCO2Tree Nextflow workflow (MVP)

This directory contains an initial Nextflow (DSL2) workflow for BUSCO2Tree.  
The goal is a minimal but robust implementation that runs Steps 1-4 using the BUSCO2Tree Docker image.

This workflow is intentionally simple:
1. It executes the step scripts located at `/opt/busco2tree/steps/` inside the container.
2. It avoids running the monolithic `BUSCO2Tree.py` orchestrator from Nextflow.
3. It produces the same directory structure as the classic run, but under an output folder chosen by the user.
4. It uses a metadata file produced by Step 3 to drive Step 4 without heuristic file searching.

## Requirements

1. Nextflow (DSL2).
2. Docker.
3. The BUSCO2Tree Docker image built locally (or available by tag).
4. BUSCO v5 results directory as input.

## Build the Docker image

From the repository root:

```bash
docker build -t busco2tree:latest -f docker/Dockerfile .
```

Important note: the Docker image must NOT define a fixed ENTRYPOINT. Nextflow needs to execute arbitrary commands inside the container.

## Quick sanity checks (recommended)

```bash
docker run --rm busco2tree:latest python /opt/busco2tree/steps/find_singlecopy.py -h | head
docker run --rm busco2tree:latest python /opt/busco2tree/steps/align_BUSCOs.py -h | head
docker run --rm busco2tree:latest python /opt/busco2tree/steps/create_matrix.py -h | head
docker run --rm busco2tree:latest python /opt/busco2tree/steps/phylogenetic_analysis.py -h | head
```

## Running the workflow

All commands below must be run from the repository root, or provide explicit paths.

### Run all steps (1–4)

```bash
nextflow run nextflow/main.nf -profile docker \
  --steps 1,2,3,4 \
  --buscodir /path/to/BUSCO_results_dir \
  --outdir B2T_output_nf
```

### Run only Step 1

```bash
nextflow run nextflow/main.nf -profile docker \
  --steps 1 \
  --buscodir /path/to/BUSCO_results_dir \
  --outdir B2T_output_nf
```

### Run Step 1 and Step 2

```bash
nextflow run nextflow/main.nf -profile docker \
  --steps 1,2 \
  --buscodir /path/to/BUSCO_results_dir \
  --outdir B2T_output_nf
```

### Step 2 with MAFFT YAML config

```bash
nextflow run nextflow/main.nf -profile docker \
  --steps 1,2 \
  --buscodir /path/to/BUSCO_results_dir \
  --outdir B2T_output_nf \
  --config /path/to/config_mafft.yaml
```

### Step 2 trimming

```bash
nextflow run nextflow/main.nf -profile docker \
  --steps 1,2 \
  --buscodir /path/to/BUSCO_results_dir \
  --outdir B2T_output_nf \
  --trim \
  --trimparams "-gt 0.3"
```

### Full run with trimming and MAFFT config

```bash
nextflow run nextflow/main.nf -profile docker \
  --steps 1,2,3,4 \
  --buscodir /path/to/BUSCO_results_dir \
  --outdir B2T_output_nf \
  --config /path/to/config_mafft.yaml \
  --trim \
  --trimparams "-gt 0.3"
```

## Parameters

This workflow uses these key parameters.

### Input/output
1. `--buscodir`  
   Directory containing multiple BUSCO outputs (one per genome/proteome).

2. `--outdir`  
   Output directory where Nextflow publishes step results. Default is `B2T_output_nf`.

3. `--steps`  
   Comma-separated list of steps to run (example: `1,2,3,4`).

### Step 1
1. `--lineage`  
   BUSCO lineage label used for naming/reporting.
2. `--odb`  
   BUSCO odb version number.
3. `--seqtype`  
   `AA` or `DNA`.

### Step 2 (alignments)
1. `--config`  
   YAML config file for MAFFT.
2. `--mafft_cmd`  
   Manual MAFFT command options string passed to the step script.
3. `--trim`  
   Enable trimming.
4. `--trimparams`  
   trimAl parameters, for example `"-gt 0.3"`.

### Step 3 (matrix)
1. `--format`  
   `phylip` or `nexus`.
2. `--seqtype`  
   `AA` or `DNA` (mapped to Step3 datatype internally).

### Step 4 (IQ-TREE)
1. `--prefix`  
   Output prefix for IQ-TREE results.
2. `--bootstrap`  
   Number of bootstrap replicates.
3. `--threads`  
   Threads for IQ-TREE.
4. `--genetrees`  
   Enable per-locus gene trees.
5. `--concordance`  
   Compute concordance factors (requires `--genetrees`).

## Output structure

Under `--outdir`:

1. `01_single-copy/`  
   Common single-copy BUSCO FASTAs and lists.

2. `02_alignments/`  
   Subdirectories:
   - `00_raw_aln/` for raw MAFFT alignments
   - `01_trim_aln/` for trimmed alignments when trimming is enabled

3. `03_matrix/`  
   Concatenated matrix and partitions, plus:
   - `metadata.yaml`

4. `04_phylogenetic_tree/`  
   IQ-TREE outputs.

## Step 3 metadata contract

Step 3 writes `03_matrix/metadata.yaml`. Step 4 reads it to locate the correct matrix file, partitions file, alignment subdir, and sequence type.

Example:

```yaml
matrix_file: matrix.phylip
partitions_file: partitions.nex
format: phylip
datatype: PROTEIN
aln_subdir: 00_raw_aln
seqtype: AA
```

## Sanity checks (expected outputs)

After a successful run (`--steps 1,2,3,4`), the output directory must contain:

- `01_single-copy/` with:
  - `common_busco_sequences/`
  - `list_common_busco.txt`
  - `genomes_names.txt`

- `02_alignments/` with:
  - `00_raw_aln/` (always)
  - `01_trim_aln/` only when `--trim` or `--trimparams` are provided

- `03_matrix/` with:
  - `matrix.<format>` (e.g. `matrix.phylip` or `matrix.nex`)
  - `partitions.nex`
  - `metadata.yaml`

- `04_phylogenetic_tree/` with IQ-TREE outputs (e.g. `*.treefile`, `*.log`).

Quick commands:

```bash
ls -lah B2T_output_nf
ls -lah B2T_output_nf/02_alignments/00_raw_aln | head
ls -lah B2T_output_nf/03_matrix
ls -lah B2T_output_nf/04_phylogenetic_tree | head
