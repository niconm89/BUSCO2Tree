#!/usr/bin/env bash
set -euo pipefail

BUSCODIR="${1:-}"
OUTDIR="${2:-B2T_output_nf}"

if [[ -z "${BUSCODIR}" ]]; then
  echo "Usage: $0 <BUSCO_results_dir> [outdir]"
  exit 1
fi

nextflow run nextflow/main.nf -profile docker \
  --steps 1,2,3,4 \
  --buscodir "${BUSCODIR}" \
  --outdir "${OUTDIR}"


