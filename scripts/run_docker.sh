#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# BUSCO2Tree Docker runner
# -----------------------------
# Usage examples:
#   ./run_docker.sh -b /home/user/BUSCO_outputs/Drosophila -o B2T_output -t 8
#   ./run_docker.sh -b /data/BUSCO -s "1 2 3 4" -l eukaryota -d 10 -fmt phylip -B 1000
#
# Notes:
# - The BUSCO dir is mounted read-only to prevent accidental modification.
# - Outputs are written to the current directory (mounted as /work).
# -----------------------------

IMAGE="${IMAGE:-busco2tree:latest}"

# Defaults (match your typical run)
STEPS="1 2 3 4"
LINEAGE="eukaryota"
ODB="10"
THREADS="8"
FORMAT="phylip"
BOOTSTRAP="1000"
SEQTYPE="AA"
OUTDIR="B2T_output"

BUSCODIR=""

die() { echo "ERROR: $*" >&2; exit 1; }

usage() {
  cat >&2 <<EOF
Usage: $0 -b <BUSCO_DIR_ON_HOST> [options]

Required:
  -b, --busco-dir   Path on HOST to BUSCO outputs directory (mounted read-only)

Optional:
  -s, --steps       Steps to run (default: "$STEPS")
  -o, --outdir      Output dir name inside /work (default: "$OUTDIR")
  -l, --lineage     BUSCO lineage (default: "$LINEAGE")
  -d, --odb         OrthoDB version (default: "$ODB")
  -t, --threads     Threads (default: "$THREADS")
  -f, --format      Matrix format: phylip|nexus (default: "$FORMAT")
  -B, --bootstrap   Bootstrap replicates (default: "$BOOTSTRAP")
  --seqtype         AA|DNA (default: "$SEQTYPE")
  --image           Docker image (default: $IMAGE)

Examples:
  $0 -b /home/user/BUSCO_outputs/Drosophila -t 8
  $0 -b /home/user/BUSCO_outputs/Drosophila -s "4" -o B2T_output -B 1000
EOF
}

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -b|--busco-dir) BUSCODIR="$2"; shift 2 ;;
    -s|--steps) STEPS="$2"; shift 2 ;;
    -o|--outdir) OUTDIR="$2"; shift 2 ;;
    -l|--lineage) LINEAGE="$2"; shift 2 ;;
    -d|--odb) ODB="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    -f|--format) FORMAT="$2"; shift 2 ;;
    -B|--bootstrap) BOOTSTRAP="$2"; shift 2 ;;
    --seqtype) SEQTYPE="$2"; shift 2 ;;
    --image) IMAGE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1 (use -h for help)" ;;
  esac
done

[[ -n "$BUSCODIR" ]] || { usage; die "Missing -b/--busco-dir"; }
[[ -d "$BUSCODIR" ]] || die "BUSCO dir not found on host: $BUSCODIR"

# Run
echo "Running BUSCO2Tree in Docker"
echo "  Image:      $IMAGE"
echo "  Host BUSCO: $BUSCODIR (mounted read-only)"
echo "  Workdir:    $PWD (mounted as /work)"
echo "  Steps:      $STEPS"
echo "  Outdir:     $OUTDIR"
echo "  Lineage:    $LINEAGE"
echo "  ODB:        $ODB"
echo "  Threads:    $THREADS"
echo "  Format:     $FORMAT"
echo "  Bootstrap:  $BOOTSTRAP"
echo "  Seqtype:    $SEQTYPE"
echo

docker run --rm \
  -v "$BUSCODIR":/data:ro \
  -v "$PWD":/work \
  -w /work \
  "$IMAGE" \
  -s $STEPS \
  -b /data \
  -o "$OUTDIR" \
  -l "$LINEAGE" \
  -d "$ODB" \
  -t "$THREADS" \
  -fmt "$FORMAT" \
  -B "$BOOTSTRAP" \
  -st "$SEQTYPE"
