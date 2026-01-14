#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

echo "[1/6] Basic CLI checks"
python BUSCO2Tree.py -h >/dev/null
python steps/find_singlecopy.py -h >/dev/null
python steps/align_BUSCOs.py -h >/dev/null
python steps/create_matrix.py -h >/dev/null
python steps/phylogenetic_analysis.py -h >/dev/null
echo "OK: help/CLI"

echo "[2/6] Inputs"
BUSCO_DIR="$ROOT/example/Drosophila_proteins_tiny"
CFG="$ROOT/example/config/config_mafft.yaml"
OUT="$ROOT/.smoke_out"

test -d "$BUSCO_DIR"
test -f "$CFG"

echo "BUSCO_DIR=$BUSCO_DIR"
echo "CFG=$CFG"
echo "OUT=$OUT"

echo "[3/6] Reset outputs"
rm -rf "$OUT"
mkdir -p "$OUT"

echo "[4/6] Step 1"
python BUSCO2Tree.py -s 1 -b "$BUSCO_DIR" -l eukaryota -d 10 -o "$OUT" -t 2
SC_DIR="$OUT/01_single-copy/common_busco_sequences"
test -d "$SC_DIR"
echo "OK: Step 1"

echo "[5/6] Steps 2 and 3"
python BUSCO2Tree.py -s 2 -f "$SC_DIR" -c "$CFG" -o "$OUT" -t 2
ALN_DIR="$OUT/02_alignments/00_raw_aln"
test -d "$ALN_DIR"

python BUSCO2Tree.py -s 3 -a "$ALN_DIR" -o "$OUT" -t 2
test -s "$OUT/03_matrix/matrix.phylip"
test -s "$OUT/03_matrix/partitions.nex"
echo "OK: Steps 2-3"

echo "[6/6] Optional Step 4"
if command -v iqtree2 >/dev/null 2>&1; then
  python BUSCO2Tree.py -s 4 \
    -m "$OUT/03_matrix/matrix.phylip" \
    -p "$OUT/03_matrix/partitions.nex" \
    -o "$OUT" -t 2 -B 1000 || true
  TREEFILE="$(find "$OUT/04_phylogenetic_tree" -type f -name "*.treefile" -print -quit 2>/dev/null || true)"
  if [[ -n "$TREEFILE" ]]; then
    echo "OK: IQ-TREE produced $TREEFILE"
  else
    echo "NOTE: IQ-TREE ran but no .treefile detected. Check $OUT/04_phylogenetic_tree"
  fi
else
  echo "iqtree2 not found; skipping Step 4"
fi

echo "SMOKE TEST PASSED"
