nextflow.enable.dsl=2

params.buscodir   = params.buscodir
params.outdir     = params.outdir ?: 'B2T_output_nf'

params.lineage    = params.lineage ?: 'eukaryota'
params.odb        = params.odb ?: 10
params.seqtype    = params.seqtype ?: 'AA'

params.config     = params.config ?: null      // YAML para Step2 (opcional)
params.mafft_cmd  = params.mafft_cmd ?: null   // string con opciones MAFFT (opcional)
params.trim       = params.trim ?: false
params.trimparams = params.trimparams ?: null

params.container  = params.container ?: 'busco2tree:latest'


process STEP1_FIND_SINGLECOPY {
  tag "step1"

  container params.container
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
    path buscodir

  output:
    path "01_single-copy/common_busco_sequences", emit: fastadir

  script:
    """
    set -euo pipefail
    set -x

    mkdir -p 01_single-copy

    python /opt/busco2tree/steps/find_singlecopy.py \
      -b "${buscodir}" \
      -o 01_single-copy \
      -l "${params.lineage}" \
      -d "${params.odb}" \
      --seqtype "${params.seqtype}"

    ls -lah 01_single-copy | head
    """
}

process STEP2_ALIGN_BUSCOS_NO_CONFIG {
  tag "step2:no_config"

  container params.container
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
    path fastadir

  output:
    path "02_alignments", emit: alndir

  script:
    def mcmd = params.mafft_cmd ? "-m \"${params.mafft_cmd}\"" : ""
    def doTrim = (params.trim instanceof Boolean) ? params.trim : params.trim.toString().toLowerCase() == 'true'
    def trimFlag = doTrim ? "-t" : ""
    def tparams = params.trimparams ? "-p \"${params.trimparams}\"" : ""

    """
    set -euo pipefail
    set -x

    python /opt/busco2tree/steps/align_BUSCOs.py \
      -f "${fastadir}" \
      -o 02_alignments \
      ${mcmd} \
      ${trimFlag} \
      ${tparams}

    ls -lah 02_alignments | head
    """
}

process STEP2_ALIGN_BUSCOS_WITH_CONFIG {
  tag "step2:with_config"

  container params.container
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
    path fastadir
    path config_file

  output:
    path "02_alignments", emit: alndir

  script:
    def mcmd = params.mafft_cmd ? "-m \"${params.mafft_cmd}\"" : ""
    def doTrim = (params.trim instanceof Boolean) ? params.trim : params.trim.toString().toLowerCase() == 'true'
    def trimFlag = doTrim ? "-t" : ""
    def tparams = params.trimparams ? "-p \"${params.trimparams}\"" : ""

    """
    set -euo pipefail
    set -x

    python /opt/busco2tree/steps/align_BUSCOs.py \
      -f "${fastadir}" \
      -o 02_alignments \
      -c "${config_file}" \
      ${mcmd} \
      ${trimFlag} \
      ${tparams}

    ls -lah 02_alignments | head
    """
}

process STEP3_CREATE_MATRIX {
  tag "step3"

  container params.container
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
    path alndir

  output:
    path "03_matrix", emit: matrixdir
    path "03_matrix/metadata.yaml", emit: meta

  script:
    def doTrim = (params.trim instanceof Boolean) ? params.trim : params.trim.toString().toLowerCase() == 'true'
    def aln_subdir = doTrim ? '01_trim_aln' : '00_raw_aln'
    // mapear params.seqtype (AA/DNA) a datatype (PROTEIN/DNA)
    def dtype = (params.seqtype == 'DNA') ? 'DNA' : 'PROTEIN'

    """
    set -euo pipefail
    set -x

    test -d "${alndir}/${aln_subdir}"
    n_aln=\$(find "${alndir}/${aln_subdir}" -maxdepth 1 -type f | wc -l)
    echo "Using: ${alndir}/${aln_subdir}  (files: \$n_aln)"
    test "\$n_aln" -gt 0

    mkdir -p 03_matrix

    python /opt/busco2tree/steps/create_matrix.py \
      -a "${alndir}/${aln_subdir}" \
      --datatype "${dtype}" \
      -o 03_matrix \
      -f "${params.format}"

    ls -lah 03_matrix | head

    # Detectar la matrix producida (convención: matrix.*)
    MATRIX_BASENAME=\$(find 03_matrix -maxdepth 1 -type f -name "matrix.*" -printf "%f\\n" | head -n 1)
    test -n "\$MATRIX_BASENAME"
    test -s "03_matrix/\$MATRIX_BASENAME"
    test -s "03_matrix/partitions.nex"

    cat > 03_matrix/metadata.yaml << EOF
matrix_file: \$MATRIX_BASENAME
partitions_file: partitions.nex
format: ${params.format}
datatype: ${dtype}
aln_subdir: ${aln_subdir}
seqtype: ${params.seqtype}
EOF
    """
}

process STEP4_PHYLOGENETIC_ANALYSIS {
  tag "step4"

  container params.container
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
    path matrixdir
    path meta
    path alndir_root

  output:
    path "04_phylogenetic_tree", emit: treedir

  script:
    def gtFlag = ((params.genetrees instanceof Boolean) ? params.genetrees : params.genetrees.toString().toLowerCase() == 'true') ? "-gt" : ""
    def cfFlag = ((params.concordance instanceof Boolean) ? params.concordance : params.concordance.toString().toLowerCase() == 'true') ? "-cf" : ""

    """
    set -euo pipefail
    set -x

    mkdir -p 04_phylogenetic_tree

    # Parse metadata -> __meta.tsv (tab-separated)
    python - "${meta}" << 'PY' > __meta.tsv
import yaml
from pathlib import Path
import sys

meta_path = Path(sys.argv[1])
d = yaml.safe_load(meta_path.read_text())

required = ["matrix_file", "partitions_file", "aln_subdir", "seqtype"]
missing = [k for k in required if k not in d]
if missing:
    raise SystemExit(f"metadata missing keys: {missing}")

print(f"{d['matrix_file']}\t{d['partitions_file']}\t{d['aln_subdir']}\t{d['seqtype']}")
PY

read -r MATRIX_BASENAME PART_BASENAME ALN_SUBDIR SEQTYPE < __meta.tsv
echo "Parsed: \$MATRIX_BASENAME \$PART_BASENAME \$ALN_SUBDIR \$SEQTYPE"

    MATRIX_FILE="${matrixdir}/\${MATRIX_BASENAME}"
    PART_FILE="${matrixdir}/\${PART_BASENAME}"
    ALN_DIR="${alndir_root}/\${ALN_SUBDIR}"

    test -s "\$MATRIX_FILE"
    test -s "\$PART_FILE"
    test -d "\$ALN_DIR"

    python /opt/busco2tree/steps/phylogenetic_analysis.py \
      -m "\$MATRIX_FILE" \
      -p "\$PART_FILE" \
      -s "\$SEQTYPE" \
      -o 04_phylogenetic_tree \
      -P "${params.prefix}" \
      -b "${params.bootstrap}" \
      -t "${params.threads}" \
      ${gtFlag} \
      ${cfFlag}

    ls -lah 04_phylogenetic_tree | head
    """
}

workflow {
  if (!params.buscodir) error "Falta --buscodir"

  def runSteps = params.steps.toString().split(',')*.trim() as Set
  def ch_buscodir = Channel.fromPath(params.buscodir, checkIfExists: true)

  // STEP 1
  if (runSteps.contains('1')) {
    STEP1_FIND_SINGLECOPY(ch_buscodir)
  } else {
    error "Step1 es requerido para este pipeline"
  }

  // STEP 2
  if (runSteps.contains('2')) {
    if (params.config) {
      def ch_cfg = Channel.fromPath(params.config, checkIfExists: true)
      STEP2_ALIGN_BUSCOS_WITH_CONFIG(STEP1_FIND_SINGLECOPY.out.fastadir, ch_cfg)
    } else {
      STEP2_ALIGN_BUSCOS_NO_CONFIG(STEP1_FIND_SINGLECOPY.out.fastadir)
    }
  }

  // STEP 3
  if (runSteps.contains('3')) {
    if (runSteps.contains('2')) {
      if (params.config) {
        STEP3_CREATE_MATRIX(STEP2_ALIGN_BUSCOS_WITH_CONFIG.out.alndir)
      } else {
        STEP3_CREATE_MATRIX(STEP2_ALIGN_BUSCOS_NO_CONFIG.out.alndir)
      }
    } else {
      error "Step3 requiere Step2"
    }
  }
  // STEP 4
  if (runSteps.contains('4')) {
    if (!runSteps.contains('3')) error "Step4 requiere Step3"
    if (!runSteps.contains('2')) error "Step4 requiere Step2 (aligndir para partitions automáticas)"

    if (params.config) {
      STEP4_PHYLOGENETIC_ANALYSIS(
        STEP3_CREATE_MATRIX.out.matrixdir,
        STEP3_CREATE_MATRIX.out.meta,
        STEP2_ALIGN_BUSCOS_WITH_CONFIG.out.alndir
      )
    } else {
      STEP4_PHYLOGENETIC_ANALYSIS(
        STEP3_CREATE_MATRIX.out.matrixdir,
        STEP3_CREATE_MATRIX.out.meta,
        STEP2_ALIGN_BUSCOS_NO_CONFIG.out.alndir
      )
    }
  }
}
