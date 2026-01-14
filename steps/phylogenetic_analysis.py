#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 00:37:16 2022

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
"""
#%% Imports
import argparse
from pathlib import Path
import subprocess
from time import time

#%% Function definition
def run_cmd(cmd, cwd=None):
    """
    Run a command robustly.
    cmd: list[str] (preferred) or str (if shell=True is really needed).
    cwd: directory where to run the command.
    """
    try:
        res = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
            cwd=cwd
        )
        return res
    except subprocess.CalledProcessError as e:
        cmd_str = cmd if isinstance(cmd, str) else " ".join(cmd)
        msg = (
            "Command failed.\n"
            f"Command: {cmd_str}\n"
            f"Return code: {e.returncode}\n"
            f"STDOUT:\n{e.stdout}\n"
            f"STDERR:\n{e.stderr}\n"
        )
        raise RuntimeError(msg) from e
#end

def model_partitions(MATRIXFILE, PARTITIONFILE, OUTDIR, SEQTYPE, PREFIX, BOOTSTRAP, THREADS):
    cwd = Path.cwd()

    outdir = Path(OUTDIR)
    outdir.mkdir(parents=True, exist_ok=True)

    # Resolver a paths absolutos
    matrixfile = Path(MATRIXFILE)
    if not matrixfile.is_absolute():
        matrixfile = (cwd / matrixfile).resolve()
    else:
        matrixfile = matrixfile.resolve()

    partitionfile = Path(PARTITIONFILE)
    if not partitionfile.is_absolute():
        partitionfile = (cwd / partitionfile).resolve()
    else:
        partitionfile = partitionfile.resolve()

    # Validar inputs
    if not matrixfile.exists():
        raise FileNotFoundError(f"Matrix file not found: {matrixfile}")
    if not partitionfile.exists():
        raise FileNotFoundError(f"Partition file not found: {partitionfile}")

    if int(BOOTSTRAP) < 1000:
        print("WARNING: bootstrap < 1000. OK for quick tests, not recommended for final inference.")


    # Comando (lista, sin shell)
    cmd = [
        "iqtree2",  # si tu binario es "iqtree", cambiá esto
        "-s", str(matrixfile),
        "-p", str(partitionfile),
        "-m", "MFP",
        "--seqtype", str(SEQTYPE),
        "--prefix", str(PREFIX),
        "-B", str(int(BOOTSTRAP)),
        "-T", str(int(THREADS)),
    ]

    # Ejecutar dentro del OUTDIR (sin chdir global)
    run_cmd(cmd, cwd=str(outdir))
#end
def _aligned_fasta_length(fasta_path: Path) -> int:
	"""
	Return alignment length reading only the first sequence in a multi-FASTA alignment.
	Assumes alignment is valid (all taxa have same length).
	"""
	length = 0
	in_seq = False
	with open(fasta_path, "rt") as fh:
		for line in fh:
			line = line.strip()
			if not line:
				continue
			if line.startswith(">"):
				if in_seq:
					break
				in_seq = True
				continue
			if in_seq:
				length += len(line)

	if length == 0:
		raise ValueError(f"Could not determine alignment length (empty/invalid FASTA): {fasta_path}")
	return length
#end
def _is_fasta_like(p: Path) -> bool:
	if not p.is_file():
		return False
	# extensiones típicas del pipeline
	fasta_suffixes = {".fa", ".faa", ".fna", ".ffn", ".fasta"}
	if p.suffix.lower() in fasta_suffixes:
		return True
	# casos como .aln.fa o .trimmed.fa
	name = p.name.lower()
	return name.endswith(".aln.fa") or name.endswith(".trimmed.fa") or name.endswith(".aln.fasta") or name.endswith(".trimmed.fasta")

def build_partitions_from_aligndir(ALIGNDIR, OUTDIR) -> Path:
	"""
	Build partitions.nex (Nexus sets/charset) using the same ordering used by create_matrix:
	sorted alignment files in ALIGNDIR.

	Works with .fa/.fasta and also .aln.fa/.trimmed.fa variants.
	"""
	aligndir = Path(ALIGNDIR)
	if not aligndir.is_dir():
		raise FileNotFoundError(f"Alignment directory not found: {aligndir}")

	outdir = Path(OUTDIR)
	outdir.mkdir(parents=True, exist_ok=True)

	aln_files = sorted([p for p in aligndir.iterdir() if _is_fasta_like(p)])
	if not aln_files:
		raise FileNotFoundError(f"No alignment FASTA files found in: {aligndir}")

	start = 1
	partitions_path = outdir / "partitions.nex"
	with open(partitions_path, "wt") as out:
		out.write("#nexus\nbegin sets;\n")
		for fp in aln_files:
			aln_len = _aligned_fasta_length(fp)
			end = start + aln_len - 1
			part_name = fp.stem.replace(".aln", "").replace(".trimmed", "")
			out.write(f"\tcharset {part_name} = {start}-{end};\n")
			start = end + 1
		out.write("end;\n")

	return partitions_path
#end
def gene_trees(MATRIXFILE, PREFIX, THREADS, OUTDIR):
	outdir = Path(OUTDIR)
	outdir.mkdir(parents=True, exist_ok=True)

	matrixfile = Path(MATRIXFILE)
	if not matrixfile.is_absolute():
		matrixfile = (Path.cwd() / matrixfile).resolve()
	else:
		matrixfile = matrixfile.resolve()

	best_scheme = Path(outdir / f"{PREFIX}.best_scheme.nex") # best_scheme = outdir / f"{PREFIX}.best_scheme.nex"
	if not best_scheme.is_absolute():
		best_scheme = (Path.cwd() / best_scheme).resolve()
	else:
		best_scheme = best_scheme.resolve()

	if not matrixfile.exists():
		raise FileNotFoundError(f"Matrix file not found: {matrixfile}")
	if not best_scheme.exists():
		raise FileNotFoundError(f"Best scheme file not found (run MFP first): {best_scheme}")

	# iqtree2 -s /home/nmoreyra/Soft/BUSCO2Tree/busco2tree/B2T_step4_genetrees/03_matrix/matrix.phylip 
	# -S test_sp.best_scheme.nex --prefix test_sp.loci -T 8
	print("### Estimating gene trees")
	cmd = [
        "iqtree2",
        "-s", str(matrixfile),
        "-S", str(best_scheme), #cambio -S por -p
        "--prefix", str(f"{PREFIX}.loci"),
        "-T", str(int(THREADS)),
    ]
	run_cmd(cmd, cwd=str(outdir))
#end

def concordance_factors(MATRIXFILE, PREFIX, THREADS, OUTDIR):
	outdir = Path(OUTDIR)
	outdir.mkdir(parents=True, exist_ok=True)

	matrixfile = Path(MATRIXFILE)
	if not matrixfile.is_absolute():
		matrixfile = (Path.cwd() / matrixfile).resolve()
	else:
		matrixfile = matrixfile.resolve()

	treefile = Path(outdir / f"{PREFIX}.treefile")
	if not treefile.is_absolute():
		treefile = (Path.cwd() / treefile).resolve()
	else:
		treefile = treefile.resolve()
	
	loci_treefile = Path(outdir / f"{PREFIX}.loci.treefile")
	if not loci_treefile.is_absolute():
		loci_treefile = (Path.cwd() / loci_treefile).resolve()
	else:
		loci_treefile = loci_treefile.resolve()

	best_scheme = Path(outdir / f"{PREFIX}.best_scheme.nex") # best_scheme = outdir / f"{PREFIX}.best_scheme.nex"
	if not best_scheme.is_absolute():
		best_scheme = (Path.cwd() / best_scheme).resolve()
	else:
		best_scheme = best_scheme.resolve()

	# Validaciones mínimas
	if not treefile.exists():
		raise FileNotFoundError(f"Missing main treefile: {treefile}")
	if not loci_treefile.exists():
		raise FileNotFoundError(f"Missing loci.treefile (run --genetrees first): {loci_treefile}")
	if not best_scheme.exists():
		raise FileNotFoundError(f"Missing best scheme file (run MFP first): {best_scheme}")
	if not matrixfile.exists():
		raise FileNotFoundError(f"Matrix file not found: {matrixfile}")

    # gCF
	print("### Estimating gene concordance factors")
	cmd_gcf = [
        "iqtree2",
        "-t", str(treefile),
        "--gcf", str(loci_treefile),
        "--prefix", f"{PREFIX}.loci.gCF",
        "-T", str(int(THREADS)),
    ]
	run_cmd(cmd_gcf, cwd=str(outdir))

    # sCF (usa el árbol resultante del paso gCF)
	print("### Estimating site concordance factors")
	#gcf_tree = outdir / f"{PREFIX}.loci.gCF.cf.tree"
	gcf_tree = Path(outdir / f"{PREFIX}.loci.gCF.cf.tree")
	print(gcf_tree)
	if not gcf_tree.is_absolute():
		gcf_tree = (Path.cwd() / gcf_tree).resolve()
	else:
		gcf_tree = gcf_tree.resolve()
	if not gcf_tree.exists():
		raise FileNotFoundError(f"Missing gCF tree output: {gcf_tree}")

	cmd_scf = [
        "iqtree2",
        "-te", str(gcf_tree),
        "-s", str(matrixfile),
        "-p", str(best_scheme),
        "-blfix",
        "--scf", "100",
        "--prefix", f"{PREFIX}.loci.sCF",
        "-T", str(int(THREADS)),
    ]
	run_cmd(cmd_scf, cwd=str(outdir))
#end

#%% Menu -> Usage
def usage():
	parser = argparse.ArgumentParser(
		description=(
            "phylogenetic_analysis.py runs IQ-TREE2 to infer a species tree from a concatenated matrix. "
            "Partitions are required for the intended BUSCO2Tree analysis. Provide --partitions, or provide "
			"--aligndir to auto-generate partitions from per-locus alignments."
		),
		epilog="End of the help")
	parser.add_argument('-m', '--matrix', type=str, required=True, help='Path to the matrix file.')
	parser.add_argument('-p', '--partitions', type=str, required=False, default=None, help='Path to a Nexus partitions file. If not provided, --aligndir is required to auto-generate partitions.')
	parser.add_argument('-a', '--aligndir', type=str, required=False, default=None, help="Directory with per-locus alignment files (*.fasta) used to build the concatenated matrix. "
						"Used only when --partitions is not provided.")
	parser.add_argument('-s', '--seqtype', type=str, choices=["DNA","AA"], default="AA", help='Sequence type. Default: AA')
	parser.add_argument('-o', '--outdir', type=str, required=False, default="04_phylogenetic_tree", help='Output directory. If it does not exist, it will be created. Default: 04_phylogenetic_tree.')
	parser.add_argument('-P', '--prefix', type=str, required=False, default="iqtree", help='Prefix for the output dataset and results. Default: iqtree')
	parser.add_argument('-b', '--bootstrap', type=int, required=False, default=1000, help='Number of bootstrap replicates. Default: 1000')
	parser.add_argument('-t', '--threads', type=int, required=False, default=8, help='Number of threads to run in IQ-Tree. Default: 8')
	parser.add_argument('-gt', '--genetrees', action='store_true', required=False, help='Generate loci (gene) trees using iqtree')
	parser.add_argument('-cf', '--concordance', action='store_true', required=False, help='Calculate concordance factors using IQ-TREE. Requires --genetrees.')
	return parser.parse_args()
#end

#%% Main program
if __name__ == '__main__':
	args = usage()
	start = time() #time 0
	print("Starting the phylogenetic analysis with IQ-Tree...")
	cwd = Path.cwd()
	matrix_path = Path(args.matrix)
	if not matrix_path.is_file():
			raise ValueError("The matrix file can not be found.")

	# Enforce: partitions siempre, pero permitir autogeneración desde aligndir
	partitions_path = None
	if args.partitions is not None:
			partitions_path = Path(args.partitions)
			if not partitions_path.is_file():
					raise ValueError("The partition file can not be found.")
	else:
			if args.aligndir is None:
					raise ValueError(
							"Partitions are required. Provide --partitions, or provide --aligndir to auto-generate partitions."
					)
			print("No partitions provided. Building partitions from --aligndir ...")
			partitions_path = build_partitions_from_aligndir(args.aligndir, args.outdir)
	# Ejecutar árbol particionado (siempre)
	model_partitions(args.matrix, str(partitions_path), args.outdir, args.seqtype, args.prefix, args.bootstrap, args.threads)

	if args.genetrees:
			matrixfile = Path(args.matrix).resolve()
			print("Estimating gene trees using IQ-Tree...")
			gene_trees(matrixfile, args.prefix, args.threads, args.outdir)
			if args.concordance:
					print("Calculating concordance factors using IQ-Tree...")
					concordance_factors(matrixfile, args.prefix, args.threads, args.outdir)
	print("IQ-Tree has finished.")
	print(f'Time taken to run: {time() - start} seconds.')
	
	# This is the main part of the script. It first parses the command line arguments, then starts the timer,
	# runs the model_partitions function with the parsed arguments, then prints the time taken to run.

#end
