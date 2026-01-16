#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
align_BUSCOs.py

Adaptado por ChatGPT para usar configuración desde config.yaml (YAML).
"""

#%% Imports
from pathlib import Path
import argparse
from time import time
import subprocess

#%% Funciones
def check_arguments(FASTADIR, OUTDIR, CONFIG_FILE=False):
    if not Path(FASTADIR).is_dir():
        raise ValueError(f"FASTA directory {FASTADIR} does not exist.")
    if not Path(OUTDIR).is_dir():
        try:
            Path(OUTDIR).mkdir(parents=True, exist_ok=True)
        except:
            raise ValueError(f"Output directory {OUTDIR} could not be created.")
    if CONFIG_FILE and not Path(CONFIG_FILE).is_file():
        raise ValueError(f"Config file {CONFIG_FILE} does not exist.")

def is_fasta(filename):
    with open(filename, "r") as handle:
        first_line = handle.readline()
        return first_line.startswith(">")

def parse_mafft_config(CONFIG_FILE):
    settings = {
        "mafft_bin": "mafft",
        "align_method": "auto",
        "threads": 4
    }

    if CONFIG_FILE:
        import yaml
        try:
            with open(CONFIG_FILE) as f:
                config = yaml.safe_load(f) or {}
                align = config.get("alignment", None)
                if align is None:
                    align = config.get("align", {})
                if not isinstance(align, dict):
                    raise ValueError("YAML section 'alignment'/'align' must be a dict.")
                
                settings.update({
                    "mafft_bin": align.get("mafft_path", "mafft"),
                    "align_method": align.get("align_method", "auto"),
                    "threads": int(align.get("threads", 4))
                })
        except Exception as e:
            raise ValueError(f"Error al leer archivo de configuración YAML: {CONFIG_FILE}\n{e}")
    return settings

def generate_cmdline(infile, settings):
    """
    Build MAFFT command line as a list of arguments.
    """
    mafft_cmd = [settings["mafft_bin"], "--thread", str(settings["threads"])]

    if settings.get("align_method", "auto") == "auto":
        mafft_cmd.append("--auto")

    mafft_cmd.append(str(infile))
    return mafft_cmd

def run_mafft(infile, outfile, settings):
    """
    Run MAFFT and write alignment to outfile.
    """
    mafft_cmd = generate_cmdline(infile, settings)

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    with open(outfile, "w") as out_fh:
        result = subprocess.run(
            mafft_cmd,
            stdout=out_fh,
            stderr=subprocess.PIPE,
            text=True
        )

    if result.returncode != 0:
        raise RuntimeError(
            f"MAFFT failed for {infile}\n"
            f"Command: {' '.join(mafft_cmd)}\n"
            f"STDERR:\n{result.stderr}"
        )

    return mafft_cmd

def align_config_mafft(FASTADIR, OUTDIR, CONFIG_FILE=False, THREADS=None):
    if CONFIG_FILE:
        try:
            import yaml
        except ImportError:
            raise ImportError("PyYAML is required to use --config in Step 2. "
            "Install it with: pip install pyyaml"
            )
    check_arguments(FASTADIR, OUTDIR, CONFIG_FILE)
    files = [f for f in Path(FASTADIR).iterdir() if f.is_file()]
    if not files:
        raise RuntimeError(f"No FASTA files found in {FASTADIR}")

    print(f"{len(files)} BUSCO FASTA files found.")
    settings = parse_mafft_config(CONFIG_FILE)
    if THREADS is not None:
        settings["threads"] = int(THREADS)

    for busco in files:
        infile = busco
        if not is_fasta(infile):
            continue
        #mafft_cmline = generate_cmdline(infile, settings)
        #stdout, stderr = mafft_cmline()
        outfile = Path(OUTDIR) / f"{busco.stem}.aln.fa"
        run_mafft(infile, outfile, settings)
        #with open(outfile, "wt") as alnout:
        #    alnout.write(stdout)

def align_command_mafft(FASTADIR, OUTDIR, COMMAND):
    check_arguments(FASTADIR, OUTDIR)
    for busco in Path(FASTADIR).iterdir():
        if not busco.is_file() or not is_fasta(busco):
            continue
        infile = busco
        outfile = Path(OUTDIR) / f"{busco.stem}.aln.fa"
        #cmd_mafft = f"mafft {COMMAND} {infile} > {outfile}"
        cmd_mafft = ["mafft"] + COMMAND.split() + [str(infile)]
        with open(outfile, "wt") as out:
            result = subprocess.run(cmd_mafft, stdout=out, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"MAFFT failed for {infile}\nSTDERR:\n{result.stderr}")

def trim_alns(ALIGNDIR, OUTDIR_TRIM, TRIMPARAMS=False):
    list_of_alns = [p for p in Path(ALIGNDIR).iterdir() if p.is_file()]
    outdir = Path(OUTDIR_TRIM)
    outdir.mkdir(parents=True, exist_ok=True)

    for aln in list_of_alns:
        if not is_fasta(aln):
            continue
        outfile = outdir / f"{aln.stem}.trimmed.fa"
        cmd = ["trimal", "-in", str(aln), "-out", str(outfile)]
        if TRIMPARAMS:
            # si pasás algo tipo "-gt 0.3", lo partimos en tokens
            cmd += TRIMPARAMS.split()
        subprocess.run(cmd, check=True)

#%% CLI
def usage():
    parser = argparse.ArgumentParser(
        description='Align BUSCO orthologs using MAFFT and optionally trim with trimAl.')
    parser.add_argument('-f', '--fastadir', type=str, required=True,
                        help='Directory with BUSCO multi FASTA files.')
    parser.add_argument('-o', '--outdir', type=str, default='02_alignments',
                        help='Directory to save alignments. Default: 02_alignments')
    parser.add_argument('-c', '--config', metavar='<config file>', type=str, required=False,
                        help='YAML config file. Overrides default settings.')
    parser.add_argument('-m', '--command', metavar='<command>', type=str, required=False,
                        help='Manual MAFFT command options.')
    parser.add_argument('-t', '--trim', action='store_true', required=False,
                        help='Trim alignments using trimAl in automated mode.')
    parser.add_argument('-p', '--trimparams', type=str, required=False,
                        help='Custom parameters for trimAl (in quotes, e.g. "-gt 0.3").')
    return parser.parse_args()

#%% Main
if __name__ == '__main__':
    args = usage()
    start = time()

    step2_root = Path(args.outdir)
    raw_aln_dir = step2_root / "00_raw_aln"
    trim_aln_dir = step2_root / "01_trim_aln"
    raw_aln_dir.mkdir(parents=True, exist_ok=True)

    print("Step 2.1: Running MAFFT alignments...")

    if args.command:
        print("→ Using user-defined MAFFT command.")
        align_command_mafft(args.fastadir, raw_aln_dir, args.command)
    else:
        align_config_mafft(args.fastadir, raw_aln_dir, args.config)

    print("Step 2.2: Alignments completed.")

    if args.trim or args.trimparams:
        trimparams = args.trimparams
        if args.trim and not args.trimparams:
            trimparams = "-automated1"
        trim_alns(raw_aln_dir, trim_aln_dir, trimparams)
        print("→ Poorly aligned regions removed.")
    else:
        print("→ Trimming skipped.")

    print(f'Total time: {time() - start:.2f} seconds.')
