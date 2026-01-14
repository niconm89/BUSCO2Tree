#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
find_singlecopy_BUSCOs.py

Adaptado por ChatGPT para leer configuración desde config.yaml opcionalmente.
"""

#%% Imports
#import yaml
import argparse
from pathlib import Path
from time import time  # ⚠️ Esto faltaba para usar `time()` en tu main

#%% Configuración
def load_config(config_path="config.yaml"):
    try:
        import yaml
    except ImportError:
        raise ImportError("PyYAML is required to use --config in Step 2. "
                          "Install it with: pip install pyyaml"
                          )
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"No se encontró el archivo de configuración: {config_path}")
    with open(config_path) as f:
        return yaml.safe_load(f)

#%% Functions definition
def check_paths(BUSCODIR, OUTDIR, out_seq_dir=False):
    if not Path(BUSCODIR).is_dir():
        raise ValueError("Error: BUSCO output directory %s does not exist." % BUSCODIR)
    if not Path(OUTDIR).is_dir():
        try:
            Path(OUTDIR).mkdir(parents=True, exist_ok=True)
        except:
            raise ValueError("Error: Output directory %s cannot be created." % OUTDIR)
    if out_seq_dir:
        if not Path(out_seq_dir).is_dir():
            try:
                Path(out_seq_dir).mkdir(parents=True, exist_ok=True)
            except:
                raise ValueError("Error: Output directory %s cannot be created." % out_seq_dir)

def find_singlecopy(BUSCODIR, OUTDIR, ODB, LINEAGE, seqtype="AA"):
    check_paths(BUSCODIR, OUTDIR)
    gene_dict = {}
    genomes_dirs = [f.name for f in Path(BUSCODIR).iterdir() if f.is_dir()]
    ngenomes = 0
    genomes_names = []
    common_busco_ids = []

    for genome in genomes_dirs:
        if genome not in ["batch_summary.txt", "logs", "list_common_busco.txt"]:
            genomes_names.append(genome)
            ngenomes += 1
            singlecopypath = (
                Path(BUSCODIR)
                / genome
                / f"run_{LINEAGE}_odb{ODB}"
                / "busco_sequences"
                / "single_copy_busco_sequences"
            )
            if not singlecopypath.is_dir():
                raise FileNotFoundError(
                    f"Missing BUSCO single-copy directory for genome '{genome}': {singlecopypath}"
                    )
            exts = [".faa"] if seqtype == "AA" else [".fna", ".ffn"]
            singlecopyfiles = []
            for ext in exts:
                singlecopyfiles.extend(singlecopypath.glob(f"*{ext}"))
            for busco in singlecopyfiles:
                if busco.name not in gene_dict:
                    gene_dict[busco.name] = 1
                else:
                    gene_dict[busco.name] += 1

    outfile = Path(OUTDIR) / "list_common_busco.txt"
    with open(outfile, 'wt') as BUSCO:
        for busco, n in gene_dict.items():
            if n == ngenomes:
                BUSCO.write(busco[:-4] + "\n")
                common_busco_ids.append(busco[:-4])

    namesfile = Path(OUTDIR) / "genomes_names.txt"
    with open(namesfile, 'wt') as NAMES:
        genomes_names = sorted(genomes_names)
        for name in genomes_names:
            NAMES.write(name + "\n")

    return genomes_names, common_busco_ids

def create_busco_fasta(BUSCODIR, OUTDIR, ODB, LINEAGE, genomes_names, common_busco_ids, common_busco_dir="common_busco_sequences", seqtype="AA"):
    out_seq_dir = Path(OUTDIR) / common_busco_dir
    check_paths(BUSCODIR, OUTDIR, out_seq_dir)
    for busco_id in common_busco_ids:
        out_ext = ".faa" if seqtype == "AA" else ".fna"
        busco_multifasta = out_seq_dir / f"{busco_id}{out_ext}"
        with open(busco_multifasta, 'wt') as MULTIFASTA:
            for genome in genomes_names:
                exts = [".faa"] if seqtype == "AA" else [".fna", ".ffn"]
                base = (
                Path(BUSCODIR)
                / genome
                / f"run_{LINEAGE}_odb{ODB}"
                / "busco_sequences"
                / "single_copy_busco_sequences"
                )

                busco_fasta_file = None
                for ext in exts:
                    candidate = base / f"{busco_id}{ext}"
                    if candidate.is_file():
                        busco_fasta_file = candidate
                        break
                if busco_fasta_file is None:
                    raise FileNotFoundError(
                        f"Missing BUSCO sequence for {genome} / {busco_id} in {base} (seqtype={seqtype})"
                        )
                with open(busco_fasta_file, 'rt') as FASTA:
                    for line in FASTA:
                        if line.startswith(">"):
                            MULTIFASTA.write(f">{genome}|{busco_id}\n")
                        else:
                            MULTIFASTA.write(line)

#%% Menu
def usage():
    parser = argparse.ArgumentParser(
        description="find_singlecopy.py busca genes BUSCO comunes entre varios genomas.",
        epilog="Fin de la ayuda.")
    parser.add_argument('-b', '--buscodir', type=str, required=True, help='Directorio con las salidas de BUSCO por genoma.')
    parser.add_argument('-o', '--outdir', type=str, default='01_single-copy', help='Directorio de salida. Default: 01_single-copy')
    parser.add_argument('-d', '--odb', type=int, help='Versión de OrthoDB (por ejemplo: 10)')
    parser.add_argument('-l', '--lineage', type=str, help='Linea BUSCO usada (por ejemplo: lepidoptera)')
    parser.add_argument('--seqtype', type=str, choices=['AA', 'DNA'], default='AA', help="Tipo de secuencia BUSCO a recopilar. AA usa .faa. DNA usa .fna/.ffn. Default: AA")
    parser.add_argument('--config', type=str, default='config.yaml', help='Archivo YAML de configuración')
    return parser.parse_args()

#%% Main
if __name__ == '__main__':
    args = usage()
    # Solo cargar config si hace falta (standalone debe funcionar con args explícitos)
    config = None
    if args.odb is None or args.lineage is None:
        config = load_config(args.config)
    
    ODB = args.odb if args.odb is not None else config["busco"]["odb_version"]
    LINEAGE = args.lineage if args.lineage is not None else config["busco"]["lineage"]

    start = time()
    print(f"Step 1.1: Buscando BUSCOs de copia única para el linaje {LINEAGE} (ODB v{ODB})...")

    genomes_names, common_busco_ids = find_singlecopy(args.buscodir, args.outdir, ODB, LINEAGE, seqtype=args.seqtype)

    print("Step 1.2: Creando archivos multifasta por BUSCO...")
    create_busco_fasta(args.buscodir, args.outdir, ODB, LINEAGE, genomes_names, common_busco_ids, seqtype=args.seqtype)

    print("BUSCOs comunes encontrados correctamente.")
    print(f'Tiempo total: {time() - start:.2f} segundos.')