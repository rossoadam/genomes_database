#!/usr/bin/env python3

import argparse
import csv
import json
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Iterable

import pyfaidx
import pymysql

# Path - will be good for making filesystem work easier & safer
# dataclass - is cleaner than passing many loose variables
# Optional - helps document what may be missing
# Pyfaidx - fast indexing

"""
run:
    python3 busco_sql_v6_step0.py --genomes-dir /home/red_data/dns/genomes --db-password PASS
"""

def parse_args():
    parser = argparse.ArgumentParser(
        description = "Load Compleasm/BUSCO-derived features into SQL."
    )

    parser.add_argument(
        "--genomes-dir",
        required = True,
        help = "Path to genomes directory"
    )

    parser.add_argument(
        "--lineage",
        default = "sauropsida_odb12",
        help = "Default lineage if not found in metadata"
    )

    parser.add_argument(
        "--db-host",
        default = "localhost"
    )

    parser.add_argument(
        "--db-user",
        default = "root"
    )

    parser.add_argument(
        "--db-password",
        required = True
    )

    parser.add_argument(
        "--db-name",
        default = "test"
    )

    parser.add_argument(
        "--mode",
        choices=["all","one"],
        default = "all"
    )

    parser.add_argument(
        "--accession",
        help = "Required if --mode one"
    )

    return parser.parse_args()

def connect_db(args: argparse.Namespace):
    conn = pymysql.connect(
       host = args.db_host,
       user = args.db_user,
       password = args.db_password,
       db = args.db_name,
       charset = "utf8mb4",
       cursorclass = pymysql.cursors.DictCursor,
       local_infile = True,
       autocommit = False,
    )
    return conn

# this object stores eveerything the script needs for one genome
# cleaner than what I had previously e.g. self.acc_id, self.sp_id, etc.

@dataclass
class GenomeRun:
    accession: str
    species: str
    organism_name: str
    genome_fna: Path
    lineage: Optional[str]
    full_table: Optional[Path]
    cds_fasta: Optional[Path]
    gff_path: Optional[Path]

# add path helpers for the repository structure I'm not sure if this should go in the class or before it...
def get_records_paths(genomes_dir): # Path -> dict[str,Path]
    genomes_dir = Path(genomes_dir)
    records_dir = genomes_dir / "records"
    compleasm_records_dir = records_dir / "compleasm" / "records"
    return {
        "genomes_metadata": records_dir / "genomes_metadata.csv",
        "compleasm_metadata": compleasm_records_dir / "metadata.csv",
    }

def require_file(path: Path, label: str) -> None: # input Path object and labelt that is a string, this is just checking that the file / dir exists
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")

def normalize_species_name(organism_name: str) -> str: # input string output string no spaces
    organism_name = organism_name.strip() # remove trailing line enders and spaces?
    if "_" in organism_name:
        parts = [p for p in organism_name.split("_") if p] # split it on underscores
    else:
        parts = organism_name.split()

    if len(parts) >= 2:
#         if len(parts) == 3: # I think what I need to do is figure out a good way of storing subspecies - maybe it's own column in the sql table?
#             return f"{parts[0]}_{parts[1]}_{parts[2]}"
#         else: # if I unindent these lines I need to indent the following line
        return f"{parts[0]}_{parts[1]}"
    elif len(parts) == 1:
        return parts[0]
    return "unknown species"

def load_genomes_metadata(genomes_metadata_csv: Path) -> dict(str, dict): # this dictionary allows me to look up, organisms_name, path_to_fna for any accession
    rows = {}
    with open(genomes_metadata_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            accession = row["accession"].strip()
            rows[accession] = row
    return rows

def build_manifest(genomes_dir: Path, lineage_default: str) -> list[GenomeRun]:
    paths = get_records_paths(genomes_dir)

    require_file(paths["genomes_metadata"], "genomes metadata")
    genomes_rows = load_genomes_metadata(paths["genomes_metadata"])

    runs: list[GenomeRun] = []

    if paths["compleasm_metadata"].exists():
        with open(paths["compleasm_metadata"], newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                accession = row["accession"].strip()
                genome_row = genomes_rows.get(accession)

                if genome_row is None:
                    print(f"WARNING: accession missing from genomes metadata: {accession}")
                    continue

                organism_name = genome_row["organism_name"].strip()
                species = normalize_species_name(organism_name)

                genome_path_str = (
                    genome_row.get("path_to_fna")
                    or genome_row.get("pathway_to_genome")
                )

                if not genome_path_str:
                    print(f"WARNING: no genome FASTA path for accession: {accession}")
                    continue

                genome_fna = Path(genome_path_str).expanduser()

                lineage = row.get("lineage", lineage_default)
                full_table = Path(row["full_table"]).expanduser() if row.get("full_table") else None
                cds_fasta = Path(row["cds_fasta"]).expanduser() if row.get("cds_fasta") else None

                gff_path = None
                if full_table is not None:
                    gff_candidate = full_table.parent / "miniprot_output.gff"
                    if gff_candidate.exists():
                        gff_path = gff_candidate
                runs.append(
                    GenomeRun(
                        accession=accession,
                        species=species,
                        organism_name=organism_name,
                        genome_fna=genome_fna,
                        lineage=lineage,
                        full_table=full_table,
                        cds_fasta=cds_fasta,
                        gff_path=gff_path,
                    )
                )
    else:
        for accession, row in genomes_rows.items():
            organism_name = row["organism_name"].strip()
            species = normalize_species_name(organism_name)

            genome_path_str = row.get("path_to_fna") or row.get("pathway_to_genome")
            if not genome_path_str:
                print(f"WARNING: no genome FASTA path for accession: {accession}")
                continue

            genome_fna = Path(genome_path_str).expanduser()

            safe_org = re.sub(r"[^A-Za-z0-9_.-]+", "_", organism_name)
            outdir = genomes_dir / "records" / "compleasm" / f"{accession}__{safe_org}"
            lin_dir = outdir /lineage_default

            full_table = None
            for candidate in [lin_dir / "full_table.csv", lin_dir / "full_table.tsv"]:
                if candidate.exists():
                    full_table = candidate
                    break

            cds_fasta = lin_dir / f"{species.lower()}_cds_compleasm.fasta"
            if not cds_fasta.exists():
                cds_fasta = None

            gff_path = lin_dir / "miniprot_output.gff"
            if not gff_path.exists():
                gff_path = None

            runs.append(
                GenomeRun(
                    accession=accession,
                    species=species,
                    organism_name=organism_name,
                    genome_fna=genome_fna,
                    lineage=lineage_default,
                    full_table=full_table,
                    cds_fasta=cds_fasta,
                    gff_path=gff_path,
                    )
                )
    return runs

if __name__ == "__main__":
    args = parse_args()
    # print(args)
    connect = connect_db(args)
    # print(connect)
    # get_records_paths(args.genomes_dir)
    genomes_dir = Path(args.genomes_dir)
    lineage = args.lineage
    print(lineage)
    print(build_manifest(genomes_dir, lineage))

















