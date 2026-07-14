#!/usr/bin/env python3
"""
Build starter SQL TSVs for the GC3 / genome dynamics database.

This script creates four load-ready TSVs:
  1. natural_history.tsv
  2. species_name_audit.tsv
  3. genomes.tsv
  4. sequences.tsv

Default project behavior:
  - Uses mass_predicts_dna_dynamics_with_s_punctatus_manifest.csv
  - Selects these five genomes by default:
      GCA_039797435.1
      GCA_051312515.2
      GCF_035046505.1
      GCA_053572275.1
      GCF_028583425.1
  - Writes outputs to genomes/records/sql_tsvs/
  - Leaves sequence.gc as SQL NULL (\\N) for now.

Example:
  python 00_build_starter_sql_tsvs.py \
    --genomes-dir /Users/rossoaa/projects/genomes \
    --manifest /Users/rossoaa/projects/genomes/records/project_manifests/mass_predicts_dna_dynamics_with_s_punctatus_manifest.csv \
    --natural-history /Users/rossoaa/projects/genomes/records/natural_history/natural_history.tsv \
    --species-name-audit /Users/rossoaa/projects/genomes/records/natural_history/species_name_audit.tsv

Optional MySQL loading:
  python 00_build_starter_sql_tsvs.py ... --load-sql --mysql-db gc3_dynamics --mysql-user root
"""

from __future__ import annotations

import argparse
import csv
import getpass
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    import pyfaidx
except ImportError:  # handled later with a clear error if sequences are requested
    pyfaidx = None


DEFAULT_ACCESSIONS = [
    "GCA_039797435.1",
    "GCA_051312515.2",
    "GCF_035046505.1",
    "GCA_053572275.1",
    "GCF_028583425.1",
]

NATURAL_HISTORY_COLUMNS = [
    "species_pk",
    "species_normalized",
    "mass_meiri",
    "mass_title",
    "mass_ji",
    "genome_size",
    "ct_min",
    "ct_max",
]

SPECIES_NAME_AUDIT_COLUMNS = [
    "species_name_audit_pk",
    "species_pk",
    "source_dataset",
    "source_species_name",
    "species_normalized",
    "match_status",
]

GENOMES_COLUMNS = [
    "genome_pk",
    "species_ncbi",
    "species_pk",
    "accession_id",
    "is_current",
]

SEQUENCES_COLUMNS = [
    "sequence_pk",
    "genome_pk",
    "sequence_id",
    "sequence_length",
    "sequence_type",
    "gc",
]

FASTA_SUFFIXES = (".fna", ".fa", ".fasta")


def normalize_species_name(value: str) -> str:
    """Normalize species names to lower-case genus_species style."""
    if value is None:
        return ""
    value = str(value).strip()
    value = re.sub(r"\([^)]*\)", " ", value)  # remove parenthetical annotations
    value = value.replace("-", " ").replace(".", " ")
    value = re.sub(r"[^A-Za-z0-9_\s]", " ", value)
    value = re.sub(r"\s+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value.lower()


def species_level_key(value: str) -> str:
    """Return genus_species from genus_species_subspecies when possible."""
    normalized = normalize_species_name(value)
    parts = normalized.split("_")
    if len(parts) >= 2:
        return "_".join(parts[:2])
    return normalized


def sql_null_if_blank(value: object) -> str:
    """Convert blank/None/NaN-ish values to MySQL LOAD DATA compatible NULL."""
    if value is None:
        return r"\N"
    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "null", "na"}:
        return r"\N"
    return text


def read_table(path: Path, delimiter: Optional[str] = None) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    if delimiter is None:
        delimiter = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"No header found in {path}")
        return [dict(row) for row in reader]


def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Dict[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: sql_null_if_blank(row.get(column, r"\N")) for column in columns})
            count += 1
    return count


def build_species_pk_lookup(natural_history_rows: Sequence[Dict[str, str]]) -> Dict[str, str]:
    lookup: Dict[str, str] = {}
    for row in natural_history_rows:
        key = normalize_species_name(row.get("species_normalized", ""))
        species_pk = sql_null_if_blank(row.get("species_pk"))
        if key and species_pk != r"\N":
            lookup[key] = species_pk
    return lookup


def filter_manifest_rows(manifest_rows: Sequence[Dict[str, str]], accessions: Sequence[str]) -> List[Dict[str, str]]:
    wanted = set(accessions)
    selected = [row for row in manifest_rows if row.get("accession") in wanted]
    found = {row.get("accession") for row in selected}
    missing = [acc for acc in accessions if acc not in found]
    if missing:
        raise ValueError("These accessions were not found in the manifest: " + ", ".join(missing))
    selected_by_accession = {row["accession"]: row for row in selected}
    return [selected_by_accession[acc] for acc in accessions]


def build_genomes_rows(
    manifest_rows: Sequence[Dict[str, str]],
    species_pk_lookup: Dict[str, str],
    start_pk: int = 1,
) -> Tuple[List[Dict[str, object]], Dict[str, int]]:
    rows: List[Dict[str, object]] = []
    accession_to_genome_pk: Dict[str, int] = {}

    for offset, row in enumerate(manifest_rows):
        genome_pk = start_pk + offset
        accession = row.get("accession", "").strip()
        species_key = row.get("species_key") or row.get("organism_name", "")
        normalized = normalize_species_name(species_key)
        fallback = species_level_key(species_key)
        species_pk = species_pk_lookup.get(normalized) or species_pk_lookup.get(fallback) or r"\N"

        rows.append(
            {
                "genome_pk": genome_pk,
                "species_ncbi": row.get("organism_name", r"\N"),
                "species_pk": species_pk,
                "accession_id": accession,
                "is_current": normalize_bool(row.get("is_current_for_species_at_freeze", "True")),
            }
        )
        accession_to_genome_pk[accession] = genome_pk
    return rows, accession_to_genome_pk


def normalize_bool(value: object) -> str:
    text = str(value).strip().lower()
    if text in {"true", "t", "1", "yes", "y"}:
        return "1"
    if text in {"false", "f", "0", "no", "n"}:
        return "0"
    return r"\N"


def candidate_fasta_paths(accession: str, manifest_row: Dict[str, str], genomes_dir: Path) -> List[Path]:
    candidates: List[Path] = []
    manifest_root = Path(manifest_row.get("accession_root", "")) if manifest_row.get("accession_root") else None

    roots = [
        genomes_dir / accession / "ncbi_dataset" / "data" / accession,
        genomes_dir / accession,
    ]
    if manifest_root:
        roots.extend([
            manifest_root / "ncbi_dataset" / "data" / accession,
            manifest_root,
        ])

    seen = set()
    for root in roots:
        if not root or root in seen:
            continue
        seen.add(root)
        if root.exists():
            for suffix in FASTA_SUFFIXES:
                candidates.extend(sorted(root.glob(f"*{suffix}")))
                candidates.extend(sorted(root.glob(f"*{suffix}.gz")))
            # recursive fallback for slightly different NCBI layouts
            for suffix in FASTA_SUFFIXES:
                candidates.extend(sorted(root.rglob(f"*{suffix}")))
                candidates.extend(sorted(root.rglob(f"*{suffix}.gz")))
    # pyfaidx does not index gzip FASTA directly in the same way; keep uncompressed first.
    unique: List[Path] = []
    seen_paths = set()
    for path in candidates:
        if path not in seen_paths:
            unique.append(path)
            seen_paths.add(path)
    unique.sort(key=lambda p: (str(p).endswith(".gz"), len(str(p))))
    return unique


def find_fasta(accession: str, manifest_row: Dict[str, str], genomes_dir: Path) -> Path:
    candidates = candidate_fasta_paths(accession, manifest_row, genomes_dir)
    uncompressed = [path for path in candidates if not str(path).endswith(".gz")]
    if uncompressed:
        return uncompressed[0]
    if candidates:
        raise ValueError(
            f"Only gzipped FASTA files were found for {accession}. "
            "Please decompress the FASTA before using pyfaidx. Candidate: " + str(candidates[0])
        )
    raise FileNotFoundError(
        f"No FASTA found for {accession}. Searched under {genomes_dir / accession} "
        "and the accession_root recorded in the manifest."
    )


def load_genome(input_fasta: Path):
    """Index/load a genome FASTA with pyfaidx, matching load_genome_pyfaidx.py behavior."""
    if pyfaidx is None:
        raise ImportError("pyfaidx is required. Install with: conda install -c bioconda pyfaidx")
    try:
        return pyfaidx.Fasta(str(input_fasta), as_raw=True, build_index=True)
    except Exception as exc:
        raise RuntimeError(f"Failed to load/index FASTA with pyfaidx: {input_fasta}\n{exc}") from exc


def infer_sequence_type(seq_name: str, description: str = "") -> str:
    """
    Infer sequence type from both FASTA sequence ID and full FASTA header.

    Priority:
      1. Explicit mitochondrial/plastid/organelle labels
      2. Explicit scaffold/contig/unplaced/unlocalized labels
      3. Explicit chromosome/chrom/lg/linkage group labels
      4. RefSeq/GenBank accession prefix fallback
      5. unclassified
    """
    text = f"{seq_name} {description}".lower()
    text = text.replace("|", " ")
    text = re.sub(r"[_\-:;,\[\]\(\)=]+", " ", text)
    text = re.sub(r"\s+", " ", text).strip()

    mitochondrial_patterns = [
        r"\bmitochondrion\b",
        r"\bmitochondrial\b",
        r"\bmitogenome\b",
        r"\bmt genome\b",
        r"\bmtDNA\b".lower(),
        r"\bchrM\b".lower(),
        r"\bchromosome m\b",
        r"\bmitochondrial chromosome\b",
    ]

    scaffold_patterns = [
        r"\bunplaced\b",
        r"\bunlocalized\b",
        r"\bscaffold\b",
        r"\bscaf\b",
        r"\bcontig\b",
        r"\bctg\b",
    ]

    wgs_scaffold_patterns = [
        r"\bwhole genome shotgun sequence\b",
        r"\bwgs sequence\b",
        r"\bshotgun sequence\b",
        r"\bgenomic scaffold\b",
        r"\bgenomic contig\b",
    ]

    chromosome_patterns = [
        r"\bchromosome\b",
        r"\bchrom\b",
        r"\bchr\b",
        r"\bmicrochromosome\b",
        r"\bmacrochromosome\b",
        r"\blinkage group\b",
        r"\blg\s*[0-9a-z]+\b",
        r"\bsex chromosome\b",
        r"\bchromosome\s*[0-9ivxlcdmzw]+\b",
        r"\bchr\s*[0-9ivxlcdmzw]+\b",
    ]

    for pattern in mitochondrial_patterns:
        if re.search(pattern, text):
            return "mitochondrion"

    for pattern in scaffold_patterns:
        if re.search(pattern, text):
            return "scaffold"
            
    for pattern in wgs_scaffold_patterns:
        if re.search(pattern, text):
            return "scaffold"
            
    for pattern in chromosome_patterns:
        if re.search(pattern, text):
            return "chromosome"

    seq_upper = str(seq_name).upper()
    
    

    # NCBI accession-style fallback:
    # NC_/CM_ are often chromosome-level, but mitochondrial records also use NC_,
    # so this fallback must come after mitochondrial text checks.
    
    if re.match(r"^[A-Z]{4,6}\d{6,}\.\d+$", seq_upper):
        return "scaffold"
    
    if seq_upper.startswith(("NC_", "CM_", "OX")):
        return "chromosome"

    if seq_upper.startswith(("NW_", "NT_", "NZ_")):
        return "scaffold"

    return "unclassified"


def build_sequences_rows(
    manifest_rows: Sequence[Dict[str, str]],
    accession_to_genome_pk: Dict[str, int],
    genomes_dir: Path,
    start_pk: int = 1,
    max_sequences_per_genome: Optional[int] = None,
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    sequence_pk = start_pk

    for manifest_row in manifest_rows:
        accession = manifest_row["accession"]
        genome_pk = accession_to_genome_pk[accession]
        fasta_path = find_fasta(accession, manifest_row, genomes_dir)
        print(f"[INFO] Loading {accession}: {fasta_path}", file=sys.stderr)
        genome = load_genome(fasta_path)

        sequence_ids = list(genome.keys())
        if max_sequences_per_genome is not None:
            sequence_ids = sequence_ids[:max_sequences_per_genome]

        for sequence_id in sequence_ids:
            rows.append(
                {
                    "sequence_pk": sequence_pk,
                    "genome_pk": genome_pk,
                    "sequence_id": sequence_id,
                    "sequence_length": len(genome[sequence_id]),
                    "sequence_type": infer_sequence_type(sequence_id, getattr(genome[sequence_id], "long_name", ""),),
                    "gc": r"\N",}
            )
            sequence_pk += 1
        genome.close()
    return rows


def load_tsvs_to_mysql(tsv_dir: Path, mysql_db: str, mysql_user: str, mysql_host: str, mysql_port: int) -> None:
    try:
        import pymysql
    except ImportError as exc:
        raise ImportError("pymysql is required for --load-sql. Install with: pip install pymysql") from exc

    password = getpass.getpass(f"MySQL password for {mysql_user}@{mysql_host}: ")
    conn = pymysql.connect(
        host=mysql_host,
        port=mysql_port,
        user=mysql_user,
        password=password,
        db=mysql_db,
        charset="utf8mb4",
        cursorclass=pymysql.cursors.DictCursor,
        local_infile=True,
    )
    load_order = ["natural_history", "species_name_audit", "genomes", "sequences"]
    columns_by_table = {
        "natural_history": NATURAL_HISTORY_COLUMNS,
        "species_name_audit": SPECIES_NAME_AUDIT_COLUMNS,
        "genomes": GENOMES_COLUMNS,
        "sequences": SEQUENCES_COLUMNS,
    }
    try:
        with conn.cursor() as cursor:
            for table in load_order:
                path = (tsv_dir / f"{table}.tsv").resolve()
                columns = ", ".join(f"`{col}`" for col in columns_by_table[table])
                sql = (
                    f"LOAD DATA LOCAL INFILE %s INTO TABLE `{table}` "
                    "FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' "
                    f"IGNORE 1 LINES ({columns});"
                )
                print(f"[INFO] Loading {path} into {mysql_db}.{table}", file=sys.stderr)
                cursor.execute(sql, (str(path),))
        conn.commit()
    finally:
        conn.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build starter SQL TSVs for natural_history, species_name_audit, genomes, and sequences."
    )
    parser.add_argument("--genomes-dir", default="genomes", help="Root genomes directory. Default: genomes")
    parser.add_argument(
        "--manifest",
        default="genomes/records/project_manifests/mass_predicts_dna_dynamics_with_s_punctatus_manifest.csv",
        help="Project manifest CSV.",
    )
    parser.add_argument(
        "--natural-history",
        default="genomes/records/natural_history/natural_history.tsv",
        help="Input natural_history TSV.",
    )
    parser.add_argument(
        "--species-name-audit",
        default="genomes/records/natural_history/species_name_audit.tsv",
        help="Input species_name_audit TSV.",
    )
    parser.add_argument(
        "--output-dir",
        default="genomes/records/sql_tsvs",
        help="Output directory for SQL TSVs. Default: genomes/records/sql_tsvs",
    )
    parser.add_argument(
        "--accessions",
        nargs="+",
        default=DEFAULT_ACCESSIONS,
        help="Accessions to include in genomes/sequences TSVs.",
    )
    parser.add_argument("--genome-pk-start", type=int, default=1, help="Starting genome_pk. Default: 1")
    parser.add_argument("--sequence-pk-start", type=int, default=1, help="Starting sequence_pk. Default: 1")
    parser.add_argument(
        "--max-sequences-per-genome",
        type=int,
        default=None,
        help="Optional testing/debug limit. Omit to include every sequence in each FASTA.",
    )
    parser.add_argument("--load-sql", action="store_true", help="Load generated TSVs into MySQL after writing them.")
    parser.add_argument("--mysql-db", default="gc3_dynamics", help="MySQL database name for --load-sql.")
    parser.add_argument("--mysql-user", default="root", help="MySQL user for --load-sql.")
    parser.add_argument("--mysql-host", default="localhost", help="MySQL host for --load-sql.")
    parser.add_argument("--mysql-port", type=int, default=3306, help="MySQL port for --load-sql.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    genomes_dir = Path(args.genomes_dir)
    manifest_path = Path(args.manifest)
    natural_history_path = Path(args.natural_history)
    species_name_audit_path = Path(args.species_name_audit)
    output_dir = Path(args.output_dir)
    
    # Resolve relative paths against genomes_dir
    if not manifest_path.is_absolute():
        manifest_path = genomes_dir / "records/project_manifests" / manifest_path.name
    
    if not natural_history_path.is_absolute():
        natural_history_path = genomes_dir / "records/sql_tsvs" / natural_history_path.name
    
    if not species_name_audit_path.is_absolute():
        species_name_audit_path = genomes_dir / "records/sql_tsvs" / species_name_audit_path.name
    
    if not output_dir.is_absolute():
        output_dir = genomes_dir / "records/sql_tsvs"
    
    manifest_path = manifest_path.resolve()
    natural_history_path = natural_history_path.resolve()
    species_name_audit_path = species_name_audit_path.resolve()
    output_dir = output_dir.resolve()
    manifest_rows = read_table(manifest_path, delimiter=",")
    natural_history_rows = read_table(natural_history_path, delimiter="\t")
    species_name_audit_rows = read_table(species_name_audit_path, delimiter="\t")

    selected_manifest_rows = filter_manifest_rows(manifest_rows, args.accessions)
    species_pk_lookup = build_species_pk_lookup(natural_history_rows)

    genomes_rows, accession_to_genome_pk = build_genomes_rows(
        selected_manifest_rows,
        species_pk_lookup,
        start_pk=args.genome_pk_start,
    )
    sequences_rows = build_sequences_rows(
        selected_manifest_rows,
        accession_to_genome_pk,
        genomes_dir=genomes_dir,
        start_pk=args.sequence_pk_start,
        max_sequences_per_genome=args.max_sequences_per_genome,
    )

    written = {
        "natural_history.tsv": write_tsv(output_dir / "natural_history.tsv", NATURAL_HISTORY_COLUMNS, natural_history_rows),
        "species_name_audit.tsv": write_tsv(output_dir / "species_name_audit.tsv", SPECIES_NAME_AUDIT_COLUMNS, species_name_audit_rows),
        "genomes.tsv": write_tsv(output_dir / "genomes.tsv", GENOMES_COLUMNS, genomes_rows),
        "sequences.tsv": write_tsv(output_dir / "sequences.tsv", SEQUENCES_COLUMNS, sequences_rows),
    }

    print("\n[OK] Wrote SQL TSVs:")
    for filename, n_rows in written.items():
        print(f"  {output_dir / filename}\t{n_rows} data rows")

    unmatched = [row for row in genomes_rows if row["species_pk"] == r"\N"]
    if unmatched:
        print("\n[WARN] These selected genomes did not map to natural_history.species_pk:", file=sys.stderr)
        for row in unmatched:
            print(f"  {row['accession_id']}\t{row['species_ncbi']}", file=sys.stderr)

    if args.load_sql:
        load_tsvs_to_mysql(output_dir, args.mysql_db, args.mysql_user, args.mysql_host, args.mysql_port)
        print(f"\n[OK] Loaded TSVs into MySQL database: {args.mysql_db}")


if __name__ == "__main__":
    main()
