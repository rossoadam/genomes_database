#!/usr/bin/env python3
"""
Build starter SQL TSVs for the GC3 / genome dynamics database.

This script creates load-ready TSVs:
  1. natural_history.tsv
  2. species_name_audit.tsv
  3. genomes.tsv
  4. sequences.tsv
  5. sequence_type_audit.tsv
  6. analysis_run.tsv
  7. window_set.tsv
  8. genomic_windows/*.tsv.gz
  9. gc_window_stats/*.tsv.gz
  10. sequence_summary/*.tsv.gz
  11. genome_summary/*.tsv.gz

Default project behavior:
  - Uses mass_predicts_dna_dynamics_with_s_punctatus_manifest.csv
  - Selects these five genomes by default:
      GCA_039797435.1
      GCA_051312515.2
      GCF_035046505.1
      GCA_053572275.1
      GCF_028583425.1
  - Use --n-accessions 10, 15, 20, etc. to keep the default five
    and add more accessions from the manifest in manifest order.
  - Use --all-accessions to include every accession in the manifest.
  - Uses genomes/records/genomes_metadata.csv to get the exact genomic FASTA path for each accession.
  - Writes small/reference outputs to genomes/records/sql_tsvs/
  - Writes high-volume window outputs as compressed chunks in table-specific directories:
      genomes/records/sql_tsvs/genomic_windows/*.tsv.gz
      genomes/records/sql_tsvs/gc_window_stats/*.tsv.gz
      genomes/records/sql_tsvs/sequence_summary/*.tsv.gz
      genomes/records/sql_tsvs/genome_summary/*.tsv.gz
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
import gzip
import getpass
import json
import math
import re
import statistics
import sys
from datetime import datetime
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

SEQUENCE_TYPE_AUDIT_COLUMNS = [
    "accession_id",
    "genome_pk",
    "n_fasta_sequences",
    "n_chromosome_inferred",
    "n_scaffold_inferred",
    "n_mitochondrion_inferred",
    "n_unclassified_inferred",
    "ncbi_total_chromosomes",
    "ncbi_number_scaffolds",
    "ncbi_number_contigs",
    "ncbi_component_sequences",
    "chromosome_count_match",
    "component_sequence_count_match",
    "audit_status",
    "assembly_report_path",
]


ANALYSIS_RUN_COLUMNS = [
    "run_pk",
    "analysis_name",
    "software",
    "software_version",
    "mask_mode",
    "gap_break_bp",
    "min_callable_frac",
    "created_at",
    "notes",
]

WINDOW_SET_COLUMNS = [
    "window_set_pk",
    "run_pk",
    "genome_pk",
    "tiling_type",
    "standard_window_size_bp",
    "step_size_bp",
    "start_offset_bp",
    "seq_scope",
    "notes",
]

GENOMIC_WINDOWS_COLUMNS = [
    "window_pk",
    "window_set_pk",
    "sequence_pk",
    "window_rank",
    "start_bp",
    "end_bp",
    "mid_bp",
    "standard_width_bp",
    "width_actual_bp",
    "callable_frac",
    "keep_flag",
]

GC_WINDOW_STATS_COLUMNS = [
    "window_pk",
    "a_count",
    "c_count",
    "g_count",
    "t_count",
    "n_count",
    "other_count",
    "callable_bp",
    "gc_bp",
    "gc_prop",
    "callable_frac",
    "masked_bp",
    "gap_bp",
]

SEQUENCE_SUMMARY_COLUMNS = [
    "sequence_summary_pk",
    "run_pk",
    "sequence_pk",
    "genome_pk",
    "standard_window_size_bp",
    "step_size_bp",
    "tiling_type",
    "mask_mode",
    "n_windows_total",
    "n_windows_kept",
    "n_windows_excluded_missing",
    "n_windows_excluded_short",
    "mean_gc",
    "weighted_mean_gc",
    "sd_gc",
    "var_gc",
    "median_gc",
    "mad_gc",
    "iqr_gc",
    "q05_gc",
    "q25_gc",
    "q75_gc",
    "q95_gc",
    "mean_callable_fraction",
    "median_callable_fraction",
    "callable_bp_total",
    "gap_bp_total",
]

GENOME_SUMMARY_COLUMNS = [
    "genome_summary_pk",
    "run_pk",
    "genome_pk",
    "species_pk",
    "standard_window_size_bp",
    "step_size_bp",
    "tiling_type",
    "mask_mode",
    "n_windows_total",
    "n_windows_kept",
    "mean_gc",
    "weighted_mean_gc",
    "sd_gc",
    "var_gc",
    "median_gc",
    "mad_gc",
    "iqr_gc",
    "q05_gc",
    "q25_gc",
    "q75_gc",
    "q95_gc",
    "mean_callable_fraction",
    "callable_bp_total",
    "gap_bp_total",
    "seq_count_used",
    "largest_seq_fraction",
]

WINDOW_CHUNK_MANIFEST_COLUMNS = [
    "table_name",
    "chunk_path",
    "run_pk",
    "genome_pk",
    "accession_id",
    "window_set_pk",
    "standard_window_size_bp",
    "step_size_bp",
    "n_rows",
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

def sql_pk(value: object) -> str:
    value = sql_null_if_blank(value)
    if value == r"\N":
        return value
    try:
        return str(int(float(value)))
    except ValueError:
        return str(value).strip()

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


PK_LIKE_COLUMNS = {
    "species_pk",
    "species_name_audit_pk",
    "genome_pk",
    "sequence_pk",
    "run_pk",
    "window_set_pk",
    "window_pk",
    "sequence_summary_pk",
    "genome_summary_pk",
}

def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Dict[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            out = {}
            for column in columns:
                value = row.get(column, r"\N")
                out[column] = sql_pk(value) if column in PK_LIKE_COLUMNS else sql_null_if_blank(value)
            writer.writerow(out)
            count += 1
    return count



def sanitize_filename_token(value: object) -> str:
    """Return a filesystem-safe token for chunked TSV filenames."""
    text = str(value).strip()
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "unknown"


def write_tsv_gz(path: Path, columns: Sequence[str], rows: Iterable[Dict[str, object]]) -> int:
    """Write a compressed TSV with the same formatting rules as write_tsv()."""
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            out = {}
            for column in columns:
                value = row.get(column, r"\N")
                out[column] = sql_pk(value) if column in PK_LIKE_COLUMNS else sql_null_if_blank(value)
            writer.writerow(out)
            count += 1
    return count


def chunk_filename(table: str, run_pk: int, genome_pk: int, accession: str, window_size_bp: int, step_size_bp: int) -> str:
    """Standard chunk filename for one table × accession × window-size combination."""
    accession_token = sanitize_filename_token(accession)
    return (
        f"{table}__run_{run_pk}__genome_pk_{genome_pk}__"
        f"accession_{accession_token}__window_{window_size_bp}bp__step_{step_size_bp}bp.tsv.gz"
    )


def write_window_chunk(
    output_dir: Path,
    table: str,
    columns: Sequence[str],
    rows: Sequence[Dict[str, object]],
    run_pk: int,
    genome_pk: int,
    accession: str,
    window_size_bp: int,
    step_size_bp: int,
) -> Tuple[str, int]:
    """Write one compressed chunk into output_dir/<table>/ and return relative path/count."""
    filename = chunk_filename(table, run_pk, genome_pk, accession, window_size_bp, step_size_bp)
    path = output_dir / table / filename
    n_rows = write_tsv_gz(path, columns, rows)
    return str(path.relative_to(output_dir)), n_rows

def build_species_pk_lookup(natural_history_rows: Sequence[Dict[str, str]]) -> Dict[str, str]:
    lookup: Dict[str, str] = {}
    for row in natural_history_rows:
        key = normalize_species_name(row.get("species_normalized", ""))
        species_pk = sql_null_if_blank(row.get("species_pk"))
        if key and species_pk != r"\N":
            lookup[key] = species_pk
    return lookup


def filter_manifest_rows(manifest_rows: Sequence[Dict[str, str]], accessions: Sequence[str]) -> List[Dict[str, str]]:
    """Return manifest rows for a requested accession list, preserving request order."""
    wanted = set(accessions)
    selected = [row for row in manifest_rows if row.get("accession") in wanted]
    found = {row.get("accession") for row in selected}
    missing = [acc for acc in accessions if acc not in found]
    if missing:
        raise ValueError("These accessions were not found in the manifest: " + ", ".join(missing))
    selected_by_accession = {row["accession"]: row for row in selected}
    return [selected_by_accession[acc] for acc in accessions]


def manifest_accessions(manifest_rows: Sequence[Dict[str, str]]) -> List[str]:
    """Return unique, nonblank accessions from the manifest in manifest order."""
    accessions: List[str] = []
    seen = set()
    for row in manifest_rows:
        accession = row.get("accession", "").strip()
        if accession and accession not in seen:
            accessions.append(accession)
            seen.add(accession)
    return accessions


def choose_accessions(
    manifest_rows: Sequence[Dict[str, str]],
    requested_accessions: Optional[Sequence[str]] = None,
    n_accessions: Optional[int] = None,
    all_accessions: bool = False,
) -> List[str]:
    """
    Decide which accessions to process.

    Default behavior is the original five DEFAULT_ACCESSIONS.
    --n-accessions keeps those defaults first, then appends additional manifest
    accessions in manifest order until the requested total is reached.
    --all-accessions returns every accession in the manifest in manifest order.
    """
    available = manifest_accessions(manifest_rows)

    if not available:
        raise ValueError("No nonblank accession values were found in the manifest.")

    if requested_accessions and (all_accessions or n_accessions is not None):
        raise ValueError("Use only one accession-selection mode: --accessions, --n-accessions, or --all-accessions.")

    if all_accessions:
        return available

    if requested_accessions:
        # Preserve user-provided order and remove accidental duplicates.
        selected: List[str] = []
        seen = set()
        for accession in requested_accessions:
            accession = str(accession).strip()
            if accession and accession not in seen:
                selected.append(accession)
                seen.add(accession)
        return selected

    selected = list(DEFAULT_ACCESSIONS)

    if n_accessions is None:
        return selected

    if n_accessions < len(DEFAULT_ACCESSIONS):
        raise ValueError(
            f"--n-accessions must be at least {len(DEFAULT_ACCESSIONS)} because "
            "the default five accessions are always included."
        )

    selected_set = set(selected)
    for accession in available:
        if len(selected) >= n_accessions:
            break
        if accession not in selected_set:
            selected.append(accession)
            selected_set.add(accession)

    if len(selected) < n_accessions:
        raise ValueError(
            f"Requested --n-accessions {n_accessions}, but only {len(selected)} unique "
            "accessions were available after combining the defaults with the manifest."
        )

    return selected


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


def accession_aliases(accession: str) -> List[str]:
    """Return versioned and unversioned accession aliases for robust metadata joins."""
    accession = str(accession).strip()
    if not accession:
        return []
    aliases = [accession]
    unversioned = re.sub(r"\.\d+$", "", accession)
    if unversioned and unversioned != accession:
        aliases.append(unversioned)
    return aliases


def resolve_project_path(path_text: object, genomes_dir: Path) -> Path:
    """Resolve paths from metadata; relative paths are treated as relative to genomes_dir."""
    path = Path(str(path_text).strip()).expanduser()
    if path.is_absolute():
        return path
    return genomes_dir / path


def looks_like_fasta_path(path: Path) -> bool:
    """Return True for common FASTA suffixes, including gzipped FASTA files."""
    name = path.name.lower()
    return name.endswith(FASTA_SUFFIXES) or any(name.endswith(f"{suffix}.gz") for suffix in FASTA_SUFFIXES)


def looks_like_genomic_fasta(path: Path) -> bool:
    """Reject obvious CDS/proteome/transcript FASTAs and accept genome FASTAs.

    Important: this function must not reject arbitrary substrings like "rna"
    inside assembly names. For example, rNatHel and rNatMau are valid assembly
    names, but both contain the letters "rna" when lower-cased. Therefore,
    rejection is based on path parts and tokenized filename words, not raw
    substring matching across the full path.
    """
    path_text = str(path).lower()
    name = path.name.lower()

    if not looks_like_fasta_path(path):
        return False

    # Always reject archived/old paths. These are where previous proteome/CDS
    # files were being found by recursive directory walking.
    path_parts = [part.lower() for part in path.parts]
    if any(part in {"old", "archive", "archived"} for part in path_parts):
        return False

    # Tokenize the filename so terms like RNA/CDS/protein are only rejected as
    # meaningful tokens, not as accidental substrings in assembly names.
    name_tokens = set(re.split(r"[^a-z0-9]+", name))
    reject_tokens = {
        "proteome",
        "protein",
        "proteins",
        "pep",
        "peptide",
        "peptides",
        "cds",
        "rna",
        "mrna",
        "transcript",
        "transcripts",
        "transcriptome",
    }
    if name_tokens.intersection(reject_tokens):
        return False

    # Strong positive signal for NCBI Datasets genome FASTAs. We still accept
    # other FASTA names from genomes_metadata.csv, because that file is the
    # project authority for exact genome paths.
    return True


def find_genomes_metadata(path: Optional[Path], genomes_dir: Path) -> Path:
    """Locate genomes_metadata.csv, defaulting to genomes/records/genomes_metadata.csv."""
    if path is None:
        return (genomes_dir / "records" / "genomes_metadata.csv").resolve()
    if not path.is_absolute():
        return (genomes_dir / "records" / path.name).resolve()
    return path.resolve()


def build_fasta_path_lookup(genomes_metadata_rows: Sequence[Dict[str, str]], genomes_dir: Path) -> Dict[str, Path]:
    """Build accession -> exact genomic FASTA path lookup from genomes_metadata.csv.

    The lookup stores both versioned and unversioned accession keys. This handles
    cases where the manifest has GCA_123.1 but genomes_metadata.csv stores
    GCA_123, or vice versa.
    """
    if not genomes_metadata_rows:
        raise ValueError("genomes_metadata.csv was empty; cannot build accession -> FASTA path lookup.")

    fieldnames = set(genomes_metadata_rows[0].keys())
    accession_columns = ["accession", "accession_id", "assembly_accession"]
    accession_column = next((col for col in accession_columns if col in fieldnames), None)
    if accession_column is None:
        raise ValueError(
            "Could not find an accession column in genomes_metadata.csv. "
            f"Tried: {', '.join(accession_columns)}"
        )

    preferred_path_columns = [
        "path_to_fna",
        "pathway_to_genome",
        "path_to_genome",
        "genome_fna_path",
        "genome_fasta_path",
        "genomic_fna_path",
        "genomic_fasta_path",
        "genome_path",
        "genomic_fna",
        "genomic_fasta",
        "fna_path",
        "fasta_path",
        "path",
    ]
    path_columns = [col for col in preferred_path_columns if col in fieldnames]

    if not path_columns:
        path_columns = [
            col for col in genomes_metadata_rows[0].keys()
            if any(token in col.lower() for token in ["path", "fna", "fasta"])
        ]

    if not path_columns:
        raise ValueError(
            "Could not find a FASTA path column in genomes_metadata.csv. "
            "Expected pathway_to_genome, genome_path, genomic_fna, or fasta_path."
        )

    lookup: Dict[str, Path] = {}
    rejected_examples: List[str] = []

    for row in genomes_metadata_rows:
        accession = row.get(accession_column, "").strip()
        if not accession:
            continue

        chosen: Optional[Path] = None
        considered: List[Path] = []
        for col in path_columns:
            raw = row.get(col, "")
            if not str(raw).strip() or str(raw).strip() == r"\N":
                continue
            path = resolve_project_path(raw, genomes_dir)
            considered.append(path)
            if looks_like_genomic_fasta(path):
                chosen = path
                break

        if chosen is None:
            if considered:
                rejected_examples.append(f"{accession}: " + "; ".join(str(p) for p in considered[:3]))
            continue

        for alias in accession_aliases(accession):
            lookup.setdefault(alias, chosen)

    if not lookup:
        example = "\n  ".join(rejected_examples[:10]) if rejected_examples else "No candidate paths were present."
        raise ValueError(
            "No usable genomic FASTA paths were found in genomes_metadata.csv. "
            "This script intentionally rejects CDS/proteome/transcript/old paths.\n"
            f"Examples considered:\n  {example}"
        )

    return lookup


def find_fasta(accession: str, fasta_path_lookup: Dict[str, Path]) -> Path:
    """Return exact genomic FASTA path for accession from genomes_metadata.csv."""
    fasta_path = None
    for alias in accession_aliases(accession):
        fasta_path = fasta_path_lookup.get(alias)
        if fasta_path is not None:
            break
    if fasta_path is None:
        raise FileNotFoundError(
            f"No genomic FASTA path was found for {accession} in genomes_metadata.csv. "
            "Update genomes_metadata.csv or remove this accession from the manifest/selection."
        )

    if str(fasta_path).endswith(".gz"):
        raise ValueError(
            f"The FASTA path for {accession} is gzipped: {fasta_path}. "
            "Please decompress it before using pyfaidx."
        )

    if not looks_like_genomic_fasta(fasta_path):
        raise ValueError(
            f"The metadata path for {accession} does not look like a genomic FASTA: {fasta_path}. "
            "This protects against accidentally loading CDS/proteome/transcript FASTAs."
        )

    if not fasta_path.exists():
        raise FileNotFoundError(f"Metadata FASTA path for {accession} does not exist: {fasta_path}")

    return fasta_path


def load_genome(input_fasta: Path):
    """Index/load a genome FASTA with pyfaidx, matching load_genome_pyfaidx.py behavior."""
    if pyfaidx is None:
        raise ImportError("pyfaidx is required. Install with: conda install -c bioconda pyfaidx")
    try:
        return pyfaidx.Fasta(str(input_fasta), as_raw=True, build_index=True)
    except Exception as exc:
        raise RuntimeError(f"Failed to load/index FASTA with pyfaidx: {input_fasta}\n{exc}") from exc

def get_sequence(genome_object, sequence_id: str, start_1_based: int, end_1_based: int) -> str:
    start_0_based = int(start_1_based) - 1
    end_0_based = int(end_1_based)
    return genome_object[sequence_id][start_0_based:end_0_based]

def infer_sequence_type(seq_name: str, description: str = "") -> str:
    """
    Infer sequence type from both FASTA sequence ID and full FASTA header.

    This function is intentionally conservative and NCBI-aware. A common
    NCBI chromosome-level header can contain both chromosome evidence and
    WGS wording, for example:

      CM012345.1 Species chromosome 1, whole genome shotgun sequence

    In v3, broad scaffold/WGS rules could be evaluated before chromosome
    evidence, causing chromosome records to be labeled as scaffolds. In v4,
    mitochondrial labels are still highest priority, but chromosome evidence
    and chromosome-level accession prefixes are evaluated before broad WGS
    wording. Explicit unplaced/unlocalized/scaffold/contig labels still win
    over generic chromosome wording when they are present.
    """
    raw_text = f"{seq_name} {description}"
    text = raw_text.lower()
    text = text.replace("|", " ")
    text = re.sub(r"[_\-:;,\[\]\(\)=]+", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    seq_upper = str(seq_name).upper()

    def has_any(patterns: Sequence[str]) -> bool:
        return any(re.search(pattern, text) for pattern in patterns)

    mitochondrial_patterns = [
        r"\bmitochondrion\b",
        r"\bmitochondrial\b",
        r"\bmitogenome\b",
        r"\bmt\s+genome\b",
        r"\bmtdna\b",
        r"\bchrm\b",
        r"\bchr\s*m\b",
        r"\bchromosome\s+m\b",
        r"\bmitochondrial\s+chromosome\b",
    ]

    # These labels are strong evidence that a record is not a named nuclear
    # chromosome, even if the organism name or other prose contains "chromosome".
    explicit_scaffold_patterns = [
        r"\bunplaced\b",
        r"\bunlocalized\b",
        r"\bscaffold\b",
        r"\bscaf\b",
        r"\bcontig\b",
        r"\bctg\b",
    ]

    chromosome_patterns = [
        r"\bchromosome\s+[0-9ivxlcdmxyzw]+\b",
        r"\bchromosome\b",
        r"\bchrom\s+[0-9ivxlcdmxyzw]+\b",
        r"\bchr\s*[0-9ivxlcdmxyzw]+\b",
        r"\bmicrochromosome\b",
        r"\bmacrochromosome\b",
        r"\blinkage\s+group\b",
        r"\blg\s*[0-9a-z]+\b",
        r"\bsex\s+chromosome\b",
    ]

    broad_wgs_scaffold_patterns = [
        r"\bwhole\s+genome\s+shotgun\s+sequence\b",
        r"\bwgs\s+sequence\b",
        r"\bshotgun\s+sequence\b",
        r"\bgenomic\s+scaffold\b",
        r"\bgenomic\s+contig\b",
    ]

    if has_any(mitochondrial_patterns):
        return "mitochondrion"

    # Accession-prefix fallbacks are useful because pyfaidx sequence IDs often
    # only contain the accession, not the full NCBI description.
    # NC_ can be nuclear chromosome or mitochondrial; mitochondrial text above
    # must be checked first.
    if seq_upper.startswith(("NC_", "CM_", "OX_", "OW_", "CP_")):
        return "chromosome"

    if seq_upper.startswith(("NW_", "NT_", "NZ_")):
        return "scaffold"

    # Some GenBank chromosome-scale records have project-style accessions rather
    # than CM_/NC_ prefixes. If the header says chromosome, trust that before
    # broad WGS wording such as "whole genome shotgun sequence".
    has_explicit_scaffold = has_any(explicit_scaffold_patterns)
    has_chromosome = has_any(chromosome_patterns)

    if has_explicit_scaffold:
        return "scaffold"

    if has_chromosome:
        return "chromosome"

    if has_any(broad_wgs_scaffold_patterns):
        return "scaffold"

    # Generic WGS contig/scaffold accessions such as JAAAAA010000001.1 are not
    # chromosome-level unless the header supplied chromosome evidence above.
    if re.match(r"^[A-Z]{4,6}\d{6,}\.\d+$", seq_upper):
        return "scaffold"

    return "unclassified"


def build_sequences_rows(
    manifest_rows: Sequence[Dict[str, str]],
    accession_to_genome_pk: Dict[str, int],
    fasta_path_lookup: Dict[str, Path],
    start_pk: int = 1,
    max_sequences_per_genome: Optional[int] = None,
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    sequence_pk = start_pk

    for manifest_row in manifest_rows:
        accession = manifest_row["accession"]
        genome_pk = accession_to_genome_pk[accession]
        fasta_path = find_fasta(accession, fasta_path_lookup)
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


def safe_int(value: object) -> Optional[int]:
    """Convert NCBI JSON numeric values to int when possible."""
    if value is None:
        return None
    try:
        text = str(value).strip()
        if text == "" or text.lower() in {"nan", "none", "null", "na"}:
            return None
        return int(float(text))
    except (TypeError, ValueError):
        return None


def find_assembly_data_report(accession: str, genomes_dir: Path) -> Optional[Path]:
    """
    Locate the NCBI Datasets assembly_data_report.jsonl for one accession.

    Expected layout:
      genomes/<accession>/ncbi_dataset/data/assembly_data_report.jsonl

    A recursive fallback is included for slightly different package layouts.
    """
    direct_path = genomes_dir / accession / "ncbi_dataset" / "data" / "assembly_data_report.jsonl"
    if direct_path.exists():
        return direct_path

    accession_root = genomes_dir / accession
    if accession_root.exists():
        matches = sorted(accession_root.rglob("assembly_data_report.jsonl"))
        if matches:
            return matches[0]

    return None


def read_ncbi_assembly_stats(accession: str, genomes_dir: Path) -> Dict[str, object]:
    """
    Read NCBI Datasets assemblyStats for one accession.

    Returns NULL-like values when the report is missing or fields are absent.
    """
    report_path = find_assembly_data_report(accession, genomes_dir)

    empty = {
        "ncbi_total_chromosomes": r"\N",
        "ncbi_number_scaffolds": r"\N",
        "ncbi_number_contigs": r"\N",
        "ncbi_component_sequences": r"\N",
        "assembly_report_path": r"\N",
    }

    if report_path is None:
        return empty

    with report_path.open("r") as handle:
        first_line = handle.readline().strip()

    if not first_line:
        empty["assembly_report_path"] = str(report_path)
        return empty

    record = json.loads(first_line)
    stats = record.get("assemblyStats", {})

    return {
        "ncbi_total_chromosomes": safe_int(stats.get("totalNumberOfChromosomes")) or r"\N",
        "ncbi_number_scaffolds": safe_int(stats.get("numberOfScaffolds")) or r"\N",
        "ncbi_number_contigs": safe_int(stats.get("numberOfContigs")) or r"\N",
        "ncbi_component_sequences": safe_int(stats.get("numberOfComponentSequences")) or r"\N",
        "assembly_report_path": str(report_path),
    }


def yes_no_null(observed: int, expected: object) -> str:
    """Return 1/0/NULL-style comparison for audit columns."""
    expected_int = safe_int(expected)
    if expected_int is None:
        return r"\N"
    return "1" if observed == expected_int else "0"


def build_sequence_type_audit_rows(
    sequences_rows: Sequence[Dict[str, object]],
    genomes_rows: Sequence[Dict[str, object]],
    genomes_dir: Path,
) -> List[Dict[str, object]]:
    """
    Summarize inferred sequence types by genome and compare to NCBI assemblyStats.

    Notes:
      - ncbi_component_sequences should usually match total FASTA sequences.
      - ncbi_total_chromosomes should usually match n_chromosome_inferred.
      - ncbi_number_scaffolds is reported for reference, but may include
        chromosome-scale scaffolds, so it is not used as a hard pass/fail check.
    """
    genome_pk_to_accession = {
        int(row["genome_pk"]): str(row["accession_id"]) for row in genomes_rows
    }

    counts: Dict[int, Dict[str, int]] = {}
    for row in sequences_rows:
        genome_pk = int(row["genome_pk"])
        seq_type = str(row.get("sequence_type", "unclassified"))
        if genome_pk not in counts:
            counts[genome_pk] = {
                "n_fasta_sequences": 0,
                "n_chromosome_inferred": 0,
                "n_scaffold_inferred": 0,
                "n_mitochondrion_inferred": 0,
                "n_unclassified_inferred": 0,
            }

        counts[genome_pk]["n_fasta_sequences"] += 1
        if seq_type == "chromosome":
            counts[genome_pk]["n_chromosome_inferred"] += 1
        elif seq_type == "scaffold":
            counts[genome_pk]["n_scaffold_inferred"] += 1
        elif seq_type == "mitochondrion":
            counts[genome_pk]["n_mitochondrion_inferred"] += 1
        else:
            counts[genome_pk]["n_unclassified_inferred"] += 1

    audit_rows: List[Dict[str, object]] = []
    for genome_pk in sorted(genome_pk_to_accession):
        accession = genome_pk_to_accession[genome_pk]
        observed = counts.get(
            genome_pk,
            {
                "n_fasta_sequences": 0,
                "n_chromosome_inferred": 0,
                "n_scaffold_inferred": 0,
                "n_mitochondrion_inferred": 0,
                "n_unclassified_inferred": 0,
            },
        )
        ncbi = read_ncbi_assembly_stats(accession, genomes_dir)

        chromosome_match = yes_no_null(
            observed["n_chromosome_inferred"],
            ncbi["ncbi_total_chromosomes"],
        )
        component_match = yes_no_null(
            observed["n_fasta_sequences"],
            ncbi["ncbi_component_sequences"],
        )

        if chromosome_match == "1" and component_match == "1":
            audit_status = "pass"
        elif chromosome_match == r"\N" or component_match == r"\N":
            audit_status = "missing_ncbi_stats"
        else:
            audit_status = "review"

        audit_rows.append(
            {
                "accession_id": accession,
                "genome_pk": genome_pk,
                **observed,
                **ncbi,
                "chromosome_count_match": chromosome_match,
                "component_sequence_count_match": component_match,
                "audit_status": audit_status,
            }
        )

    return audit_rows



def parse_positive_int_list(values: Sequence[int]) -> List[int]:
    """Return sorted unique positive integers, preserving the user's broad scale order descending."""
    clean: List[int] = []
    seen = set()
    for value in values:
        ivalue = int(value)
        if ivalue <= 0:
            raise ValueError("Window sizes and step sizes must be positive integers.")
        if ivalue not in seen:
            clean.append(ivalue)
            seen.add(ivalue)
    return clean


def base_counts(sequence: str) -> Dict[str, int]:
    """Count bases for a window; callable bases are A/C/G/T only."""
    seq = str(sequence).upper()
    a = seq.count("A")
    c = seq.count("C")
    g = seq.count("G")
    t = seq.count("T")
    n = seq.count("N")
    gap = seq.count("-")
    other = len(seq) - (a + c + g + t + n + gap)
    callable_bp = a + c + g + t
    gc_bp = g + c
    gc_prop = (gc_bp / callable_bp) if callable_bp > 0 else None
    callable_frac = (callable_bp / len(seq)) if len(seq) > 0 else None
    return {
        "a_count": a,
        "c_count": c,
        "g_count": g,
        "t_count": t,
        "n_count": n,
        "other_count": other,
        "callable_bp": callable_bp,
        "gc_bp": gc_bp,
        "gc_prop": gc_prop if gc_prop is not None else r"\N",
        "callable_frac": callable_frac if callable_frac is not None else r"\N",
        "masked_bp": n,
        "gap_bp": gap,
    }


def quantile(values: Sequence[float], q: float) -> object:
    """Linear-interpolated quantile with SQL NULL for empty input."""
    vals = sorted(float(v) for v in values if v is not None)
    if not vals:
        return r"\N"
    if len(vals) == 1:
        return vals[0]
    pos = (len(vals) - 1) * q
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return vals[lo]
    return vals[lo] * (hi - pos) + vals[hi] * (pos - lo)


def median_abs_deviation(values: Sequence[float]) -> object:
    vals = [float(v) for v in values if v is not None]
    if not vals:
        return r"\N"
    med = statistics.median(vals)
    return statistics.median([abs(v - med) for v in vals])


def summarize_gc_values(gc_values: Sequence[float], callable_fracs: Sequence[float], callable_bps: Sequence[int], gap_bps: Sequence[int]) -> Dict[str, object]:
    """Summary stats used by sequence_summary and genome_summary."""
    vals = [float(v) for v in gc_values if v is not None]
    call_fracs = [float(v) for v in callable_fracs if v is not None]
    total_callable = sum(int(v) for v in callable_bps)
    weighted_mean = r"\N"
    if total_callable > 0 and vals:
        weighted_mean = sum(float(gc) * int(bp) for gc, bp in zip(gc_values, callable_bps) if gc is not None) / total_callable

    if vals:
        mean_gc = statistics.mean(vals)
        var_gc = statistics.variance(vals) if len(vals) > 1 else 0.0
        sd_gc = math.sqrt(var_gc)
        median_gc = statistics.median(vals)
        iqr_gc = quantile(vals, 0.75) - quantile(vals, 0.25) if len(vals) > 0 else r"\N"
    else:
        mean_gc = weighted_mean_gc = sd_gc = var_gc = median_gc = iqr_gc = r"\N"

    return {
        "mean_gc": mean_gc if vals else r"\N",
        "weighted_mean_gc": weighted_mean,
        "sd_gc": sd_gc if vals else r"\N",
        "var_gc": var_gc if vals else r"\N",
        "median_gc": median_gc if vals else r"\N",
        "mad_gc": median_abs_deviation(vals),
        "iqr_gc": iqr_gc if vals else r"\N",
        "q05_gc": quantile(vals, 0.05),
        "q25_gc": quantile(vals, 0.25),
        "q75_gc": quantile(vals, 0.75),
        "q95_gc": quantile(vals, 0.95),
        "mean_callable_fraction": statistics.mean(call_fracs) if call_fracs else r"\N",
        "median_callable_fraction": statistics.median(call_fracs) if call_fracs else r"\N",
        "callable_bp_total": total_callable,
        "gap_bp_total": sum(int(v) for v in gap_bps),
    }


def build_analysis_run_rows(run_pk: int, analysis_name: str, mask_mode: str, min_callable_frac: float, notes: str) -> List[Dict[str, object]]:
    return [{
        "run_pk": run_pk,
        "analysis_name": analysis_name,
        "software": "00_build_starter_sql_tsvs_v12.py",
        "software_version": "v12",
        "mask_mode": mask_mode,
        "gap_break_bp": r"\N",
        "min_callable_frac": min_callable_frac,
        "created_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "notes": notes,
    }]


def build_window_tables(
    manifest_rows: Sequence[Dict[str, str]],
    genomes_rows: Sequence[Dict[str, object]],
    sequences_rows: Sequence[Dict[str, object]],
    fasta_path_lookup: Dict[str, Path],
    run_pk: int,
    window_sizes_bp: Sequence[int],
    step_sizes_bp: Optional[Sequence[int]],
    sequence_types: Sequence[str],
    min_callable_frac: float,
    tiling_type: str,
    start_offset_bp: int,
    seq_scope: str,
    mask_mode: str,
    output_dir: Path,
    start_window_set_pk: int = 1,
    start_window_pk: int = 1,
    start_sequence_summary_pk: int = 1,
    start_genome_summary_pk: int = 1,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    """Build window_set rows and write window-derived tables as compressed chunks.

    The high-volume tables are written directly to disk as accession × window-size
    chunks to avoid holding millions of rows in memory and to make R exploration
    easier. Output layout:

      output_dir/genomic_windows/*.tsv.gz
      output_dir/gc_window_stats/*.tsv.gz
      output_dir/sequence_summary/*.tsv.gz
      output_dir/genome_summary/*.tsv.gz

    Coordinates are 1-based inclusive in the SQL tables. pyfaidx slicing is 0-based,
    so start/end are converted at slice time only.
    """
    window_sizes = parse_positive_int_list(window_sizes_bp)
    if step_sizes_bp is None:
        step_sizes = list(window_sizes)
    else:
        step_sizes = parse_positive_int_list(step_sizes_bp)
        if len(step_sizes) != len(window_sizes):
            raise ValueError("--step-sizes-bp must have the same number of values as --window-sizes-bp.")

    allowed_types = {str(t).strip().lower() for t in sequence_types if str(t).strip()}
    genome_pk_to_accession = {int(row["genome_pk"]): str(row["accession_id"]) for row in genomes_rows}
    genome_pk_to_species_pk = {int(row["genome_pk"]): row.get("species_pk", r"\N") for row in genomes_rows}
    sequence_by_genome: Dict[int, List[Dict[str, object]]] = {}
    for row in sequences_rows:
        seq_type = str(row.get("sequence_type", "")).lower()
        if allowed_types and seq_type not in allowed_types:
            continue
        sequence_by_genome.setdefault(int(row["genome_pk"]), []).append(row)

    window_set_rows: List[Dict[str, object]] = []
    chunk_manifest_rows: List[Dict[str, object]] = []

    window_set_pk = start_window_set_pk
    window_pk = start_window_pk
    sequence_summary_pk = start_sequence_summary_pk
    genome_summary_pk = start_genome_summary_pk

    for genome_pk in sorted(genome_pk_to_accession):
        accession = genome_pk_to_accession[genome_pk]
        fasta_path = find_fasta(accession, fasta_path_lookup)
        print(f"[INFO] Building window chunks for {accession}: {fasta_path}", file=sys.stderr)
        genome = load_genome(fasta_path)
        selected_sequences = sequence_by_genome.get(genome_pk, [])

        for size_bp, step_bp in zip(window_sizes, step_sizes):
            current_window_set_pk = window_set_pk
            window_set_rows.append({
                "window_set_pk": current_window_set_pk,
                "run_pk": run_pk,
                "genome_pk": genome_pk,
                "tiling_type": tiling_type,
                "standard_window_size_bp": size_bp,
                "step_size_bp": step_bp,
                "start_offset_bp": start_offset_bp,
                "seq_scope": seq_scope,
                "notes": f"Window chunks written by accession/window size for sequence_type in {sorted(allowed_types)}; created for GC variance decay analyses.",
            })
            window_set_pk += 1

            genomic_window_rows: List[Dict[str, object]] = []
            gc_window_rows: List[Dict[str, object]] = []
            sequence_summary_rows: List[Dict[str, object]] = []
            genome_summary_rows: List[Dict[str, object]] = []

            genome_gc_values: List[float] = []
            genome_callable_fracs: List[float] = []
            genome_callable_bps: List[int] = []
            genome_gap_bps: List[int] = []
            genome_total_windows = 0
            genome_kept_windows = 0
            genome_sequence_bp_used = 0
            largest_sequence_bp = 0
            seq_count_used = 0

            for seq_row in selected_sequences:
                sequence_pk = int(seq_row["sequence_pk"])
                sequence_id = str(seq_row["sequence_id"])
                seq_len = int(seq_row["sequence_length"])
                if sequence_id not in genome:
                    print(f"[WARN] Sequence {sequence_id} not found in FASTA for {accession}; skipping.", file=sys.stderr)
                    continue
                seq_count_used += 1
                genome_sequence_bp_used += seq_len
                largest_sequence_bp = max(largest_sequence_bp, seq_len)

                seq_gc_values: List[float] = []
                seq_callable_fracs: List[float] = []
                seq_callable_bps: List[int] = []
                seq_gap_bps: List[int] = []
                seq_total_windows = 0
                seq_kept_windows = 0
                seq_excluded_missing = 0
                seq_excluded_short = 0

                window_rank = 0
                start0 = max(0, int(start_offset_bp))
                for window_start_0 in range(start0, seq_len, step_bp):
                    window_end_0 = min(window_start_0 + size_bp, seq_len)
                    if window_end_0 <= window_start_0:
                        continue
                    width_actual = window_end_0 - window_start_0
                    window_rank += 1
                    seq_total_windows += 1
                    genome_total_windows += 1
                    start_1_based = window_start_0 + 1
                    end_1_based = window_end_0

                    seq_fragment = get_sequence(genome, sequence_id, start_1_based, end_1_based)
                    counts = base_counts(str(seq_fragment))
                    callable_frac = counts["callable_frac"] if counts["callable_frac"] != r"\N" else 0.0
                    gc_prop = counts["gc_prop"] if counts["gc_prop"] != r"\N" else None
                    keep_flag = 1 if (width_actual == size_bp and callable_frac >= min_callable_frac and gc_prop is not None) else 0
                    if keep_flag:
                        seq_kept_windows += 1
                        genome_kept_windows += 1
                        seq_gc_values.append(float(gc_prop))
                        seq_callable_fracs.append(float(callable_frac))
                        seq_callable_bps.append(int(counts["callable_bp"]))
                        seq_gap_bps.append(int(counts["gap_bp"]))
                        genome_gc_values.append(float(gc_prop))
                        genome_callable_fracs.append(float(callable_frac))
                        genome_callable_bps.append(int(counts["callable_bp"]))
                        genome_gap_bps.append(int(counts["gap_bp"]))
                    else:
                        if width_actual < size_bp:
                            seq_excluded_short += 1
                        else:
                            seq_excluded_missing += 1

                    start_bp = window_start_0 + 1
                    end_bp = window_end_0
                    genomic_window_rows.append({
                        "window_pk": window_pk,
                        "window_set_pk": current_window_set_pk,
                        "sequence_pk": sequence_pk,
                        "window_rank": window_rank,
                        "start_bp": start_bp,
                        "end_bp": end_bp,
                        "mid_bp": (start_bp + end_bp) // 2,
                        "standard_width_bp": size_bp,
                        "width_actual_bp": width_actual,
                        "callable_frac": callable_frac,
                        "keep_flag": keep_flag,
                    })
                    counts["window_pk"] = window_pk
                    gc_window_rows.append(counts)
                    window_pk += 1

                seq_summary = summarize_gc_values(seq_gc_values, seq_callable_fracs, seq_callable_bps, seq_gap_bps)
                sequence_summary_rows.append({
                    "sequence_summary_pk": sequence_summary_pk,
                    "run_pk": run_pk,
                    "sequence_pk": sequence_pk,
                    "genome_pk": genome_pk,
                    "standard_window_size_bp": size_bp,
                    "step_size_bp": step_bp,
                    "tiling_type": tiling_type,
                    "mask_mode": mask_mode,
                    "n_windows_total": seq_total_windows,
                    "n_windows_kept": seq_kept_windows,
                    "n_windows_excluded_missing": seq_excluded_missing,
                    "n_windows_excluded_short": seq_excluded_short,
                    **seq_summary,
                })
                sequence_summary_pk += 1

            genome_summary = summarize_gc_values(genome_gc_values, genome_callable_fracs, genome_callable_bps, genome_gap_bps)
            genome_summary.pop("median_callable_fraction", None)
            genome_summary_rows.append({
                "genome_summary_pk": genome_summary_pk,
                "run_pk": run_pk,
                "genome_pk": genome_pk,
                "species_pk": genome_pk_to_species_pk.get(genome_pk, r"\N"),
                "standard_window_size_bp": size_bp,
                "step_size_bp": step_bp,
                "tiling_type": tiling_type,
                "mask_mode": mask_mode,
                "n_windows_total": genome_total_windows,
                "n_windows_kept": genome_kept_windows,
                **genome_summary,
                "seq_count_used": seq_count_used,
                "largest_seq_fraction": (largest_sequence_bp / genome_sequence_bp_used) if genome_sequence_bp_used > 0 else r"\N",
            })
            genome_summary_pk += 1

            for table, columns, rows_for_table in [
                ("genomic_windows", GENOMIC_WINDOWS_COLUMNS, genomic_window_rows),
                ("gc_window_stats", GC_WINDOW_STATS_COLUMNS, gc_window_rows),
                ("sequence_summary", SEQUENCE_SUMMARY_COLUMNS, sequence_summary_rows),
                ("genome_summary", GENOME_SUMMARY_COLUMNS, genome_summary_rows),
            ]:
                rel_path, n_rows = write_window_chunk(
                    output_dir=output_dir,
                    table=table,
                    columns=columns,
                    rows=rows_for_table,
                    run_pk=run_pk,
                    genome_pk=genome_pk,
                    accession=accession,
                    window_size_bp=size_bp,
                    step_size_bp=step_bp,
                )
                chunk_manifest_rows.append({
                    "table_name": table,
                    "chunk_path": rel_path,
                    "run_pk": run_pk,
                    "genome_pk": genome_pk,
                    "accession_id": accession,
                    "window_set_pk": current_window_set_pk,
                    "standard_window_size_bp": size_bp,
                    "step_size_bp": step_bp,
                    "n_rows": n_rows,
                })
                print(f"[INFO] Wrote {rel_path} ({n_rows} rows)", file=sys.stderr)
        genome.close()

    return window_set_rows, chunk_manifest_rows


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
    load_order = ["natural_history", "species_name_audit", "genomes", "sequences", "analysis_run", "window_set", "genomic_windows", "gc_window_stats", "sequence_summary", "genome_summary"]
    columns_by_table = {
        "natural_history": NATURAL_HISTORY_COLUMNS,
        "species_name_audit": SPECIES_NAME_AUDIT_COLUMNS,
        "genomes": GENOMES_COLUMNS,
        "sequences": SEQUENCES_COLUMNS,
        "analysis_run": ANALYSIS_RUN_COLUMNS,
        "window_set": WINDOW_SET_COLUMNS,
        "genomic_windows": GENOMIC_WINDOWS_COLUMNS,
        "gc_window_stats": GC_WINDOW_STATS_COLUMNS,
        "sequence_summary": SEQUENCE_SUMMARY_COLUMNS,
        "genome_summary": GENOME_SUMMARY_COLUMNS,
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
        description="Build starter SQL TSVs and compressed chunked window TSVs for GC3 dynamics."
    )
    parser.add_argument("--genomes-dir", default="genomes", help="Root genomes directory. Default: genomes")
    parser.add_argument(
        "--manifest",
        default="genomes/records/project_manifests/mass_predicts_dna_dynamics_with_s_punctatus_manifest.csv",
        help="Project manifest CSV.",
    )
    parser.add_argument(
        "--natural-history",
        default="genomes/records/sql_tsvs/natural_history.tsv",
        help="Input natural_history TSV.",
    )
    parser.add_argument(
        "--species-name-audit",
        default="genomes/records/sql_tsvs/species_name_audit.tsv",
        help="Input species_name_audit TSV.",
    )
    parser.add_argument(
        "--output-dir",
        default="genomes/records/sql_tsvs",
        help="Output directory for SQL TSVs. Default: genomes/records/sql_tsvs",
    )
    parser.add_argument(
        "--genomes-metadata",
        default="genomes/records/genomes_metadata.csv",
        help=(
            "CSV containing accession and exact genomic FASTA path. "
            "Default: genomes/records/genomes_metadata.csv. "
            "The FASTA path column is usually pathway_to_genome."
        ),
    )
    accession_group = parser.add_mutually_exclusive_group()
    accession_group.add_argument(
        "--accessions",
        nargs="+",
        default=None,
        help="Explicit accession IDs to include in genomes/sequences TSVs. Default: the built-in five accessions.",
    )
    accession_group.add_argument(
        "--n-accessions",
        type=int,
        default=None,
        help=(
            "Total number of accessions to include. Keeps the default five first, "
            "then adds additional accessions from the manifest in manifest order. "
            "Examples: 10 adds 5 more; 15 adds 10 more; 20 adds 15 more."
        ),
    )
    accession_group.add_argument(
        "--all-accessions",
        action="store_true",
        help="Include every accession found in the manifest, in manifest order.",
    )
    parser.add_argument("--genome-pk-start", type=int, default=1, help="Starting genome_pk. Default: 1")
    parser.add_argument("--sequence-pk-start", type=int, default=1, help="Starting sequence_pk. Default: 1")
    parser.add_argument(
        "--max-sequences-per-genome",
        type=int,
        default=None,
        help="Optional testing/debug limit. Omit to include every sequence in each FASTA.",
    )
    parser.add_argument(
        "--skip-windows",
        action="store_true",
        help="Skip creation of analysis_run/window_set/genomic_windows/gc_window_stats/sequence_summary/genome_summary TSVs.",
    )
    parser.add_argument(
        "--window-sizes-bp",
        nargs="+",
        type=int,
        default=[300000, 100000, 20000, 5000],
        help="Window sizes to build, in bp. Default: 300000, 100000, 20000, 5000.",
    )
    parser.add_argument(
        "--step-sizes-bp",
        nargs="+",
        type=int,
        default=None,
        help="Optional step sizes in bp. Must match --window-sizes-bp length. Default: same as window size for non-overlapping windows.",
    )
    parser.add_argument(
        "--window-sequence-types",
        nargs="+",
        default=["chromosome"],
        help="Sequence types to include in window tables. Default: chromosome.",
    )
    parser.add_argument("--min-callable-frac", type=float, default=0.8, help="Minimum A/C/G/T fraction for keep_flag=1. Default: 0.8")
    parser.add_argument("--run-pk", type=int, default=1, help="analysis_run.run_pk for generated window tables. Default: 1")
    parser.add_argument("--analysis-name", default="gc_decay_by_window_size", help="analysis_run.analysis_name. Default: gc_decay_by_window_size")
    parser.add_argument("--mask-mode", default="raw_fasta_N_excluded", help="analysis_run/window summary mask_mode. Default: raw_fasta_N_excluded")
    parser.add_argument("--tiling-type", default="non_overlapping", help="window_set tiling_type. Default: non_overlapping")
    parser.add_argument("--start-offset-bp", type=int, default=0, help="0-based offset for first window start. Default: 0")
    parser.add_argument("--seq-scope", default="chromosome", help="window_set seq_scope label. Default: chromosome")
    parser.add_argument("--load-sql", action="store_true", help="Load generated TSVs into MySQL after writing them.")
    parser.add_argument("--mysql-db", default="gc3_dynamics", help="MySQL database name for --load-sql.")
    parser.add_argument("--mysql-user", default="root", help="MySQL user for --load-sql.")
    parser.add_argument("--mysql-host", default="localhost", help="MySQL host for --load-sql.")
    parser.add_argument("--mysql-port", type=int, default=3306, help="MySQL port for --load-sql.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    genomes_dir = Path(args.genomes_dir).resolve()
    manifest_path = Path(args.manifest)
    natural_history_path = Path(args.natural_history)
    species_name_audit_path = Path(args.species_name_audit)
    output_dir = Path(args.output_dir)
    genomes_metadata_path = Path(args.genomes_metadata) if args.genomes_metadata else None

    # Resolve relative paths against genomes_dir, not the shell working directory.
    if not manifest_path.is_absolute():
        manifest_path = genomes_dir / "records/project_manifests" / manifest_path.name

    if not natural_history_path.is_absolute():
        natural_history_path = genomes_dir / "records/sql_tsvs" / natural_history_path.name

    if not species_name_audit_path.is_absolute():
        species_name_audit_path = genomes_dir / "records/sql_tsvs" / species_name_audit_path.name

    if not output_dir.is_absolute():
        output_dir = genomes_dir / "records/sql_tsvs"

    genomes_metadata_path = find_genomes_metadata(genomes_metadata_path, genomes_dir)

    manifest_path = manifest_path.resolve()
    natural_history_path = natural_history_path.resolve()
    species_name_audit_path = species_name_audit_path.resolve()
    output_dir = output_dir.resolve()

    manifest_rows = read_table(manifest_path, delimiter=",")
    genomes_metadata_rows = read_table(genomes_metadata_path, delimiter=",")
    natural_history_rows = read_table(natural_history_path, delimiter="\t")
    species_name_audit_rows = read_table(species_name_audit_path, delimiter="\t")

    selected_accessions = choose_accessions(
        manifest_rows,
        requested_accessions=args.accessions,
        n_accessions=args.n_accessions,
        all_accessions=args.all_accessions,
    )
    selected_manifest_rows = filter_manifest_rows(manifest_rows, selected_accessions)
    print(f"[INFO] Selected {len(selected_accessions)} accession(s) for genomes/sequences TSVs.", file=sys.stderr)

    species_pk_lookup = build_species_pk_lookup(natural_history_rows)

    genomes_rows, accession_to_genome_pk = build_genomes_rows(
        selected_manifest_rows,
        species_pk_lookup,
        start_pk=args.genome_pk_start,
    )

    accession_to_fasta_path = build_fasta_path_lookup(genomes_metadata_rows, genomes_dir)

    missing_fasta_accessions = [
        accession for accession in selected_accessions
        if not any(alias in accession_to_fasta_path for alias in accession_aliases(accession))
    ]
    if missing_fasta_accessions:
        raise FileNotFoundError(
            "No usable genomic FASTA path was found in genomes_metadata.csv for these selected accession(s): "
            + ", ".join(missing_fasta_accessions[:25])
            + (" ..." if len(missing_fasta_accessions) > 25 else "")
        )

    sequences_rows = build_sequences_rows(
        selected_manifest_rows,
        accession_to_genome_pk,
        fasta_path_lookup=accession_to_fasta_path,
        start_pk=args.sequence_pk_start,
        max_sequences_per_genome=args.max_sequences_per_genome,
    )
    sequence_type_audit_rows = build_sequence_type_audit_rows(
        sequences_rows,
        genomes_rows,
        genomes_dir=genomes_dir,
    )

    analysis_run_rows: List[Dict[str, object]] = []
    window_set_rows: List[Dict[str, object]] = []
    window_chunk_manifest_rows: List[Dict[str, object]] = []

    if not args.skip_windows:
        analysis_run_rows = build_analysis_run_rows(
            run_pk=args.run_pk,
            analysis_name=args.analysis_name,
            mask_mode=args.mask_mode,
            min_callable_frac=args.min_callable_frac,
            notes=(
                "Starter GC window analysis for visualizing and modeling decay in GC variation "
                "as window size increases; intended for Zuur-style exploration and mixed models "
                "with mass, chromosome size, and other predictors."
            ),
        )
        (
            window_set_rows,
            window_chunk_manifest_rows,
        ) = build_window_tables(
            selected_manifest_rows,
            genomes_rows,
            sequences_rows,
            fasta_path_lookup=accession_to_fasta_path,
            run_pk=args.run_pk,
            window_sizes_bp=args.window_sizes_bp,
            step_sizes_bp=args.step_sizes_bp,
            sequence_types=args.window_sequence_types,
            min_callable_frac=args.min_callable_frac,
            tiling_type=args.tiling_type,
            start_offset_bp=args.start_offset_bp,
            seq_scope=args.seq_scope,
            mask_mode=args.mask_mode,
            output_dir=output_dir,
        )

    written = {
        "natural_history.tsv": write_tsv(output_dir / "natural_history.tsv", NATURAL_HISTORY_COLUMNS, natural_history_rows),
        "species_name_audit.tsv": write_tsv(output_dir / "species_name_audit.tsv", SPECIES_NAME_AUDIT_COLUMNS, species_name_audit_rows),
        "genomes.tsv": write_tsv(output_dir / "genomes.tsv", GENOMES_COLUMNS, genomes_rows),
        "sequences.tsv": write_tsv(output_dir / "sequences.tsv", SEQUENCES_COLUMNS, sequences_rows),
        "sequence_type_audit.tsv": write_tsv(
            output_dir / "sequence_type_audit.tsv",
            SEQUENCE_TYPE_AUDIT_COLUMNS,
            sequence_type_audit_rows,
        ),
    }

    if not args.skip_windows:
        written.update({
            "analysis_run.tsv": write_tsv(output_dir / "analysis_run.tsv", ANALYSIS_RUN_COLUMNS, analysis_run_rows),
            "window_set.tsv": write_tsv(output_dir / "window_set.tsv", WINDOW_SET_COLUMNS, window_set_rows),
            "window_chunk_manifest.tsv": write_tsv(
                output_dir / "window_chunk_manifest.tsv",
                WINDOW_CHUNK_MANIFEST_COLUMNS,
                window_chunk_manifest_rows,
            ),
        })

    print("\n[OK] Wrote SQL TSVs:")
    for filename, n_rows in written.items():
        print(f"  {output_dir / filename}\t{n_rows} data rows")

    review_rows = [row for row in sequence_type_audit_rows if row["audit_status"] != "pass"]
    if review_rows:
        print("\n[WARN] Sequence type audit rows needing review:", file=sys.stderr)
        for row in review_rows:
            print(
                f"  {row['accession_id']}\tstatus={row['audit_status']}\t"
                f"chromosomes={row['n_chromosome_inferred']}/{row['ncbi_total_chromosomes']}\t"
                f"components={row['n_fasta_sequences']}/{row['ncbi_component_sequences']}",
                file=sys.stderr,
            )

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
