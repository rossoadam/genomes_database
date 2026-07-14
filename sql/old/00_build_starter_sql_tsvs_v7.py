#!/usr/bin/env python3
"""
Build starter SQL TSVs for the GC3 / genome dynamics database.

This script creates four load-ready TSVs:
  1. natural_history.tsv
  2. species_name_audit.tsv
  3. genomes.tsv
  4. sequences.tsv
  5. sequence_type_audit.tsv

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
import json
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
    """Reject obvious CDS/proteome/transcript FASTAs and accept genome FASTAs."""
    path_text = str(path).lower()
    name = path.name.lower()

    reject_terms = [
        "proteome",
        "protein",
        "pep",
        "peptide",
        "cds",
        "rna",
        "transcript",
        "mrna",
        "/old/",
        "/old",
    ]
    if any(term in path_text for term in reject_terms):
        return False

    if not looks_like_fasta_path(path):
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
    """Build accession -> exact genomic FASTA path lookup from genomes_metadata.csv."""
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
        "pathway_to_genome",
        "path_to_genome",
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
            if any(token in col.lower() for token in ["path", "genome", "fna", "fasta"])
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
        if not accession or accession in lookup:
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

        lookup[accession] = chosen

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
    fasta_path = fasta_path_lookup.get(accession)
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
        if accession not in accession_to_fasta_path
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
