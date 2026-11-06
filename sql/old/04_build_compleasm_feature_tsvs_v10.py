#!/usr/bin/env python3
"""
04_build_compleasm_feature_tsvs_v10.py

Build Compleasm-derived feature TSVs compatible with gc3_dynamics_v13_reduced.txt.

Required inputs:
    --genomes       Path to the genomes directory.
    --manifest      Project manifest CSV/TSV with an accession column.
    --orthodb       Path to sauropsida_odb12 *.tar.gz archive or extracted directory.
    --sequences-tsv Base sequences lookup table containing sequence_pk and sequence_id.
    --genomes-tsv   Base genomes lookup table containing genome_pk/accession_id/species_pk.

Default supporting inputs expected below --genomes:
    records/compleasm/records/metadata.csv
    records/genomes_metadata.csv

Outputs:
    sauropsida_odb12.tsv
    orthologs.tsv
    ortholog_summary.tsv
    intron_compleasm.tsv
    intron_compleasm_summary.tsv
    flanks_compleasm.tsv
    flank_sets_compleasm.tsv
    flank_compleasm_summary.tsv
    ortholog_source_audit.tsv (helper file; not a SQL-schema table)

Notes:
    - Genomic feature coordinates are 1-based inclusive. The source audit
      preserves Compleasm coordinate fields as parsed, including zero starts.
    - Introns are inferred from Compleasm full_table Codons/CDS coordinate tokens.
    - Flanks are clipped at sequence ends and may overlap neighboring genes.
    - GC ignores N/ambiguous bases in the denominator.

Version v5 updates:
    - Include QC-passing Duplicated (or Duplicate) orthologs in summaries,
      introns, and flanks for the exact accession GWHFILM00000000.1
      (triploid Indotyphlops braminus). Other accessions remain Single-only.
    - Preserve raw CDS QC, the reported Compleasm status, and TSV schemas.
      Per-ortholog GC/GC3/GC4 calculations already accept any QC-passing CDS.
    - Compatible with existing validation TSVs and the additional accession
      policy column written by 03_raw_ortholog_validity_batch_v3.py.
    - Retain the existing one-record-per-OrthoDB-ID input handling. Multiple
      copies are not expanded or selected by this exception.

Version v6 updates:
    - Retain ortholog rows and their QC-passing CDS-based GC/GC3/GC4 values
      when their Compleasm input sequence cannot be linked to sequences.tsv.
      Genomic linkage is not a requirement for calculating CDS composition.
    - Leave sequence_pk, strand, start, and end blank for unlinked orthologs.
      Transcript/input-sequence coordinates are not genomic coordinates.
    - Write ortholog_source_audit.tsv for every ortholog, preserving accession,
      genome_pk, original sequence name, original strand/coordinates/Codons,
      mapping reason, and inclusion in CDS summaries and genomic features.
      This audit is a helper file, not an extra SQL table. Version v6 retained
      the v8 schema columns; empty genomic fields must import as SQL NULL.
    - Continue CDS summaries for QC-passing status-eligible orthologs even
      when genomic mapping is unavailable. Skip genomic introns/flanks for
      those records and report CDS-summary and mapping counts separately.
    - No new GFF input, transcript remapping, or artificial sequence IDs.
      sequences.tsv is not modified. The v5 accession exception is retained.
    - In the v8 schema, unlinked orthologs cannot join to genomes through
      sequences; their accession/genome association is preserved in the audit.

Version v8 updates:
    - Preserve v7 feature calculations while separating initial and run QC.
    - Generate ortholog_initial_qc.tsv and ortholog_run_qc.tsv with SQL NULLs.
    - Require strict exact QC mappings and explicit per-sequence final run QC.
    - Audit actual nhPhyML input sequences without inferring membership from QC.
    - Stage outputs and atomically switch output directories on success.

Version v9 updates:
    - Represent all missing values in newly generated TSV files as the literal
      SQL import marker \\N, rather than empty tab-separated cells.
    - Use one scalar serializer for regular and streaming TSV writers, audit
      records, and failure diagnostics, retaining nonmissing values unchanged.
    - Validate the SQL NULL convention and TSV widths/headers on the staging
      directory before its atomic replacement of previous output files.
    - Preserve v8 record enumeration, all ortholog IDs, QC eligibility, GC
      calculations, input rules, output schemas and relational validations.

Version v10 updates:
    - Collapse identical repeated Compleasm full_table records into one ortholog.
    - Retain distinct, QC-ineligible genomic copies as a single unresolved
      ortholog identity without copying any genomic location or CDS GC values.
    - Fail closed for ambiguous eligible orthologs, conflicting source statuses,
      and excluded duplicates contradicted by initial or run-specific QC.
    - Preserve each duplicate source line in a separate copy-level audit and
      retain that audit outside staging on duplicate-related failures.
    - Preserve v9 calculations, SQL schemas, SQL NULLs, and staged replacement.

Version v7 updates:
    - Add genome_pk immediately after sequence_pk in orthologs.tsv, matching
      the column order in gc3_dynamics_9_reduced.txt.
    - Populate genome_pk for every ortholog row from the exact accession
      match in --genomes-tsv, independently of sequence_pk. This provides a
      direct genome/species link even when the genomic location is unknown.
    - Retain v6 CDS composition, QC, summaries, and source auditing. Unmapped
      orthologs still have blank sequence_pk/strand/start/end and retain their
      QC-passing CDS GC/GC3/GC4 values. Existing CLI arguments are unchanged.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import tarfile
from dataclasses import dataclass
from collections import Counter, defaultdict
import os
import shutil
import tempfile
from datetime import datetime
from pathlib import Path
from statistics import mean, median
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pyfaidx
from Bio.SeqIO.FastaIO import SimpleFastaParser


STOP_CODONS = {"TAA", "TAG", "TGA"}
DUPLICATED_STATUS_ALLOWED_ACCESSIONS = frozenset({"GWHFILM00000000.1"})


def status_allowed_for_features(accession: str, status: str) -> bool:
    """Preserve Single-only filtering except for the specified triploid assembly."""
    return status == "Single" or (
        str(accession).strip() in DUPLICATED_STATUS_ALLOWED_ACCESSIONS
        and str(status).strip().lower() in {"duplicated", "duplicate"}
    )


def describe_sequence_mapping(record, sequence_pk):
    """Describe exact mapping, distinguishing deliberate copy ambiguity."""
    if record is None:
        return (
            "missing_full_table_record", "unknown",
            "No matching Compleasm full_table record; genomic location unavailable.",
        )
    if record.mapping_state != "resolved":
        return (
            record.mapping_state, "unresolved_genomic_copies",
            record.ambiguity_reason or
            "Multiple distinct QC-ineligible genomic copies; no copy selected. "
            "Source-copy coordinates are recorded only in the duplicate audit.",
        )
    if sequence_pk is not None:
        return "mapped_genomic_sequence", "genomic_sequence", "Exact genome_pk/sequence_id match in sequences.tsv."
    if "|rna-" in (record.sequence_id or "").lower():
        return (
            "unmapped_transcript_sequence", "transcript_sequence",
            "Compleasm names a transcript/CDS input sequence; no genomic sequence_pk match. "
            "CDS composition is retained when QC passes; genomic introns/flanks are unavailable.",
        )
    return (
        "unmapped_sequence_id", "compleasm_input_sequence_unresolved",
        "Compleasm sequence_id has no exact match for this genome_pk in sequences.tsv. "
        "Its coordinates are preserved only in this audit until genomic mapping is established.",
    )


FLANK_SPECS = [
    (1, "50_up_50_down", 50, 50),
    (2, "100_up_100_down", 100, 100),
    (3, "200_up_200_down", 200, 200),
]


@dataclass
class FullRecord:
    odb12_id: str
    status: str
    sequence_id: Optional[str]
    strand: Optional[str]
    start: Optional[int]
    end: Optional[int]
    exons: List[Dict[str, object]]
    mapping_state: str = "resolved"
    ambiguity_reason: Optional[str] = None


def root_acc(accession: str) -> str:
    accession = str(accession).strip()
    return accession.split(".")[0] if accession else accession


def accession_version(accession: str) -> int:
    accession = str(accession).strip()
    if "." not in accession:
        return -1
    suffix = accession.rsplit(".", 1)[-1]
    return int(suffix) if suffix.isdigit() else -1


def read_delimited(path: Path) -> pd.DataFrame:
    """Read CSV/TSV by sniffing the delimiter from the suffix and first line."""
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    if suffix == ".csv":
        return pd.read_csv(path)
    with open(path, "r", newline="") as handle:
        sample = handle.readline()
    sep = "\t" if "\t" in sample else ","
    return pd.read_csv(path, sep=sep)

def parse_bool(value) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes", "y", "pass", "passed"}

def repath(path_to_fix, genomes_dir: Path) -> Path:
    path_str = str(path_to_fix).strip()
    if not path_str or path_str.lower() == "nan":
        return Path("")
    path_obj = Path(path_str)
    if path_obj.exists():
        return path_obj
    root_name = genomes_dir.name
    if root_name in path_str:
        suffix = path_str.split(root_name, 1)[-1].lstrip("/\\")
        candidate = genomes_dir / suffix
        if candidate.exists():
            return candidate
    return path_obj

def now_str() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def load_genome(input_fasta: Path) -> pyfaidx.Fasta:
    try:
        return pyfaidx.Fasta(str(input_fasta), as_raw=True, build_index=True)
    except pyfaidx.FaidxException as e:
        raise RuntimeError(f"Failed to load genome FASTA with pyfaidx: {input_fasta}; {e}") from e


def reverse_complement(seq: str) -> str:
    complement = str.maketrans(
        "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
        "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
    )
    return str(seq).translate(complement)[::-1]


def calculate_gc(sequence: str) -> Tuple[Optional[float], int, int, int, int, int, int]:
    sequence = str(sequence).upper()
    g_count = sequence.count("G")
    c_count = sequence.count("C")
    a_count = sequence.count("A")
    t_count = sequence.count("T")
    n_count = sequence.count("N")
    valid = a_count + t_count + g_count + c_count
    gc = None if valid == 0 else (g_count + c_count) / valid
    return gc, g_count, c_count, a_count, t_count, n_count, valid


def remove_terminal_stop(sequence: str) -> str:
    sequence = str(sequence).upper().replace("-", "")
    if len(sequence) >= 3 and sequence[-3:] in STOP_CODONS:
        return sequence[:-3]
    return sequence


def check_cds_sequence(sequence: str) -> Tuple[bool, int, int, bool]:
    sequence = str(sequence).upper().replace("-", "")
    length_mod_3 = len(sequence) % 3
    if length_mod_3 != 0:
        return False, length_mod_3, 0, False
    codons = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
    terminal_stop = bool(codons and codons[-1] in STOP_CODONS)
    internal_stops = sum(1 for codon in codons[:-1] if codon in STOP_CODONS)
    return length_mod_3 == 0 and internal_stops == 0, length_mod_3, internal_stops, terminal_stop


def calculate_gc3(sequence: str) -> Optional[float]:
    original = str(sequence).upper().replace("-", "")
    passes, length_mod_3, internal_stops, _terminal_stop = check_cds_sequence(original)
    if not passes:
        if length_mod_3 != 0:
            raise ValueError(f"Sequence length ({len(original)}) is not divisible by 3")
        if internal_stops > 0:
            raise ValueError(f"Sequence contains {internal_stops} internal stop codon(s)")
    sequence = remove_terminal_stop(original)
    arr = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
    third = arr[2::3]
    A, C, G, T = map(ord, "ACGT")
    valid = (third == A) | (third == C) | (third == G) | (third == T)
    gc3_seq = third[valid].tobytes().decode("ascii")
    return calculate_gc(gc3_seq)[0]


def calculate_gc4(sequence: str) -> Optional[float]:
    sequence = remove_terminal_stop(str(sequence).upper().replace("-", ""))
    if len(sequence) % 3 != 0:
        raise ValueError(f"Sequence length ({len(sequence)}) is not divisible by 3")
    arr = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
    n1, n2, n3 = arr[0::3], arr[1::3], arr[2::3]
    A, C, G, T = map(ord, "ACGT")
    valid_third = (n3 == A) | (n3 == C) | (n3 == G) | (n3 == T)
    candidate = (
        ((n1 == G) & (n2 == T)) |  # Val
        ((n1 == C) & (n2 == C)) |  # Pro
        ((n1 == A) & (n2 == C)) |  # Thr
        ((n1 == G) & (n2 == C)) |  # Ala
        ((n1 == G) & (n2 == G)) |  # Gly
        ((n1 == C) & (n2 == T)) |  # Leu CTN
        ((n1 == C) & (n2 == G)) |  # Arg CGN
        ((n1 == T) & (n2 == C))    # Ser TCN
    )
    gc4_seq = n3[candidate & valid_third].tobytes().decode("ascii")
    return calculate_gc(gc4_seq)[0]


def write_tsv(rows: List[Dict[str, object]], path: Path, columns: List[str]) -> None:
    """Write rows through the same SQL-NULL serializer as streaming outputs."""
    handle, writer = tsv_writer(path, columns)
    try:
        for row in rows:
            writer.writerow(row)
    finally:
        handle.close()


def find_col(columns: Iterable[str], candidates: Sequence[str], required: bool = True) -> Optional[str]:
    lower = {c.lower().strip(): c for c in columns}
    for candidate in candidates:
        key = candidate.lower().strip()
        if key in lower:
            return lower[key]
    compact = {c.lower().replace(" ", "").replace("_", ""): c for c in columns}
    for candidate in candidates:
        key = candidate.lower().replace(" ", "").replace("_", "")
        if key in compact:
            return compact[key]
    if required:
        raise ValueError(f"Missing required column. Tried: {', '.join(candidates)}")
    return None


def load_manifest_accessions(manifest_path: Path) -> List[str]:
    df = read_delimited(manifest_path)
    col = find_col(df.columns, ["accession", "accession_id", "assembly_accession"])
    accessions = df[col].dropna().astype(str).str.strip()
    return sorted(set(x for x in accessions if x))

def load_genomes_lookup(genomes_tsv: Path) -> Dict[str, Dict[str, object]]:
    df = read_delimited(genomes_tsv)
    accession_col = find_col(df.columns, ["accession_id", "accession", "assembly_accession"])
    genome_pk_col = find_col(df.columns, ["genome_pk"])
    species_pk_col = find_col(df.columns, ["species_pk"], required=False)

    lookup = {}

    for _, row in df.iterrows():
        accession = str(row[accession_col]).strip()
        if not accession or accession.lower() == "nan":
            continue

        lookup[accession] = {
            "genome_pk": int(row[genome_pk_col]),
            "accession_id": accession,
            "species_pk": int(row[species_pk_col]) if species_pk_col and not pd.isna(row[species_pk_col]) else None,
        }

    return lookup

def load_sequences_lookup(sequences_tsv: Path) -> Dict[Tuple[int, str], int]:
    df = read_delimited(sequences_tsv)
    sequence_pk_col = find_col(df.columns, ["sequence_pk"])
    genome_pk_col = find_col(df.columns, ["genome_pk"])
    sequence_id_col = find_col(df.columns, ["sequence_id", "seq_id", "chromosome", "scaffold"])
    lookup: Dict[Tuple[int, str], int] = {}
    for _, row in df.iterrows():
        if pd.isna(row[genome_pk_col]) or pd.isna(row[sequence_pk_col]) or pd.isna(row[sequence_id_col]):
            continue
        lookup[(int(row[genome_pk_col]), str(row[sequence_id_col]).strip())] = int(row[sequence_pk_col])
    return lookup

def resolve_metadata_rows(metadata_path: Path, manifest_accessions: List[str]) -> Tuple[pd.DataFrame, List[str]]:
    metadata_df = read_delimited(metadata_path)
    accession_col = find_col(metadata_df.columns, ["accession", "accession_id", "assembly_accession"])
    cds_col = find_col(metadata_df.columns, ["cds_fasta"])
    full_col = find_col(metadata_df.columns, ["full_table", "full"])

    metadata_df = metadata_df.rename(
        columns={accession_col: "accession", cds_col: "cds_fasta", full_col: "full_table"}
    ).copy()

    metadata_df["accession"] = metadata_df["accession"].astype(str).str.strip()
    allowed = set(str(a).strip() for a in manifest_accessions)

    resolved = metadata_df[metadata_df["accession"].isin(allowed)].copy().reset_index(drop=True)
    missing = sorted(allowed - set(resolved["accession"]))

    return resolved, missing

def attach_genome_fna(metadata_rows: pd.DataFrame, genomes_metadata_path: Path, genomes_dir: Path) -> pd.DataFrame:
    gm = read_delimited(genomes_metadata_path)
    accession_col = find_col(gm.columns, ["accession", "accession_id", "assembly_accession"])
    fna_col = find_col(gm.columns, ["path_to_fna", "genome_fasta", "fna_path", "genomic_fna", "path_to_genome"])

    gm = gm.rename(columns={accession_col: "genome_metadata_accession", fna_col: "path_to_fna"}).copy()
    gm["genome_metadata_accession"] = gm["genome_metadata_accession"].astype(str).str.strip()

    out = metadata_rows.merge(
        gm,
        left_on="accession",
        right_on="genome_metadata_accession",
        how="left",
        suffixes=("", "_genome_metadata"),
    )

    if out["path_to_fna"].isna().any():
        missing = out.loc[out["path_to_fna"].isna(), "accession"].tolist()
        raise ValueError(
            "Missing exact accession matches in genomes_metadata.csv for: "
            + ", ".join(missing[:20])
        )

    out["path_to_fna"] = out["path_to_fna"].map(lambda x: repath(x, genomes_dir))
    return out

def parse_orthodb_archive(orthodb_path: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    if orthodb_path.is_dir():
        candidates = list(orthodb_path.rglob("links_to_ODB12.txt"))
        if not candidates:
            raise FileNotFoundError(f"links_to_ODB12.txt not found below {orthodb_path}")
        with open(candidates[0], "r") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    rows.append({"odb12_id": parts[0], "name": parts[1], "orthodb": parts[2]})
        return rows
    if not tarfile.is_tarfile(orthodb_path):
        raise ValueError(f"--orthodb must be a tar archive or extracted directory: {orthodb_path}")
    with tarfile.open(orthodb_path, "r:*") as tar:
        member = next((m for m in tar.getmembers() if m.name.endswith("links_to_ODB12.txt")), None)
        if member is None:
            raise FileNotFoundError(f"links_to_ODB12.txt not found inside {orthodb_path}")
        extracted = tar.extractfile(member)
        if extracted is None:
            raise FileNotFoundError(f"Could not read {member.name} inside {orthodb_path}")
        for raw in extracted:
            line = raw.decode("utf-8").rstrip("\n")
            parts = line.split("\t")
            if len(parts) >= 3:
                rows.append({"odb12_id": parts[0], "name": parts[1], "orthodb": parts[2]})
    return rows


def find_validity_tsv(cds_fasta: Path) -> Path:
    parent = cds_fasta.parent
    candidates = sorted(parent.glob("*ortholog_validity.tsv")) or sorted(parent.glob("*raw_ortholog_validity.tsv"))
    if not candidates:
        raise FileNotFoundError(f"No *ortholog_validity.tsv found in {parent}")
    nonsummary = [p for p in candidates if "summary" not in p.name.lower()]
    return nonsummary[0] if nonsummary else candidates[0]


def load_validity_rows(validity_tsv: Path) -> Dict[str, bool]:
    """
    Legacy helper for one per-genome validity TSV. Kept for compatibility,
    but v4 uses --ortholog-validation and load_ortholog_validation_lookup().
    """
    df = read_delimited(validity_tsv)
    seq_col = find_col(df.columns, ["sequence_id", "Gene", "odb12_id", "busco_id"])
    qc_col = find_col(df.columns, ["passes_raw_cds_qc"])
    return {str(row[seq_col]).strip(): parse_bool(row[qc_col]) for _, row in df.iterrows() if str(row[seq_col]).strip()}


def load_ortholog_validation_lookup(ortholog_validation_tsv: Path) -> Dict[str, Dict[str, bool]]:
    """
    Load the project-wide ortholog validation table.

    Expected columns:
        accession
        sequence_id
        passes_raw_cds_qc

    Returns:
        {
            accession: {
                sequence_id: passes_qc
            }
        }

    Matching is performed using exact accession versions only.
    """
    df = read_delimited(ortholog_validation_tsv)

    accession_col = find_col(df.columns, ["accession", "accession_id"])
    seq_col = find_col(df.columns, ["sequence_id", "Gene", "odb12_id", "busco_id"])
    qc_col = find_col(df.columns, ["passes_raw_cds_qc"])

    lookup: Dict[str, Dict[str, bool]] = {}

    for _, row in df.iterrows():

        accession = str(row[accession_col]).strip()
        if not accession or accession.lower() == "nan":
            continue

        sequence_id = str(row[seq_col]).strip()
        if not sequence_id or sequence_id.lower() == "nan":
            continue

        passes_qc = parse_bool(row[qc_col])

        lookup.setdefault(accession, {})[sequence_id] = passes_qc

    return lookup


def load_cds_sequences(cds_fasta: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    with open(cds_fasta, "r") as handle:
        for header, seq in SimpleFastaParser(handle):
            odb12_id = header.split()[0]
            if not odb12_id:
                raise ValueError(f'Empty CDS FASTA identifier in {cds_fasta}')
            if odb12_id in seqs:
                raise ValueError(f'Duplicate CDS FASTA identifier {odb12_id} in {cds_fasta}; cannot select one copy silently')
            seqs[odb12_id] = seq.upper().replace("-", "")
    return seqs


def parse_coordinate_token(token: str, odb12_id: str) -> Optional[Dict[str, object]]:
    token = str(token).strip()
    if not token:
        return None
    parts = token.split("_")
    if len(parts) != 3:
        return None
    try:
        start = int(parts[0])
        end = int(parts[1])
    except ValueError:
        return None
    if start > end:
        start, end = end, start
    return {"odb12_id": odb12_id, "start": start, "end": end, "strand": parts[2]}


# Auxiliary report; deliberately NOT part of the gc3_dynamics_v13 SQL schema.
DUPLICATE_AUDIT_COLUMNS = [
    'accession', 'genome_pk', 'odb12_id', 'source_file', 'source_line',
    'compleasm_status', 'source_sequence_id', 'source_start', 'source_end',
    'source_strand', 'source_score', 'source_length', 'source_identity',
    'source_codon_intervals', 'n_records_for_ortholog', 'qc_eligible',
    'resolution_status', 'resolution_reason', 'genomic_copy_selected',
]
# Kept in memory so the failure handler can persist rejected groups even when
# the staged directory is subsequently removed. The audit only contains
# repeated groups, not all ~500,000 single-source orthologs.
DUPLICATE_AUDIT_ROWS = []
DUPLICATE_COUNTS = Counter()


def parse_full_table(full_table: Path, accession: str = None,
                     genome_pk: int = None, initial_qc_keys=None,
                     duplicate_audit_rows=None, duplicate_counts=None) -> Dict[str, FullRecord]:
    """Resolve repeated Compleasm IDs without selecting a distinct genomic copy.

    A normal one-record group retains v9 parsing semantics. Identical parsed
    records are collapsed. Distinct genomic copies can be represented by one
    unresolved ortholog *only* if they are QC-ineligible and have no matching
    initial QC record. All ambiguous eligible groups fail closed.

    Source line numbers are physical 1-based lines, including comments/header.
    For compatibility with v9, rows too short to hold its essential positional
    fields are ignored; no fields of a valid parsed record are silently used
    to break ties among distinct copies.
    """
    if accession is None:
        raise ValueError('parse_full_table requires an exact accession to evaluate QC eligibility')
    if genome_pk is None:
        raise ValueError('parse_full_table requires genome_pk for the audit')
    if initial_qc_keys is None:
        initial_qc_keys = set()
    if duplicate_audit_rows is None:
        duplicate_audit_rows = DUPLICATE_AUDIT_ROWS
    if duplicate_counts is None:
        duplicate_counts = DUPLICATE_COUNTS

    with open(full_table, 'r', encoding='utf-8-sig') as handle:
        header: Optional[List[str]] = None
        raw_rows = []
        for line_num, line in enumerate(handle, 1):
            line = line.rstrip('\r\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if header is None and parts[0].lower().strip() == 'gene':
                header = parts
                continue
            raw_rows.append((line_num, parts))
    if header is None:
        # Preserve the positional fallback of v9 for older Compleasm versions.
        header = ["Gene", "Status", "Sequence", "Score", "Length", "Strand",
                  "Col7", "Col8", "Col9", "Gene Start", "Gene End", "Col12", "Codons"]
    colmap = {name.lower().strip(): i for i, name in enumerate(header)}

    def idx(*names: str, default: Optional[int] = None) -> Optional[int]:
        for name in names:
            key = name.lower().strip()
            if key in colmap:
                return colmap[key]
        compact = {h.lower().replace(' ', '').replace('_', ''): i
                   for i, h in enumerate(header)}
        for name in names:
            key = name.lower().replace(' ', '').replace('_', '')
            if key in compact:
                return compact[key]
        return default

    gene_i = idx('Gene', default=0)
    status_i = idx('Status', default=1)
    sequence_i = idx('Sequence', 'Chromosome', 'Contig', default=2)
    strand_i = idx('Strand', default=5)
    start_i = idx('Gene Start', 'Start', default=9)
    end_i = idx('Gene End', 'End', default=10)
    codons_i = idx('Codons', 'Coordinates', 'Exons', default=12)
    score_i = idx('Score', default=3)
    length_i = idx('Length', default=4)
    identity_i = idx('Identity', 'Percent Identity', 'Sequence Identity',
                     'Identity Percentage', 'Identity %')

    def cell(parts, i):
        if i is None or i >= len(parts):
            return None
        v = parts[i].strip()
        return v if v else None

    groups = defaultdict(list)
    for line_num, parts in raw_rows:
        # Retain v9's treatment of short non-genomic entries (e.g., Missing).
        if len(parts) <= max(gene_i or 0, status_i or 0, sequence_i or 0, strand_i or 0):
            continue
        odb12_id = (cell(parts, gene_i) or '')
        if not odb12_id:
            continue
        status = cell(parts, status_i) or ''
        sequence_id = cell(parts, sequence_i) or ''
        strand = cell(parts, strand_i) or ''
        exons = []
        codon_text = cell(parts, codons_i)
        if codon_text:
            for token in codon_text.split('|'):
                exon = parse_coordinate_token(token, odb12_id)
                if exon and (not strand or exon['strand'] == strand):
                    exons.append(exon)
        exons.sort(key=lambda e: int(e['start']))
        start = end = None
        try:
            if cell(parts, start_i):
                start = int(float(cell(parts, start_i)))
            if cell(parts, end_i):
                end = int(float(cell(parts, end_i)))
        except ValueError:
            start = end = None
        if start is None and exons:
            start = min(int(e['start']) for e in exons)
        if end is None and exons:
            end = max(int(e['end']) for e in exons)
        parsed = FullRecord(odb12_id, status, sequence_id, strand, start, end, exons)
        copy_audit = {
            'accession': accession, 'genome_pk': genome_pk, 'odb12_id': odb12_id,
            'source_file': str(full_table), 'source_line': line_num,
            'compleasm_status': status, 'source_sequence_id': sequence_id,
            'source_start': cell(parts, start_i),
            'source_end': cell(parts, end_i), 'source_strand': strand,
            'source_score': cell(parts, score_i),
            'source_length': cell(parts, length_i),
            'source_identity': cell(parts, identity_i),
            'source_codon_intervals': codon_text,
        }
        groups[odb12_id].append((parsed, copy_audit))
        duplicate_counts['source_records'] += 1

    records = {}
    valid_statuses = {'Single', 'Duplicated', 'Duplicate', 'Fragmented', 'Missing'}
    for odb12_id, copies in groups.items():
        if len(copies) == 1:
            records[odb12_id] = copies[0][0]
            continue

        duplicate_counts['repeated_groups'] += 1
        n = len(copies)
        statuses = {r.status for r, _ in copies}
        identical = all(r == copies[0][0] for r, _ in copies[1:])
        has_qc = (accession, odb12_id) in initial_qc_keys
        if not statuses.issubset(valid_statuses):
            resolution = 'invalid_status_error'
            reason = f'Missing or uninterpretable Compleasm statuses: {sorted(statuses)!r}'
        elif len(statuses) != 1:
            resolution = 'conflicting_status_error'
            reason = f'Conflicting Compleasm statuses: {sorted(statuses)!r}'
        else:
            status = copies[0][0].status
            eligible = status_allowed_for_features(accession, status)
            if not eligible and has_qc:
                resolution = 'contradictory_qc_error'
                reason = ('QC-ineligible ortholog has matching initial QC evidence; '
                          'cannot suppress contradictory observations')
            elif identical:
                resolution = 'identical_records_collapsed'
                reason = ('All parsed genomic records are identical; biological '
                          'identity and v9 feature calculations preserved')
            elif eligible:
                resolution = 'ambiguous_eligible_error'
                reason = ('Multiple distinct genomic copies of a QC-eligible '
                          f'{status} ortholog; no copy may be selected')
            else:
                resolution = 'ambiguous_ineligible_retained'
                reason = ('Multiple distinct genomic copies with the same '
                          f'QC-ineligible status {status}; no genomic copy selected')
        eligible_any = any(status_allowed_for_features(accession, r.status)
                           for r, _ in copies)
        for parsed, row in copies:
            duplicate_audit_rows.append({
                **row, 'n_records_for_ortholog': n,
                'qc_eligible': 'TRUE' if eligible_any else 'FALSE',
                'resolution_status': resolution,
                'resolution_reason': reason,
                'genomic_copy_selected': 'FALSE',
            })
        duplicate_counts['audit_source_copy_rows'] += n
        if resolution == 'identical_records_collapsed':
            duplicate_counts['identical_groups_collapsed'] += 1
            records[odb12_id] = copies[0][0]
        elif resolution == 'ambiguous_ineligible_retained':
            duplicate_counts['ambiguous_ineligible_retained'] += 1
            status = copies[0][0].status
            mapping_state = ('ambiguous_duplicated_ineligible'
                             if status in ('Duplicated', 'Duplicate')
                             else 'ambiguous_ineligible_genomic_copies')
            records[odb12_id] = FullRecord(
                odb12_id=odb12_id, status=status, sequence_id=None, strand=None,
                start=None, end=None, exons=[], mapping_state=mapping_state,
                ambiguity_reason=(
                    f'{n} distinct Compleasm genomic records exist for this '
                    f'QC-ineligible {status} ortholog. No genomic copy was '
                    'selected; see ortholog_duplicate_full_table_audit.tsv '
                    f'for {full_table} source lines '
                    + ', '.join(str(row['source_line']) for _, row in copies) + '.'
                ),
            )
        else:
            if resolution == 'ambiguous_eligible_error':
                duplicate_counts['ambiguous_eligible_rejected'] += 1
            duplicate_counts['rejected_groups'] += 1
            lines = [row['source_line'] for _, row in copies]
            raise ValueError(
                f'{full_table}: accession={accession}, genome_pk={genome_pk}, '
                f'odb12_id={odb12_id}, source_lines={lines}: '
                f'{resolution}: {reason}. '
                'Individual source records are saved in the v10 duplicate failure audit.'
            )
    duplicate_counts['unique_ortholog_identities'] += len(groups)
    duplicate_counts['deliberately_unresolved'] += sum(
        r.mapping_state != 'resolved' for r in records.values())
    return records


def slice_genome(genome: pyfaidx.Fasta, sequence_id: str, start_1: int, end_1: int) -> str:
    if sequence_id not in genome:
        raise KeyError(f"Sequence '{sequence_id}' not found in FASTA")
    chrom_len = len(genome[sequence_id])
    start_1 = max(1, int(start_1))
    end_1 = min(chrom_len, int(end_1))
    if start_1 > end_1:
        return ""
    return str(genome[sequence_id][start_1 - 1:end_1])


def infer_introns(exons: List[Dict[str, object]], strand: str) -> List[Dict[str, object]]:
    if len(exons) <= 1:
        return []
    genomic_introns: List[Dict[str, object]] = []
    for i in range(len(exons) - 1):
        left, right = exons[i], exons[i + 1]
        intron_start = int(left["end"]) + 1
        intron_end = int(right["start"]) - 1
        if intron_start > intron_end:
            continue
        genomic_introns.append({"start": intron_start, "end": intron_end, "length": intron_end - intron_start + 1})
    if strand == "-":
        genomic_introns.reverse()
    return genomic_introns


def extract_introns(genome: pyfaidx.Fasta, record: FullRecord) -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    for i, intron in enumerate(infer_introns(record.exons, record.strand), start=1):
        seq = slice_genome(genome, record.sequence_id, int(intron["start"]), int(intron["end"]))
        if record.strand == "-":
            seq = reverse_complement(seq)
        gc, *_counts, valid = calculate_gc(seq)
        out.append({**intron, "intron_id": f"i_{record.odb12_id}_{i}", "gc": gc, "valid_bases": valid})
    return out


def extract_flank(genome: pyfaidx.Fasta, record: FullRecord, up_bp: int, down_bp: int) -> Dict[str, object]:
    if record.start is None or record.end is None:
        raise ValueError(f"Cannot extract flanks for {record.odb12_id}; missing start/end")
    gene_start = int(record.start)
    gene_end = int(record.end)
    if record.strand == "+":
        up_start, up_end = max(1, gene_start - up_bp), gene_start - 1
        down_start, down_end = gene_end + 1, gene_end + down_bp
        up_seq = slice_genome(genome, record.sequence_id, up_start, up_end)
        down_seq = slice_genome(genome, record.sequence_id, down_start, down_end)
    else:
        up_start, up_end = gene_end + 1, gene_end + up_bp
        down_start, down_end = max(1, gene_start - down_bp), gene_start - 1
        up_seq = reverse_complement(slice_genome(genome, record.sequence_id, up_start, up_end))
        down_seq = reverse_complement(slice_genome(genome, record.sequence_id, down_start, down_end))
    combined = up_seq + down_seq
    gc, *_counts, valid = calculate_gc(combined)
    up_gc = calculate_gc(up_seq)[0]
    down_gc = calculate_gc(down_seq)[0]
    return {
        "up_start": up_start, "up_end": up_end,
        "down_start": down_start, "down_end": down_end,
        "gc": gc,
        "length": valid,
        "upstream_gc": up_gc,
        "downstream_gc": down_gc,
    }


def numeric_values(values: Sequence[object]) -> List[float]:
    out = []
    for v in values:
        if v is None or v == "":
            continue
        try:
            f = float(v)
        except (TypeError, ValueError):
            continue
        if not math.isnan(f):
            out.append(f)
    return out


def weighted_mean(values: Sequence[object], weights: Sequence[object]) -> Optional[float]:
    pairs = []
    for v, w in zip(values, weights):
        try:
            if v is None or w is None:
                continue
            vf, wi = float(v), int(w)
            if math.isnan(vf) or wi <= 0:
                continue
            pairs.append((vf, wi))
        except (TypeError, ValueError):
            continue
    if not pairs:
        return None
    denom = sum(w for _, w in pairs)
    return sum(v * w for v, w in pairs) / denom if denom else None


def q(values: Sequence[object], quantile: float) -> Optional[float]:
    vals = numeric_values(values)
    return float(np.quantile(vals, quantile)) if vals else None


def sd(values: Sequence[object]) -> Optional[float]:
    vals = numeric_values(values)
    return float(np.std(vals, ddof=1)) if len(vals) >= 2 else None


def var(values: Sequence[object]) -> Optional[float]:
    vals = numeric_values(values)
    return float(np.var(vals, ddof=1)) if len(vals) >= 2 else None


def mad(values: Sequence[object]) -> Optional[float]:
    vals = numeric_values(values)
    if not vals:
        return None
    med = median(vals)
    return float(median([abs(v - med) for v in vals]))


def stats(values: Sequence[object], lengths: Sequence[object], prefix: str) -> Dict[str, object]:
    vals = numeric_values(values)
    return {
        f"mean_{prefix}": mean(vals) if vals else None,
        f"weighted_mean_{prefix}": weighted_mean(values, lengths),
        f"sd_{prefix}": sd(values),
        f"var_{prefix}": var(values),
        f"median_{prefix}": median(vals) if vals else None,
        f"mad_{prefix}": mad(values),
        f"iqr_{prefix}": (q(values, 0.75) - q(values, 0.25)) if vals else None,
    }


# ---------------------------------------------------------------------------
# v8: integrity / provenance / run-specific QC infrastructure.
# The v7 composition and genomic-feature calculation functions above remain
# unchanged; this section only organizes existing measurement sources.
# ---------------------------------------------------------------------------

INITIAL_QC_COLUMNS = [
    'accession', 'odb12_id', 'final_qc', 'macse_qc', 'raw_cds_qc',
    'raw_length', 'raw_length_mod_3', 'raw_internal_stop_count',
    'raw_terminal_stop', 'raw_ambiguous_fraction', 'raw_invalid_fraction',
    'raw_failure_reasons', 'macse_frameshift', 'macse_internal_stop',
    'macse_terminal_stop', 'macse_failure_reasons', 'qc_version',
]
ORTHOLOG_V8_COLUMNS = [
    'ortholog_pk', 'genome_pk', 'sequence_pk', 'odb12_id', 'status',
    'strand', 'start', 'end', 'gc', 'gc3', 'gc4',
]
INITIAL_OUTPUT_COLUMNS = [
    'ortholog_pk', 'genome_pk', 'sequence_pk', 'odb12_id',
    *INITIAL_QC_COLUMNS[2:],
]
RUN_OUTPUT_COLUMNS = [
    'ortholog_run_qc_pk', 'ortholog_pk', 'nhphyml_run_pk',
    'final_alignment_qc', 'macse_failure_reasons',
    'included_in_nhphyml_alignment', 'qc_version',
]
RUN_CONFIG_COLUMNS = [
    'nhphyml_run_pk', 'run_name', 'final_macse_validation',
    'included_accessions', 'shared_orthologs', 'nhphyml_input_dir',
    'qc_version',
]
AUDIT_COLUMNS = [
    'ortholog_run_qc_pk', 'ortholog_pk', 'genome_pk', 'accession',
    'odb12_id', 'nhphyml_run_pk', 'run_name', 'initial_raw_cds_qc',
    'initial_final_qc', 'genome_selected_for_run',
    'ortholog_selected_for_run', 'evaluated_in_final_macse_validation',
    'final_alignment_qc', 'macse_failure_reasons',
    'included_in_nhphyml_alignment', 'selection_status',
    'exclusion_reason', 'validation_source', 'nhphyml_input_path',
    'mapping_status',
]
RUN_SUMMARY_COLUMNS = [
    'nhphyml_run_pk', 'run_name', 'selected_genomes',
    'selected_distinct_odb12_ids', 'selected_accession_ortholog_pairs',
    'evaluated_accession_ortholog_pairs', 'evaluated_distinct_odb12_ids',
    'final_qc_pass', 'final_qc_fail', 'membership_true',
    'membership_false', 'membership_unknown', 'pass_but_absent',
    'fail_but_present', 'run_qc_primary_keys', 'pk_unique',
    'composite_key_unique',
]
NULL_SQL = r'\N'

# Registered when each v9 TSV is generated. Carried-forward legacy files are
# deliberately excluded because they are not rewritten by this version.
GENERATED_TSV_COLUMNS = {}
REQUIRED_V9_TSVS = frozenset({
    'sauropsida_odb12.tsv', 'orthologs.tsv', 'ortholog_initial_qc.tsv',
    'ortholog_run_qc.tsv', 'ortholog_summary.tsv', 'intron_compleasm.tsv',
    'intron_compleasm_summary.tsv', 'flanks_compleasm.tsv',
    'flank_sets_compleasm.tsv', 'flank_compleasm_summary.tsv',
    'ortholog_source_audit.tsv', 'ortholog_run_selection_audit.tsv',
    'ortholog_run_qc_summary.tsv', 'feature_summary_count_comparison.tsv',
    'flank_summary.tsv', 'ortholog_duplicate_full_table_audit.tsv',
})


def sql_cell(value):
    """Serialize a scalar for SQL-import TSVs, preserving nonmissing values.

    Only actual scalar missing values and blank/whitespace-only strings become
    the literal backslash-N marker. Strings such as 'NA' and 'None' remain
    legitimate strings; source-file parsing rules are not changed.
    """
    if value is None:
        return NULL_SQL
    if isinstance(value, str):
        return NULL_SQL if not value.strip() else value
    # pd.isna recognizes Python/NumPy NaN, pandas.NA, pandas.NaT, NaT-like
    # scalars. Never evaluate a potentially array-valued result as a boolean.
    try:
        is_missing = pd.isna(value)
    except (TypeError, ValueError):
        is_missing = False
    if isinstance(is_missing, (bool, np.bool_)) and is_missing:
        return NULL_SQL
    return value


class SQLNullDictWriter(csv.DictWriter):
    """DictWriter that consistently serializes missing cells for all paths."""

    def writerow(self, rowdict):
        return super().writerow({col: sql_cell(rowdict.get(col))
                                 for col in self.fieldnames})


def clean_cell(value):
    """Retain 0/FALSE, but normalize source nulls to Python None."""
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    text_value = str(value).strip()
    if text_value in ('', r'\N') or text_value.lower() in ('nan', 'none', 'null', 'na'):
        return None
    return text_value


def nonempty(value, where):
    val = clean_cell(value)
    if val is None:
        raise ValueError(f'{where} must be non-NULL and non-empty')
    return val


def strict_rows(path, required, delimiter=None):
    if not path.is_file():
        raise FileNotFoundError(f'Required source file is absent: {path}')
    with path.open(newline='', encoding='utf-8-sig') as handle:
        if delimiter is None:
            delimiter = ',' if path.suffix.lower() == '.csv' else '\t'
        rdr = csv.DictReader(handle, delimiter=delimiter)
        columns = rdr.fieldnames or []
        if len(columns) != len(set(columns)):
            raise ValueError(f'Duplicate column names in {path}: {columns}')
        missing = sorted(set(required).difference(columns))
        if missing:
            raise ValueError(f'{path}: missing columns {missing}; found {columns}')
        for line, row in enumerate(rdr, 2):
            if None in row:
                raise ValueError(f'{path}:{line}: too many delimited fields')
            yield line, {c: clean_cell(v) for c, v in row.items()}


def required_header(path, columns):
    # Eagerly verify even empty source files (strict_rows is lazy).
    if not path.is_file():
        raise FileNotFoundError(f'Required source not found: {path}')
    delimiter = ',' if path.suffix.lower() == '.csv' else '\t'
    with path.open(newline='', encoding='utf-8-sig') as h:
        rd = csv.reader(h, delimiter=delimiter)
        header = next(rd, [])
    if len(header) != len(set(header)) or not set(columns).issubset(header):
        raise ValueError(f'{path}: expected columns {columns}; observed {header}')
    return header


def read_initial_qc(path, scope):
    required_header(path, INITIAL_QC_COLUMNS)
    qc = {}
    # Scope is current processing set, not entire manifest when --test.
    for line, row in strict_rows(path, INITIAL_QC_COLUMNS):
        accession = nonempty(row['accession'], f'{path}:{line} accession')
        if accession not in scope:
            continue
        odb = nonempty(row['odb12_id'], f'{path}:{line} odb12_id')
        key = (accession, odb)
        if key in qc:
            raise ValueError(f'{path}:{line}: duplicate initial-QC key accession={accession}, odb12_id={odb}; previous line {qc[key][0]}')
        for col in ('raw_cds_qc', 'macse_qc'):
            if row[col] is not None and row[col] not in ('pass', 'fail'):
                raise ValueError(f'{path}:{line}: invalid {col}={row[col]!r} for {key}')
        if row['final_qc'] is not None and row['final_qc'] not in ('passing', 'rescued', 'unsuitable', 'alignment-problematic'):
            raise ValueError(f'{path}:{line}: unrecognized final_qc={row["final_qc"]!r} for {key}')
        if row['raw_cds_qc'] is None:
            raise ValueError(f'{path}:{line}: QC record lacks raw_cds_qc for {key}')
        if row['macse_qc'] is None:
            if row['final_qc'] is not None or any(row[k] is not None for k in (
                'macse_frameshift', 'macse_internal_stop', 'macse_terminal_stop', 'macse_failure_reasons'
            )):
                raise ValueError(f'{path}:{line}: raw-only record contains MACSE/final QC fields for {key}')
        qc[key] = (line, row)
    return qc


def enforce_one_to_one_genome_mappings(metadata, genome_lookup, scope):
    if metadata['accession'].duplicated().any():
        dup = metadata.loc[metadata['accession'].duplicated(False), 'accession'].tolist()
        raise ValueError(f'Duplicate Compleasm metadata accessions: {dup[:20]}')
    if set(metadata['accession']) != set(scope):
        raise ValueError(f'Missing Compleasm metadata for manifest-scope genomes: {sorted(set(scope)-set(metadata["accession"]))[:20]}')
    for accession in scope:
        if accession not in genome_lookup:
            raise ValueError(f'No exact genome_pk match in genomes.tsv for {accession}')
    pks = [genome_lookup[a]['genome_pk'] for a in scope]
    if len(pks) != len(set(pks)):
        raise ValueError('Two manifest accessions resolve to the same genome_pk; mapping ambiguous')


def check_lookup_table_uniqueness(path, columns, key_cols, id_col=None):
    required_header(path, columns)
    keys, ids = set(), set()
    for line, row in strict_rows(path, columns):
        key = tuple(row.get(k) for k in key_cols)
        if any(v is None for v in key):
            raise ValueError(f'{path}:{line}: NULL lookup key {key_cols}')
        if key in keys:
            raise ValueError(f'{path}:{line}: duplicate lookup key {key_cols}={key}')
        keys.add(key)
        if id_col:
            identity = nonempty(row[id_col], f'{path}:{line} {id_col}')
            if identity in ids:
                raise ValueError(f'{path}:{line}: duplicate {id_col}={identity}')
            ids.add(identity)


def read_manifest_taxa(path, scope):
    delimiter = ',' if path.suffix.lower() == '.csv' else '\t'
    with path.open(encoding='utf-8-sig', newline='') as h:
        reader = csv.DictReader(h, delimiter=delimiter)
        if not reader.fieldnames or not {'accession','species_key'}.issubset(reader.fieldnames):
            raise ValueError(f'{path}: manifest must contain accession and species_key for run-label validation')
        taxa = {}
        for line, row in enumerate(reader, 2):
            a = clean_cell(row.get('accession'))
            if a is None:
                continue
            name = nonempty(row.get('species_key'), f'{path}:{line} species_key')
            if a in taxa:
                raise ValueError(f'{path}:{line}: duplicate manifest accession {a}')
            parts = name.split('_')
            if len(parts) < 2 or any(not x for x in parts[:2]):
                raise ValueError(f'{path}:{line}: invalid species_key {name!r}')
            taxa[a] = (name, '_'.join(parts[:2]))
    if not scope.issubset(taxa):
        raise ValueError(f'Missing manifest taxon labels: {sorted(scope-set(taxa))[:20]}')
    return taxa


def config_path(config_file, value, is_dir=False):
    value = nonempty(value, f'configuration path in {config_file}')
    p = Path(value).expanduser()
    if not p.is_absolute():
        p = config_file.parent / p
    p = p.resolve()
    if not (p.is_dir() if is_dir else p.is_file()):
        raise FileNotFoundError(f'Configured {"directory" if is_dir else "file"} missing: {p}')
    return p


def read_nhphyml_runs(path):
    runs = {}
    for line, row in strict_rows(path, ['nhphyml_run_pk', 'run_name', 'alignment_name']):
        pk_str = nonempty(row['nhphyml_run_pk'], f'{path}:{line} nhphyml_run_pk')
        try:
            pk = int(pk_str)
        except ValueError:
            raise ValueError(f'{path}:{line}: invalid nhphyml_run_pk {pk_str}')
        if pk < 1 or pk in runs:
            raise ValueError(f'{path}:{line}: repeated/invalid nhphyml_run_pk {pk}')
        runs[pk] = (nonempty(row['run_name'], f'{path}:{line} run_name'),
                    nonempty(row['alignment_name'], f'{path}:{line} alignment_name'))
    return runs


def read_run_config(path, runs_metadata):
    required_header(path, RUN_CONFIG_COLUMNS)
    configs = []
    seen = set()
    for line, row in strict_rows(path, RUN_CONFIG_COLUMNS):
        try:
            pk = int(nonempty(row['nhphyml_run_pk'], f'{path}:{line} nhphyml_run_pk'))
        except ValueError:
            raise ValueError(f'{path}:{line}: invalid nhphyml_run_pk')
        if pk not in runs_metadata:
            raise ValueError(f'{path}:{line}: unknown nhphyml_run_pk {pk} in nhphyml_runs.tsv')
        if pk in seen:
            raise ValueError(f'{path}:{line}: repeated run configuration {pk}')
        seen.add(pk)
        run_name = nonempty(row['run_name'], f'{path}:{line} run_name')
        metadata_run_name, alignment_name = runs_metadata[pk]
        if run_name not in {metadata_run_name, alignment_name}:
            raise ValueError(f'{path}:{line}: run_name={run_name!r} incompatible with nhphyml_runs.tsv run_name={metadata_run_name!r} alignment_name={alignment_name!r}')
        cfg = {'nhphyml_run_pk': pk, 'run_name': run_name, 'alignment_name': alignment_name}
        for key in ('final_macse_validation','included_accessions','shared_orthologs'):
            cfg[key] = config_path(path, row[key])
        cfg['nhphyml_input_dir'] = config_path(path, row['nhphyml_input_dir'], is_dir=True)
        cfg['qc_version'] = row['qc_version']
        configs.append(cfg)
    if not configs:
        raise ValueError(f'{path}: no run configurations')
    return sorted(configs, key=lambda c: c['nhphyml_run_pk'])


def read_one_column_set(path, label):
    result = set()
    with path.open(encoding='utf-8-sig') as h:
        for line, raw in enumerate(h, 1):
            s = raw.strip()
            if not s or s.startswith('#'):
                continue
            if '\t' in s or ' ' in s:
                raise ValueError(f'{path}:{line}: {label} must be one identifier per line; observed {s!r}')
            if s in result:
                raise ValueError(f'{path}:{line}: duplicate {label}: {s}')
            result.add(s)
    if not result:
        raise ValueError(f'{path}: empty {label}')
    return result


def get_run_validation(path, known_ortholog_keys, manifest_taxa, cfg_version, test=False, scope=None):
    """Use ONLY individual-sequence observations; alignment-level records fail closed.

    Accepted exact source identity is (accession, odb12_id). Alternatively
    (species_key/species_key_short/taxon/sequence_id, odb12_id) is allowed only
    when the given taxon label uniquely determines an accession in the manifest.
    An 'alignment'-only row is *not* per-sequence evidence.
    """
    header = required_header(path, ['final_alignment_qc'])
    odb_col = next((c for c in ('odb12_id', 'ortholog_id', 'Gene') if c in header), None)
    if odb_col is None and 'alignment' in header:
        odb_col = 'alignment'
    acc_col = next((c for c in ('accession','accession_id') if c in header), None)
    taxon_col = next((c for c in ('species_key','species_key_short','taxon','sequence_id') if c in header), None)
    if odb_col is None or (acc_col is None and taxon_col is None):
        raise ValueError(
            f'{path}: cannot assign whole-alignment QC to ortholog_pk. Observed {header}; '
            'need per-sequence rows with accession + odb12_id (preferred), or unique '
            'manifest taxon label + odb12_id; otherwise generate a genuine per-sequence QC file.'
        )
    if 'macse_failure_reasons' not in header:
        raise ValueError(f'{path}: missing macse_failure_reasons; source columns: {header}')
    version_col = 'qc_version' if 'qc_version' in header else None
    taxa_to_accessions = defaultdict(set)
    for accession, (full, short) in manifest_taxa.items():
        for alias in {accession, full, short}:
            taxa_to_accessions[alias].add(accession)
    outcomes = {}
    for line, row in strict_rows(path, [odb_col, 'final_alignment_qc','macse_failure_reasons']):
        odb = nonempty(row[odb_col], f'{path}:{line} {odb_col}')
        for ext in ('.fasta', '.fas', '.fa', '.phylip', '.phy'):
            if odb_col == 'alignment' and odb.endswith(ext):
                odb = odb[:-len(ext)]
                break
        if acc_col:
            accession = nonempty(row[acc_col], f'{path}:{line} {acc_col}')
        else:
            taxon = nonempty(row[taxon_col], f'{path}:{line} {taxon_col}')
            matches = taxa_to_accessions.get(taxon, set())
            if len(matches) != 1:
                raise ValueError(f'{path}:{line}: ambiguous/unmapped {taxon_col}={taxon!r} -> {sorted(matches)}')
            accession = next(iter(matches))
        key = (accession, odb)
        if test and accession not in scope:
            continue
        if key not in known_ortholog_keys:
            raise ValueError(f'{path}:{line}: run QC identity {key} missing from current orthologs')
        if key in outcomes:
            raise ValueError(f'{path}:{line}: repeated sequence QC record {key}; previous line={outcomes[key]["line"]}')
        status = nonempty(row['final_alignment_qc'], f'{path}:{line} final_alignment_qc')
        if status not in ('pass','fail'):
            raise ValueError(f'{path}:{line}: unrecognized final_alignment_qc={status!r} for {key}')
        version = row.get(version_col) if version_col else None
        if version is None:
            version = cfg_version
        if version is None:
            raise ValueError(f'{path}:{line}: no qc_version in source or run config for {key}')
        outcomes[key] = {'line':line,'final_alignment_qc':status,
                         'macse_failure_reasons':row['macse_failure_reasons'],
                         'qc_version':version}
    return outcomes


def read_nhphyml_input(path):
    """Return exact, complete taxon labels; reject ambiguous or truncated ids."""
    from Bio import AlignIO
    if path.suffix.lower() in ('.phylip', '.phy'):
        fmt = 'phylip-relaxed'
    elif path.suffix.lower() in ('.fa', '.fasta', '.fas', '.fna'):
        fmt = 'fasta'
    else:
        raise ValueError(f'Unknown nhPhyML input alignment format: {path}')
    try:
        aln = AlignIO.read(str(path), fmt)
    except Exception as exc:
        raise ValueError(f'Could not parse nhPhyML input alignment {path}: {exc}') from exc
    labels = [seq.id for seq in aln]
    if len(labels) != len(set(labels)):
        raise ValueError(f'{path}: duplicated taxon identifiers in nhPhyML input')
    if not labels:
        raise ValueError(f'{path}: empty nhPhyML alignment')
    return set(labels)


def discover_input_paths(cfg):
    """Require exact <odb12_id>.<extension> names, without guessing paths."""
    directory = cfg['nhphyml_input_dir']
    extensions = ('.phylip','.phy','.fasta','.fas','.fa','.fna')
    result = {}
    for p in sorted(directory.iterdir()):
        if not p.is_file() or p.suffix.lower() not in extensions:
            continue
        odb = p.name[:-len(p.suffix)]
        if odb in result:
            raise ValueError(f'Run {cfg["run_name"]}: multiple possible inputs for odb12_id={odb}: {result[odb]}, {p}')
        result[odb] = p
    if not result:
        raise ValueError(f'{directory}: no exact <odb12_id>.<phylip/phy/fasta> inputs found; configure the actual nhPhyML input directory')
    return result


def run_label_resolver(accessions, manifest_taxa):
    alias = defaultdict(set)
    for accession in accessions:
        full, short = manifest_taxa[accession]
        for label in {accession, full, short}:
            alias[label].add(accession)
    return alias


def membership_for_alignment(path, resolver, expected_accessions, odb):
    """Parse once per ODB alignment, return full membership by accession."""
    labels = read_nhphyml_input(path)
    present = set()
    for label in labels:
        possible = resolver.get(label, set())
        if len(possible) != 1:
            raise ValueError(f'{path}: unknown/ambiguous/truncated taxon label {label!r}; matching accessions={sorted(possible)}. Exact untruncated labels are required')
        accession = next(iter(possible))
        if accession not in expected_accessions:
            raise ValueError(f'{path}: taxon {label!r} maps to excluded accession={accession}')
        if accession in present:
            raise ValueError(f'{path}: multiple taxon labels resolve to accession={accession}')
        present.add(accession)
    return present


def sql_boolean(value):
    if sql_cell(value) == NULL_SQL:
        return NULL_SQL
    return 'TRUE' if value else 'FALSE'


def tsv_writer(path, columns, sql_null=False):
    """Open a streaming TSV and register the precise expected output header.

    The legacy sql_null argument remains accepted for compatibility; all v9
    outputs use the same serialization regardless of its value.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    f = path.open('w', newline='', encoding='utf-8')
    writer = SQLNullDictWriter(f, fieldnames=columns, delimiter='\t',
                               extrasaction='ignore')
    try:
        writer.writeheader()
    except BaseException:
        f.close()
        raise
    GENERATED_TSV_COLUMNS[path.resolve()] = list(columns)
    return f, writer


def audit_row(writer, row):
    writer.writerow(row)


def write_run_tables(stage, configs, ortholog_rows, ortholog_qc, manifest_taxa, strict_inclusion=False, test=False):
    """Streaming audit (potentially 2M+ rows); stream run QC in sorted order.

    Ortholog rows are in accession, odb12 order; sorted index by ortholog_pk.
    Run QC PKs increase with (nhphyml_run_pk, ortholog_pk).
    """
    index = {(row['accession'],row['odb12_id']):row for row in ortholog_rows}
    if len(index) != len(ortholog_rows):
        raise ValueError('Ortholog identity is not unique in v8 output')
    all_accessions = {r['accession'] for r in ortholog_rows}
    qcf, qcw = tsv_writer(stage/'ortholog_run_qc.tsv', RUN_OUTPUT_COLUMNS)
    audf, audw = tsv_writer(stage/'ortholog_run_selection_audit.tsv', AUDIT_COLUMNS)
    totals = []
    pk_count = 0
    audit_count = 0
    composite_seen = set()
    try:
        for cfg in configs:
            pk = cfg['nhphyml_run_pk']
            included = read_one_column_set(cfg['included_accessions'], 'included_accession')
            shared = read_one_column_set(cfg['shared_orthologs'], 'shared_ortholog')
            full_included = set(included)
            if not included.issubset(all_accessions):
                unknown = sorted(included-all_accessions)
                if test:
                    print(f'TEST: run {pk}, {len(unknown)} selected accessions outside two-genome scope; ignoring those for this execution')
                    included &= all_accessions
                else:
                    raise ValueError(f'Run {pk}: included_accessions not in current ortholog scope: {unknown[:20]}')
            if not included:
                if test:
                    print(f'TEST: run {pk}: no selected genomes in test scope; skipping run QC')
                    continue
                raise ValueError(f'Run {pk}: no included accessions')
            observed = get_run_validation(cfg['final_macse_validation'], index, manifest_taxa, cfg['qc_version'], test=test, scope=all_accessions)
            # Exclude validation records outside test scope only if they could
            # not have been loaded; full runs require all records resolvable.
            evaluated_outside = [key for key in observed if key[0] not in included or key[1] not in shared]
            if evaluated_outside:
                print(f'WARNING: run {pk}: {len(evaluated_outside)} sequences explicitly evaluated in final validation were not in selected genome/ortholog sets; retained in run QC, see audit. Examples: {evaluated_outside[:5]}')
            input_paths = discover_input_paths(cfg)
            unknown_input_ids = set(input_paths)-shared
            if unknown_input_ids:
                raise ValueError(f'Run {pk}: nhPhyML input filenames contain orthologs outside shared_orthologs: {sorted(unknown_input_ids)[:10]}')
            if not full_included.issubset(manifest_taxa):
                raise ValueError(f'Run {pk}: selected accessions missing from manifest: {sorted(full_included-set(manifest_taxa))[:10]}')
            resolver = run_label_resolver(full_included, manifest_taxa)
            membership = {}
            # Check all input labels once and reuse across ortholog records.
            for odb, path in input_paths.items():
                membership[odb] = membership_for_alignment(path, resolver, full_included, odb)
            summary = {'nhphyml_run_pk':pk,'run_name':cfg['run_name'],
                'selected_genomes':len(included),'selected_distinct_odb12_ids':len(shared),
                'selected_accession_ortholog_pairs':sum((r['accession'] in included and r['odb12_id'] in shared) for r in ortholog_rows),
                'evaluated_accession_ortholog_pairs':len(observed),
                'evaluated_distinct_odb12_ids':len({o for _,o in observed}),
                'final_qc_pass':0,'final_qc_fail':0,'membership_true':0,
                'membership_false':0,'membership_unknown':0,'pass_but_absent':0,
                'fail_but_present':0,'run_qc_primary_keys':0,
                'pk_unique':'TRUE','composite_key_unique':'TRUE'}
            for ortho in ortholog_rows:
                key = (ortho['accession'], ortho['odb12_id'])
                selected_genome = ortho['accession'] in included
                selected_ortholog = ortho['odb12_id'] in shared
                selection = bool(selected_genome and selected_ortholog)
                observed_qc = observed.get(key)
                input_path = input_paths.get(key[1])
                present = (ortho['accession'] in membership[key[1]]) if key[1] in membership else None
                if ortho['mapping_status'].startswith('ambiguous_') and (observed_qc is not None or present is True):
                    raise ValueError(
                        f'Run {pk}: QC-ineligible unresolved genomic copies contradict '
                        f'run evidence: accession={key[0]}, odb12_id={key[1]}, '
                        f'ortholog_pk={ortho["ortholog_pk"]}, '
                        f'final_validation_record={observed_qc is not None}, '
                        f'nhphyml_input_member={present}, input={input_path}; '
                        'see ortholog_duplicate_full_table_audit.tsv'
                    )
                status, reason = None, None
                this_pk = None
                if observed_qc:
                    if key[1] not in input_paths:
                        present = None
                    if not selected_genome and present is True:
                        raise ValueError(f'Run {pk}: excluded genome {key[0]} found in nhPhyML input {input_path}')
                    if not selected_ortholog and present is True:
                        raise ValueError(f'Run {pk}: ortholog {key[1]} not in shared_orthologs found in nhPhyML input {input_path}')
                    if present is True and observed_qc['final_alignment_qc'] == 'fail':
                        raise ValueError(f'FATAL final QC fail but included in nhPhyML: accession={key[0]}, odb12_id={key[1]}, ortholog_pk={ortho["ortholog_pk"]}, ortholog_run_qc_pk={pk_count + 1}, nhphyml_run_pk={pk}, qc={observed_qc["final_alignment_qc"]}, input={input_path}')
                    if observed_qc['final_alignment_qc']=='pass' and present is False:
                        if selection:
                            summary['pass_but_absent'] += 1
                            status,reason='pass_not_in_input','unexpected_downstream_exclusion'
                            if strict_inclusion:
                                raise ValueError(f'Run {pk}: QC pass but absent from nhPhyML input: {key}; input={input_path}')
                        else:
                            status,reason='evaluated_not_selected','evaluated_but_excluded_from_subset' 
                    elif observed_qc['final_alignment_qc']=='fail':
                        status,reason='failed_final_alignment_qc','final_alignment_qc_fail'
                    elif not selection:
                        status,reason='evaluated_not_selected','evaluated_but_excluded_from_subset'
                    elif present is True:
                        status,reason='included_in_nhphyml_input',None
                    elif present is None:
                        status,reason='membership_unverifiable','nhphyml_input_unavailable'
                    else:
                        status,reason='passed_final_alignment_qc',None
                    summary['final_qc_'+observed_qc['final_alignment_qc']] += 1
                    summary['membership_true' if present is True else 'membership_false' if present is False else 'membership_unknown'] += 1
                    pk_count += 1
                    this_pk = pk_count
                    alt = (ortho['ortholog_pk'],pk)
                    if alt in composite_seen:
                        raise ValueError(f'Duplicate run QC alternate key {alt}')
                    composite_seen.add(alt)
                    qcw.writerow({
                        'ortholog_run_qc_pk':this_pk, 'ortholog_pk':ortho['ortholog_pk'],
                        'nhphyml_run_pk':pk,
                        'final_alignment_qc':observed_qc['final_alignment_qc'],
                        'macse_failure_reasons':sql_cell(observed_qc['macse_failure_reasons']),
                        'included_in_nhphyml_alignment':sql_boolean(present),
                        'qc_version':observed_qc['qc_version'],
                    })
                    summary['run_qc_primary_keys'] += 1
                else:
                    if not selected_genome:
                        status,reason='genome_excluded','genome_excluded_from_subset'
                    elif not selected_ortholog:
                        status,reason='ortholog_not_selected','ortholog_not_in_shared_subset'
                    elif present is True:
                        raise ValueError(f'Run {pk}: ortholog {key} appears in nhPhyML input but is absent from final MACSE validation')
                    elif present is False:
                        status,reason='selected_not_evaluated','no_final_alignment_qc_record'
                    else:
                        status,reason='membership_unverifiable','no_final_qc_or_input_evidence'
                audit_row(audw, {
                    'ortholog_run_qc_pk':this_pk,'ortholog_pk':ortho['ortholog_pk'],
                    'genome_pk':ortho['genome_pk'],'accession':ortho['accession'],
                    'odb12_id':ortho['odb12_id'],'nhphyml_run_pk':pk,
                    'run_name':cfg['run_name'],
                    'initial_raw_cds_qc':ortholog_qc[key].get('raw_cds_qc'),
                    'initial_final_qc':ortholog_qc[key].get('final_qc'),
                    'genome_selected_for_run':sql_boolean(selected_genome),
                    'ortholog_selected_for_run':sql_boolean(selected_ortholog),
                    'evaluated_in_final_macse_validation':sql_boolean(bool(observed_qc)),
                    'final_alignment_qc':observed_qc['final_alignment_qc'] if observed_qc else None,
                    'macse_failure_reasons':observed_qc['macse_failure_reasons'] if observed_qc else None,
                    'included_in_nhphyml_alignment':sql_boolean(present),
                    'selection_status':status,'exclusion_reason':reason,
                    'validation_source':str(cfg['final_macse_validation']) if observed_qc else None,
                    'nhphyml_input_path':str(input_path) if input_path else None,
                    'mapping_status':ortho['mapping_status'],
                })
                audit_count += 1
            totals.append(summary)
            if summary['pass_but_absent']:
                print(f'WARNING: Run {pk}: {summary["pass_but_absent"]} QC-passing sequences were CONFIRMED ABSENT from the actual nhPhyML input; see ortholog_run_selection_audit.tsv')
            if summary['membership_unknown']:
                print(f'WARNING: Run {pk}: {summary["membership_unknown"]} evaluated sequences have unverifiable nhPhyML input membership; see ortholog_run_selection_audit.tsv')
            print(f'  Run {pk} ({cfg["run_name"]}): selected genomes={len(included)}, distinct ODB12={len(shared)}, evaluated sequences={len(observed)}, pass={summary["final_qc_pass"]}, fail={summary["final_qc_fail"]}, present={summary["membership_true"]}, absent={summary["membership_false"]}, unknown={summary["membership_unknown"]}, pass-but-absent={summary["pass_but_absent"]}')
    finally:
        qcf.close()
        audf.close()
    write_tsv(totals, stage/'ortholog_run_qc_summary.tsv', RUN_SUMMARY_COLUMNS)
    print(f'Run QC rows: {pk_count}; run-selection audit rows: {audit_count}')
    return totals


def validate_staged_relations(stage, expected_orthologs, expected_run_qc):
    """Read staged TSVs and verify PK, FKs, ordering and audit correspondence."""
    from itertools import zip_longest
    with (stage/'orthologs.tsv').open(newline='',encoding='utf-8') as o, \
         (stage/'ortholog_initial_qc.tsv').open(newline='',encoding='utf-8') as i:
        ordr = csv.DictReader(o,delimiter='\t')
        iqdr = csv.DictReader(i,delimiter='\t')
        count = 0
        for ortho, qc in zip_longest(ordr,iqdr):
            count += 1
            if ortho is None or qc is None:
                raise ValueError('Staged orthologs and initial QC row counts differ')
            if int(ortho['ortholog_pk']) != count or int(qc['ortholog_pk']) != count:
                raise ValueError(f'Staged ortholog/initial-QC PK out of order at row {count}')
            for field in ('genome_pk','sequence_pk','odb12_id'):
                if ortho[field] != qc[field]:
                    raise ValueError(f'Staged initial QC FK mismatch at ortholog_pk={count}, field={field}')
        if count != expected_orthologs:
            raise ValueError(f'Staged ortholog row count mismatch: expected {expected_orthologs}, found {count}')
    with (stage/'ortholog_run_qc.tsv').open(newline='',encoding='utf-8') as q, \
         (stage/'ortholog_run_selection_audit.tsv').open(newline='',encoding='utf-8') as a:
        runr = csv.DictReader(q,delimiter='\t')
        audr = csv.DictReader(a,delimiter='\t')
        next_qc = next(runr,None)
        count, previous_pair = 0, (0,0)
        for ar in audr:
            if ar['ortholog_run_qc_pk']==NULL_SQL:
                continue
            count += 1
            if next_qc is None:
                raise ValueError(f'Audit references extra run QC record at position {count}')
            match_fields = ('ortholog_run_qc_pk','ortholog_pk','nhphyml_run_pk')
            for field in match_fields:
                if ar[field] != next_qc[field]:
                    raise ValueError(f'Audit/run QC FK discrepancy at row {count}, field {field}: {ar[field]!r} != {next_qc[field]!r}')
            pk = int(next_qc['ortholog_run_qc_pk'])
            o_pk = int(next_qc['ortholog_pk'])
            r_pk = int(next_qc['nhphyml_run_pk'])
            if pk != count or not (1<=o_pk<=expected_orthologs):
                raise ValueError(f'Invalid staged ortholog_run_qc_pk={pk} or ortholog_pk={o_pk}')
            pair = (r_pk,o_pk)
            if pair <= previous_pair:
                raise ValueError(f'Run QC composite key not strictly ordered/unique: {pair} <= {previous_pair}')
            previous_pair = pair
            next_qc = next(runr,None)
        if next_qc is not None or count != expected_run_qc:
            raise ValueError(f'Staged run QC audit references incomplete: {count} vs expected {expected_run_qc}')
    print(f'Staged relational checks passed: {expected_orthologs} orthologs/initial-QC pairs, {expected_run_qc} run-QC records, audit FK associations')



def validate_v10_duplicate_integrity(stage, ortholog_rows, initial_qc_rows,
                                     intron_rows, flank_rows, source_audit_rows):
    """Enforce v13 per-genome uniqueness and absence of copy-derived features."""
    seen_pk = set()
    seen_identity = set()
    initial = {row['ortholog_pk']: row for row in initial_qc_rows}
    source = {row['ortholog_pk']: row for row in source_audit_rows}
    intron_ids = {row['ortholog_pk'] for row in intron_rows}
    flank_ids = {row['ortholog_pk'] for row in flank_rows}
    if len(initial) != len(ortholog_rows) or len(source) != len(ortholog_rows):
        raise ValueError('v13: duplicate/missing ortholog initial-QC or source-audit PK')
    for ortho in ortholog_rows:
        pk = ortho['ortholog_pk']
        key = (ortho['genome_pk'], ortho['odb12_id'])
        if pk in seen_pk or key in seen_identity:
            raise ValueError(f'v13: duplicate ortholog primary key or genome/ODB12 identity {key}')
        seen_pk.add(pk)
        seen_identity.add(key)
        if ortho['genome_pk'] is None or not ortho['odb12_id']:
            raise ValueError(f'v13: ortholog_pk={pk} has NULL genome_pk or odb12_id')
        qc = initial[pk]
        for field in ('genome_pk', 'sequence_pk', 'odb12_id'):
            if qc[field] != ortho[field]:
                raise ValueError(f'v13: ortholog_pk={pk} initial-QC mismatch in {field}')
        if source[pk]['mapping_status'] not in {
                'ambiguous_duplicated_ineligible', 'ambiguous_ineligible_genomic_copies'}:
            continue
        for field in ('sequence_pk', 'strand', 'start', 'end', 'gc', 'gc3', 'gc4'):
            if ortho[field] is not None:
                raise ValueError(f'Unresolved copy ortholog_pk={pk}: {field} must be SQL NULL')
        for field in ('source_sequence_id','source_strand','source_start',
                      'source_end','source_codon_intervals','sequence_pk'):
            if source[pk][field] is not None:
                raise ValueError(f'Unresolved copy ortholog_pk={pk}: source {field} must be SQL NULL')
        if any(qc[field] is not None for field in INITIAL_QC_COLUMNS[2:]):
            raise ValueError(f'Unresolved copy ortholog_pk={pk}: must have NULL initial QC')
        if pk in intron_ids or pk in flank_ids:
            raise ValueError(f'Unresolved copy ortholog_pk={pk}: generated genomic intron/flank')
    recorded = Counter((row['accession'], row['odb12_id']) for row in DUPLICATE_AUDIT_ROWS)
    for row in DUPLICATE_AUDIT_ROWS:
        key = (row['accession'], row['odb12_id'])
        if row['n_records_for_ortholog'] != recorded[key]:
            raise ValueError(f'Duplicate audit count discrepancy for {key}')
        if row['genomic_copy_selected'] != 'FALSE':
            raise ValueError(f'Duplicate audit falsely selects a genomic copy for {key}')
    print(f'v13 ortholog uniqueness and unresolved-copy feature exclusion passed: {len(ortholog_rows)} identities, {DUPLICATE_COUNTS["deliberately_unresolved"]} unresolved')

def compare_prior_summary_counts(stage, outdir):
    """Audit observed count differences vs previous build; not causal claims."""
    names = ['ortholog_summary.tsv', 'intron_compleasm_summary.tsv',
             'flank_compleasm_summary.tsv', 'orthologs.tsv',
             'intron_compleasm.tsv', 'flanks_compleasm.tsv']
    rows = []
    for name in names:
        old = outdir/name
        new = stage/name
        def describe(path):
            if not path.is_file():
                return None, None
            with path.open(encoding='utf-8', newline='') as f:
                rd = csv.DictReader(f,delimiter='\t')
                n, total = 0, 0
                has_n_orthologs = bool(rd.fieldnames and 'n_orthologs' in rd.fieldnames)
                for row in rd:
                    n += 1
                    if has_n_orthologs and clean_cell(row['n_orthologs']) is not None:
                        total += int(row['n_orthologs'])
            return n, total if has_n_orthologs else None
        old_n, old_summary = describe(old)
        new_n, new_summary = describe(new)
        changed = (old_n is not None and (new_n != old_n or new_summary != old_summary))
        rows.append({'table':name, 'previous_rows':old_n, 'new_rows':new_n,
                     'previous_sum_n_orthologs':old_summary,
                     'new_sum_n_orthologs':new_summary,
                     'changed_from_previous': 'TRUE' if changed else 'FALSE' if old_n is not None else None})
        if changed:
            print(f'WARNING: {name}: previous rows={old_n}, new rows={new_n}; previous sum n_orthologs={old_summary}, new sum={new_summary}; investigate changes in input/source validation')
    write_sql_tsv(rows,stage/'feature_summary_count_comparison.tsv',
                  ['table','previous_rows','new_rows','previous_sum_n_orthologs',
                   'new_sum_n_orthologs','changed_from_previous'])


def commit_staging(stage, outdir):
    """Atomic directory switch on one filesystem; restore old outputs on error."""
    import os
    import shutil
    backup = outdir.with_name(outdir.name + '.v10_previous_build')
    if backup.exists():
        raise FileExistsError(f'Stale backup found; inspect before continuing: {backup}')
    had_output = outdir.exists()
    if had_output:
        # Preserve files v7 may have left that v8 does not explicitly own.
        for old in outdir.iterdir():
            new = stage/old.name
            if not new.exists():
                if old.is_dir():
                    shutil.copytree(old, new)
                elif old.is_file():
                    shutil.copy2(old, new)
        os.replace(outdir, backup)
    try:
        os.replace(stage, outdir)
    except BaseException:
        if had_output and backup.exists():
            os.replace(backup, outdir)
        raise
    if had_output:
        shutil.rmtree(backup)


def write_failure_details(outdir, exc):
    import traceback
    path = outdir.parent/(outdir.name + '.v10_failure_report.tsv')
    row = {'time':datetime.now().isoformat(timespec='seconds'),
           'error_type':type(exc).__name__, 'message':str(exc),
           'traceback':traceback.format_exc()}
    write_tsv([row], path, ['time','error_type','message','traceback'])
    print(f'ERROR details: {path}', file=sys.stderr)
    if DUPLICATE_AUDIT_ROWS:
        duplicate_path = outdir.parent / (outdir.name + '.v10_duplicate_failure_audit.tsv')
        write_tsv(DUPLICATE_AUDIT_ROWS, duplicate_path, DUPLICATE_AUDIT_COLUMNS)
        print(f'Duplicate copy diagnostics: {duplicate_path}', file=sys.stderr)


def write_sql_tsv(rows, path, columns):
    """Compatibility entry point: SQL and feature TSVs share one serializer."""
    write_tsv(rows, path, columns)


def is_typed_null_field(column):
    """Flag NULL-like *text* only where it cannot be a legitimate text value.

    Do not classify arbitrary 'NA'/'None' strings in free-text fields as NULL;
    the output serializer intentionally preserves them as valid strings.
    """
    numeric_exact = {
        'start', 'end', 'up_start', 'up_end', 'gc', 'gc3', 'gc4',
        'source_start', 'source_end', 'source_score',
        'source_length', 'source_identity', 'source_line',
        'n_records_for_ortholog', 'upstream_bp', 'downstream_bp',
        'raw_length', 'raw_length_mod_3', 'raw_terminal_stop',
        'raw_ambiguous_fraction', 'raw_invalid_fraction',
        'raw_internal_stop_count', 'macse_frameshift',
        'macse_internal_stop', 'macse_terminal_stop',
        'n_orthologs', 'n_flanks', 'n_introns',
        'callable_bp_total', 'mean_upstream_gc', 'mean_downstream_gc',
        'selected_genomes', 'selected_distinct_odb12_ids',
        'selected_accession_ortholog_pairs', 'evaluated_accession_ortholog_pairs',
        'evaluated_distinct_odb12_ids', 'final_qc_pass', 'final_qc_fail',
        'membership_true', 'membership_false', 'membership_unknown',
        'pass_but_absent', 'fail_but_present', 'run_qc_primary_keys',
        'previous_rows', 'new_rows', 'previous_sum_n_orthologs',
        'new_sum_n_orthologs',
    }
    boolean_exact = {
        'passes_raw_cds_qc', 'cds_present', 'gc_calculated',
        'included_in_ortholog_summary', 'genomic_features_eligible',
        'included_in_nhphyml_alignment', 'genome_selected_for_run',
        'ortholog_selected_for_run', 'evaluated_in_final_macse_validation',
        'changed_from_previous', 'pk_unique', 'composite_key_unique',
        'qc_eligible', 'genomic_copy_selected',
    }
    return (column in numeric_exact or column in boolean_exact
            or column.endswith('_pk') or column.startswith((
                'mean_', 'weighted_mean_', 'sd_', 'var_', 'median_', 'mad_',
                'iqr_', 'q05_', 'q25_', 'q75_', 'q95_')))


def validate_sql_null_serialization(stage):
    """Stream-check generated staged TSVs for consistent SQL NULL encoding.

    This intentionally ignores legacy TSVs copied forward by commit_staging.
    It runs only on files created and registered by the v9 writers, checking
    exact registered headers, every data-row width and forbidden empty cells.
    """
    generated = {p.name: columns for p, columns in GENERATED_TSV_COLUMNS.items()
                 if p.parent == stage.resolve()}
    missing_files = sorted(REQUIRED_V9_TSVS - set(generated))
    if missing_files:
        raise ValueError(f'v9 did not generate expected TSVs: {missing_files}')
    staged_files = {p.name for p in stage.glob('*.tsv')}
    unregistered = sorted(staged_files - set(generated))
    if unregistered:
        raise ValueError(f'Staged TSVs lack registered expected headers: {unregistered}')
    missing_on_disk = sorted(set(generated) - staged_files)
    if missing_on_disk:
        raise ValueError(f'Registered TSVs missing on disk: {missing_on_disk}')

    file_count = null_count = row_count = blank_count = 0
    for name, expected_columns in sorted(generated.items()):
        file_count += 1
        path = stage/name
        with path.open('r', newline='', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader, None)
            if header != expected_columns:
                raise ValueError(f'{name}: output header/order mismatch: '
                                 f'{header!r} != {expected_columns!r}')
            if len(header) != len(set(header)):
                raise ValueError(f'{name}: duplicate header column names')
            typed_indices = [i for i, col in enumerate(header)
                             if is_typed_null_field(col)]
            for row in reader:
                row_count += 1
                if len(row) != len(header):
                    raise ValueError(f'{name}: line {reader.line_num} has '
                                     f'{len(row)} fields, expected {len(header)}')
                for index, value in enumerate(row):
                    if not value.strip():
                        blank_count += 1
                        raise ValueError(f'{name}: line {reader.line_num}, '
                                         f'column {header[index]}: empty field; '
                                         r'expected \N for SQL NULL')
                    if value == NULL_SQL:
                        null_count += 1
                for index in typed_indices:
                    if row[index].strip().lower() in ('nan', 'none', 'null'):
                        raise ValueError(f'{name}: line {reader.line_num}, '
                                         f'column {header[index]}: unexpected '
                                         f'non-SQL NULL marker {row[index]!r}')
    print('SQL NULL serialization validation passed: '
          f'{file_count} generated TSVs; {row_count} data rows; '
          f'{null_count} \\N cells; {blank_count} unexpected empty fields')
    return {'files':file_count, 'rows':row_count,
            'null_cells':null_count, 'unexpected_blank_fields':blank_count}


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Compleasm feature TSVs with initial and run-specific QC, preserving v9 calculations and SQL NULLs (v10 duplicate resolution).")
    parser.add_argument("--genomes", required=True, help="Path to genomes directory")
    parser.add_argument("--manifest", required=True, help="Manifest CSV/TSV with accession column")
    parser.add_argument("--orthodb", required=True, help="sauropsida_odb12 tar.gz archive or extracted directory containing links_to_ODB12.txt")
    parser.add_argument("--sequences-tsv", required=True, help="Base sequences.tsv lookup table")
    parser.add_argument("--genomes-tsv", required=True, help="Base genomes.tsv lookup table")
    parser.add_argument("--ortholog-validation", required=True, help="Expanded 17-column ortholog_validity.tsv keyed by accession + odb12_id")
    parser.add_argument("--compleasm-metadata", default=None, help="Optional Compleasm metadata.csv path")
    parser.add_argument("--genomes-metadata", default=None, help="Optional genomes_metadata.csv path with path_to_fna")
    parser.add_argument("--outdir", default=None, help="Output directory; default: <genomes>/records/sql_tsvs/compleasm_features")
    parser.add_argument("--test", action="store_true", help="Restrict execution to two genomes; default output moves to compleasm_features_test")
    parser.add_argument("--run-config", default=None, help="Run-specific QC config TSV; optional for legacy-only builds")
    parser.add_argument("--nhphyml-runs", default=None, help="nhphyml_runs.tsv FK reference; defaults to <genomes>/records/sql_tsvs/nhphyml_runs.tsv when run config supplied")
    parser.add_argument("--strict-inclusion", action="store_true", help="Fail when a final QC-passing sequence is confirmed absent from nhPhyML input")
    args = parser.parse_args()
    # Reset v10 diagnostics when main() is invoked more than once in a process.
    DUPLICATE_AUDIT_ROWS.clear()
    DUPLICATE_COUNTS.clear()
    GENERATED_TSV_COLUMNS.clear()

    genomes_dir = Path(args.genomes).resolve()
    manifest_path = Path(args.manifest).resolve()
    orthodb_path = Path(args.orthodb).resolve()
    sequences_tsv = Path(args.sequences_tsv).resolve()
    genomes_tsv = Path(args.genomes_tsv).resolve()
    ortholog_validation_tsv = Path(args.ortholog_validation).resolve()
    compleasm_metadata = Path(args.compleasm_metadata).resolve() if args.compleasm_metadata else genomes_dir / "records/compleasm/records/metadata.csv"
    genomes_metadata = Path(args.genomes_metadata).resolve() if args.genomes_metadata else genomes_dir / "records/genomes_metadata.csv"
    outdir = Path(args.outdir).resolve() if args.outdir else genomes_dir / ('records/sql_tsvs/compleasm_features_test' if args.test else 'records/sql_tsvs/compleasm_features')
    outdir.parent.mkdir(parents=True, exist_ok=True)
    # Old outputs remain unchanged throughout source loading, calculations,
    # validation, and staging. Commit only after all tables pass checks.
    stage = Path(tempfile.mkdtemp(prefix='.compleasm_features_v10_stage_', dir=outdir.parent))
    global CURRENT_STAGING
    CURRENT_STAGING = stage

    manifest_accessions = load_manifest_accessions(manifest_path)
    metadata_rows, missing_roots = resolve_metadata_rows(compleasm_metadata, manifest_accessions)
    if missing_roots:
        raise ValueError(f'Compleasm metadata lacks {len(missing_roots)} manifest accessions: {missing_roots[:20]}')
    metadata_rows = attach_genome_fna(metadata_rows, genomes_metadata, genomes_dir)
    metadata_rows = metadata_rows.sort_values('accession', kind='stable').reset_index(drop=True)
    if args.test:
        metadata_rows = metadata_rows.head(2).copy()
    scope = set(metadata_rows['accession'])
    check_lookup_table_uniqueness(genomes_tsv, ['genome_pk','accession_id'], ['accession_id'], 'genome_pk')
    check_lookup_table_uniqueness(sequences_tsv, ['genome_pk','sequence_pk','sequence_id'], ['genome_pk','sequence_id'], 'sequence_pk')
    genomes_lookup = load_genomes_lookup(genomes_tsv)
    sequences_lookup = load_sequences_lookup(sequences_tsv)
    enforce_one_to_one_genome_mappings(metadata_rows, genomes_lookup, scope)
    ortholog_validation_lookup = read_initial_qc(ortholog_validation_tsv, scope)
    manifest_taxa = read_manifest_taxa(manifest_path, scope)
    configs = []
    if args.run_config:
        runs_file = Path(args.nhphyml_runs).resolve() if args.nhphyml_runs else genomes_dir / 'records/sql_tsvs/nhphyml_runs.tsv'
        runs_metadata = read_nhphyml_runs(runs_file)
        configs = read_run_config(Path(args.run_config).resolve(), runs_metadata)
    else:
        print('WARNING: --run-config not supplied; ortholog_run_qc.tsv and selection audit will be empty (headers only).')
    sauropsida_rows = parse_orthodb_archive(orthodb_path)

    print(f"Loaded {len(manifest_accessions)} accessions from manifest")
    print(f"Resolved {len(metadata_rows)} Compleasm metadata rows")
    print(f"Loaded {len(ortholog_validation_lookup)} initial QC records scoped to processed accessions")
    print(f"Missing accession roots in Compleasm metadata: {len(missing_roots)}")
    if args.test:
        print("TEST MODE: processing first 2 resolved genomes")

    ortholog_rows: List[Dict[str, object]] = []
    ortholog_source_audit_rows: List[Dict[str, object]] = []
    ortholog_summary_rows: List[Dict[str, object]] = []
    intron_rows: List[Dict[str, object]] = []
    intron_summary_rows: List[Dict[str, object]] = []
    flank_rows: List[Dict[str, object]] = []
    flank_summary_rows: List[Dict[str, object]] = []
    flank_set_rows = [
        {"flank_set_pk": pk, "name": name, "upstream_bp": up, "downstream_bp": down, "created_at": now_str(), "notes": "Compleasm flank set; upstream and downstream concatenated for GC"}
        for pk, name, up, down in FLANK_SPECS
    ]

    ortholog_pk = intron_pk = flank_pk = 0
    ortholog_summary_pk = intron_summary_pk = flank_summary_pk = 0
    ortholog_key_to_pk: Dict[Tuple[int, str], int] = {}
    used_initial_qc_keys = set()
    initial_qc_rows = []
    ortholog_run_index = []
    ortholog_qc_by_key = {}
    diagnostics = Counter()
    final_qc_counts = Counter()

    for _, meta in metadata_rows.iterrows():
        accession = str(meta["accession"]).strip()
        lookup = genomes_lookup.get(accession)
        if lookup is None:
            raise ValueError(f'{accession}: missing genome_pk in genomes.tsv')
        genome_pk = int(lookup["genome_pk"])
        species_pk = lookup.get("species_pk")

        cds_fasta = repath(meta["cds_fasta"], genomes_dir)
        full_table = repath(meta["full_table"], genomes_dir)
        genome_fasta = Path(meta["path_to_fna"])
        if not cds_fasta.is_file() or not full_table.is_file() or not genome_fasta.is_file():
            raise FileNotFoundError(f'{accession}: required Compleasm/genomic FASTA missing; cds={cds_fasta}, full={full_table}, genome={genome_fasta}')

        try:
            cds_sequences = load_cds_sequences(cds_fasta)
            full_records = parse_full_table(
                full_table, accession=accession, genome_pk=genome_pk,
                initial_qc_keys=ortholog_validation_lookup,
            )
            genome = load_genome(genome_fasta)
        except Exception as e:
            raise RuntimeError(f'{accession}: failed loading Compleasm sources: {e}') from e

        genome_valid_orthologs: List[Dict[str, object]] = []
        genome_introns: List[Dict[str, object]] = []
        genome_flanks_by_set: Dict[int, List[Dict[str, object]]] = {pk: [] for pk, _, _, _ in FLANK_SPECS}

        all_ortholog_ids = sorted(set(cds_sequences) | set(full_records))
        missing_sequence_debug_count = 0
        genome_unmapped_with_cds_gc = 0
        for odb12_id in all_ortholog_ids:
            record = full_records.get(odb12_id)
            seq = cds_sequences.get(odb12_id)
            key = (accession, odb12_id)
            if record is None or not record.status:
                raise ValueError(f'{accession} / {odb12_id}: no unambiguous Compleasm full_table status; cannot determine QC eligibility')
            eligible = status_allowed_for_features(accession, record.status)
            qc_tuple = ortholog_validation_lookup.get(key)
            if eligible and qc_tuple is None:
                raise ValueError(f'{accession} / {odb12_id}: QC-eligible {record.status} ortholog missing required initial QC record')
            if not eligible and qc_tuple is not None:
                raise ValueError(f'{accession} / {odb12_id}: QC-ineligible status {record.status!r} has contradictory initial QC record at line {qc_tuple[0]}')
            qc = qc_tuple[1] if qc_tuple else {c:None for c in INITIAL_QC_COLUMNS}
            if qc_tuple:
                used_initial_qc_keys.add(key)
            diagnostics['qc_eligible' if eligible else 'qc_ineligible'] += 1
            if qc_tuple:
                diagnostics['qc_matched'] += 1
                diagnostics['raw_'+qc['raw_cds_qc']] += 1
                if qc['macse_qc'] is None:
                    diagnostics['without_macse'] += 1
                if qc['final_qc']:
                    final_qc_counts[qc['final_qc']] += 1
            passes_qc = qc['raw_cds_qc'] == 'pass'
            if record.mapping_state != 'resolved' and (eligible or qc_tuple is not None):
                raise ValueError(f'{accession} / {odb12_id}: unresolved genomic '
                                 'copies must be ineligible and have no initial QC record')
            sequence_pk = None
            if record is not None and record.mapping_state == 'resolved':
                sequence_pk = sequences_lookup.get((genome_pk, record.sequence_id))
            
                if passes_qc and sequence_pk is None and missing_sequence_debug_count < 10:
                    print(
                        f"INFO {accession}: "
                        f"genome_pk={genome_pk} "
                        f"sequence_id='{record.sequence_id}' "
                        f"not found in sequences.tsv; retaining the ortholog "
                        "and QC-passing CDS GC values with SQL NULL genomic fields"
                    )
                    missing_sequence_debug_count += 1
            gc = gc3 = gc4 = None
            ortholog_length = None
            if passes_qc and seq:
                ortholog_length = len(seq)
                gc = calculate_gc(seq)[0]
                try:
                    gc3 = calculate_gc3(seq)
                except ValueError:
                    gc3 = None
                try:
                    gc4 = calculate_gc4(seq)
                except ValueError:
                    gc4 = None

            ortholog_pk += 1
            this_pk = ortholog_pk
            if (genome_pk, odb12_id) in ortholog_key_to_pk:
                raise ValueError(f'{accession} / {odb12_id}: duplicate ortholog identity in construction')
            ortholog_key_to_pk[(genome_pk, odb12_id)] = this_pk
            mapping_status, coordinate_basis, mapping_notes = describe_sequence_mapping(record, sequence_pk)
            genomic_mapped = sequence_pk is not None
            included_in_summary = bool(
                passes_qc and seq and record
                and status_allowed_for_features(accession, record.status)
            )
            genomic_features_eligible = bool(
                included_in_summary and genomic_mapped and record.exons
            )
            if not genomic_mapped and gc is not None:
                genome_unmapped_with_cds_gc += 1
            ortholog_rows.append({
                "ortholog_pk": this_pk,
                "sequence_pk": sequence_pk,
                "genome_pk": genome_pk,
                "odb12_id": odb12_id,
                "status": record.status if record else None,
                "strand": record.strand if record and genomic_mapped else None,
                "start": record.start if record and genomic_mapped else None,
                "end": record.end if record and genomic_mapped else None,
                "gc": gc,
                "gc3": gc3,
                "gc4": gc4,
            })
            initial_qc_rows.append({
                'ortholog_pk':this_pk,'genome_pk':genome_pk,
                'sequence_pk':sequence_pk,'odb12_id':odb12_id,
                **{c:qc.get(c) for c in INITIAL_QC_COLUMNS[2:]},
            })
            ortholog_run_index.append({
                'ortholog_pk':this_pk, 'genome_pk':genome_pk,
                'accession':accession,'odb12_id':odb12_id,
                'mapping_status':mapping_status,
            })
            ortholog_qc_by_key[key] = qc
            if sequence_pk is None:
                diagnostics['unmapped_sequence_pk'] += 1
            ortholog_source_audit_rows.append({
                "ortholog_pk": this_pk,
                "genome_pk": genome_pk,
                "accession": accession,
                "odb12_id": odb12_id,
                "sequence_pk": sequence_pk,
                "mapping_status": mapping_status,
                "source_sequence_id": record.sequence_id if record else None,
                "source_strand": record.strand if record else None,
                "source_start": record.start if record else None,
                "source_end": record.end if record else None,
                "source_coordinate_basis": coordinate_basis,
                # Gene Start/End and Codons are copied as parsed. They are
                # deliberately not normalized to one shared coordinate origin.
                "source_codon_intervals": "|".join(
                    f"{exon['start']}_{exon['end']}_{exon['strand']}"
                    for exon in record.exons
                ) if record and record.mapping_state == 'resolved' else None,
                "passes_raw_cds_qc": passes_qc,
                "cds_present": bool(seq),
                "gc_calculated": gc is not None,
                "included_in_ortholog_summary": included_in_summary,
                "genomic_features_eligible": genomic_features_eligible,
                "mapping_notes": mapping_notes,
            })

            if not included_in_summary:
                continue
            
            valid_summary_row = {
                "ortholog_pk": this_pk,
                "gc": gc,
                "gc3": gc3,
                "gc4": gc4,
                "length": ortholog_length,
                "genomic_mapped": genomic_mapped,
            }
            genome_valid_orthologs.append(valid_summary_row)

            if not genomic_features_eligible:
                continue

            try:
                for intron in extract_introns(genome, record):
                    intron_pk += 1
                    row = {
                        "intron_pk": intron_pk,
                        "ortholog_pk": this_pk,
                        "parent_id": odb12_id,
                        "intron_id": intron["intron_id"],
                        "status": record.status,
                        "strand": record.strand,
                        "start": intron["start"],
                        "end": intron["end"],
                        "length": intron["length"],
                        "gc": intron["gc"],
                    }
                    intron_rows.append(row)
                    genome_introns.append(row)
            except Exception as e:
                print(f"WARNING: {accession} {odb12_id}: intron extraction failed: {e}")

            for flank_set_pk, _name, up_bp, down_bp in FLANK_SPECS:
                try:
                    flank = extract_flank(genome, record, up_bp, down_bp)
                except Exception as e:
                    print(f"WARNING: {accession} {odb12_id}: flank extraction failed: {e}")
                    continue
                flank_pk += 1
                row = {
                    "flank_pk": flank_pk,
                    "ortholog_pk": this_pk,
                    "flank_set_pk": flank_set_pk,
                    "status": record.status,
                    "strand": record.strand,
                    "up_start": flank["up_start"],
                    "up_end": flank["up_end"],
                    "down_start": flank["down_start"],
                    "down_end": flank["down_end"],
                    "gc": flank["gc"],
                    "length": flank["length"],
                    "upstream_gc": flank["upstream_gc"],
                    "downstream_gc": flank["downstream_gc"],
                }
                flank_rows.append(row)
                genome_flanks_by_set[flank_set_pk].append(row)

        if genome_valid_orthologs:
            ortholog_summary_pk += 1
            gcs = [r["gc"] for r in genome_valid_orthologs]
            gc3s = [r["gc3"] for r in genome_valid_orthologs]
            gc4s = [r["gc4"] for r in genome_valid_orthologs]
            lengths = [r["length"] for r in genome_valid_orthologs]
            ortholog_summary_rows.append({
                "ortholog_summary_pk": ortholog_summary_pk,
                "genome_pk": genome_pk,
                "species_pk": species_pk,
                "n_orthologs": len(genome_valid_orthologs),
                "callable_bp_total": sum(int(x) for x in lengths if x),
                **stats(gcs, lengths, "gc"),
                **stats(gc3s, lengths, "gc3"),
                **stats(gc4s, lengths, "gc4"),
                "q05_gc3": q(gc3s, 0.05),
                "q25_gc3": q(gc3s, 0.25),
                "q75_gc3": q(gc3s, 0.75),
                "q95_gc3": q(gc3s, 0.95),
                "mean_ortholog_length": mean([x for x in lengths if x]) if lengths else None,
                "median_ortholog_length": median([x for x in lengths if x]) if lengths else None,
                "created_at": now_str(),
            })

        n_summary_mapped = sum(bool(row["genomic_mapped"]) for row in genome_valid_orthologs)
        mapping_summary_note = (
            f"{n_summary_mapped}/{len(genome_valid_orthologs)} CDS-summary orthologs "
            "linked to genomic sequences. Unlinked records retain CDS GC values "
            "but do not contribute genomic introns/flanks."
        )
        if not genome_valid_orthologs:
            no_features_note = "No QC-passing status-eligible CDS orthologs for genomic features."
        elif n_summary_mapped == 0:
            no_features_note = "Genomic sequence_pk could not be mapped; CDS GC summaries retained."
        else:
            no_features_note = (
                "No feature records recovered from the available mapped Compleasm intervals; "
                "see ortholog_source_audit.tsv and extraction messages."
            )

        if genome_introns:
            intron_summary_pk += 1
            vals = [r["gc"] for r in genome_introns]
            lengths = [r["length"] for r in genome_introns]
            intron_summary_rows.append({
                "intron_summary_pk": intron_summary_pk,
                "genome_pk": genome_pk,
                "species_pk": species_pk,
                "n_introns": len(genome_introns),
                "n_orthologs": len(set(r["ortholog_pk"] for r in genome_introns)),
                "callable_bp_total": sum(int(x) for x in lengths if x),
                **stats(vals, lengths, "gc"),
                "q05_gc": q(vals, 0.05),
                "q25_gc": q(vals, 0.25),
                "q75_gc": q(vals, 0.75),
                "q95_gc": q(vals, 0.95),
                "mean_intron_length": mean(lengths),
                "median_intron_length": median(lengths),
                "recovery_status": "recovered",
                "recovery_notes": mapping_summary_note,
                "created_at": now_str(),
            })
        else:
            intron_summary_pk += 1
            intron_summary_rows.append({
                "intron_summary_pk": intron_summary_pk,
                "genome_pk": genome_pk,
                "species_pk": species_pk,
                "n_introns": 0,
                "n_orthologs": 0,
                "callable_bp_total": 0,
                "recovery_status": "not_recoverable_from_compleasm",
                "recovery_notes": no_features_note,
                "created_at": now_str(),
            })

        for flank_set_pk, rows_for_set in genome_flanks_by_set.items():
            if not rows_for_set:
                flank_summary_pk += 1
                flank_summary_rows.append({
                    "flank_summary_pk": flank_summary_pk,
                    "genome_pk": genome_pk,
                    "species_pk": species_pk,
                    "flank_set_pk": flank_set_pk,
                    "n_flanks": 0,
                    "n_orthologs": 0,
                    "callable_bp_total": 0,
                    "recovery_status": "not_recoverable_from_compleasm",
                    "recovery_notes": no_features_note,
                    "created_at": now_str(),
                })
                continue
            flank_summary_pk += 1
            vals = [r["gc"] for r in rows_for_set]
            lengths = [r["length"] for r in rows_for_set]
            flank_summary_rows.append({
                "flank_summary_pk": flank_summary_pk,
                "genome_pk": genome_pk,
                "species_pk": species_pk,
                "flank_set_pk": flank_set_pk,
                "n_flanks": len(rows_for_set),
                "n_orthologs": len(set(r["ortholog_pk"] for r in rows_for_set)),
                "callable_bp_total": sum(int(x) for x in lengths if x),
                **stats(vals, lengths, "gc"),
                "q05_gc": q(vals, 0.05),
                "q25_gc": q(vals, 0.25),
                "q75_gc": q(vals, 0.75),
                "q95_gc": q(vals, 0.95),
                "mean_flank_length": mean(lengths),
                "median_flank_length": median(lengths),
                "mean_upstream_gc": mean(numeric_values([r["upstream_gc"] for r in rows_for_set])) if numeric_values([r["upstream_gc"] for r in rows_for_set]) else None,
                "mean_downstream_gc": mean(numeric_values([r["downstream_gc"] for r in rows_for_set])) if numeric_values([r["downstream_gc"] for r in rows_for_set]) else None,
                "recovery_status": "recovered",
                "recovery_notes": mapping_summary_note,
                "created_at": now_str(),
            })

        allowed_statuses = (
            "Single or Duplicated/Duplicate"
            if accession in DUPLICATED_STATUS_ALLOWED_ACCESSIONS else "Single"
        )
        print(
            f"Processed {accession}: {len(genome_valid_orthologs)} valid {allowed_statuses} "
            f"orthologs in CDS summaries; {n_summary_mapped} linked to genomic sequences; "
            f"{genome_unmapped_with_cds_gc} unlinked ortholog rows retain CDS GC; "
            f"{len(genome_introns)} genomic introns; "
            f"{sum(len(rows) for rows in genome_flanks_by_set.values())} genomic flank rows"
        )

    unused_qc = set(ortholog_validation_lookup)-used_initial_qc_keys
    if unused_qc:
        raise ValueError(f'{len(unused_qc)} manifest/test-scoped initial QC records unmatched to constructed orthologs; examples: {sorted(unused_qc)[:20]}')
    if len(ortholog_rows) != len(initial_qc_rows) or len(ortholog_key_to_pk) != len(ortholog_rows):
        raise ValueError('Initial QC/ortholog one-to-one relational integrity failure')
    if [r['ortholog_pk'] for r in ortholog_rows] != list(range(1,len(ortholog_rows)+1)):
        raise ValueError('Non-consecutive ortholog primary keys')
    ortholog_ids = {r['ortholog_pk'] for r in ortholog_rows}
    for table_name, rows in [('intron',intron_rows),('flank',flank_rows)]:
        for row in rows:
            if row['ortholog_pk'] not in ortholog_ids:
                raise ValueError(f'{table_name}: unknown ortholog_pk {row["ortholog_pk"]}')
    print('Initial QC: total=',len(ortholog_rows), 'eligible=',diagnostics['qc_eligible'],
          'ineligible=',diagnostics['qc_ineligible'], 'matched=',diagnostics['qc_matched'],
          'raw pass=',diagnostics['raw_pass'], 'raw fail=',diagnostics['raw_fail'],
          'without MACSE=',diagnostics['without_macse'],
          'unmapped genomic sequence_pk=',diagnostics['unmapped_sequence_pk'],
          'initial final QC categories=',dict(final_qc_counts))
    write_tsv(sauropsida_rows, stage / "sauropsida_odb12.tsv", ["odb12_id", "name", "orthodb"])
    write_tsv(ortholog_rows, stage / "orthologs.tsv", ORTHOLOG_V8_COLUMNS)
    write_tsv(DUPLICATE_AUDIT_ROWS, stage / 'ortholog_duplicate_full_table_audit.tsv',
              DUPLICATE_AUDIT_COLUMNS)
    write_sql_tsv(initial_qc_rows, stage/'ortholog_initial_qc.tsv', INITIAL_OUTPUT_COLUMNS)
    write_tsv(ortholog_source_audit_rows, stage / "ortholog_source_audit.tsv", [
        "ortholog_pk", "genome_pk", "accession", "odb12_id", "sequence_pk", "mapping_status",
        "source_sequence_id", "source_strand", "source_start", "source_end",
        "source_coordinate_basis", "source_codon_intervals", "passes_raw_cds_qc",
        "cds_present", "gc_calculated", "included_in_ortholog_summary",
        "genomic_features_eligible", "mapping_notes",
    ])
    write_tsv(ortholog_summary_rows, stage / "ortholog_summary.tsv", [
        "ortholog_summary_pk", "genome_pk", "species_pk", "n_orthologs", "callable_bp_total",
        "mean_gc", "weighted_mean_gc", "mean_gc3", "weighted_mean_gc3", "mean_gc4", "weighted_mean_gc4",
        "sd_gc", "sd_gc3", "sd_gc4", "var_gc", "var_gc3", "var_gc4", "median_gc", "median_gc3", "median_gc4",
        "mad_gc", "mad_gc3", "mad_gc4", "iqr_gc", "iqr_gc3", "iqr_gc4",
        "q05_gc3", "q25_gc3", "q75_gc3", "q95_gc3", "mean_ortholog_length", "median_ortholog_length", "created_at"
    ])
    write_tsv(intron_rows, stage / "intron_compleasm.tsv", ["intron_pk", "ortholog_pk", "parent_id", "intron_id", "status", "strand", "start", "end", "length", "gc"])
    write_tsv(intron_summary_rows, stage / "intron_compleasm_summary.tsv", [
        "intron_summary_pk", "genome_pk", "species_pk", "n_introns", "n_orthologs", "callable_bp_total",
        "mean_gc", "weighted_mean_gc", "sd_gc", "var_gc", "median_gc", "mad_gc", "iqr_gc",
        "q05_gc", "q25_gc", "q75_gc", "q95_gc", "mean_intron_length", "median_intron_length", "recovery_status", "recovery_notes", "created_at"
    ])
    write_tsv(flank_rows, stage / "flanks_compleasm.tsv", ["flank_pk", "ortholog_pk", "flank_set_pk", "status", "strand", "up_start", "up_end", "down_start", "down_end", "gc"])
    write_tsv(flank_set_rows, stage / "flank_sets_compleasm.tsv", ["flank_set_pk", "name", "upstream_bp", "downstream_bp", "created_at", "notes"])
    write_tsv(flank_summary_rows, stage / "flank_compleasm_summary.tsv", [
        "flank_summary_pk", "genome_pk", "species_pk", "flank_set_pk", "n_flanks", "n_orthologs", "callable_bp_total",
        "mean_gc", "weighted_mean_gc", "sd_gc", "var_gc", "median_gc", "mad_gc", "iqr_gc",
        "q05_gc", "q25_gc", "q75_gc", "q95_gc", "mean_flank_length", "median_flank_length",
        "mean_upstream_gc", "mean_downstream_gc", "recovery_status", "recovery_notes", "created_at"
    ])
    # Compatibility alias requested in the prompt wording.
    write_tsv(flank_summary_rows, stage / "flank_summary.tsv", [
        "flank_summary_pk", "genome_pk", "species_pk", "flank_set_pk", "n_flanks", "n_orthologs", "callable_bp_total",
        "mean_gc", "weighted_mean_gc", "sd_gc", "var_gc", "median_gc", "mad_gc", "iqr_gc",
        "q05_gc", "q25_gc", "q75_gc", "q95_gc", "mean_flank_length", "median_flank_length",
        "mean_upstream_gc", "mean_downstream_gc", "created_at"
    ])

    print('Processing run-specific QC ...')
    run_summaries = write_run_tables(stage, configs, ortholog_run_index,
                                    ortholog_qc_by_key, manifest_taxa,
                                    strict_inclusion=args.strict_inclusion, test=args.test)
    # Ensure new table schema / keys and v7 retained table schemas are intact.
    expected_output_headers = {
        'orthologs.tsv': ORTHOLOG_V8_COLUMNS,
        'ortholog_initial_qc.tsv': INITIAL_OUTPUT_COLUMNS,
        'ortholog_run_qc.tsv': RUN_OUTPUT_COLUMNS,
        'ortholog_run_selection_audit.tsv': AUDIT_COLUMNS,
    }
    for filename, cols in expected_output_headers.items():
        required_header(stage/filename, cols)
        with (stage/filename).open(encoding='utf-8') as handle:
            actual_header = next(csv.reader(handle,delimiter='\t'))
        if actual_header != cols:
            raise ValueError(f'{filename}: output column order mismatch: {actual_header}')
    assert sum(s['run_qc_primary_keys'] for s in run_summaries) == sum(s['evaluated_accession_ortholog_pairs'] for s in run_summaries)
    validate_staged_relations(stage, len(ortholog_rows),
                              sum(s['run_qc_primary_keys'] for s in run_summaries))
    validate_v10_duplicate_integrity(stage, ortholog_rows, initial_qc_rows,
                                     intron_rows, flank_rows,
                                     ortholog_source_audit_rows)
    compare_prior_summary_counts(stage, outdir)
    null_validation_summary = validate_sql_null_serialization(stage)
    commit_staging(stage, outdir)
    CURRENT_STAGING = None
    print('v10 serialization checks passed: '
          f'{null_validation_summary["files"]} TSVs validated, '
          f'{null_validation_summary["null_cells"]} SQL NULL markers, '
          f'{null_validation_summary["unexpected_blank_fields"]} unexpected empty cells; '
          'relational integrity passed')
    print('Compleasm duplicate-resolution summary:')
    print(f'  Parsed full_table records: {DUPLICATE_COUNTS["source_records"]}')
    print(f'  Unique accession/ODB12 identities: {DUPLICATE_COUNTS["unique_ortholog_identities"]}')
    print(f'  Identities with repeated records: {DUPLICATE_COUNTS["repeated_groups"]}')
    print(f'  Identical groups collapsed: {DUPLICATE_COUNTS["identical_groups_collapsed"]}')
    print(f'  Ambiguous ineligible groups retained: {DUPLICATE_COUNTS["ambiguous_ineligible_retained"]}')
    print(f'  Ambiguous eligible groups rejected: {DUPLICATE_COUNTS["ambiguous_eligible_rejected"]} (successful build requires 0)')
    print(f'  Duplicate-audit source-copy records: {len(DUPLICATE_AUDIT_ROWS)}')
    print(f'  Orthologs with deliberately unresolved coordinates: {DUPLICATE_COUNTS["deliberately_unresolved"]}')
    print('  v13 relational integrity: passed')
    print('  SQL NULL serialization: passed')
    print("\nFinished building Compleasm feature TSVs.")
    print(f"Output directory: {outdir}")
    print(f"sauropsida_odb12 rows: {len(sauropsida_rows)}")
    print(f"orthologs rows: {len(ortholog_rows)}")
    print(f"ortholog_source_audit rows: {len(ortholog_source_audit_rows)}")
    print(
        "Ortholog rows retaining CDS GC without sequence_pk: "
        f"{sum(row['sequence_pk'] is None and row['gc'] is not None for row in ortholog_rows)}"
    )
    print("ortholog_source_audit.tsv is a helper file; do not load it as a SQL table.")
    print(f"ortholog_summary rows: {len(ortholog_summary_rows)}")
    print(f"intron_compleasm rows: {len(intron_rows)}")
    print(f"intron_compleasm_summary rows: {len(intron_summary_rows)}")
    print(f"flanks_compleasm rows: {len(flank_rows)}")
    print(f"flank_sets_compleasm rows: {len(flank_set_rows)}")
    print(f"flank_compleasm_summary rows: {len(flank_summary_rows)}")


CURRENT_STAGING = None

if __name__ == "__main__":
    try:
        main()
    except BaseException as error:
        # Preserve existing output directory; diagnostics live alongside it.
        import traceback
        print(f'ERROR: {type(error).__name__}: {error}',file=sys.stderr)
        if CURRENT_STAGING is not None and CURRENT_STAGING.exists():
            shutil.rmtree(CURRENT_STAGING)
        try:
            argv = sys.argv
            if '--outdir' in argv:
                dest = Path(argv[argv.index('--outdir')+1]).resolve()
            elif '--genomes' in argv:
                dest = Path(argv[argv.index('--genomes')+1]).resolve()/('records/sql_tsvs/compleasm_features_test' if '--test' in argv else 'records/sql_tsvs/compleasm_features')
            else:
                dest = Path.cwd()/'compleasm_features'
            write_failure_details(dest, error)
        except Exception as reporting_error:
            print(f'Could not write failure diagnostics: {reporting_error}', file=sys.stderr)
        raise SystemExit(1)

