#!/usr/bin/env python3
"""
04_build_compleasm_feature_tsvs_v3.py

Build Compleasm-derived feature TSVs for the gc3_dynamics_v6 schema.

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

Notes:
    - Coordinates written to TSV are 1-based inclusive.
    - Introns are inferred from Compleasm full_table Codons/CDS coordinate tokens.
    - Flanks are clipped at sequence ends and may overlap neighboring genes.
    - GC ignores N/ambiguous bases in the denominator.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import tarfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from statistics import mean, median
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pyfaidx
from Bio.SeqIO.FastaIO import SimpleFastaParser


STOP_CODONS = {"TAA", "TAG", "TGA"}
FLANK_SPECS = [
    (1, "50_up_50_down", 50, 50),
    (2, "100_up_100_down", 100, 100),
    (3, "200_up_200_down", 200, 200),
]


@dataclass
class FullRecord:
    odb12_id: str
    status: str
    sequence_id: str
    strand: str
    start: Optional[int]
    end: Optional[int]
    exons: List[Dict[str, object]]


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

def blank_if_none(value):
    if value is None:
        return ""
    try:
        if isinstance(value, float) and math.isnan(value):
            return ""
    except TypeError:
        pass
    return value


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
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: blank_if_none(row.get(col)) for col in columns})


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
            seqs[header.split()[0]] = seq.upper().replace("-", "")
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


def parse_full_table(full_table: Path) -> Dict[str, FullRecord]:
    with open(full_table, "r") as handle:
        header: Optional[List[str]] = None
        raw_rows: List[List[str]] = []
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None and parts[0].lower() == "gene":
                header = parts
                continue
            raw_rows.append(parts)
    if header is None:
        # Common Compleasm full_table positional layout used by previous script.
        header = ["Gene", "Status", "Sequence", "Score", "Length", "Strand", "Col7", "Col8", "Col9", "Gene Start", "Gene End", "Col12", "Codons"]
    records: Dict[str, FullRecord] = {}
    colmap = {name.lower().strip(): i for i, name in enumerate(header)}

    def idx(*names: str, default: Optional[int] = None) -> Optional[int]:
        for name in names:
            key = name.lower().strip()
            if key in colmap:
                return colmap[key]
        compact = {h.lower().replace(" ", "").replace("_", ""): i for i, h in enumerate(header)}
        for name in names:
            key = name.lower().replace(" ", "").replace("_", "")
            if key in compact:
                return compact[key]
        return default

    gene_i = idx("Gene", default=0)
    status_i = idx("Status", default=1)
    sequence_i = idx("Sequence", "Chromosome", "Contig", default=2)
    strand_i = idx("Strand", default=5)
    start_i = idx("Gene Start", "Start", default=9)
    end_i = idx("Gene End", "End", default=10)
    codons_i = idx("Codons", "Coordinates", "Exons", default=12)

    for parts in raw_rows:
        if len(parts) <= max(gene_i or 0, status_i or 0, sequence_i or 0, strand_i or 0):
            continue
        odb12_id = parts[gene_i].strip()
        if not odb12_id:
            continue
        status = parts[status_i].strip() if status_i is not None and status_i < len(parts) else ""
        sequence_id = parts[sequence_i].strip() if sequence_i is not None and sequence_i < len(parts) else ""
        strand = parts[strand_i].strip() if strand_i is not None and strand_i < len(parts) else ""
        exons: List[Dict[str, object]] = []
        if codons_i is not None and codons_i < len(parts):
            for token in str(parts[codons_i]).split("|"):
                exon = parse_coordinate_token(token, odb12_id)
                if exon and (not strand or exon["strand"] == strand):
                    exons.append(exon)
        exons.sort(key=lambda e: int(e["start"]))
        start: Optional[int] = None
        end: Optional[int] = None
        for source_i, attr in [(start_i, "start"), (end_i, "end")]:
            pass
        try:
            if start_i is not None and start_i < len(parts) and str(parts[start_i]).strip():
                start = int(float(parts[start_i]))
            if end_i is not None and end_i < len(parts) and str(parts[end_i]).strip():
                end = int(float(parts[end_i]))
        except ValueError:
            start = None
            end = None
        if start is None and exons:
            start = min(int(e["start"]) for e in exons)
        if end is None and exons:
            end = max(int(e["end"]) for e in exons)
        records[odb12_id] = FullRecord(odb12_id, status, sequence_id, strand, start, end, exons)
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


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Compleasm feature TSVs for gc3_dynamics_v6.")
    parser.add_argument("--genomes", required=True, help="Path to genomes directory")
    parser.add_argument("--manifest", required=True, help="Manifest CSV/TSV with accession column")
    parser.add_argument("--orthodb", required=True, help="sauropsida_odb12 tar.gz archive or extracted directory containing links_to_ODB12.txt")
    parser.add_argument("--sequences-tsv", required=True, help="Base sequences.tsv lookup table")
    parser.add_argument("--genomes-tsv", required=True, help="Base genomes.tsv lookup table")
    parser.add_argument("--ortholog-validation", required=True, help="Project-wide ortholog validation TSV with accession, sequence_id, and passes_raw_cds_qc columns")
    parser.add_argument("--compleasm-metadata", default=None, help="Optional Compleasm metadata.csv path")
    parser.add_argument("--genomes-metadata", default=None, help="Optional genomes_metadata.csv path with path_to_fna")
    parser.add_argument("--outdir", default=None, help="Output directory; default: <genomes>/records/sql_tsvs/compleasm_features")
    parser.add_argument("--test", action="store_true", help="Restrict execution to two genomes")
    args = parser.parse_args()

    genomes_dir = Path(args.genomes).resolve()
    manifest_path = Path(args.manifest).resolve()
    orthodb_path = Path(args.orthodb).resolve()
    sequences_tsv = Path(args.sequences_tsv).resolve()
    genomes_tsv = Path(args.genomes_tsv).resolve()
    ortholog_validation_tsv = Path(args.ortholog_validation).resolve()
    compleasm_metadata = Path(args.compleasm_metadata).resolve() if args.compleasm_metadata else genomes_dir / "records/compleasm/records/metadata.csv"
    genomes_metadata = Path(args.genomes_metadata).resolve() if args.genomes_metadata else genomes_dir / "records/genomes_metadata.csv"
    outdir = Path(args.outdir).resolve() if args.outdir else genomes_dir / "records/sql_tsvs/compleasm_features"
    outdir.mkdir(parents=True, exist_ok=True)

    manifest_accessions = load_manifest_accessions(manifest_path)
    metadata_rows, missing_roots = resolve_metadata_rows(compleasm_metadata, manifest_accessions)
    metadata_rows = attach_genome_fna(metadata_rows, genomes_metadata, genomes_dir)
    if args.test:
        metadata_rows = metadata_rows.head(2).copy()

    genomes_lookup = load_genomes_lookup(genomes_tsv)
    sequences_lookup = load_sequences_lookup(sequences_tsv)
    ortholog_validation_lookup = load_ortholog_validation_lookup(ortholog_validation_tsv)
    sauropsida_rows = parse_orthodb_archive(orthodb_path)

    print(f"Loaded {len(manifest_accessions)} accessions from manifest")
    print(f"Resolved {len(metadata_rows)} Compleasm metadata rows")
    print(f"Loaded ortholog-validation entries for {len(ortholog_validation_lookup)} accession/accession-root keys")
    print(f"Missing accession roots in Compleasm metadata: {len(missing_roots)}")
    if args.test:
        print("TEST MODE: processing first 2 resolved genomes")

    ortholog_rows: List[Dict[str, object]] = []
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

    for _, meta in metadata_rows.iterrows():
        accession = str(meta["accession"]).strip()
        lookup = genomes_lookup.get(accession)
        if lookup is None:
            print(f"WARNING: skipping {accession}; not found in --genomes-tsv")
            continue
        genome_pk = int(lookup["genome_pk"])
        species_pk = lookup.get("species_pk")

        cds_fasta = repath(meta["cds_fasta"], genomes_dir)
        full_table = repath(meta["full_table"], genomes_dir)
        genome_fasta = Path(meta["path_to_fna"])
        if not cds_fasta.exists() or not full_table.exists() or not genome_fasta.exists():
            print(f"WARNING: skipping {accession}; missing cds_fasta/full_table/genome_fasta")
            print(f"  cds_fasta: {cds_fasta}")
            print(f"  full_table: {full_table}")
            print(f"  genome_fasta: {genome_fasta}")
            continue
        validity = ortholog_validation_lookup.get(accession)
        if validity is None:
            print(f"WARNING: skipping {accession}; not found in --ortholog-validation")
            continue

        try:
            cds_sequences = load_cds_sequences(cds_fasta)
            full_records = parse_full_table(full_table)
            genome = load_genome(genome_fasta)
        except Exception as e:
            print(f"WARNING: skipping {accession}; failed input parsing/loading: {e}")
            continue

        genome_valid_orthologs: List[Dict[str, object]] = []
        genome_introns: List[Dict[str, object]] = []
        genome_flanks_by_set: Dict[int, List[Dict[str, object]]] = {pk: [] for pk, _, _, _ in FLANK_SPECS}

        all_ortholog_ids = sorted(set(validity) | set(cds_sequences) | set(full_records))
        missing_sequence_debug_count = 0
        for odb12_id in all_ortholog_ids:
            passes_qc = bool(validity.get(odb12_id, False))
            record = full_records.get(odb12_id)
            seq = cds_sequences.get(odb12_id)
            sequence_pk = None
            if record is not None:
                sequence_pk = sequences_lookup.get((genome_pk, record.sequence_id))
            
                if passes_qc and sequence_pk is None and missing_sequence_debug_count < 10:
                    print(
                        f"DEBUG {accession}: "
                        f"genome_pk={genome_pk} "
                        f"sequence_id='{record.sequence_id}' "
                        f"not found in sequences.tsv"
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
            ortholog_key_to_pk[(genome_pk, odb12_id)] = this_pk
            ortholog_rows.append({
                "ortholog_pk": this_pk,
                "sequence_pk": sequence_pk,
                "odb12_id": odb12_id,
                "passes_raw_cds_qc": "TRUE" if passes_qc else "FALSE",
                "status": record.status if record else None,
                "strand": record.strand if record else None,
                "start": record.start if record else None,
                "end": record.end if record else None,
                "gc": gc,
                "gc3": gc3,
                "gc4": gc4,
            })

            if not (passes_qc and seq and record and record.status == "Single"):
                continue
            
            valid_summary_row = {
                "ortholog_pk": this_pk,
                "gc": gc,
                "gc3": gc3,
                "gc4": gc4,
                "length": ortholog_length,
            }
            genome_valid_orthologs.append(valid_summary_row)

            if not (record.exons and sequence_pk is not None):
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
                "recovery_notes": "",
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
                "recovery_notes": "Compleasm coordinates are not genomic sequence IDs; sequence_pk could not be mapped.",
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
                    "recovery_notes": "Compleasm coordinates are not genomic sequence IDs; sequence_pk could not be mapped.",
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
                "recovery_notes": "",
                "created_at": now_str(),
            })

        print(f"Processed {accession}: {len(genome_valid_orthologs)} valid Single orthologs for summaries/features")

    write_tsv(sauropsida_rows, outdir / "sauropsida_odb12.tsv", ["odb12_id", "name", "orthodb"])
    write_tsv(ortholog_rows, outdir / "orthologs.tsv", ["ortholog_pk", "sequence_pk", "odb12_id", "passes_raw_cds_qc", "status", "strand", "start", "end", "gc", "gc3", "gc4"])
    write_tsv(ortholog_summary_rows, outdir / "ortholog_summary.tsv", [
        "ortholog_summary_pk", "genome_pk", "species_pk", "n_orthologs", "callable_bp_total",
        "mean_gc", "weighted_mean_gc", "mean_gc3", "weighted_mean_gc3", "mean_gc4", "weighted_mean_gc4",
        "sd_gc", "sd_gc3", "sd_gc4", "var_gc", "var_gc3", "var_gc4", "median_gc", "median_gc3", "median_gc4",
        "mad_gc", "mad_gc3", "mad_gc4", "iqr_gc", "iqr_gc3", "iqr_gc4",
        "q05_gc3", "q25_gc3", "q75_gc3", "q95_gc3", "mean_ortholog_length", "median_ortholog_length", "created_at"
    ])
    write_tsv(intron_rows, outdir / "intron_compleasm.tsv", ["intron_pk", "ortholog_pk", "parent_id", "intron_id", "status", "strand", "start", "end", "length", "gc"])
    write_tsv(intron_summary_rows, outdir / "intron_compleasm_summary.tsv", [
        "intron_summary_pk", "genome_pk", "species_pk", "n_introns", "n_orthologs", "callable_bp_total",
        "mean_gc", "weighted_mean_gc", "sd_gc", "var_gc", "median_gc", "mad_gc", "iqr_gc",
        "q05_gc", "q25_gc", "q75_gc", "q95_gc", "mean_intron_length", "median_intron_length", "recovery_status", "recovery_notes", "created_at"
    ])
    write_tsv(flank_rows, outdir / "flanks_compleasm.tsv", ["flank_pk", "ortholog_pk", "flank_set_pk", "status", "strand", "up_start", "up_end", "down_start", "down_end", "gc"])
    write_tsv(flank_set_rows, outdir / "flank_sets_compleasm.tsv", ["flank_set_pk", "name", "upstream_bp", "downstream_bp", "created_at", "notes"])
    write_tsv(flank_summary_rows, outdir / "flank_compleasm_summary.tsv", [
        "flank_summary_pk", "genome_pk", "species_pk", "flank_set_pk", "n_flanks", "n_orthologs", "callable_bp_total",
        "mean_gc", "weighted_mean_gc", "sd_gc", "var_gc", "median_gc", "mad_gc", "iqr_gc",
        "q05_gc", "q25_gc", "q75_gc", "q95_gc", "mean_flank_length", "median_flank_length",
        "mean_upstream_gc", "mean_downstream_gc", "recovery_status", "recovery_notes", "created_at"
    ])
    # Compatibility alias requested in the prompt wording.
    write_tsv(flank_summary_rows, outdir / "flank_summary.tsv", [
        "flank_summary_pk", "genome_pk", "species_pk", "flank_set_pk", "n_flanks", "n_orthologs", "callable_bp_total",
        "mean_gc", "weighted_mean_gc", "sd_gc", "var_gc", "median_gc", "mad_gc", "iqr_gc",
        "q05_gc", "q25_gc", "q75_gc", "q95_gc", "mean_flank_length", "median_flank_length",
        "mean_upstream_gc", "mean_downstream_gc", "created_at"
    ])

    print("\nFinished building Compleasm feature TSVs.")
    print(f"Output directory: {outdir}")
    print(f"sauropsida_odb12 rows: {len(sauropsida_rows)}")
    print(f"orthologs rows: {len(ortholog_rows)}")
    print(f"ortholog_summary rows: {len(ortholog_summary_rows)}")
    print(f"intron_compleasm rows: {len(intron_rows)}")
    print(f"intron_compleasm_summary rows: {len(intron_summary_rows)}")
    print(f"flanks_compleasm rows: {len(flank_rows)}")
    print(f"flank_sets_compleasm rows: {len(flank_set_rows)}")
    print(f"flank_compleasm_summary rows: {len(flank_summary_rows)}")


if __name__ == "__main__":
    main()
