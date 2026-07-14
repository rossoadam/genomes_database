#!/usr/bin/env python3

"""
04_build_compleasm_feature_tsvs.py

Build Compleasm-derived feature TSVs for the gc3_dynamics_v5 schema.

This script uses:
    - genomes/records/compleasm/records/metadata.csv
    - genomes/records/genomes_metadata.csv
    - a project manifest with an accession column

The script obtains:
    - cds_fasta and full_table paths from Compleasm metadata
    - genome FASTA (.fna) paths from genomes_metadata.csv

For each genome in the project manifest with Compleasm output, the script:
    1. Uses cds_fasta to find the parent Compleasm output directory.
    2. Finds *ortholog_validity.tsv in that parent directory.
    3. Builds a set of valid ortholog IDs using passes_raw_cds_qc == True.
    4. Opens genus_species_cds_compleasm.fasta and calculates GC, GC3, GC4
       for each valid ortholog.
    5. Opens full_table.tsv and extracts introns for each valid ortholog.
    6. Extracts three flank sets for each valid ortholog:
          50 bp upstream + 50 bp downstream
          100 bp upstream + 100 bp downstream
          200 bp upstream + 200 bp downstream
    7. Writes TSVs matching the intended Compleasm feature tables:
          sauropsida_odb12.tsv
          orthologs.tsv
          ortholog_summary.tsv
          intron_compleasm.tsv
          intron_compleasm_summary.tsv
          flanks_compleasm.tsv
          flank_sets_compleasm.tsv
          flank_compleasm_summary.tsv

Example usage:
    python3 04_build_compleasm_feature_tsvs.py \
        /Users/rossoaa/projects/genomes \
        /Users/rossoaa/projects/genomes/records/project_manifests/mass_predicts_dna_dynamics_with_s_punctatus_manifest.csv

Test on two genomes:
    python3 04_build_compleasm_feature_tsvs.py \
        /Users/rossoaa/projects/genomes \
        /Users/rossoaa/projects/genomes/records/project_manifests/mass_predicts_dna_dynamics_with_s_punctatus_manifest.csv \
        --test
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
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
    (1, "flank_50_up_50_down", 50, 50),
    (2, "flank_100_up_100_down", 100, 100),
    (3, "flank_200_up_200_down", 200, 200),
]


def root_acc(accession: str) -> str:
    accession = str(accession).strip()
    return accession.split(".")[0] if accession else accession


def accession_version(accession: str) -> int:
    accession = str(accession).strip()
    if "." not in accession:
        return -1
    suffix = accession.rsplit(".", 1)[-1]
    return int(suffix) if suffix.isdigit() else -1


def normalize_species_name(value: str) -> str:
    value = str(value).strip()
    if not value or value.lower() == "nan":
        return "unknown_species"
    if "_" in value:
        parts = [p for p in value.split("_") if p]
    else:
        parts = [p for p in value.split() if p]
    if len(parts) >= 2:
        return f"{parts[0].lower()}_{parts[1].lower()}"
    return value.replace(" ", "_").lower()


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
        return genomes_dir / suffix
    return path_obj


def load_genome(input_fasta: Path) -> pyfaidx.Fasta:
    # 1. Index the FASTA file with pyfaidx for fast access
    try:
        genome = pyfaidx.Fasta(
            str(input_fasta),
            as_raw=True,        # slices return plain strings
            build_index=True     # makes .fai if missing
        )
        print(f"Status: Genome FASTA indexed successfully: {input_fasta}")
        return genome
    except pyfaidx.FaidxException as e:
        print(
            f"CRITICAL ERROR: Failed to load genome FASTA with pyfaidx. "
            f"Ensure the .fai file exists for {input_fasta}. Error: {e}"
        )
        sys.exit(1)
    except Exception as e:
        print(f"CRITICAL ERROR: Failed to load genome FASTA. Error: {e}")
        sys.exit(1)


def reverse_complement(seq):
    complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
        'D': 'H', 'H': 'D', 'N': 'N',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
        'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
        'k': 'm', 'm': 'k', 'b': 'v', 'v': 'b',
        'd': 'h', 'h': 'd', 'n': 'n'
    }
    return ''.join(complement.get(base, base) for base in reversed(seq))


def calculate_gc(sequence: str) -> Tuple[float, int, int, int, int, int, int]:
    """
    Calculate GC content and nucleotide counts for a DNA sequence.

    GC content is calculated as:
        (G + C) / (A + T + G + C)

    N bases are counted separately but excluded from the denominator.
    """
    sequence = str(sequence).upper()

    g_count = sequence.count("G")
    c_count = sequence.count("C")
    a_count = sequence.count("A")
    t_count = sequence.count("T")
    n_count = sequence.count("N")

    total_valid_bases = a_count + t_count + g_count + c_count

    gc_content_float = 0.0
    if total_valid_bases > 0:
        gc_content_float = (g_count + c_count) / total_valid_bases

    return (
        gc_content_float,
        g_count,
        c_count,
        a_count,
        t_count,
        n_count,
        total_valid_bases,
    )


def check_cds_sequence(sequence: str) -> Tuple[bool, int, int, bool]:
    """
    Check whether a coding DNA sequence passes basic frame checks.

    The sequence passes if:
        - sequence length is divisible by 3
        - sequence has no internal stop codons

    Terminal stop codons are allowed but reported separately.
    """
    sequence = str(sequence).upper().replace("-", "")
    length_mod_3 = len(sequence) % 3

    if length_mod_3 != 0:
        return False, length_mod_3, 0, False

    codons = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]

    terminal_stop = False
    if codons:
        terminal_stop = codons[-1] in STOP_CODONS

    internal_stops = sum(
        1 for codon in codons[:-1]
        if codon in STOP_CODONS
    )

    passes_cds_check = (
        length_mod_3 == 0
        and internal_stops == 0
    )

    return (
        passes_cds_check,
        length_mod_3,
        internal_stops,
        terminal_stop,
    )


def remove_terminal_stop(sequence: str) -> str:
    sequence = str(sequence).upper().replace("-", "")
    if len(sequence) >= 3 and sequence[-3:] in STOP_CODONS:
        return sequence[:-3]
    return sequence


def mask_for_gc3(sequence: str) -> np.ndarray:
    """
    Return a nucleotide-level boolean mask for GC3 sites.

    True values mark the third codon position of each codon.
    The input sequence must be in frame and have a length divisible by 3.
    """
    sequence = str(sequence).upper().replace("-", "")
    sequence = remove_terminal_stop(sequence)

    (
        passes_cds_check,
        length_mod_3,
        internal_stops,
        terminal_stop,
    ) = check_cds_sequence(sequence)

    if length_mod_3 != 0:
        raise ValueError(
            f"Sequence length ({len(sequence)}) is not divisible by 3"
        )

    if internal_stops > 0:
        raise ValueError(
            f"Sequence contains {internal_stops} internal stop codon(s)"
        )

    integer_sequence = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)

    nuc3 = integer_sequence[2::3]

    A, C, G, T = map(ord, "ACGT")

    valid_third_position = (
        (nuc3 == A)
        | (nuc3 == C)
        | (nuc3 == G)
        | (nuc3 == T)
    )

    # all third-codon positions regardless of third base validity
    candidate_gc3_mask = np.zeros(len(integer_sequence), dtype=bool)
    candidate_gc3_mask[2::3] = True

    # only third-codon positions with valid A/T/G/C bases
    valid_gc3_mask = np.zeros(len(integer_sequence), dtype=bool)
    valid_gc3_mask[2::3] = valid_third_position

    return candidate_gc3_mask, valid_gc3_mask


def calculate_gc3(sequence: str) -> Tuple[float, int, int, int, int, int, int, int, int, int, int, bool]:
    """
    Calculate GC content at GC3 sites only.
    """
    original_sequence = str(sequence).upper().replace("-", "")

    (
        passes_cds_check,
        length_mod_3,
        internal_stops,
        terminal_stop,
    ) = check_cds_sequence(original_sequence)

    if not passes_cds_check:
        if length_mod_3 != 0:
            raise ValueError(
                f"Sequence length ({len(original_sequence)}) is not divisible by 3"
            )
        if internal_stops > 0:
            raise ValueError(
                f"Sequence contains {internal_stops} internal stop codon(s)"
            )

    sequence = remove_terminal_stop(original_sequence)
    candidate_gc3_mask, valid_gc3_mask = mask_for_gc3(sequence)

    sequence_array = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)

    gc3_sequence = (
        sequence_array[valid_gc3_mask]
        .tobytes()
        .decode("ascii")
    )

    (
        gc3_content_float,
        g_count,
        c_count,
        a_count,
        t_count,
        n_count,
        total_valid_bases,
    ) = calculate_gc(gc3_sequence)

    candidate_gc3_sites = int(candidate_gc3_mask.sum())
    valid_gc3_sites = int(valid_gc3_mask.sum())
    invalid_gc3_sites = (
        candidate_gc3_sites - valid_gc3_sites
    )

    return (
        gc3_content_float,
        g_count,
        c_count,
        a_count,
        t_count,
        n_count,
        valid_gc3_sites,
        invalid_gc3_sites,
        candidate_gc3_sites,
        length_mod_3,
        internal_stops,
        terminal_stop,
    )


def mask_for_gc4(sequence: str) -> np.ndarray:
    """
    Return a nucleotide-level boolean mask for GC4 sites.

    True values mark the third codon position of fourfold-degenerate codons.
    The input sequence must be in frame and have a length divisible by 3.
    """
    sequence = str(sequence).upper().replace("-", "")
    sequence = remove_terminal_stop(sequence)

    if len(sequence) % 3 != 0:
        raise ValueError(
            f"Sequence length ({len(sequence)}) is not divisible by 3"
        )

    integer_sequence = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)

    nuc1 = integer_sequence[0::3]
    nuc2 = integer_sequence[1::3]
    nuc3 = integer_sequence[2::3]

    A, C, G, T = map(ord, "ACGT")

    valid_third_position = (
        (nuc3 == A)
        | (nuc3 == C)
        | (nuc3 == G)
        | (nuc3 == T)
    )

    is_candidate_codon = (
        ((nuc1 == G) & (nuc2 == T))  # Val: GTN
        | ((nuc1 == C) & (nuc2 == C))  # Pro: CCN
        | ((nuc1 == A) & (nuc2 == C))  # Thr: ACN
        | ((nuc1 == G) & (nuc2 == C))  # Ala: GCN
        | ((nuc1 == G) & (nuc2 == G))  # Gly: GGN
        | ((nuc1 == C) & (nuc2 == T))  # Leu: CTN
        | ((nuc1 == C) & (nuc2 == G))  # Arg: CGN
        | ((nuc1 == T) & (nuc2 == C))  # Ser: TCN
    )

    # all GC4-family codons regardless of third base validity
    candidate_gc4_mask = np.zeros(len(integer_sequence), dtype=bool)
    candidate_gc4_mask[2::3] = is_candidate_codon

    # only GC4 codons with valid A/T/G/C third bases
    valid_gc4_mask = np.zeros(len(integer_sequence), dtype=bool)
    valid_gc4_mask[2::3] = (
        is_candidate_codon & valid_third_position
    )

    return candidate_gc4_mask, valid_gc4_mask


def calculate_gc4(sequence: str) -> Tuple[float, int, int, int, int, int, int, int, int]:
    """
    Calculate GC content at GC4 sites only.
    """
    sequence = str(sequence).upper().replace("-", "")
    sequence = remove_terminal_stop(sequence)
    candidate_gc4_mask, valid_gc4_mask = mask_for_gc4(sequence)

    sequence_array = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)

    gc4_sequence = (
        sequence_array[valid_gc4_mask]
        .tobytes()
        .decode("ascii")
    )
    (
        gc4_content_float,
        g_count,
        c_count,
        a_count,
        t_count,
        n_count,
        total_valid_bases,
    ) = calculate_gc(gc4_sequence)

    candidate_gc4_sites = int(candidate_gc4_mask.sum())
    valid_gc4_sites = int(valid_gc4_mask.sum())
    invalid_gc4_sites = (
        candidate_gc4_sites - valid_gc4_sites
    )

    return (
        gc4_content_float,
        g_count,
        c_count,
        a_count,
        t_count,
        n_count,
        valid_gc4_sites,
        invalid_gc4_sites,
        candidate_gc4_sites,
    )


def safe_float(value) -> Optional[float]:
    if value is None:
        return None
    try:
        if math.isnan(value):
            return None
    except TypeError:
        pass
    return float(value)


def parse_bool(value) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def parse_coordinate_token(coordinates: str, busco_id: str) -> Dict[str, object]:
    parts = coordinates.split("_")
    if len(parts) != 3:
        raise ValueError(f"WARNING: {busco_id}: unexpected exon token '{coordinates}'")
    return {
        "busco_id": busco_id,
        "start": int(parts[0]),
        "end": int(parts[1]),
        "strand": parts[2],
    }


def parse_full_table(full_table: Path) -> Dict[str, Dict[str, object]]:
    records: Dict[str, Dict[str, object]] = {}
    with open(full_table, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("Gene"):
                continue
            data = line.split("\t")
            if len(data) < 13:
                continue
            busco_id = data[0]
            status = data[1]
            chromosome = data[2]
            strand = data[5]
            coordinate_tokens = [x for x in data[12].split("|") if x]
            exons = []
            for token in coordinate_tokens:
                exon = parse_coordinate_token(token, busco_id)
                if exon["strand"] != strand:
                    continue
                exons.append(exon)
            exons.sort(key=lambda x: int(x["start"]))
            if exons:
                start = min(int(e["start"]) for e in exons)
                end = max(int(e["end"]) for e in exons)
            else:
                start = None
                end = None
            records[busco_id] = {
                "odb12_id": busco_id,
                "status": status,
                "chromosome": chromosome,
                "strand": strand,
                "start": start,
                "end": end,
                "exons": exons,
            }
    return records


def infer_introns(sorted_exons: List[Dict[str, object]], strand: str) -> List[Dict[str, object]]:
    introns = []
    if len(sorted_exons) <= 1:
        return introns
    genomic_introns = []
    for i in range(len(sorted_exons) - 1):
        left_exon = sorted_exons[i]
        right_exon = sorted_exons[i + 1]
        intron_start = int(left_exon["end"]) + 1
        intron_end = int(right_exon["start"]) - 1
        if intron_start > intron_end:
            continue
        genomic_introns.append({
            "busco_id": sorted_exons[0]["busco_id"],
            "intron_start": intron_start,
            "intron_end": intron_end,
            "intron_length": intron_end - intron_start + 1,
            "strand": strand,
        })
    if strand == "-":
        genomic_introns = list(reversed(genomic_introns))
    for ij, intron in enumerate(genomic_introns, start=1):
        intron["intron_number"] = ij
        intron["intron_id"] = f"i_{intron['busco_id']}_{ij}"
    return genomic_introns


def slice_genome(genome: pyfaidx.Fasta, chromosome: str, start_1: int, end_1: int) -> str:
    if start_1 < 1:
        start_1 = 1
    chrom_len = len(genome[chromosome])
    if end_1 > chrom_len:
        end_1 = chrom_len
    if start_1 > end_1:
        return ""
    return genome[chromosome][start_1 - 1:end_1]


def extract_introns(genome: pyfaidx.Fasta, record: Dict[str, object]) -> List[Dict[str, object]]:
    chromosome = str(record["chromosome"])
    strand = str(record["strand"])
    exons = record["exons"]
    introns = infer_introns(exons, strand)
    rows = []
    for intron in introns:
        seq = slice_genome(genome, chromosome, int(intron["intron_start"]), int(intron["intron_end"]))
        if strand == "-":
            seq = reverse_complement(seq)
        gc, g, c, a, t, n, valid = calculate_gc(seq)
        rows.append({
            **intron,
            "sequence": seq,
            "gc": gc,
            "valid_bases": valid,
        })
    return rows


def extract_flanks(genome: pyfaidx.Fasta, record: Dict[str, object], up_bp: int, down_bp: int) -> Dict[str, object]:
    chromosome = str(record["chromosome"])
    strand = str(record["strand"])
    gene_start = int(record["start"])
    gene_end = int(record["end"])

    if strand == "+":
        up_start = max(1, gene_start - up_bp)
        up_end = gene_start - 1
        down_start = gene_end + 1
        down_end = gene_end + down_bp
        up_seq = slice_genome(genome, chromosome, up_start, up_end)
        down_seq = slice_genome(genome, chromosome, down_start, down_end)
    else:
        # Biological upstream for a minus-strand gene lies after the gene in genomic coordinates.
        up_start = gene_end + 1
        up_end = gene_end + up_bp
        down_start = max(1, gene_start - down_bp)
        down_end = gene_start - 1
        up_seq = reverse_complement(slice_genome(genome, chromosome, up_start, up_end))
        down_seq = reverse_complement(slice_genome(genome, chromosome, down_start, down_end))

    combined = up_seq + down_seq
    gc, g, c, a, t, n, valid = calculate_gc(combined)
    up_gc = calculate_gc(up_seq)[0]
    down_gc = calculate_gc(down_seq)[0]

    return {
        "up_start": up_start,
        "up_end": up_end,
        "down_start": down_start,
        "down_end": down_end,
        "up_sequence": up_seq,
        "down_sequence": down_seq,
        "combined_sequence": combined,
        "up_length": len(up_seq),
        "down_length": len(down_seq),
        "combined_length": len(combined),
        "gc": gc,
        "upstream_gc": up_gc,
        "downstream_gc": down_gc,
        "valid_bases": valid,
    }


def find_genome_fasta(row: pd.Series, genomes_dir: Path) -> Path:
    candidate_cols = [
        "genome_fasta", "path_to_fna", "fasta", "input_fasta", "path_to_genome",
        "genomic_fna", "genome_path", "fna_path"
    ]
    for col in candidate_cols:
        if col in row.index:
            path = repath(row[col], genomes_dir)
            if str(path) and path.exists():
                return path
    raise FileNotFoundError(
        "Could not find genome FASTA path in metadata row. Tried columns: "
        + ", ".join(candidate_cols)
    )


def find_validity_tsv(cds_fasta: Path) -> Path:
    parent = cds_fasta.parent
    candidates = sorted(parent.glob("*ortholog_validity.tsv"))
    if not candidates:
        candidates = sorted(parent.glob("*raw_ortholog_validity.tsv"))
    if not candidates:
        raise FileNotFoundError(f"No *ortholog_validity.tsv found in {parent}")
    if len(candidates) > 1:
        # Prefer the file that is not a summary.
        nonsummary = [p for p in candidates if "summary" not in p.name]
        if nonsummary:
            return nonsummary[0]
    return candidates[0]


def load_valid_orthologs(validity_tsv: Path) -> set[str]:
    df = pd.read_csv(validity_tsv, sep="\t")
    if "sequence_id" not in df.columns or "passes_raw_cds_qc" not in df.columns:
        raise ValueError(f"{validity_tsv} must contain sequence_id and passes_raw_cds_qc columns")
    return set(df.loc[df["passes_raw_cds_qc"].map(parse_bool), "sequence_id"].astype(str))


def load_cds_sequences(cds_fasta: Path) -> Dict[str, str]:
    seqs = {}
    with open(cds_fasta, "r") as handle:
        for header, seq in SimpleFastaParser(handle):
            seq_id = header.split()[0]
            seqs[seq_id] = seq.upper().replace("-", "")
    return seqs


def read_manifest_accessions(manifest_path: Path) -> List[str]:
    df = pd.read_csv(manifest_path)
    if "accession" not in df.columns:
        raise ValueError(f"Manifest is missing required accession column: {manifest_path}")
    accessions = df["accession"].dropna().astype(str).str.strip()
    return sorted(set(accessions[accessions != ""]))


def resolve_metadata_rows(metadata_path: Path, manifest_accessions: List[str]) -> Tuple[pd.DataFrame, List[str]]:
    metadata_df = pd.read_csv(metadata_path)
    required = {"accession", "cds_fasta", "full_table"}
    missing = required.difference(metadata_df.columns)
    if missing:
        raise ValueError("Compleasm metadata missing required columns: " + ", ".join(sorted(missing)))

    allowed_accessions = set(str(a).strip() for a in manifest_accessions)
    allowed_roots = {root_acc(a) for a in allowed_accessions}
    metadata_df = metadata_df.copy()
    metadata_df["accession"] = metadata_df["accession"].astype(str).str.strip()
    metadata_df["accession_root"] = metadata_df["accession"].map(root_acc)
    metadata_df["accession_version"] = metadata_df["accession"].map(accession_version)

    exact = metadata_df[metadata_df["accession"].isin(allowed_accessions)].copy()
    matched_roots = set(exact["accession_root"])
    unresolved_roots = allowed_roots - matched_roots
    fallback = metadata_df[
        (~metadata_df["accession"].isin(allowed_accessions))
        & (metadata_df["accession_root"].isin(unresolved_roots))
    ].copy()
    if not fallback.empty:
        fallback = (
            fallback.sort_values(["accession_root", "accession_version", "accession"], ascending=[True, False, True])
            .drop_duplicates(subset=["accession_root"], keep="first")
        )
    resolved = pd.concat([exact, fallback], ignore_index=True)
    if resolved.empty:
        return resolved, sorted(allowed_roots)
    resolved = (
        resolved.sort_values(["accession_root", "accession_version", "accession"], ascending=[True, False, True])
        .drop_duplicates(subset=["accession_root"], keep="first")
        .reset_index(drop=True)
    )
    missing_roots = sorted(allowed_roots - set(resolved["accession_root"]))
    return resolved, missing_roots


def std(values: Sequence[float]) -> Optional[float]:
    values = [float(v) for v in values if v is not None and not math.isnan(float(v))]
    if len(values) < 2:
        return None
    return float(np.std(values, ddof=1))


def var(values: Sequence[float]) -> Optional[float]:
    values = [float(v) for v in values if v is not None and not math.isnan(float(v))]
    if len(values) < 2:
        return None
    return float(np.var(values, ddof=1))


def mad(values: Sequence[float]) -> Optional[float]:
    values = [float(v) for v in values if v is not None and not math.isnan(float(v))]
    if not values:
        return None
    med = median(values)
    return float(median([abs(v - med) for v in values]))


def q(values: Sequence[float], quantile: float) -> Optional[float]:
    values = [float(v) for v in values if v is not None and not math.isnan(float(v))]
    if not values:
        return None
    return float(np.quantile(values, quantile))


def weighted_mean(values: Sequence[float], weights: Sequence[int]) -> Optional[float]:
    pairs = [(float(v), int(w)) for v, w in zip(values, weights) if v is not None and w is not None and int(w) > 0]
    if not pairs:
        return None
    denom = sum(w for _, w in pairs)
    if denom == 0:
        return None
    return sum(v * w for v, w in pairs) / denom


def summary_stats(values: Sequence[float], lengths: Sequence[int], prefix: str) -> Dict[str, object]:
    vals = [float(v) for v in values if v is not None and not math.isnan(float(v))]
    return {
        f"mean_{prefix}": mean(vals) if vals else None,
        f"weighted_mean_{prefix}": weighted_mean(values, lengths),
        f"sd_{prefix}": std(values),
        f"var_{prefix}": var(values),
        f"median_{prefix}": median(vals) if vals else None,
        f"mad_{prefix}": mad(values),
        f"iqr_{prefix}": (q(values, 0.75) - q(values, 0.25)) if vals else None,
    }


def write_tsv(rows: List[Dict[str, object]], path: Path, columns: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            clean = {col: ("" if row.get(col) is None else row.get(col)) for col in columns}
            writer.writerow(clean)


def now_str() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build Compleasm ortholog, intron, and flank TSVs for gc3_dynamics_v5."
    )
    parser.add_argument("genomes_dir", help="Path to genomes directory")
    parser.add_argument("project_manifest", help="Project manifest CSV with accession column")
    parser.add_argument("--metadata", help="Optional Compleasm metadata.csv path")
    parser.add_argument("--outdir", help="Optional output directory")
    parser.add_argument("--test", action="store_true", help="Process only the first two resolved genomes")
    parser.add_argument("--test_n", type=int, default=2, help="Number of genomes to process in --test mode")
    args = parser.parse_args()

    genomes_dir = Path(args.genomes_dir).resolve()
    manifest_path = Path(args.project_manifest).resolve()
    metadata_path = Path(args.metadata).resolve() if args.metadata else genomes_dir / "records/compleasm/records/metadata.csv"
    outdir = Path(args.outdir).resolve() if args.outdir else genomes_dir / "records/sql_tsvs"
    outdir.mkdir(parents=True, exist_ok=True)

    manifest_accessions = read_manifest_accessions(manifest_path)
    
    metadata_rows, missing_roots = resolve_metadata_rows(
    	metadata_path,
    	manifest_accessions
    )
    
    genomes_metadata_path = genomes_dir / "records/genomes_metadata.csv"
    genomes_metadata_df = pd.read_csv(genomes_metadata_path)
    
    genomes_metadata_df["accession"] = genomes_metadata_df["accession"].astype(str).str.strip()
    genomes_metadata_df["accession_root"] = genomes_metadata_df["accession"].map(root_acc)
    
    metadata_rows = metadata_rows.merge(
    	genomes_metadata_df,
    	on="accession_root",
    	how="left",
    	suffixes=("", "_genome")
    )
    if args.test:
        metadata_rows = metadata_rows.head(args.test_n).copy()

    print(f"Loaded {len(manifest_accessions)} accessions from manifest")
    print(f"Resolved {len(metadata_rows)} metadata rows")
    print(f"Missing accession roots in metadata: {len(missing_roots)}")
    if args.test:
        print(f"TEST MODE: processing {len(metadata_rows)} genome(s)")

    sauropsida_ids = set()
    ortholog_rows = []
    ortholog_summary_rows = []
    intron_rows = []
    intron_summary_rows = []
    flank_rows = []
    flank_summary_rows = []

    flank_set_rows = [
        {
            "flank_set_pk": pk,
            "name": name,
            "upstream_bp": up,
            "downstream_bp": down,
            "created_at": now_str(),
            "notes": "Compleasm flank set generated from valid raw CDS orthologs",
        }
        for pk, name, up, down in FLANK_SPECS
    ]

    ortholog_pk = 0
    intron_pk = 0
    flank_pk = 0
    ortholog_summary_pk = 0
    intron_summary_pk = 0
    flank_summary_pk = 0

    for genome_index, (_, row) in enumerate(metadata_rows.iterrows(), start=1):
        accession = str(row["accession"]).strip()
        accession_root = str(row["accession_root"]).strip()
        genome_pk = genome_index
        species_pk = genome_index
        species = normalize_species_name(row.get("species", row.get("organism_name", accession)))

        cds_fasta = repath(row["cds_fasta"], genomes_dir)
        full_table = repath(row["full_table"], genomes_dir)
        if not cds_fasta.exists():
            print(f"WARNING: missing cds_fasta for {accession}: {cds_fasta}")
            continue
        if not full_table.exists():
            print(f"WARNING: missing full_table for {accession}: {full_table}")
            continue

        try:
            genome_fasta = repath(row["path_to_fna"], genomes_dir)
            if not genome_fasta.exists():
                print(f"WARNING: missing genome FASTA for {accession}: {genome_fasta}")
                continue
        except FileNotFoundError as e:
            print(f"WARNING: {accession}: {e}")
            continue

        try:
            validity_tsv = find_validity_tsv(cds_fasta)
            valid_ids = load_valid_orthologs(validity_tsv)
        except Exception as e:
            print(f"WARNING: {accession}: could not load ortholog validity TSV: {e}")
            continue

        cds_sequences = load_cds_sequences(cds_fasta)
        full_records = parse_full_table(full_table)
        genome = load_genome(genome_fasta)

        genome_ortholog_pks = []
        genome_gc = []
        genome_gc3 = []
        genome_gc4 = []
        genome_ortholog_lengths = []

        genome_introns = []
        genome_flanks_by_set: Dict[int, List[Dict[str, object]]] = {pk: [] for pk, _, _, _ in FLANK_SPECS}

        for odb12_id in sorted(valid_ids):
            if odb12_id not in cds_sequences:
                continue
            if odb12_id not in full_records:
                continue

            seq = cds_sequences[odb12_id]
            record = full_records[odb12_id]
            if record["status"] != "Single" or not record["exons"] or record["start"] is None or record["end"] is None:
                continue

            sauropsida_ids.add(odb12_id)
            gc = calculate_gc(seq)[0]
            try:
                gc3 = calculate_gc3(seq)[0]
            except ValueError:
                gc3 = None
            try:
                gc4 = calculate_gc4(seq)[0]
            except ValueError:
                gc4 = None

            ortholog_pk += 1
            this_ortholog_pk = ortholog_pk
            genome_ortholog_pks.append(this_ortholog_pk)
            genome_gc.append(gc)
            genome_gc3.append(gc3)
            genome_gc4.append(gc4)
            genome_ortholog_lengths.append(len(seq))

            ortholog_rows.append({
                "ortholog_pk": this_ortholog_pk,
                "sequence_pk": genome_pk,
                "odb12_id": odb12_id,
                "status": record["status"],
                "strand": record["strand"],
                "start": record["start"],
                "end": record["end"],
                "gc": gc,
                "gc3": gc3,
                "gc4": gc4,
            })

            for intron in extract_introns(genome, record):
                intron_pk += 1
                intron_row = {
                    "intron_pk": intron_pk,
                    "ortholog_pk": this_ortholog_pk,
                    "parent_id": odb12_id,
                    "intron_id": intron["intron_id"],
                    "status": record["status"],
                    "strand": record["strand"],
                    "start": intron["intron_start"],
                    "end": intron["intron_end"],
                    "length": intron["intron_length"],
                    "gc": intron["gc"],
                }
                intron_rows.append(intron_row)
                genome_introns.append(intron_row)

            for flank_set_pk, flank_name, up_bp, down_bp in FLANK_SPECS:
                flank = extract_flanks(genome, record, up_bp, down_bp)
                flank_pk += 1
                flank_row = {
                    "flank_pk": flank_pk,
                    "ortholog_pk": this_ortholog_pk,
                    "flank_set_pk": flank_set_pk,
                    "status": record["status"],
                    "strand": record["strand"],
                    "up_start": flank["up_start"],
                    "up_end": flank["up_end"],
                    "down_start": flank["down_start"],
                    "down_end": flank["down_end"],
                    "gc": flank["gc"],
                    # extra values used internally for summary only, not written to flanks_compleasm.tsv
                    "length": flank["combined_length"],
                    "upstream_gc": flank["upstream_gc"],
                    "downstream_gc": flank["downstream_gc"],
                }
                flank_rows.append(flank_row)
                genome_flanks_by_set[flank_set_pk].append(flank_row)

        # Ortholog summary
        if genome_ortholog_pks:
            ortholog_summary_pk += 1
            gc_stats = summary_stats(genome_gc, genome_ortholog_lengths, "gc")
            gc3_stats = summary_stats(genome_gc3, genome_ortholog_lengths, "gc3")
            gc4_stats = summary_stats(genome_gc4, genome_ortholog_lengths, "gc4")
            ortholog_summary_rows.append({
                "ortholog_summary_pk": ortholog_summary_pk,
                "genome_pk": genome_pk,
                "species_pk": species_pk,
                "n_orthologs": len(genome_ortholog_pks),
                "callable_bp_total": sum(genome_ortholog_lengths),
                **gc_stats,
                **gc3_stats,
                **gc4_stats,
                "q05_gc3": q(genome_gc3, 0.05),
                "q25_gc3": q(genome_gc3, 0.25),
                "q75_gc3": q(genome_gc3, 0.75),
                "q95_gc3": q(genome_gc3, 0.95),
                "mean_ortholog_length": mean(genome_ortholog_lengths),
                "median_ortholog_length": median(genome_ortholog_lengths),
                "created_at": now_str(),
            })

        if genome_introns:
            intron_summary_pk += 1
            intron_gc = [r["gc"] for r in genome_introns]
            intron_lengths = [r["length"] for r in genome_introns]
            intron_stats = summary_stats(intron_gc, intron_lengths, "gc")
            intron_summary_rows.append({
                "intron_summary_pk": intron_summary_pk,
                "genome_pk": genome_pk,
                "species_pk": species_pk,
                "n_introns": len(genome_introns),
                "n_orthologs": len(set(r["ortholog_pk"] for r in genome_introns)),
                "callable_bp_total": sum(intron_lengths),
                **intron_stats,
                "q05_gc": q(intron_gc, 0.05),
                "q25_gc": q(intron_gc, 0.25),
                "q75_gc": q(intron_gc, 0.75),
                "q95_gc": q(intron_gc, 0.95),
                "mean_intron_length": mean(intron_lengths),
                "median_intron_length": median(intron_lengths),
                "created_at": now_str(),
            })

        for flank_set_pk, rows_for_set in genome_flanks_by_set.items():
            if not rows_for_set:
                continue
            flank_summary_pk += 1
            flank_gc = [r["gc"] for r in rows_for_set]
            flank_lengths = [r["length"] for r in rows_for_set]
            flank_stats = summary_stats(flank_gc, flank_lengths, "gc")
            flank_summary_rows.append({
                "flank_summary_pk": flank_summary_pk,
                "genome_pk": genome_pk,
                "species_pk": species_pk,
                "flank_set_pk": flank_set_pk,
                "n_flanks": len(rows_for_set),
                "n_orthologs": len(set(r["ortholog_pk"] for r in rows_for_set)),
                "callable_bp_total": sum(flank_lengths),
                **flank_stats,
                "q05_gc": q(flank_gc, 0.05),
                "q25_gc": q(flank_gc, 0.25),
                "q75_gc": q(flank_gc, 0.75),
                "q95_gc": q(flank_gc, 0.95),
                "mean_flank_length": mean(flank_lengths),
                "median_flank_length": median(flank_lengths),
                "mean_upstream_gc": mean([r["upstream_gc"] for r in rows_for_set]),
                "mean_downstream_gc": mean([r["downstream_gc"] for r in rows_for_set]),
                "created_at": now_str(),
            })

        print(f"Processed {accession}: {len(genome_ortholog_pks)} valid orthologs")

    sauropsida_rows = [{"odb12_id": x, "name": "", "orthodb": ""} for x in sorted(sauropsida_ids)]

    write_tsv(sauropsida_rows, outdir / "sauropsida_odb12.tsv", ["odb12_id", "name", "orthodb"])
    write_tsv(ortholog_rows, outdir / "orthologs.tsv", ["ortholog_pk", "sequence_pk", "odb12_id", "status", "strand", "start", "end", "gc", "gc3", "gc4"])
    write_tsv(ortholog_summary_rows, outdir / "ortholog_summary.tsv", [
        "ortholog_summary_pk", "genome_pk", "species_pk", "n_orthologs", "callable_bp_total",
        "mean_gc", "weighted_mean_gc", "mean_gc3", "weighted_mean_gc3", "mean_gc4", "weighted_mean_gc4",
        "sd_gc", "sd_gc3", "sd_gc4", "var_gc", "var_gc3", "var_gc4", "median_gc", "median_gc3",
        "median_gc4", "mad_gc", "mad_gc3", "mad_gc4", "iqr_gc", "iqr_gc3", "iqr_gc4",
        "q05_gc3", "q25_gc3", "q75_gc3", "q95_gc3", "mean_ortholog_length", "median_ortholog_length", "created_at"
    ])
    write_tsv(intron_rows, outdir / "intron_compleasm.tsv", ["intron_pk", "ortholog_pk", "parent_id", "intron_id", "status", "strand", "start", "end", "length", "gc"])
    write_tsv(intron_summary_rows, outdir / "intron_compleasm_summary.tsv", [
        "intron_summary_pk", "genome_pk", "species_pk", "n_introns", "n_orthologs", "callable_bp_total",
        "mean_gc", "weighted_mean_gc", "sd_gc", "var_gc", "median_gc", "mad_gc", "iqr_gc",
        "q05_gc", "q25_gc", "q75_gc", "q95_gc", "mean_intron_length", "median_intron_length", "created_at"
    ])
    write_tsv(flank_rows, outdir / "flanks_compleasm.tsv", ["flank_pk", "ortholog_pk", "flank_set_pk", "status", "strand", "up_start", "up_end", "down_start", "down_end", "gc"])
    write_tsv(flank_set_rows, outdir / "flank_sets_compleasm.tsv", ["flank_set_pk", "name", "upstream_bp", "downstream_bp", "created_at", "notes"])
    write_tsv(flank_summary_rows, outdir / "flank_compleasm_summary.tsv", [
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
