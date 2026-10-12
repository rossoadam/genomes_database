#!/usr/bin/env python3

"""
03_raw_ortholog_validity.py

Evaluate raw ortholog CDS sequences created by 01a_get_cds_from_compleasm_v6.py.

Each input sequence is expected to be one ortholog sequence composed of multiple
CDS fragments extracted from Compleasm output.

This script checks each sequence for:
    - length divisibility by 3
    - internal stop codons
    - terminal stop codon
    - ambiguous bases
    - valid nucleotide content
    - translated protein length using translate_dna()

This batch version takes:
    - genomes_dir
    - project_manifest

It finds:
    genomes/records/compleasm/records/metadata.csv

and uses the cds_fasta column in that metadata file to run raw ortholog
validity QC for all genomes in the project manifest that have Compleasm CDS
FASTA output.

Outputs are written to:
    genomes/records/compleasm/records/raw_ortholog_validity/

Example usage:
    python3 03_raw_ortholog_validity.py \
        /Users/rossoaa/projects/genomes \
        /Users/rossoaa/projects/genomes/records/project_manifests/mass_predicts_dna_dynamics_with_s_punctatus_manifest.csv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
from Bio import SeqIO


STOP_CODONS = {"TAA", "TAG", "TGA"}
VALID_BASES = {"A", "T", "G", "C"}
AMBIGUOUS_BASES = {"N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}


def root_acc(accession: str) -> str:
    """
    Return accession root without version suffix.
    """
    accession = str(accession).strip()
    return accession.split(".")[0] if accession else accession


def accession_version(accession: str) -> int:
    """
    Return accession version as an integer.

    Accessions without a numeric version return -1.
    """
    accession = str(accession).strip()

    if "." not in accession:
        return -1

    suffix = accession.rsplit(".", 1)[-1]

    return int(suffix) if suffix.isdigit() else -1


def translate_dna(dna: str) -> str:
    """
    Translate DNA into amino acids using the standard genetic code.

    Unknown or ambiguous codons are translated as X.
    Stop codons are translated as *.
    """
    table = {
        "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
        "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
        "TCT":"S","TCC":"S","TCA":"S","TCG":"S","AGT":"S","AGC":"S",
        "CCT":"P","CCC":"P","CCA":"P","CCG":"P","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
        "GCT":"A","GCC":"A","GCA":"A","GCG":"A","TAT":"Y","TAC":"Y",
        "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","AAT":"N","AAC":"N","AAA":"K","AAG":"K",
        "GAT":"D","GAC":"D","GAA":"E","GAG":"E","TGT":"C","TGC":"C","TGG":"W",
        "CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
        "TAA":"*","TAG":"*","TGA":"*"
    }
    up = dna.upper()
    stop = len(up) - (len(up) % 3)
    return "".join(table.get(up[i:i+3], "X") for i in range(0, stop, 3))


def repath(path_to_fix, genomes_dir: Path) -> Path:
    """
    Re-map stored paths to the current genomes_dir when needed.

    This helps when metadata.csv contains paths from another machine or an older
    location but still includes the genomes directory name.
    """
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


def read_manifest_accessions(manifest_path: Path) -> List[str]:
    """
    Read accession values from the project manifest.
    """
    df_manifest = pd.read_csv(manifest_path)

    if "accession" not in df_manifest.columns:
        raise ValueError(
            f"Manifest is missing required 'accession' column: {manifest_path}"
        )

    accessions = (
        df_manifest["accession"]
        .dropna()
        .astype(str)
        .str.strip()
    )

    accessions = accessions[accessions != ""]

    return sorted(set(accessions))


def resolve_metadata_rows(
    metadata_path: Path,
    manifest_accessions: List[str],
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Resolve metadata rows for accessions in the project manifest.

    Exact accession matches are preferred. If exact versions are not present,
    the highest available version for the same accession root is used.
    """
    metadata_df = pd.read_csv(metadata_path)

    required_columns = {"accession", "cds_fasta"}

    missing_columns = required_columns.difference(metadata_df.columns)

    if missing_columns:
        raise ValueError(
            f"Compleasm metadata is missing required columns: "
            f"{', '.join(sorted(missing_columns))}"
        )

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
            fallback
            .sort_values(
                ["accession_root", "accession_version", "accession"],
                ascending=[True, False, True],
            )
            .drop_duplicates(subset=["accession_root"], keep="first")
        )

    resolved = pd.concat([exact, fallback], ignore_index=True)

    if resolved.empty:
        return resolved, sorted(allowed_roots)

    resolved = (
        resolved
        .sort_values(
            ["accession_root", "accession_version", "accession"],
            ascending=[True, False, True],
        )
        .drop_duplicates(subset=["accession_root"], keep="first")
        .reset_index(drop=True)
    )

    missing_roots = sorted(allowed_roots - set(resolved["accession_root"]))

    return resolved, missing_roots


def count_bases(sequence: str) -> Dict[str, int]:
    """
    Count canonical, ambiguous, and invalid characters in a DNA sequence.
    """
    sequence = str(sequence).upper().replace("-", "")

    a_count = sequence.count("A")
    t_count = sequence.count("T")
    g_count = sequence.count("G")
    c_count = sequence.count("C")

    n_count = sequence.count("N")

    ambiguous_count = sum(
        1 for base in sequence
        if base in AMBIGUOUS_BASES
    )

    invalid_count = sum(
        1 for base in sequence
        if base not in VALID_BASES and base not in AMBIGUOUS_BASES
    )

    total_valid_bases = a_count + t_count + g_count + c_count

    return {
        "a_count": a_count,
        "t_count": t_count,
        "g_count": g_count,
        "c_count": c_count,
        "n_count": n_count,
        "ambiguous_count": ambiguous_count,
        "invalid_count": invalid_count,
        "total_valid_bases": total_valid_bases,
    }


def check_stop_codons(sequence: str) -> Tuple[int, bool]:
    """
    Count internal stop codons and determine whether a terminal stop is present.

    The sequence is expected to already be in frame. If length is not divisible
    by 3, stop-codon results should be interpreted cautiously.
    """
    sequence = str(sequence).upper().replace("-", "")

    codons = [
        sequence[i:i + 3]
        for i in range(0, len(sequence) - 2, 3)
    ]

    if not codons:
        return 0, False

    internal_stops = sum(
        1 for codon in codons[:-1]
        if codon in STOP_CODONS
    )

    terminal_stop = codons[-1] in STOP_CODONS

    return internal_stops, terminal_stop


def translate_sequence(sequence: str, terminal_stop: bool) -> str:
    """
    Translate a DNA sequence.

    If a terminal stop codon is present, remove it before translation length
    summaries so protein_length reflects amino acids rather than the stop.
    """
    sequence = str(sequence).upper().replace("-", "")

    if terminal_stop:
        sequence = sequence[:-3]

    trimmed_length = len(sequence) - (len(sequence) % 3)
    sequence = sequence[:trimmed_length]

    if not sequence:
        return ""

    return translate_dna(sequence)


def evaluate_sequence(record_id: str, sequence: str) -> Dict[str, object]:
    """
    Evaluate one raw ortholog CDS sequence.
    """
    sequence = str(sequence).upper().replace("-", "")

    seq_length = len(sequence)
    length_mod_3 = seq_length % 3

    base_counts = count_bases(sequence)

    ambiguous_fraction = 0.0
    invalid_fraction = 0.0

    if seq_length > 0:
        ambiguous_fraction = base_counts["ambiguous_count"] / seq_length
        invalid_fraction = base_counts["invalid_count"] / seq_length

    internal_stop_count, terminal_stop = check_stop_codons(sequence)

    protein = translate_sequence(sequence, terminal_stop=terminal_stop)
    protein_length = len(protein)

    passes_length_mod_3 = length_mod_3 == 0
    passes_internal_stop_check = internal_stop_count == 0
    passes_invalid_base_check = base_counts["invalid_count"] == 0

    passes_raw_cds_qc = (
        passes_length_mod_3
        and passes_internal_stop_check
        and passes_invalid_base_check
    )

    failure_reasons: List[str] = []

    if not passes_length_mod_3:
        failure_reasons.append("length_not_divisible_by_3")

    if not passes_internal_stop_check:
        failure_reasons.append("internal_stop_codons")

    if not passes_invalid_base_check:
        failure_reasons.append("invalid_characters")

    if not failure_reasons:
        failure_reasons.append("pass")

    return {
        "sequence_id": record_id,
        "sequence_length": seq_length,
        "length_mod_3": length_mod_3,
        "internal_stop_count": internal_stop_count,
        "terminal_stop": terminal_stop,
        "protein_length": protein_length,
        "a_count": base_counts["a_count"],
        "t_count": base_counts["t_count"],
        "g_count": base_counts["g_count"],
        "c_count": base_counts["c_count"],
        "n_count": base_counts["n_count"],
        "ambiguous_count": base_counts["ambiguous_count"],
        "ambiguous_fraction": ambiguous_fraction,
        "invalid_count": base_counts["invalid_count"],
        "invalid_fraction": invalid_fraction,
        "total_valid_bases": base_counts["total_valid_bases"],
        "passes_length_mod_3": passes_length_mod_3,
        "passes_internal_stop_check": passes_internal_stop_check,
        "passes_invalid_base_check": passes_invalid_base_check,
        "passes_raw_cds_qc": passes_raw_cds_qc,
        "failure_reasons": ";".join(failure_reasons),
    }


def evaluate_fasta(
    fasta_path: Path,
    accession: str,
    accession_root: str,
) -> List[Dict[str, object]]:
    """
    Evaluate all raw ortholog CDS sequences in one FASTA file.
    """
    rows = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        row = evaluate_sequence(
            record_id=record.id,
            sequence=str(record.seq),
        )

        row = {
            "accession": accession,
            "accession_root": accession_root,
            "cds_fasta": str(fasta_path),
            **row,
        }

        rows.append(row)

    return rows


def write_tsv(rows: List[Dict[str, object]], output_path: Path) -> None:
    """
    Write QC results to a TSV file.
    """
    if not rows:
        raise ValueError("No rows to write.")

    columns = list(rows[0].keys())

    with open(output_path, "w") as out_handle:
        out_handle.write("\t".join(columns) + "\n")

        for row in rows:
            values = [
                str(row[column])
                for column in columns
            ]
            out_handle.write("\t".join(values) + "\n")


def summarize(rows: List[Dict[str, object]]) -> Dict[str, int]:
    """
    Return a concise summary dictionary.
    """
    total_sequences = len(rows)

    passed = sum(
        1 for row in rows
        if row["passes_raw_cds_qc"]
    )

    failed = total_sequences - passed

    mod3_failed = sum(
        1 for row in rows
        if not row["passes_length_mod_3"]
    )

    internal_stop_failed = sum(
        1 for row in rows
        if not row["passes_internal_stop_check"]
    )

    invalid_base_failed = sum(
        1 for row in rows
        if not row["passes_invalid_base_check"]
    )

    terminal_stop_count = sum(
        1 for row in rows
        if row["terminal_stop"]
    )

    return {
        "total_sequences": total_sequences,
        "passed_raw_cds_qc": passed,
        "failed_raw_cds_qc": failed,
        "failed_length_mod_3_check": mod3_failed,
        "failed_internal_stop_check": internal_stop_failed,
        "failed_invalid_base_check": invalid_base_failed,
        "sequences_with_terminal_stop_codon": terminal_stop_count,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate raw ortholog CDS sequences produced by "
            "01a_get_cds_from_compleasm_v6.py for all genomes in a project "
            "manifest with Compleasm CDS FASTA output."
        )
    )
    parser.add_argument(
        "genomes_dir",
        help="Path to genomes directory.",
    )
    parser.add_argument(
        "project_manifest",
        help="Project manifest CSV with an accession column.",
    )
    parser.add_argument(
        "--metadata",
        required=False,
        help=(
            "Optional path to Compleasm metadata.csv. "
            "Default: genomes/records/compleasm/records/metadata.csv"
        ),
    )
    parser.add_argument(
        "--outdir",
        required=False,
        help=(
            "Optional output directory. "
            "Default: genomes/records/compleasm/records/raw_ortholog_validity"
        ),
    )

    args = parser.parse_args()

    genomes_dir = Path(args.genomes_dir).resolve()
    manifest_path = Path(args.project_manifest).resolve()

    if not genomes_dir.exists():
        raise FileNotFoundError(f"genomes_dir not found: {genomes_dir}")

    if not manifest_path.exists():
        raise FileNotFoundError(f"Project manifest not found: {manifest_path}")

    if args.metadata:
        metadata_path = Path(args.metadata).resolve()
    else:
        metadata_path = genomes_dir / "records/compleasm/records/metadata.csv"

    if not metadata_path.exists():
        raise FileNotFoundError(f"Compleasm metadata not found: {metadata_path}")

    if args.outdir:
        outdir = Path(args.outdir).resolve()
    else:
        outdir = genomes_dir / "records/compleasm/records/raw_ortholog_validity"

    outdir.mkdir(parents=True, exist_ok=True)

    manifest_accessions = read_manifest_accessions(manifest_path)

    metadata_rows, missing_roots = resolve_metadata_rows(
        metadata_path=metadata_path,
        manifest_accessions=manifest_accessions,
    )

    all_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []

    skipped_missing_cds_fasta = 0
    skipped_empty_fastas = 0

    print(f"Loaded {len(manifest_accessions)} accessions from manifest: {manifest_path}")
    print(f"Resolved {len(metadata_rows)} Compleasm metadata rows.")
    print(f"Missing accession roots in Compleasm metadata: {len(missing_roots)}")

    for _, metadata_row in metadata_rows.iterrows():
        accession = str(metadata_row["accession"]).strip()
        accession_root = str(metadata_row["accession_root"]).strip()

        cds_fasta = repath(metadata_row["cds_fasta"], genomes_dir=genomes_dir)

        if not cds_fasta.exists():
            print(f"WARNING: cds_fasta missing for {accession}: {cds_fasta}")
            skipped_missing_cds_fasta += 1
            continue

        rows = evaluate_fasta(
            fasta_path=cds_fasta,
            accession=accession,
            accession_root=accession_root,
        )

        if not rows:
            print(f"WARNING: no FASTA records found for {accession}: {cds_fasta}")
            skipped_empty_fastas += 1
            continue

        all_rows.extend(rows)

        summary = summarize(rows)
        summary_rows.append({
            "accession": accession,
            "accession_root": accession_root,
            "cds_fasta": str(cds_fasta),
            **summary,
        })

        print(
            f"Processed {accession}: "
            f"{summary['passed_raw_cds_qc']}/{summary['total_sequences']} passed"
        )

    if not all_rows:
        raise ValueError(
            "No raw ortholog CDS records were evaluated. "
            "Check manifest, metadata.csv, and cds_fasta paths."
        )

    cohort_label = manifest_path.stem

    detailed_output = outdir / f"{cohort_label}_raw_ortholog_validity.tsv"
    summary_output = outdir / f"{cohort_label}_raw_ortholog_validity_summary.tsv"

    write_tsv(all_rows, detailed_output)
    write_tsv(summary_rows, summary_output)

    total_summary = summarize(all_rows)

    print("\nFinished raw ortholog CDS QC.")
    print(f"Accessions requested: {len(manifest_accessions)}")
    print(f"Compleasm rows resolved: {len(metadata_rows)}")
    print(f"Missing accession roots in metadata: {len(missing_roots)}")
    print(f"Rows skipped because cds_fasta was missing: {skipped_missing_cds_fasta}")
    print(f"Rows skipped because FASTA was empty: {skipped_empty_fastas}")
    print(f"Total sequences: {total_summary['total_sequences']}")
    print(f"Passed raw CDS QC: {total_summary['passed_raw_cds_qc']}")
    print(f"Failed raw CDS QC: {total_summary['failed_raw_cds_qc']}")
    print(f"Failed length % 3 check: {total_summary['failed_length_mod_3_check']}")
    print(f"Failed internal stop check: {total_summary['failed_internal_stop_check']}")
    print(f"Failed invalid base check: {total_summary['failed_invalid_base_check']}")
    print(
        "Sequences with terminal stop codon: "
        f"{total_summary['sequences_with_terminal_stop_codon']}"
    )
    print(f"Detailed QC table written to: {detailed_output}")
    print(f"Summary QC table written to: {summary_output}")


if __name__ == "__main__":
    main()
