#!/usr/bin/env python3
"""
Compare compleasm CDS FASTA records to BUSCO single-copy FASTA files.

Assumptions (based on your message):
- The filename stem of each BUSCO fasta in single_copy_busco_sequences/
  matches the FASTA header ID (first token) in the compleasm fasta.
  Example:
    BUSCO file:   EOG0ABC123.fna
    compleasm hdr: >EOG0ABC123 some other stuff...

Outputs:
- A CSV report summarizing matches and sequence similarity.
"""

from __future__ import annotations

import argparse
import csv
import difflib
from pathlib import Path
from typing import Dict, List, Tuple


def read_fasta_as_dict(path: Path) -> Dict[str, str]:
    """
    Read a (multi)fasta into dict: {record_id: sequence}
    record_id = first token after '>' (up to whitespace).
    """
    records: Dict[str, List[str]] = {}
    current_id: str | None = None

    with path.open("r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                if current_id in records:
                    # If duplicates, append; later we will join and keep first encountered
                    # but duplicates usually indicate an upstream issue.
                    pass
                records.setdefault(current_id, [])
            else:
                if current_id is None:
                    raise ValueError(f"FASTA parsing error: sequence before header in {path}")
                records[current_id].append(line.strip())

    return {rid: "".join(seq_lines).upper() for rid, seq_lines in records.items()}


def read_single_fasta(path: Path) -> Tuple[str, str]:
    """
    Read a fasta file that should contain 1 record.
    Returns (record_id, sequence).
    If multiple records exist, it reads the first and ignores the rest.
    """
    rec_id = None
    seq_lines: List[str] = []

    with path.open("r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if rec_id is None:
                    rec_id = line[1:].split()[0]
                else:
                    # second header -> stop (we only take first record)
                    break
            else:
                if rec_id is None:
                    raise ValueError(f"FASTA parsing error: sequence before header in {path}")
                seq_lines.append(line.strip())

    if rec_id is None:
        raise ValueError(f"No FASTA records found in {path}")
    return rec_id, "".join(seq_lines).upper()


def mismatch_stats_if_same_length(a: str, b: str) -> Tuple[int, float]:
    """
    If same length, return (mismatches, identity_fraction).
    """
    assert len(a) == len(b)
    mismatches = 0
    matches = 0
    for ca, cb in zip(a, b):
        if ca == cb:
            matches += 1
        else:
            mismatches += 1
    identity = matches / len(a) if a else 0.0
    return mismatches, identity


def approx_similarity(a: str, b: str) -> float:
    """
    Approximate similarity using difflib SequenceMatcher.
    Not a biological alignment; useful as a quick diagnostic when lengths differ.
    Returns a value in [0, 1].
    """
    sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
    # ratio is 2*M / (len(a)+len(b)) where M is total matches in matching blocks
    return sm.ratio()


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Compare compleasm CDS fasta to BUSCO single-copy fasta files."
    )
    ap.add_argument(
        "--compleasm-fasta",
        required=True,
        type=Path,
        help="Path to podarcis_muralis_cds_compleasm.fasta (multi-fasta).",
    )
    ap.add_argument(
        "--busco-dir",
        required=True,
        type=Path,
        help="Path to BUSCO single_copy_busco_sequences directory (contains many fasta files).",
    )
    ap.add_argument(
        "--out-csv",
        required=True,
        type=Path,
        help="Path to write comparison report CSV.",
    )
    args = ap.parse_args()

    compleasm_fa: Path = args.compleasm_fasta
    busco_dir: Path = args.busco_dir
    out_csv: Path = args.out_csv

    if not compleasm_fa.exists():
        raise FileNotFoundError(compleasm_fa)
    if not busco_dir.exists():
        raise FileNotFoundError(busco_dir)
    if not busco_dir.is_dir():
        raise NotADirectoryError(busco_dir)

    compleasm = read_fasta_as_dict(compleasm_fa)

    # BUSCO single-copy files can be .fna, .fa, .fasta, etc.
    fasta_exts = {".fa", ".fasta", ".fna", ".fas", ".faa"}
    busco_files = sorted([p for p in busco_dir.iterdir() if p.is_file() and p.suffix.lower() in fasta_exts])

    rows = []
    matched = 0
    exact = 0
    missing = 0

    for bf in busco_files:
        busco_key = bf.stem  # filename without extension, per your assumption
        busco_id, busco_seq = read_single_fasta(bf)

        # Primary match: filename stem
        comp_seq = compleasm.get(busco_key)

        # Fallback match: sometimes BUSCO header itself might be the key
        if comp_seq is None:
            comp_seq = compleasm.get(busco_id)

        if comp_seq is None:
            missing += 1
            rows.append(
                {
                    "busco_file": bf.name,
                    "busco_key": busco_key,
                    "busco_header_id": busco_id,
                    "found_in_compleasm": "NO",
                    "exact_match": "",
                    "len_busco": len(busco_seq),
                    "len_compleasm": "",
                    "mismatches_same_len": "",
                    "identity_same_len": "",
                    "approx_similarity_if_diff_len": "",
                }
            )
            continue

        matched += 1

        is_exact = "YES" if comp_seq == busco_seq else "NO"
        if is_exact == "YES":
            exact += 1

        if len(comp_seq) == len(busco_seq):
            mism, ident = mismatch_stats_if_same_length(comp_seq, busco_seq)
            approx = ""
        else:
            mism = ""
            ident = ""
            approx = f"{approx_similarity(comp_seq, busco_seq):.6f}"

        rows.append(
            {
                "busco_file": bf.name,
                "busco_key": busco_key,
                "busco_header_id": busco_id,
                "found_in_compleasm": "YES",
                "exact_match": is_exact,
                "len_busco": len(busco_seq),
                "len_compleasm": len(comp_seq),
                "mismatches_same_len": mism,
                "identity_same_len": f"{ident:.6f}" if isinstance(ident, float) else "",
                "approx_similarity_if_diff_len": approx,
            }
        )

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "busco_file",
        "busco_key",
        "busco_header_id",
        "found_in_compleasm",
        "exact_match",
        "len_busco",
        "len_compleasm",
        "mismatches_same_len",
        "identity_same_len",
        "approx_similarity_if_diff_len",
    ]
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print("Done.")
    print(f"BUSCO files scanned: {len(busco_files)}")
    print(f"Matched in compleasm: {matched}")
    print(f"Exact sequence matches: {exact}")
    print(f"Missing in compleasm: {missing}")
    print(f"CSV report: {out_csv}")


if __name__ == "__main__":
    main()
