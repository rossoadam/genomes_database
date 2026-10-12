#!/usr/bin/env python3
"""Compress genome assemblies that are not active in a project manifest.

Usage
-----
python compress_temp.py /path/to/genomes \
    --active-manifest /path/to/project_manifest.csv [--dry-run]

The active manifest protects listed accession directories from storage changes.
Every other accession directory is swept recursively for uncompressed FASTA
(.fna, .fasta, .fa) and GFF (.gff, .gff3) files. Files are gzip-compressed
atomically through a temporary ``.part`` file, verified before the original is
removed, and metadata path references are synchronized afterward.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import gzip
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Iterable, Optional


PLAIN_SUFFIXES = (".fna", ".fasta", ".fa", ".gff", ".gff3")
GZIP_SUFFIXES = tuple(suffix + ".gz" for suffix in PLAIN_SUFFIXES)
PATH_FIELDS = (
    "path_to_fna",
    "path_to_gff",
    "path_to_lift_gff",
    "drive_to_fna",
    "drive_to_gff",
    "drive_to_lift_gff",
)
SYSTEM_DIRS = {"records", "trash"}


def now_iso() -> str:
    return datetime.datetime.now().isoformat(timespec="seconds")


def require_genomes_dir(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if resolved.name != "genomes":
        raise ValueError(
            f"Expected a directory ending with 'genomes', got: {resolved}"
        )
    if not resolved.is_dir():
        raise FileNotFoundError(f"Genome directory does not exist: {resolved}")
    return resolved


def is_nonempty_file(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def is_plain_target(path: Path) -> bool:
    lower = path.name.lower()
    return path.is_file() and any(lower.endswith(suffix) for suffix in PLAIN_SUFFIXES)


def is_gzip_target(path: Path) -> bool:
    lower = path.name.lower()
    return path.is_file() and any(lower.endswith(suffix) for suffix in GZIP_SUFFIXES)


def gzip_target(path: Path) -> Path:
    return Path(str(path) + ".gz")


def read_csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        return [], []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def write_csv_rows_atomic(
    path: Path,
    fieldnames: list[str],
    rows: Iterable[dict[str, Any]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    temporary.unlink(missing_ok=True)
    try:
        with temporary.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return payload


def write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    temporary.unlink(missing_ok=True)
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def read_json_lines(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip():
                continue
            try:
                payload = json.loads(raw_line)
            except json.JSONDecodeError as exc:
                raise ValueError(
                    f"Invalid JSON on line {line_number} of {path}: {exc}"
                ) from exc
            if isinstance(payload, dict):
                rows.append(payload)
    return rows


def write_json_lines_atomic(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    temporary.unlink(missing_ok=True)
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            for row in rows:
                handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True))
                handle.write("\n")
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def sha256_plain(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    total = 0
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
            total += len(chunk)
    return digest.hexdigest(), total


def sha256_gzip(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    total = 0
    with gzip.open(path, "rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
            total += len(chunk)
    return digest.hexdigest(), total


def gzip_matches_plain(compressed: Path, plain: Path) -> bool:
    if not is_nonempty_file(compressed) or not is_nonempty_file(plain):
        return False
    try:
        return sha256_gzip(compressed) == sha256_plain(plain)
    except (OSError, EOFError, gzip.BadGzipFile):
        return False


def atomic_compress(source: Path) -> tuple[Path, bool]:
    """Gzip *source* through .part and remove source only after verification."""

    if not is_plain_target(source):
        raise ValueError(f"Unsupported compression target: {source}")
    if not is_nonempty_file(source):
        raise ValueError(f"Cannot compress missing or empty file: {source}")

    destination = gzip_target(source)
    if is_nonempty_file(destination):
        if gzip_matches_plain(destination, source):
            source.unlink()
            return destination, True
        raise FileExistsError(
            f"Refusing to overwrite non-equivalent existing gzip file: {destination}"
        )
    if destination.exists():
        destination.unlink()

    partial = destination.with_name(destination.name + ".part")
    partial.unlink(missing_ok=True)
    source_stat = source.stat()
    try:
        source_digest = hashlib.sha256()
        source_size = 0
        with source.open("rb") as input_handle, partial.open("wb") as raw_output:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                fileobj=raw_output,
                mtime=int(source_stat.st_mtime),
            ) as output_handle:
                while True:
                    chunk = input_handle.read(1024 * 1024)
                    if not chunk:
                        break
                    output_handle.write(chunk)
                    source_digest.update(chunk)
                    source_size += len(chunk)

        if not is_nonempty_file(partial):
            raise RuntimeError(f"Compression produced an empty file: {partial}")

        compressed_hash, expanded_size = sha256_gzip(partial)
        if expanded_size != source_size or compressed_hash != source_digest.hexdigest():
            raise RuntimeError(f"Compression verification failed for {source}")

        os.replace(partial, destination)
        try:
            os.chmod(destination, source_stat.st_mode)
            os.utime(destination, (source_stat.st_atime, source_stat.st_mtime))
        except OSError:
            pass
        source.unlink()
        return destination, True
    except Exception:
        partial.unlink(missing_ok=True)
        raise


def parse_active_accessions(manifest_path: Path) -> set[str]:
    resolved = manifest_path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"Active manifest does not exist: {resolved}")
    fieldnames, rows = read_csv_rows(resolved)
    if "accession" not in fieldnames:
        raise ValueError(
            f"Active manifest lacks required 'accession' column: {resolved}"
        )
    return {
        str(row.get("accession", "") or "").strip()
        for row in rows
        if str(row.get("accession", "") or "").strip()
    }


def discover_accession_dirs(genomes_dir: Path) -> list[Path]:
    directories: list[Path] = []
    for child in sorted(genomes_dir.iterdir()):
        if (
            not child.is_dir()
            or child.name in SYSTEM_DIRS
            or child.name.startswith(".")
        ):
            continue
        has_metadata = (child / "source_metadata.json").is_file()
        has_target = any(
            is_plain_target(path) or is_gzip_target(path)
            for path in child.rglob("*")
            if path.is_file()
        )
        if has_metadata or has_target:
            directories.append(child)
    return directories


def canonical_candidates(accession_dir: Path, accession: str) -> dict[str, tuple[Path, Path]]:
    fna = accession_dir / f"{accession}_genomic.fna"
    gff = accession_dir / f"{accession}_genomic.gff"
    lift = accession_dir / f"{accession}_lift.gff"
    return {
        "fna": (fna, gzip_target(fna)),
        "gff": (gff, gzip_target(gff)),
        "lift": (lift, gzip_target(lift)),
    }


def find_accession_files(
    genomes_dir: Path,
    accession: str,
) -> tuple[str, str, str]:
    accession_dir = genomes_dir / accession
    if not accession_dir.is_dir():
        return "", "", ""
    canonical = canonical_candidates(accession_dir, accession)

    def first_nonempty(pair: tuple[Path, Path]) -> str:
        for path in pair:  # prefer uncompressed active copy
            if is_nonempty_file(path):
                return str(path.resolve())
        return ""

    found_fna = first_nonempty(canonical["fna"])
    found_gff = first_nonempty(canonical["gff"])
    found_lift = first_nonempty(canonical["lift"])

    fasta_suffixes = (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")
    gff_suffixes = (".gff", ".gff3", ".gff.gz", ".gff3.gz")
    fasta_candidates: list[Path] = []
    gff_candidates: list[Path] = []
    lift_candidates: list[Path] = []

    for path in accession_dir.rglob("*"):
        if not is_nonempty_file(path) or path.name.endswith(".part"):
            continue
        lower = path.name.lower()
        if not found_lift and "lift" in lower and lower.endswith(gff_suffixes):
            lift_candidates.append(path)
            continue
        if not found_gff and "lift" not in lower and lower.endswith(gff_suffixes):
            gff_candidates.append(path)
        if not found_fna and lower.endswith(fasta_suffixes):
            if any(token in lower for token in ("genomic", "genome")):
                fasta_candidates.append(path)

    def sort_key(path: Path) -> tuple[int, int, str]:
        return (
            1 if path.name.lower().endswith(".gz") else 0,
            len(path.parts),
            str(path),
        )

    if not found_fna and fasta_candidates:
        found_fna = str(sorted(fasta_candidates, key=sort_key)[0].resolve())
    if not found_gff and gff_candidates:
        found_gff = str(sorted(gff_candidates, key=sort_key)[0].resolve())
    if not found_lift and lift_candidates:
        found_lift = str(sorted(lift_candidates, key=sort_key)[0].resolve())
    return found_fna, found_gff, found_lift


def replace_basename(existing_path: str, local_path: str) -> str:
    """Preserve an existing path root while replacing its final filename."""

    if not existing_path or not local_path:
        return "" if not local_path else existing_path
    return str(Path(existing_path).parent / Path(local_path).name)


def path_values_for_accession(
    genomes_dir: Path,
    accession: str,
    existing: Optional[dict[str, Any]] = None,
) -> dict[str, str]:
    found_fna, found_gff, found_lift = find_accession_files(genomes_dir, accession)
    existing = existing or {}
    return {
        "path_to_fna": found_fna,
        "path_to_gff": found_gff,
        "path_to_lift_gff": found_lift,
        "drive_to_fna": replace_basename(
            str(existing.get("drive_to_fna", "") or ""), found_fna
        ),
        "drive_to_gff": replace_basename(
            str(existing.get("drive_to_gff", "") or ""), found_gff
        ),
        "drive_to_lift_gff": replace_basename(
            str(existing.get("drive_to_lift_gff", "") or ""), found_lift
        ),
    }


def catalog_accession(row: dict[str, Any]) -> str:
    accession = str(row.get("accession", "") or "").strip()
    if accession:
        return accession
    normalized = row.get("normalized_record")
    if isinstance(normalized, dict):
        return str(normalized.get("accession", "") or "").strip()
    return ""


def ensure_fields(fieldnames: list[str], additions: Iterable[str]) -> list[str]:
    result = list(fieldnames)
    for field in additions:
        if field not in result:
            result.append(field)
    return result


def synchronize_metadata(
    genomes_dir: Path,
    accession_dirs: list[Path],
    *,
    events_by_accession: dict[str, list[dict[str, str]]],
    active_manifest: Path,
    dry_run: bool,
) -> list[Path]:
    """Synchronize repository metadata after the filesystem sweep."""

    records_dir = genomes_dir / "records"
    metadata_csv = records_dir / "genomes_metadata.csv"
    catalog_jsonl = records_dir / "genomes_catalog.jsonl"
    species_index_csv = records_dir / "species_assembly_index.csv"
    accessions = {path.name for path in accession_dirs}
    timestamp = now_iso()
    would_update: list[Path] = []

    metadata_fields, metadata_rows = read_csv_rows(metadata_csv)
    metadata_by_accession = {
        str(row.get("accession", "") or ""): row
        for row in metadata_rows
        if str(row.get("accession", "") or "")
    }

    path_values: dict[str, dict[str, str]] = {}
    for accession in accessions:
        path_values[accession] = path_values_for_accession(
            genomes_dir,
            accession,
            metadata_by_accession.get(accession),
        )

    if metadata_csv.exists():
        would_update.append(metadata_csv)
        if not dry_run:
            metadata_fields = ensure_fields(
                metadata_fields,
                (*PATH_FIELDS, "metadata_last_updated"),
            )
            for row in metadata_rows:
                accession = str(row.get("accession", "") or "")
                if accession not in accessions:
                    continue
                row.update(path_values[accession])
                row["metadata_last_updated"] = timestamp
            write_csv_rows_atomic(metadata_csv, metadata_fields, metadata_rows)

    # Update/create source_metadata.json for every accession that changed.
    for accession, changes in sorted(events_by_accession.items()):
        source_metadata = genomes_dir / accession / "source_metadata.json"
        would_update.append(source_metadata)
        if dry_run:
            continue
        payload = read_json(source_metadata)
        payload.setdefault("requested_accession", accession)
        payload["metadata_last_updated"] = timestamp
        normalized = payload.get("normalized_record")
        if not isinstance(normalized, dict):
            normalized = {"accession": accession}
            payload["normalized_record"] = normalized
        normalized.update(path_values[accession])
        normalized["metadata_last_updated"] = timestamp
        events = payload.get("storage_events")
        if not isinstance(events, list):
            events = []
            payload["storage_events"] = events
        event = {
            "timestamp": timestamp,
            "action": "compress_inactive",
            "active_manifest": str(active_manifest.expanduser().resolve()),
            "path_changes": changes,
        }
        events.append(event)
        payload["storage_state"] = {
            "updated_at": timestamp,
            **path_values[accession],
        }
        write_json_atomic(source_metadata, payload)

    if catalog_jsonl.exists():
        would_update.append(catalog_jsonl)
        if not dry_run:
            catalog_rows = read_json_lines(catalog_jsonl)
            for row in catalog_rows:
                accession = catalog_accession(row)
                if accession not in accessions:
                    continue
                normalized = row.get("normalized_record")
                if isinstance(normalized, dict):
                    # Preserve any drive path root that may exist only in catalog.
                    values = path_values_for_accession(
                        genomes_dir,
                        accession,
                        normalized,
                    )
                    normalized.update(values)
                    normalized["metadata_last_updated"] = timestamp
                if accession in events_by_accession:
                    events = row.get("storage_events")
                    if not isinstance(events, list):
                        events = []
                        row["storage_events"] = events
                    events.append(
                        {
                            "timestamp": timestamp,
                            "action": "compress_inactive",
                            "active_manifest": str(active_manifest.expanduser().resolve()),
                            "path_changes": events_by_accession[accession],
                        }
                    )
            write_json_lines_atomic(catalog_jsonl, catalog_rows)

    if species_index_csv.exists():
        would_update.append(species_index_csv)
        if not dry_run:
            index_fields, index_rows = read_csv_rows(species_index_csv)
            index_fields = ensure_fields(index_fields, PATH_FIELDS)
            for row in index_rows:
                accession = str(row.get("accession", "") or "")
                if accession not in accessions:
                    continue
                row.update(path_values[accession])
            write_csv_rows_atomic(species_index_csv, index_fields, index_rows)

    # Preserve order while de-duplicating paths.
    seen: set[Path] = set()
    unique: list[Path] = []
    for path in would_update:
        if path not in seen:
            seen.add(path)
            unique.append(path)
    return unique


def build_plan(
    accession_dirs: list[Path],
    active_accessions: set[str],
) -> tuple[list[Path], dict[str, list[Path]]]:
    protected: list[Path] = []
    candidates: dict[str, list[Path]] = {}
    for accession_dir in accession_dirs:
        accession = accession_dir.name
        if accession in active_accessions:
            protected.append(accession_dir)
            continue
        files = [
            path
            for path in accession_dir.rglob("*")
            if path.is_file()
            and not path.name.endswith(".part")
            and is_plain_target(path)
        ]
        files.sort(key=lambda path: str(path))
        candidates[accession] = files
    return protected, candidates


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compress FASTA/GFF files for genome assemblies that are not active "
            "in a designated project manifest, while synchronizing repository metadata."
        )
    )
    parser.add_argument(
        "genomes_dir",
        help="Path to the genomes directory (final component must be 'genomes')",
    )
    parser.add_argument(
        "--active-manifest",
        required=True,
        help="Project manifest CSV whose accession column protects active assemblies",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report planned file and metadata changes without writing anything",
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    genomes_dir = require_genomes_dir(Path(args.genomes_dir))
    active_manifest = Path(args.active_manifest).expanduser().resolve()
    active_accessions = parse_active_accessions(active_manifest)
    accession_dirs = discover_accession_dirs(genomes_dir)
    protected, candidates = build_plan(accession_dirs, active_accessions)

    print(f"Genome repository: {genomes_dir}")
    print(f"Active manifest:   {active_manifest}")
    print(f"Manifest accessions protected: {len(active_accessions)}")
    print(f"Accession directories scanned: {len(accession_dirs)}")

    for accession_dir in protected:
        print(f"SKIP active accession: {accession_dir.name}")

    planned_file_count = sum(len(paths) for paths in candidates.values())
    print(f"Uncompressed inactive FASTA/GFF files found: {planned_file_count}")
    for accession, paths in candidates.items():
        for source in paths:
            print(f"{'WOULD COMPRESS' if args.dry_run else 'COMPRESS'}: {source} -> {gzip_target(source)}")

    if args.dry_run:
        metadata_targets = [
            path
            for path in (
                genomes_dir / "records" / "genomes_metadata.csv",
                genomes_dir / "records" / "genomes_catalog.jsonl",
                genomes_dir / "records" / "species_assembly_index.csv",
            )
            if path.exists()
        ]
        metadata_targets.extend(
            genomes_dir / accession / "source_metadata.json"
            for accession, files in candidates.items()
            if files
        )
        print("\nDRY RUN: no files were changed.")
        print("Metadata files that would be synchronized:")
        for path in metadata_targets:
            print(f"  - {path}")
        return 0

    events_by_accession: dict[str, list[dict[str, str]]] = {}
    errors: list[str] = []
    for accession, paths in candidates.items():
        changes: list[dict[str, str]] = []
        for source in paths:
            try:
                destination, changed = atomic_compress(source)
                if changed:
                    changes.append({"from": str(source), "to": str(destination)})
            except Exception as exc:
                message = f"{source}: {type(exc).__name__}: {exc}"
                errors.append(message)
                print(f"ERROR: {message}", file=sys.stderr)
        if changes:
            events_by_accession[accession] = changes

    # Synchronize all scanned accessions, including already-compressed inactive
    # assemblies and protected active assemblies, so metadata reflects disk state.
    updated_metadata = synchronize_metadata(
        genomes_dir,
        accession_dirs,
        events_by_accession=events_by_accession,
        active_manifest=active_manifest,
        dry_run=False,
    )

    changed_files = sum(len(changes) for changes in events_by_accession.values())
    print(f"\nCompressed files: {changed_files}")
    print(f"Accessions changed: {len(events_by_accession)}")
    print("Metadata synchronized:")
    for path in updated_metadata:
        print(f"  - {path}")

    if errors:
        print(
            f"Completed with {len(errors)} file error(s). Successful compression "
            "changes were retained and metadata was synchronized.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
