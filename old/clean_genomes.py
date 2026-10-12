#!/usr/bin/env python3
"""
clean_genomes.py

Conservative cleanup utility for a local genome repository.

The script identifies local NCBI assembly directories that are eligible for
cleanup only when all required conditions can be established confidently:

  * the accession is an NCBI assembly (GCA_ or GCF_)
  * NCBI Datasets confirms that it belongs to the requested taxon
  * it is not a current NCBI reference assembly
  * it is not protected by any current project manifest
  * it is not marked as locally preferred
  * its local accession directory exists and is safe to move
  * required local metadata are available and readable

Dry-run is the default.  Execute mode never permanently deletes genome data;
it moves candidate accession directories into a timestamped trash archive and
then removes successfully moved accessions from active repository metadata.

The script also includes an offline self-test suite:

    python clean_genomes.py --self-test

Typical usage:

    python clean_genomes.py /path/to/genomes --taxon squamata --dry-run

    python clean_genomes.py /path/to/genomes --taxon squamata --execute

    python clean_genomes.py /path/to/genomes --taxon squamata \
        --accessions suspicious_accessions.txt --dry-run

Notes on dry-run:
-----------------
Dry-run never moves accession directories and never changes active repository
metadata or project manifests.  It does write the requested audit TSV under
records/cleanup/ so the assessment can be reviewed and retained.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import errno
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import traceback
import unittest
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Set, Tuple
from unittest import mock


SCRIPT_VERSION = "1.0.0"
NCBI_ACCESSION_RE = re.compile(r"^(GC[AF]_\d+\.\d+)$", re.IGNORECASE)
NCBI_ACCESSION_ANY_RE = re.compile(r"\bGC[AF]_\d+(?:\.\d+)?\b", re.IGNORECASE)
GWH_ACCESSION_RE = re.compile(r"^GWH", re.IGNORECASE)
TRUE_VALUES = {"1", "true", "t", "yes", "y"}

AUDIT_COLUMNS = [
    "accession",
    "repository",
    "organism_name",
    "species_key",
    "assembly_name",
    "assembly_level",
    "source_database",
    "refseq_category",
    "local_preferred",
    "current_ncbi_reference",
    "protected_by_manifest",
    "protecting_manifests",
    "candidate_for_cleanup",
    "decision",
    "decision_reason",
    "accession_directory",
    "path_to_fna",
    "path_to_gff",
    "compressed_or_uncompressed",
    "directory_size_bytes",
    "directory_size_gb",
    "action",
    "trash_destination",
    "timestamp",
]


class SafetyAbort(RuntimeError):
    """Raised when a safety-critical prerequisite prevents cleanup."""


@dataclass
class NCBIDiscovery:
    reference_accessions: Set[str] = field(default_factory=set)
    taxon_accessions: Set[str] = field(default_factory=set)
    metadata_by_accession: Dict[str, Dict[str, str]] = field(default_factory=dict)
    reference_command: List[str] = field(default_factory=list)
    taxon_command: List[str] = field(default_factory=list)
    reference_stdout: str = ""
    taxon_stdout: str = ""
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    reference_query_ok: bool = False
    taxon_query_ok: bool = False

    @property
    def safe_for_candidate_evaluation(self) -> bool:
        return self.reference_query_ok and self.taxon_query_ok and not self.errors


@dataclass
class RepositoryContext:
    genomes_dir: Path
    metadata_rows: Dict[str, Dict[str, str]]
    metadata_file_exists: bool
    preferred_accessions: Set[str]
    protecting_manifests: Dict[str, Set[str]]
    manifest_files: List[Path]
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


@dataclass
class Assessment:
    rows: List[Dict[str, Any]]
    context: RepositoryContext
    ncbi: NCBIDiscovery
    local_existing_accessions: Set[str]
    candidate_accessions: List[str]
    fatal_errors: List[str]
    warnings: List[str]
    timestamp: str
    taxon: str
    accession_filter: Optional[Set[str]]


# ---------------------------------------------------------------------------
# Small utility helpers
# ---------------------------------------------------------------------------


def now_timestamp() -> str:
    return _dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def iso_now() -> str:
    return _dt.datetime.now().astimezone().isoformat(timespec="seconds")


def normalize_accession(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip().upper()


def normalize_text(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def bool_text(value: bool) -> str:
    return "True" if value else "False"


def is_true(value: Any) -> bool:
    return normalize_text(value).lower() in TRUE_VALUES


def relpath_for_report(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except Exception:
        return str(path)


def command_string(cmd: Sequence[str]) -> str:
    # Human-readable shell-ish command without introducing a shell dependency.
    import shlex

    return " ".join(shlex.quote(str(x)) for x in cmd)


def directory_size_bytes(path: Path) -> int:
    total = 0
    try:
        for base, dirs, files in os.walk(path, followlinks=False):
            base_path = Path(base)
            # Do not descend through symlinked directories.
            dirs[:] = [d for d in dirs if not (base_path / d).is_symlink()]
            for name in files:
                p = base_path / name
                try:
                    if p.is_symlink():
                        continue
                    total += p.stat().st_size
                except OSError:
                    # Size is for reporting only; inability to stat an individual
                    # file should not make deletion more permissive.
                    continue
    except OSError:
        return 0
    return total


def discover_sequence_files(accession_dir: Path) -> Tuple[str, str, str]:
    fna_candidates: List[Path] = []
    gff_candidates: List[Path] = []
    try:
        for p in accession_dir.rglob("*"):
            if not p.is_file():
                continue
            low = p.name.lower()
            if low.endswith(("_genomic.fna", "_genomic.fna.gz", ".fna", ".fna.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz")):
                fna_candidates.append(p)
            if low.endswith(("_genomic.gff", "_genomic.gff.gz", ".gff", ".gff.gz", ".gff3", ".gff3.gz")):
                gff_candidates.append(p)
    except OSError:
        pass

    def best(candidates: List[Path]) -> str:
        if not candidates:
            return ""
        candidates.sort(key=lambda p: (len(p.parts), str(p)))
        return str(candidates[0])

    fna = best(fna_candidates)
    gff = best(gff_candidates)
    compression_states: Set[str] = set()
    for value in (fna, gff):
        if value:
            compression_states.add("compressed" if value.lower().endswith(".gz") else "uncompressed")
    if not compression_states:
        state = "unknown"
    elif len(compression_states) == 1:
        state = next(iter(compression_states))
    else:
        state = "mixed"
    return fna, gff, state


def species_key_from_name(name: str) -> str:
    name = normalize_text(name)
    if not name:
        return ""
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    temp_path = Path(temp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_path, path)
    except Exception:
        try:
            temp_path.unlink(missing_ok=True)
        except Exception:
            pass
        raise


def write_tsv(path: Path, rows: Sequence[Mapping[str, Any]], columns: Sequence[str] = AUDIT_COLUMNS) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    temp_path = Path(temp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", extrasaction="ignore", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({col: row.get(col, "") for col in columns})
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_path, path)
    except Exception:
        try:
            temp_path.unlink(missing_ok=True)
        except Exception:
            pass
        raise


# ---------------------------------------------------------------------------
# NCBI Datasets discovery and JSON parsing
# ---------------------------------------------------------------------------


def run_command(cmd: Sequence[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(cmd),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )


def parse_json_lines(text: str, label: str) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        try:
            obj = json.loads(line)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Malformed JSON from {label} at line {lineno}: {exc}") from exc
        if not isinstance(obj, dict):
            raise ValueError(f"Unexpected non-object JSON from {label} at line {lineno}")
        records.append(obj)
    return records


def recursively_find_accessions(obj: Any) -> Set[str]:
    found: Set[str] = set()
    if isinstance(obj, Mapping):
        for value in obj.values():
            found.update(recursively_find_accessions(value))
    elif isinstance(obj, list):
        for value in obj:
            found.update(recursively_find_accessions(value))
    elif isinstance(obj, str):
        for m in NCBI_ACCESSION_ANY_RE.finditer(obj):
            acc = normalize_accession(m.group(0))
            # Cleanup operates on versioned local directories.  Unversioned
            # accessions are retained in metadata but do not match local dirs.
            found.add(acc)
    return found


def get_path(obj: Mapping[str, Any], *path: str) -> Any:
    cur: Any = obj
    for key in path:
        if not isinstance(cur, Mapping) or key not in cur:
            return None
        cur = cur[key]
    return cur


def recursive_first_scalar_by_key(obj: Any, keys: Set[str]) -> str:
    if isinstance(obj, Mapping):
        for k, v in obj.items():
            if str(k).lower() in keys and isinstance(v, (str, int, float, bool)):
                return normalize_text(v)
        for v in obj.values():
            found = recursive_first_scalar_by_key(v, keys)
            if found:
                return found
    elif isinstance(obj, list):
        for v in obj:
            found = recursive_first_scalar_by_key(v, keys)
            if found:
                return found
    return ""


def clean_enum(value: str) -> str:
    value = normalize_text(value)
    for prefix in (
        "ASSEMBLY_LEVEL_",
        "SOURCE_DATABASE_",
        "REFSEQ_CATEGORY_",
    ):
        if value.upper().startswith(prefix):
            value = value[len(prefix) :]
            break
    return value.replace("_", " ").strip().lower()


def extract_ncbi_metadata(report: Mapping[str, Any]) -> Dict[str, str]:
    organism_name = normalize_text(get_path(report, "organism", "organism_name"))
    if not organism_name:
        organism_name = recursive_first_scalar_by_key(report, {"organism_name", "scientific_name", "tax_name"})

    assembly_name = normalize_text(get_path(report, "assembly_info", "assembly_name"))
    if not assembly_name:
        assembly_name = recursive_first_scalar_by_key(report, {"assembly_name"})

    assembly_level = normalize_text(get_path(report, "assembly_info", "assembly_level"))
    if not assembly_level:
        assembly_level = recursive_first_scalar_by_key(report, {"assembly_level"})

    source_database = normalize_text(report.get("source_database"))
    if not source_database:
        source_database = recursive_first_scalar_by_key(report, {"source_database"})

    refseq_category = normalize_text(get_path(report, "assembly_info", "refseq_category"))
    if not refseq_category:
        refseq_category = recursive_first_scalar_by_key(report, {"refseq_category"})

    return {
        "organism_name": organism_name,
        "species_key": species_key_from_name(organism_name),
        "assembly_name": assembly_name,
        "assembly_level": clean_enum(assembly_level),
        "source_database": clean_enum(source_database),
        "refseq_category": clean_enum(refseq_category),
    }


def discover_ncbi(
    taxon: str,
    runner: Callable[[Sequence[str]], subprocess.CompletedProcess[str]] = run_command,
) -> NCBIDiscovery:
    result = NCBIDiscovery()
    result.reference_command = [
        "datasets",
        "summary",
        "genome",
        "taxon",
        taxon,
        "--reference",
        "--assembly-level",
        "chromosome,complete",
        "--as-json-lines",
    ]
    # A separate all-assemblies query is required to establish that a local
    # non-reference accession actually belongs to the requested taxon.  Using
    # --assembly-version all is conservative for older/superseded local copies.
    result.taxon_command = [
        "datasets",
        "summary",
        "genome",
        "taxon",
        taxon,
        "--assembly-version",
        "all",
        "--as-json-lines",
    ]

    for label, cmd in (("reference", result.reference_command), ("taxon", result.taxon_command)):
        try:
            cp = runner(cmd)
        except FileNotFoundError:
            result.errors.append(
                "NCBI Datasets CLI executable 'datasets' was not found. No cleanup candidates can be authorized."
            )
            continue
        except Exception as exc:
            result.errors.append(f"Failed to run NCBI Datasets {label} query: {exc}")
            continue

        stdout = cp.stdout or ""
        stderr = cp.stderr or ""
        if label == "reference":
            result.reference_stdout = stdout
        else:
            result.taxon_stdout = stdout

        if cp.returncode != 0:
            result.errors.append(
                f"NCBI Datasets {label} query failed with exit code {cp.returncode}: {stderr.strip() or '(no stderr)'}"
            )
            continue

        try:
            reports = parse_json_lines(stdout, f"NCBI Datasets {label} query")
        except ValueError as exc:
            result.errors.append(str(exc))
            continue

        if label == "reference":
            if not reports:
                result.warnings.append("NCBI reference query returned zero records.")
            for report in reports:
                accessions = recursively_find_accessions(report)
                meta = extract_ncbi_metadata(report)
                for acc in accessions:
                    result.reference_accessions.add(acc)
                    result.metadata_by_accession.setdefault(acc, {}).update(meta)
            result.reference_query_ok = bool(reports)
        else:
            if not reports:
                result.warnings.append("NCBI taxon query returned zero records.")
            for report in reports:
                accessions = recursively_find_accessions(report)
                meta = extract_ncbi_metadata(report)
                for acc in accessions:
                    result.taxon_accessions.add(acc)
                    existing = result.metadata_by_accession.setdefault(acc, {})
                    for key, value in meta.items():
                        if value or key not in existing:
                            existing[key] = value
                    # If the general taxon record itself identifies a current
                    # reference genome, protect it even if the filtered query's
                    # assembly-level restriction omitted it.
                    if meta.get("refseq_category", "").lower() == "reference genome":
                        result.reference_accessions.add(acc)
            result.taxon_query_ok = bool(reports)

    return result


# ---------------------------------------------------------------------------
# Local metadata and manifest discovery
# ---------------------------------------------------------------------------


def read_csv_rows(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        fields = [str(x) for x in reader.fieldnames]
        rows: List[Dict[str, str]] = []
        for row in reader:
            rows.append({str(k): ("" if v is None else str(v)) for k, v in row.items() if k is not None})
    return fields, rows


def find_column(fieldnames: Sequence[str], aliases: Iterable[str]) -> Optional[str]:
    lookup = {str(name).strip().lower(): str(name) for name in fieldnames}
    for alias in aliases:
        found = lookup.get(alias.lower())
        if found:
            return found
    return None


def load_genomes_metadata(path: Path) -> Tuple[Dict[str, Dict[str, str]], bool]:
    if not path.exists():
        return {}, False
    fields, rows = read_csv_rows(path)
    acc_col = find_column(fields, ["accession", "assembly_accession", "genome_accession"])
    if not acc_col:
        raise ValueError(f"No accession column found in {path}")
    out: Dict[str, Dict[str, str]] = {}
    for row in rows:
        acc = normalize_accession(row.get(acc_col, ""))
        if not acc:
            continue
        out[acc] = row
    return out, True


def preferred_from_csv(path: Path) -> Set[str]:
    preferred: Set[str] = set()
    if not path.exists():
        return preferred
    fields, rows = read_csv_rows(path)
    acc_col = find_column(fields, ["accession", "assembly_accession", "genome_accession"])
    if not acc_col:
        raise ValueError(f"No accession column found in {path}")
    preferred_col = find_column(fields, ["is_preferred_for_species", "preferred", "is_preferred"])
    status_col = find_column(fields, ["status", "assembly_status"])
    for row in rows:
        acc = normalize_accession(row.get(acc_col, ""))
        if not acc:
            continue
        is_pref = False
        if preferred_col and is_true(row.get(preferred_col, "")):
            is_pref = True
        if status_col and normalize_text(row.get(status_col, "")).lower() == "preferred":
            is_pref = True
        if is_pref:
            preferred.add(acc)
    return preferred


def resolve_manifest_candidate(raw: str, registry_path: Path, genomes_dir: Path) -> Optional[Path]:
    raw = normalize_text(raw)
    if not raw:
        return None
    p = Path(os.path.expanduser(raw))
    candidates: List[Path] = []
    if p.is_absolute():
        candidates.append(p)
    else:
        candidates.extend(
            [
                registry_path.parent / p,
                genomes_dir / p,
                genomes_dir / "records" / p,
                genomes_dir / "records" / "project_manifests" / p,
            ]
        )
    for candidate in candidates:
        if candidate.exists() and candidate.is_file():
            return candidate.resolve()
    return None


def registry_manifest_paths(registry_path: Path, genomes_dir: Path) -> Tuple[Set[Path], List[str]]:
    paths: Set[Path] = set()
    warnings: List[str] = []
    if not registry_path.exists():
        return paths, warnings
    fields, rows = read_csv_rows(registry_path)
    likely_columns = [
        f for f in fields if "manifest" in f.lower() or f.lower() in {"path", "file", "filepath", "csv"}
    ]
    for row in rows:
        values: List[str] = []
        if likely_columns:
            values.extend(row.get(col, "") for col in likely_columns)
        else:
            values.extend(row.values())
        for raw in values:
            text = normalize_text(raw)
            if not text:
                continue
            if not (text.lower().endswith(".csv") or "manifest" in text.lower()):
                continue
            resolved = resolve_manifest_candidate(text, registry_path, genomes_dir)
            if resolved:
                paths.add(resolved)
            elif "manifest" in text.lower():
                warnings.append(f"Project registry listed a manifest path that does not exist and was skipped: {text}")
    return paths, warnings


def scan_manifests(genomes_dir: Path) -> Tuple[Dict[str, Set[str]], List[Path], List[str], List[str]]:
    manifest_dir = genomes_dir / "records" / "project_manifests"
    manifests: Set[Path] = set()
    if manifest_dir.exists():
        manifests.update(p.resolve() for p in manifest_dir.glob("*.csv") if p.is_file())

    warnings: List[str] = []
    errors: List[str] = []
    registry_path = genomes_dir / "records" / "project_registry.csv"
    try:
        registry_paths, registry_warnings = registry_manifest_paths(registry_path, genomes_dir)
        manifests.update(registry_paths)
        warnings.extend(registry_warnings)
    except Exception as exc:
        errors.append(f"Could not read project registry {registry_path}: {exc}")

    protecting: Dict[str, Set[str]] = defaultdict(set)
    for manifest in sorted(manifests):
        try:
            fields, rows = read_csv_rows(manifest)
            acc_col = find_column(fields, ["accession"])
            if not acc_col:
                raise ValueError("missing required 'accession' column")
            display_name = relpath_for_report(manifest, genomes_dir)
            for row in rows:
                acc = normalize_accession(row.get(acc_col, ""))
                if acc:
                    protecting[acc].add(display_name)
        except Exception as exc:
            errors.append(f"Could not safely read manifest {manifest}: {exc}")

    return dict(protecting), sorted(manifests), warnings, errors


def validate_genomes_catalog(path: Path) -> None:
    """Validate JSONL syntax before any execute-mode filesystem moves."""
    if not path.exists():
        return
    with path.open("r", encoding="utf-8") as handle:
        for lineno, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"malformed JSON at line {lineno}: {exc}") from exc
            if not isinstance(obj, Mapping):
                raise ValueError(f"non-object JSON record at line {lineno}")


def build_repository_context(genomes_dir: Path) -> RepositoryContext:
    records_dir = genomes_dir / "records"
    metadata_path = records_dir / "genomes_metadata.csv"
    species_index_path = records_dir / "species_assembly_index.csv"
    catalog_path = records_dir / "genomes_catalog.jsonl"

    errors: List[str] = []
    warnings: List[str] = []
    metadata_rows: Dict[str, Dict[str, str]] = {}
    metadata_exists = metadata_path.exists()

    try:
        metadata_rows, metadata_exists = load_genomes_metadata(metadata_path)
    except Exception as exc:
        errors.append(f"Could not safely read {metadata_path}: {exc}")

    preferred: Set[str] = set()
    for path in (metadata_path, species_index_path):
        try:
            preferred.update(preferred_from_csv(path))
        except Exception as exc:
            errors.append(f"Could not safely determine locally preferred assemblies from {path}: {exc}")

    try:
        validate_genomes_catalog(catalog_path)
    except Exception as exc:
        errors.append(f"Could not safely read {catalog_path}: {exc}")

    protecting, manifest_files, manifest_warnings, manifest_errors = scan_manifests(genomes_dir)
    warnings.extend(manifest_warnings)
    errors.extend(manifest_errors)

    if not metadata_exists:
        warnings.append(
            f"{metadata_path} does not exist. Local NCBI directories cannot be automatically cleaned because the expected active metadata table is unavailable."
        )

    return RepositoryContext(
        genomes_dir=genomes_dir,
        metadata_rows=metadata_rows,
        metadata_file_exists=metadata_exists,
        preferred_accessions=preferred,
        protecting_manifests=protecting,
        manifest_files=manifest_files,
        errors=errors,
        warnings=warnings,
    )


def metadata_value(row: Mapping[str, str], aliases: Iterable[str]) -> str:
    lookup = {str(k).strip().lower(): v for k, v in row.items()}
    for alias in aliases:
        value = normalize_text(lookup.get(alias.lower(), ""))
        if value:
            return value
    return ""


def local_accession_names(genomes_dir: Path, metadata_accessions: Iterable[str]) -> Set[str]:
    names: Set[str] = set()
    metadata_set = {normalize_accession(x) for x in metadata_accessions if normalize_accession(x)}
    if not genomes_dir.exists():
        return names
    for p in genomes_dir.iterdir():
        if not p.is_dir():
            continue
        name = normalize_accession(p.name)
        if NCBI_ACCESSION_RE.match(name) or GWH_ACCESSION_RE.match(name) or name in metadata_set:
            names.add(name)
    return names


# ---------------------------------------------------------------------------
# Assessment logic
# ---------------------------------------------------------------------------


def merge_row_metadata(
    accession: str,
    local_row: Mapping[str, str],
    ncbi_meta: Mapping[str, str],
) -> Dict[str, str]:
    organism_name = (
        metadata_value(local_row, ["organism_name", "organism", "species", "species_name", "scientific_name"])
        or ncbi_meta.get("organism_name", "")
    )
    species_key = (
        metadata_value(local_row, ["species_key", "species_normalized", "genus_species"])
        or ncbi_meta.get("species_key", "")
        or species_key_from_name(organism_name)
    )
    assembly_name = metadata_value(local_row, ["assembly_name", "name"]) or ncbi_meta.get("assembly_name", "")
    assembly_level = metadata_value(local_row, ["assembly_level"]) or ncbi_meta.get("assembly_level", "")
    source_database = metadata_value(local_row, ["source_database", "database", "source"]) or ncbi_meta.get("source_database", "")
    refseq_category = metadata_value(local_row, ["refseq_category"]) or ncbi_meta.get("refseq_category", "")

    if not source_database:
        if accession.startswith("GCF_"):
            source_database = "refseq"
        elif accession.startswith("GCA_"):
            source_database = "genbank"

    return {
        "organism_name": organism_name,
        "species_key": species_key,
        "assembly_name": assembly_name,
        "assembly_level": assembly_level,
        "source_database": source_database,
        "refseq_category": refseq_category,
    }


def classify_repository(accession: str) -> str:
    if NCBI_ACCESSION_RE.match(accession):
        return "NCBI"
    if GWH_ACCESSION_RE.match(accession):
        return "GWH"
    return "local_or_unknown"


def assess_repository(
    genomes_dir: Path,
    taxon: str,
    accession_filter: Optional[Set[str]] = None,
    ncbi_data: Optional[NCBIDiscovery] = None,
    timestamp: Optional[str] = None,
) -> Assessment:
    genomes_dir = genomes_dir.resolve()
    timestamp = timestamp or now_timestamp()
    context = build_repository_context(genomes_dir)
    ncbi = ncbi_data if ncbi_data is not None else discover_ncbi(taxon)

    local_existing = local_accession_names(genomes_dir, context.metadata_rows.keys())
    metadata_only_missing = set(context.metadata_rows) - local_existing
    all_accessions = sorted(local_existing | metadata_only_missing)

    fatal_errors: List[str] = []
    warnings = list(context.warnings) + list(ncbi.warnings)
    fatal_errors.extend(context.errors)
    fatal_errors.extend(ncbi.errors)

    local_ncbi = {acc for acc in local_existing if NCBI_ACCESSION_RE.match(acc)}
    local_taxon_confirmed = local_ncbi & ncbi.taxon_accessions

    if local_ncbi and not ncbi.reference_query_ok:
        fatal_errors.append(
            "Authoritative NCBI reference discovery is unavailable or returned no usable reference records while local NCBI assemblies exist."
        )
    if local_ncbi and not ncbi.taxon_query_ok:
        fatal_errors.append(
            "NCBI taxon-membership discovery is unavailable or returned no usable taxon records while local NCBI assemblies exist."
        )
    if local_taxon_confirmed and not ncbi.reference_accessions:
        fatal_errors.append(
            "NCBI reference query produced an unexpectedly empty reference-accession set for a taxon that contains local NCBI assemblies."
        )

    if not context.metadata_file_exists and local_ncbi:
        # Not a fatal global abort for reporting, but it prevents every local
        # NCBI accession from becoming a candidate.
        warnings.append("Expected genomes_metadata.csv is missing; affected NCBI accessions will be retained for manual review.")

    rows: List[Dict[str, Any]] = []
    candidates: List[str] = []
    assessment_time = iso_now()

    filter_set = None
    if accession_filter is not None:
        filter_set = {normalize_accession(x) for x in accession_filter if normalize_accession(x)}

    fatal_state = bool(fatal_errors)

    for accession in all_accessions:
        accession_dir = genomes_dir / accession
        dir_exists = accession_dir.exists() and accession_dir.is_dir()
        repository = classify_repository(accession)
        local_row = context.metadata_rows.get(accession, {})
        ncbi_meta = ncbi.metadata_by_accession.get(accession, {})
        merged = merge_row_metadata(accession, local_row, ncbi_meta)
        local_preferred = accession in context.preferred_accessions
        manifest_names = sorted(context.protecting_manifests.get(accession, set()))
        protected_by_manifest = bool(manifest_names)
        current_reference = accession in ncbi.reference_accessions
        in_taxon = accession in ncbi.taxon_accessions
        in_filter = filter_set is None or accession in filter_set

        size_bytes = directory_size_bytes(accession_dir) if dir_exists else 0
        fna = gff = ""
        compression = "unknown"
        if dir_exists:
            fna, gff, compression = discover_sequence_files(accession_dir)

        candidate = False
        decision = ""
        reason = ""
        action = "none"

        # Most conservative decisions first.
        if not dir_exists:
            decision = "manual_review"
            reason = "metadata_lists_accession_but_local_directory_is_missing"
        elif accession_dir.is_symlink():
            decision = "manual_review"
            reason = "accession_directory_is_symlink"
        elif not in_filter:
            decision = "keep_outside_accession_filter"
            reason = "accession_not_listed_in_accessions_filter"
        elif repository != "NCBI":
            decision = "keep_non_ncbi"
            reason = "not_evaluated_non_ncbi"
        elif protected_by_manifest:
            decision = "keep_manifest"
            reason = "protected_by_project_manifest"
        elif local_preferred:
            if current_reference:
                decision = "keep_reference"
                reason = "current_ncbi_reference_and_locally_preferred"
            else:
                # The user explicitly identified this conflict as a situation
                # requiring manual review.  It remains fully protected.
                decision = "manual_review"
                reason = "locally_preferred_but_not_in_current_ncbi_reference_set"
        elif current_reference:
            decision = "keep_reference"
            reason = "current_ncbi_reference"
        elif fatal_state:
            decision = "manual_review"
            reason = "safety_prerequisite_failed_no_cleanup_authorized"
        elif not context.metadata_file_exists:
            decision = "manual_review"
            reason = "genomes_metadata_csv_missing"
        elif accession not in context.metadata_rows:
            decision = "manual_review"
            reason = "local_accession_directory_missing_from_genomes_metadata_csv"
        elif not in_taxon:
            decision = "manual_review"
            reason = "taxonomic_membership_in_requested_taxon_not_established_by_ncbi"
        elif not merged.get("organism_name"):
            decision = "manual_review"
            reason = "organism_or_species_information_could_not_be_determined"
        else:
            candidate = True
            decision = "cleanup_candidate"
            reason = "ncbi_nonreference_in_requested_taxon_not_manifest_protected_not_locally_preferred"
            action = "would_move_to_trash"
            candidates.append(accession)

        row = {
            "accession": accession,
            "repository": repository,
            "organism_name": merged.get("organism_name", ""),
            "species_key": merged.get("species_key", ""),
            "assembly_name": merged.get("assembly_name", ""),
            "assembly_level": merged.get("assembly_level", ""),
            "source_database": merged.get("source_database", ""),
            "refseq_category": merged.get("refseq_category", ""),
            "local_preferred": bool_text(local_preferred),
            "current_ncbi_reference": bool_text(current_reference),
            "protected_by_manifest": bool_text(protected_by_manifest),
            "protecting_manifests": ";".join(manifest_names),
            "candidate_for_cleanup": bool_text(candidate),
            "decision": decision,
            "decision_reason": reason,
            "accession_directory": str(accession_dir),
            "path_to_fna": fna,
            "path_to_gff": gff,
            "compressed_or_uncompressed": compression,
            "directory_size_bytes": size_bytes,
            "directory_size_gb": f"{size_bytes / (1024 ** 3):.6f}",
            "action": action,
            "trash_destination": "",
            "timestamp": assessment_time,
        }
        rows.append(row)

    return Assessment(
        rows=rows,
        context=context,
        ncbi=ncbi,
        local_existing_accessions=local_existing,
        candidate_accessions=candidates,
        fatal_errors=fatal_errors,
        warnings=warnings,
        timestamp=timestamp,
        taxon=taxon,
        accession_filter=filter_set,
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def summary_counts(assessment: Assessment) -> Dict[str, Any]:
    existing_rows = [r for r in assessment.rows if Path(str(r["accession_directory"])).is_dir()]
    ncbi_rows = [r for r in existing_rows if r["repository"] == "NCBI"]
    non_ncbi_rows = [r for r in existing_rows if r["repository"] != "NCBI"]
    candidates = [r for r in assessment.rows if r["candidate_for_cleanup"] == "True"]
    manual = [r for r in assessment.rows if r["decision"] == "manual_review"]
    reference = [r for r in assessment.rows if r["current_ncbi_reference"] == "True"]
    manifest = [r for r in assessment.rows if r["protected_by_manifest"] == "True"]
    preferred = [r for r in assessment.rows if r["local_preferred"] == "True"]
    recoverable = sum(int(r.get("directory_size_bytes", 0) or 0) for r in candidates)
    return {
        "local": len(existing_rows),
        "ncbi": len(ncbi_rows),
        "reference": len(reference),
        "manifest": len(manifest),
        "preferred": len(preferred),
        "non_ncbi": len(non_ncbi_rows),
        "candidates": len(candidates),
        "manual": len(manual),
        "recoverable_bytes": recoverable,
    }


def print_ncbi_commands(ncbi: NCBIDiscovery) -> None:
    if ncbi.reference_command:
        print("NCBI reference command:")
        print(f"  {command_string(ncbi.reference_command)}")
    if ncbi.taxon_command:
        print("NCBI taxon-membership command:")
        print(f"  {command_string(ncbi.taxon_command)}")


def print_assessment_summary(assessment: Assessment, mode: str) -> None:
    counts = summary_counts(assessment)
    print("\nGenome cleanup assessment")
    print("-------------------------")
    print(f"Taxon: {assessment.taxon}")
    print(f"Local accession directories examined: {counts['local']}")
    print(f"NCBI assemblies evaluated: {counts['ncbi']}")
    print(f"Current NCBI reference assemblies retained: {counts['reference']}")
    print(f"Protected by project manifests: {counts['manifest']}")
    print(f"Protected as locally preferred: {counts['preferred']}")
    print(f"Non-NCBI/GWH assemblies skipped: {counts['non_ncbi']}")
    print(f"Cleanup candidates: {counts['candidates']}")
    print(f"Manual review: {counts['manual']}")
    print(f"Estimated recoverable disk space: {counts['recoverable_bytes'] / (1024 ** 3):.3f} GB")

    if assessment.accession_filter is not None:
        print(f"Accession filter entries: {len(assessment.accession_filter)}")

    if assessment.warnings:
        print("\nWarnings:")
        for warning in assessment.warnings:
            print(f"  - {warning}")
    if assessment.fatal_errors:
        print("\nSAFETY BLOCKERS — no accession directory is authorized to move:")
        for error in assessment.fatal_errors:
            print(f"  - {error}")

    candidate_rows = [r for r in assessment.rows if r["candidate_for_cleanup"] == "True"]
    if candidate_rows:
        print("\nCandidate accessions:")
        for row in candidate_rows:
            print(f"  {row['accession']}\t{float(row['directory_size_gb']):.3f} GB")

    if mode == "dry-run":
        print("\nDRY RUN — no genome directories were moved and no active metadata or manifests were modified.")


def print_execution_plan(assessment: Assessment) -> None:
    counts = summary_counts(assessment)
    print("\nCleanup plan")
    print("------------")
    print(f"Taxon: {assessment.taxon}")
    print(f"Candidate accession directories: {counts['candidates']}")
    print(f"Estimated disk space recoverable: {counts['recoverable_bytes'] / (1024 ** 3):.3f} GB")
    print("\nProtected:")
    print(f"  NCBI references: {counts['reference']}")
    print(f"  Project manifests: {counts['manifest']}")
    print(f"  Local preferred assemblies: {counts['preferred']}")
    print(f"  Non-NCBI assemblies: {counts['non_ncbi']}")
    print(f"\n{counts['candidates']} directories will be MOVED TO TRASH.")
    print("No files will be permanently deleted.")


def write_run_log(path: Path, assessment: Assessment, mode: str) -> None:
    lines = [
        f"clean_genomes.py version: {SCRIPT_VERSION}",
        f"timestamp: {iso_now()}",
        f"mode: {mode}",
        f"genomes_dir: {assessment.context.genomes_dir}",
        f"taxon: {assessment.taxon}",
        f"reference_command: {command_string(assessment.ncbi.reference_command)}",
        f"taxon_command: {command_string(assessment.ncbi.taxon_command)}",
        f"reference_accessions_found: {len(assessment.ncbi.reference_accessions)}",
        f"taxon_accessions_found: {len(assessment.ncbi.taxon_accessions)}",
    ]
    for warning in assessment.warnings:
        lines.append(f"warning: {warning}")
    for error in assessment.fatal_errors:
        lines.append(f"safety_blocker: {error}")
    atomic_write_text(path, "\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Execute-mode moves and metadata synchronization
# ---------------------------------------------------------------------------


def safe_move_directory(src: Path, dst: Path) -> None:
    if not src.exists() or not src.is_dir():
        raise FileNotFoundError(f"Source accession directory does not exist: {src}")
    if src.is_symlink():
        raise SafetyAbort(f"Refusing to move symlinked accession directory: {src}")
    if dst.exists():
        raise FileExistsError(f"Trash destination already exists: {dst}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.rename(src, dst)
    except OSError as exc:
        if exc.errno == errno.EXDEV:
            # This should be uncommon because trash is inside the same repository,
            # but shutil.move is the safest portable fallback across filesystems.
            shutil.move(str(src), str(dst))
        else:
            raise


def backup_metadata_files(genomes_dir: Path, backup_dir: Path) -> List[Path]:
    records = genomes_dir / "records"
    names = ["genomes_metadata.csv", "species_assembly_index.csv", "genomes_catalog.jsonl"]
    backed_up: List[Path] = []
    backup_dir.mkdir(parents=True, exist_ok=True)
    for name in names:
        src = records / name
        if src.exists() and src.is_file():
            shutil.copy2(src, backup_dir / name)
            backed_up.append(src)
    return backed_up


def accession_column_for_file(path: Path) -> str:
    fields, _ = read_csv_rows(path)
    col = find_column(fields, ["accession", "assembly_accession", "genome_accession"])
    if not col:
        raise ValueError(f"No accession column found in active metadata file {path}")
    return col


def remove_accessions_from_csv(path: Path, remove: Set[str]) -> None:
    fields, rows = read_csv_rows(path)
    acc_col = find_column(fields, ["accession", "assembly_accession", "genome_accession"])
    if not acc_col:
        raise ValueError(f"No accession column found in {path}")
    kept = [row for row in rows if normalize_accession(row.get(acc_col, "")) not in remove]

    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    temp_path = Path(temp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            for row in kept:
                writer.writerow({field: row.get(field, "") for field in fields})
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_path, path)
    except Exception:
        try:
            temp_path.unlink(missing_ok=True)
        except Exception:
            pass
        raise


def json_record_accession(obj: Mapping[str, Any]) -> str:
    # Prefer only record-identifying fields.  Avoid recursively matching a
    # paired accession and accidentally deleting the wrong catalog record.
    for key in ("accession", "assembly_accession", "genome_accession"):
        value = obj.get(key)
        if isinstance(value, str) and value.strip():
            return normalize_accession(value)
    for container_key in ("assembly", "genome"):
        child = obj.get(container_key)
        if isinstance(child, Mapping):
            for key in ("accession", "assembly_accession", "genome_accession"):
                value = child.get(key)
                if isinstance(value, str) and value.strip():
                    return normalize_accession(value)
    return ""


def remove_accessions_from_jsonl(path: Path, remove: Set[str]) -> None:
    kept_lines: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for lineno, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Malformed JSONL in {path} at line {lineno}: {exc}") from exc
            if not isinstance(obj, Mapping):
                raise ValueError(f"Unexpected non-object JSONL record in {path} at line {lineno}")
            acc = json_record_accession(obj)
            if acc and acc in remove:
                continue
            kept_lines.append(line)
    atomic_write_text(path, ("\n".join(kept_lines) + "\n") if kept_lines else "")


def synchronize_active_metadata(genomes_dir: Path, moved_accessions: Set[str]) -> None:
    if not moved_accessions:
        return
    records = genomes_dir / "records"
    csv_files = [records / "genomes_metadata.csv", records / "species_assembly_index.csv"]
    jsonl = records / "genomes_catalog.jsonl"
    for path in csv_files:
        if path.exists():
            remove_accessions_from_csv(path, moved_accessions)
    if jsonl.exists():
        remove_accessions_from_jsonl(jsonl, moved_accessions)


def append_history(history_path: Path, entries: Sequence[Mapping[str, Any]]) -> None:
    if not entries:
        return
    history_path.parent.mkdir(parents=True, exist_ok=True)
    with history_path.open("a", encoding="utf-8") as handle:
        for entry in entries:
            handle.write(json.dumps(dict(entry), sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def execute_cleanup(assessment: Assessment) -> Tuple[Path, Set[str], Set[str]]:
    if assessment.fatal_errors:
        raise SafetyAbort("Safety prerequisites failed; execute mode is blocked and no genome directory was moved.")

    candidates = [r for r in assessment.rows if r["candidate_for_cleanup"] == "True"]
    trash_root = assessment.context.genomes_dir / "trash" / f"clean_genomes_{assessment.timestamp}"
    moved: Set[str] = set()
    failed: Set[str] = set()

    if not candidates:
        return trash_root, moved, failed

    backup_dir = trash_root / "metadata_backup"
    # Create the metadata snapshot before any move so recovery has an exact
    # pre-cleanup state even if a later step fails.
    backup_metadata_files(assessment.context.genomes_dir, backup_dir)

    row_by_accession = {r["accession"]: r for r in assessment.rows}
    history_entries: List[Dict[str, Any]] = []

    for row in candidates:
        accession = str(row["accession"])
        src = Path(str(row["accession_directory"]))
        dst = trash_root / accession
        try:
            safe_move_directory(src, dst)
            moved.add(accession)
            row["decision"] = "moved_to_trash"
            row["decision_reason"] = "candidate_directory_successfully_moved_to_trash"
            row["action"] = "moved_to_trash"
            row["trash_destination"] = str(dst)
            history_entries.append(
                {
                    "timestamp": iso_now(),
                    "accession": accession,
                    "action": "moved_to_trash",
                    "reason": "not_current_ncbi_reference_and_not_manifest_protected_and_not_locally_preferred",
                    "original_path": str(src),
                    "trash_path": str(dst),
                    "taxon_scope": assessment.taxon,
                }
            )
        except Exception as exc:
            failed.add(accession)
            row["candidate_for_cleanup"] = "False"
            row["decision"] = "move_failed"
            row["decision_reason"] = f"move_failed: {exc}"
            row["action"] = "move_failed"
            row["trash_destination"] = str(dst)
            print(f"ERROR: failed to move {accession}: {exc}", file=sys.stderr)

    # Only accessions actually moved are removed from active metadata.
    if moved:
        try:
            synchronize_active_metadata(assessment.context.genomes_dir, moved)
        except Exception as exc:
            # Data remain recoverable in trash and the pre-cleanup metadata
            # snapshot is available.  Surface the failure loudly and record it.
            for accession in moved:
                row = row_by_accession[accession]
                row["decision_reason"] += f"; metadata_sync_failed: {exc}"
            history_entries.append(
                {
                    "timestamp": iso_now(),
                    "action": "metadata_sync_failed",
                    "reason": str(exc),
                    "taxon_scope": assessment.taxon,
                    "moved_accessions": sorted(moved),
                }
            )
            append_history(assessment.context.genomes_dir / "records" / "cleanup" / "clean_genomes_history.jsonl", history_entries)
            raise SafetyAbort(
                f"Genome directories were moved to trash, but active metadata synchronization failed: {exc}. "
                f"Use the backups in {backup_dir} to restore metadata if needed."
            ) from exc

    append_history(assessment.context.genomes_dir / "records" / "cleanup" / "clean_genomes_history.jsonl", history_entries)
    return trash_root, moved, failed


# ---------------------------------------------------------------------------
# CLI helpers
# ---------------------------------------------------------------------------


def load_accession_filter(path: Path) -> Set[str]:
    accessions: Set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            token = normalize_accession(line.split()[0])
            if token:
                accessions.add(token)
    return accessions


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Conservatively identify excess local NCBI genome assemblies. "
            "Dry-run is the default; execute mode moves eligible accession directories to trash, never permanently deletes them."
        )
    )
    parser.add_argument("genomes_dir", nargs="?", help="Root of the local genome repository")
    parser.add_argument("--taxon", help="NCBI taxon name or taxonomy ID used to scope cleanup")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--dry-run", action="store_true", help="Assess and report only (default)")
    mode.add_argument("--execute", action="store_true", help="Move confirmed candidates to timestamped trash and update active metadata")
    parser.add_argument(
        "--accessions",
        type=Path,
        help="Optional file containing one accession per line; restricts evaluation without bypassing any protections",
    )
    parser.add_argument("--self-test", action="store_true", help="Run the built-in offline test suite and exit")
    parser.add_argument("--version", action="version", version=f"%(prog)s {SCRIPT_VERSION}")
    return parser


# ---------------------------------------------------------------------------
# Built-in offline tests
# ---------------------------------------------------------------------------


def _write_csv(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _fake_ncbi(
    reference: Iterable[str],
    taxon: Iterable[str],
    metadata: Optional[Mapping[str, Mapping[str, str]]] = None,
    ok: bool = True,
) -> NCBIDiscovery:
    ref = {normalize_accession(x) for x in reference}
    tax = {normalize_accession(x) for x in taxon}
    md: Dict[str, Dict[str, str]] = {}
    for acc in tax | ref:
        md[acc] = {
            "organism_name": "Testus species",
            "species_key": "testus_species",
            "assembly_name": f"assembly_{acc}",
            "assembly_level": "chromosome",
            "source_database": "refseq" if acc.startswith("GCF_") else "genbank",
            "refseq_category": "reference genome" if acc in ref else "",
        }
    if metadata:
        for acc, row in metadata.items():
            md.setdefault(normalize_accession(acc), {}).update(dict(row))
    return NCBIDiscovery(
        reference_accessions=ref,
        taxon_accessions=tax,
        metadata_by_accession=md,
        reference_command=["datasets", "summary", "genome", "taxon", "test", "--reference", "--as-json-lines"],
        taxon_command=["datasets", "summary", "genome", "taxon", "test", "--assembly-version", "all", "--as-json-lines"],
        reference_query_ok=ok,
        taxon_query_ok=ok,
        errors=[] if ok else ["mock query failure"],
    )


class CleanGenomesSelfTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name) / "genomes"
        (self.root / "records" / "project_manifests").mkdir(parents=True)
        (self.root / "records" / "cleanup").mkdir(parents=True)

    def tearDown(self) -> None:
        self.temp.cleanup()

    def add_dir(self, accession: str, compressed: bool = False, payload_size: int = 16) -> Path:
        d = self.root / accession
        d.mkdir(parents=True, exist_ok=True)
        ext = ".fna.gz" if compressed else ".fna"
        (d / f"{accession}_genomic{ext}").write_bytes(b"A" * payload_size)
        return d

    def set_metadata(self, accessions: Sequence[str], preferred: Iterable[str] = ()) -> None:
        preferred_set = {normalize_accession(x) for x in preferred}
        rows = []
        for acc in accessions:
            rows.append(
                {
                    "accession": acc,
                    "organism_name": "Testus species",
                    "is_preferred_for_species": "True" if acc in preferred_set else "False",
                    "status": "preferred" if acc in preferred_set else "available",
                }
            )
        _write_csv(
            self.root / "records" / "genomes_metadata.csv",
            ["accession", "organism_name", "is_preferred_for_species", "status"],
            rows,
        )
        _write_csv(
            self.root / "records" / "species_assembly_index.csv",
            ["accession", "status"],
            [{"accession": acc, "status": "preferred" if acc in preferred_set else "available"} for acc in accessions],
        )
        with (self.root / "records" / "genomes_catalog.jsonl").open("w", encoding="utf-8") as handle:
            for acc in accessions:
                handle.write(json.dumps({"accession": acc, "active": True}) + "\n")

    def add_manifest(self, name: str, accessions: Sequence[str]) -> Path:
        path = self.root / "records" / "project_manifests" / name
        _write_csv(path, ["accession"], [{"accession": x} for x in accessions])
        return path

    def row(self, assessment: Assessment, accession: str) -> Dict[str, Any]:
        return next(r for r in assessment.rows if r["accession"] == accession)

    def test_01_reference_is_protected(self) -> None:
        acc = "GCF_000000001.1"
        self.add_dir(acc)
        self.set_metadata([acc])
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([acc], [acc]))
        self.assertEqual(self.row(a, acc)["decision"], "keep_reference")

    def test_02_nonreference_becomes_candidate(self) -> None:
        ref = "GCF_000000001.1"
        acc = "GCA_000000002.1"
        self.add_dir(ref)
        self.add_dir(acc)
        self.set_metadata([ref, acc])
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([ref], [ref, acc]))
        self.assertEqual(self.row(a, acc)["decision"], "cleanup_candidate")

    def test_03_manifest_protects_nonreference(self) -> None:
        ref = "GCF_000000001.1"
        acc = "GCA_000000003.1"
        self.add_dir(ref); self.add_dir(acc)
        self.set_metadata([ref, acc])
        self.add_manifest("a_manifest.csv", [acc])
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([ref], [ref, acc]))
        self.assertEqual(self.row(a, acc)["decision"], "keep_manifest")

    def test_04_multiple_manifests_protect_without_duplicate_count(self) -> None:
        ref = "GCF_000000001.1"
        acc = "GCA_000000004.1"
        self.add_dir(ref); self.add_dir(acc)
        self.set_metadata([ref, acc])
        self.add_manifest("a_manifest.csv", [acc, acc])
        self.add_manifest("b_manifest.csv", [acc])
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([ref], [ref, acc]))
        row = self.row(a, acc)
        self.assertEqual(row["decision"], "keep_manifest")
        self.assertEqual(len(row["protecting_manifests"].split(";")), 2)

    def test_05_local_preferred_nonreference_is_protected(self) -> None:
        ref = "GCF_000000001.1"
        acc = "GCA_000000005.1"
        self.add_dir(ref); self.add_dir(acc)
        self.set_metadata([ref, acc], preferred=[acc])
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([ref], [ref, acc]))
        self.assertFalse(self.row(a, acc)["candidate_for_cleanup"] == "True")
        self.assertEqual(self.row(a, acc)["decision"], "manual_review")

    def test_06_gwh_never_automatically_removed(self) -> None:
        ref = "GCF_000000001.1"
        gwh = "GWHABC00000000.1"
        self.add_dir(ref); self.add_dir(gwh)
        self.set_metadata([ref, gwh])
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([ref], [ref]))
        self.assertEqual(self.row(a, gwh)["decision"], "keep_non_ncbi")

    def test_07_unknown_local_source_is_protected(self) -> None:
        ref = "GCF_000000001.1"
        local = "LOCAL_ASSEMBLY_1"
        self.add_dir(ref); self.add_dir(local)
        self.set_metadata([ref, local])
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([ref], [ref]))
        self.assertEqual(self.row(a, local)["decision"], "keep_non_ncbi")

    def test_08_compressed_is_evaluated_without_decompression(self) -> None:
        ref = "GCF_000000001.1"
        acc = "GCA_000000008.1"
        self.add_dir(ref); d = self.add_dir(acc, compressed=True)
        self.set_metadata([ref, acc])
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([ref], [ref, acc]))
        self.assertEqual(self.row(a, acc)["compressed_or_uncompressed"], "compressed")
        self.assertTrue(any(p.name.endswith(".gz") for p in d.iterdir()))

    def test_09_dry_assessment_does_not_move_or_edit_active_files(self) -> None:
        ref = "GCF_000000001.1"
        acc = "GCA_000000009.1"
        self.add_dir(ref); self.add_dir(acc)
        self.set_metadata([ref, acc])
        metadata_before = (self.root / "records" / "genomes_metadata.csv").read_bytes()
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([ref], [ref, acc]))
        self.assertTrue((self.root / acc).exists())
        self.assertEqual(metadata_before, (self.root / "records" / "genomes_metadata.csv").read_bytes())
        self.assertEqual(self.row(a, acc)["decision"], "cleanup_candidate")

    def test_10_execute_moves_complete_directory(self) -> None:
        ref = "GCF_000000001.1"
        acc = "GCA_000000010.1"
        self.add_dir(ref); d = self.add_dir(acc)
        (d / "source_metadata.json").write_text("{}", encoding="utf-8")
        self.set_metadata([ref, acc])
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([ref], [ref, acc]), timestamp="20260101_000000")
        trash, moved, failed = execute_cleanup(a)
        self.assertIn(acc, moved)
        self.assertFalse(failed)
        self.assertFalse((self.root / acc).exists())
        self.assertTrue((trash / acc / "source_metadata.json").exists())

    def test_11_metadata_updates_only_after_successful_move(self) -> None:
        ref = "GCF_000000001.1"
        good = "GCA_000000011.1"
        bad = "GCA_000000012.1"
        self.add_dir(ref); self.add_dir(good); self.add_dir(bad)
        self.set_metadata([ref, good, bad])
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([ref], [ref, good, bad]), timestamp="20260101_000001")
        original_safe_move = safe_move_directory
        def fake_move(src: Path, dst: Path) -> None:
            if src.name == bad:
                raise OSError("mock move failure")
            original_safe_move(src, dst)
        with mock.patch(__name__ + ".safe_move_directory", side_effect=fake_move):
            _, moved, failed = execute_cleanup(a)
        self.assertIn(good, moved)
        self.assertIn(bad, failed)
        _, rows = read_csv_rows(self.root / "records" / "genomes_metadata.csv")
        remaining = {normalize_accession(r["accession"]) for r in rows}
        self.assertNotIn(good, remaining)
        self.assertIn(bad, remaining)

    def test_12_metadata_backups_created(self) -> None:
        ref = "GCF_000000001.1"
        acc = "GCA_000000013.1"
        self.add_dir(ref); self.add_dir(acc)
        self.set_metadata([ref, acc])
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([ref], [ref, acc]), timestamp="20260101_000002")
        trash, _, _ = execute_cleanup(a)
        self.assertTrue((trash / "metadata_backup" / "genomes_metadata.csv").exists())
        self.assertTrue((trash / "metadata_backup" / "species_assembly_index.csv").exists())
        self.assertTrue((trash / "metadata_backup" / "genomes_catalog.jsonl").exists())

    def test_13_project_manifests_never_modified(self) -> None:
        ref = "GCF_000000001.1"
        acc = "GCA_000000014.1"
        self.add_dir(ref); self.add_dir(acc)
        self.set_metadata([ref, acc])
        manifest = self.add_manifest("keep_manifest.csv", [ref])
        before = manifest.read_bytes()
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([ref], [ref, acc]), timestamp="20260101_000003")
        execute_cleanup(a)
        self.assertEqual(before, manifest.read_bytes())

    def test_14_one_failed_move_does_not_stop_later_candidates(self) -> None:
        ref = "GCF_000000001.1"
        first = "GCA_000000015.1"
        second = "GCA_000000016.1"
        self.add_dir(ref); self.add_dir(first); self.add_dir(second)
        self.set_metadata([ref, first, second])
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([ref], [ref, first, second]), timestamp="20260101_000004")
        original_safe_move = safe_move_directory
        def fake_move(src: Path, dst: Path) -> None:
            if src.name == first:
                raise OSError("mock failure")
            original_safe_move(src, dst)
        with mock.patch(__name__ + ".safe_move_directory", side_effect=fake_move):
            _, moved, failed = execute_cleanup(a)
        self.assertIn(first, failed)
        self.assertIn(second, moved)
        self.assertFalse((self.root / second).exists())

    def test_15_ncbi_query_failure_results_in_no_candidate(self) -> None:
        acc = "GCA_000000017.1"
        self.add_dir(acc)
        self.set_metadata([acc])
        fake = _fake_ncbi([], [], ok=False)
        a = assess_repository(self.root, "test", ncbi_data=fake)
        self.assertTrue(a.fatal_errors)
        self.assertFalse(a.candidate_accessions)

    def test_16_empty_reference_query_is_suspicious(self) -> None:
        acc = "GCA_000000018.1"
        self.add_dir(acc)
        self.set_metadata([acc])
        fake = _fake_ncbi([], [acc], ok=True)
        fake.reference_query_ok = False
        a = assess_repository(self.root, "test", ncbi_data=fake)
        self.assertTrue(a.fatal_errors)
        self.assertFalse(a.candidate_accessions)

    def test_17_duplicate_manifest_accessions_handled(self) -> None:
        ref = "GCF_000000001.1"
        acc = "GCA_000000019.1"
        self.add_dir(ref); self.add_dir(acc)
        self.set_metadata([ref, acc])
        self.add_manifest("dup_manifest.csv", [acc, acc, acc])
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([ref], [ref, acc]))
        row = self.row(a, acc)
        self.assertEqual(row["decision"], "keep_manifest")
        self.assertEqual(row["protecting_manifests"].count("dup_manifest.csv"), 1)

    def test_18_accession_filter_restricts_without_bypassing_protection(self) -> None:
        ref = "GCF_000000001.1"
        a1 = "GCA_000000020.1"
        a2 = "GCA_000000021.1"
        self.add_dir(ref); self.add_dir(a1); self.add_dir(a2)
        self.set_metadata([ref, a1, a2])
        self.add_manifest("protect.csv", [a1])
        a = assess_repository(self.root, "test", accession_filter={a1}, ncbi_data=_fake_ncbi([ref], [ref, a1, a2]))
        self.assertEqual(self.row(a, a1)["decision"], "keep_manifest")
        self.assertEqual(self.row(a, a2)["decision"], "keep_outside_accession_filter")

    def test_19_audit_reason_records_keep_and_remove(self) -> None:
        ref = "GCF_000000001.1"
        acc = "GCA_000000022.1"
        self.add_dir(ref); self.add_dir(acc)
        self.set_metadata([ref, acc])
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([ref], [ref, acc]))
        out = self.root / "records" / "cleanup" / "audit.tsv"
        write_tsv(out, a.rows)
        text = out.read_text(encoding="utf-8")
        self.assertIn("current_ncbi_reference", text)
        self.assertIn("cleanup_candidate", text)
        self.assertIn("current_ncbi_reference", text)

    def test_20_directory_size_for_compressed_and_uncompressed(self) -> None:
        ref = "GCF_000000001.1"
        a1 = "GCA_000000023.1"
        a2 = "GCA_000000024.1"
        self.add_dir(ref, payload_size=1)
        self.add_dir(a1, compressed=False, payload_size=101)
        self.add_dir(a2, compressed=True, payload_size=202)
        self.set_metadata([ref, a1, a2])
        a = assess_repository(self.root, "test", ncbi_data=_fake_ncbi([ref], [ref, a1, a2]))
        self.assertGreaterEqual(int(self.row(a, a1)["directory_size_bytes"]), 101)
        self.assertGreaterEqual(int(self.row(a, a2)["directory_size_bytes"]), 202)
        self.assertEqual(self.row(a, a1)["compressed_or_uncompressed"], "uncompressed")
        self.assertEqual(self.row(a, a2)["compressed_or_uncompressed"], "compressed")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def run_self_tests() -> int:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(CleanGenomesSelfTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.self_test:
        return run_self_tests()

    if not args.genomes_dir:
        parser.error("genomes_dir is required unless --self-test is used")
    if not args.taxon:
        parser.error("--taxon is required unless --self-test is used")

    genomes_dir = Path(args.genomes_dir).expanduser().resolve()
    if not genomes_dir.exists() or not genomes_dir.is_dir():
        parser.error(f"genomes_dir does not exist or is not a directory: {genomes_dir}")

    mode = "execute" if args.execute else "dry-run"
    accession_filter: Optional[Set[str]] = None
    if args.accessions:
        try:
            accession_filter = load_accession_filter(args.accessions.expanduser().resolve())
        except Exception as exc:
            parser.error(f"could not read --accessions file: {exc}")
        if not accession_filter:
            parser.error("--accessions file did not contain any accessions")

    timestamp = now_timestamp()

    print(f"clean_genomes.py {SCRIPT_VERSION}")
    print(f"Repository: {genomes_dir}")
    print(f"Mode: {mode}")
    print(f"Taxon: {args.taxon}")

    assessment = assess_repository(
        genomes_dir=genomes_dir,
        taxon=args.taxon,
        accession_filter=accession_filter,
        timestamp=timestamp,
    )

    print_ncbi_commands(assessment.ncbi)
    print_assessment_summary(assessment, mode)

    cleanup_dir = genomes_dir / "records" / "cleanup"
    audit_path = cleanup_dir / f"clean_genomes_{timestamp}.tsv"
    log_path = cleanup_dir / f"clean_genomes_{timestamp}.log"

    # Audit files are the only intentional dry-run output files. Active genome
    # directories, metadata tables, catalogs, and project manifests are never
    # modified in dry-run mode.
    try:
        write_tsv(audit_path, assessment.rows)
        write_run_log(log_path, assessment, mode)
        print(f"\nAudit report: {audit_path}")
        print(f"Run log:      {log_path}")
    except Exception as exc:
        print(f"ERROR: could not write cleanup audit files: {exc}", file=sys.stderr)
        return 2

    if mode == "dry-run":
        return 2 if assessment.fatal_errors else 0

    print_execution_plan(assessment)
    if assessment.fatal_errors:
        print("\nEXECUTION BLOCKED — safety prerequisites failed. No genome directory was moved.", file=sys.stderr)
        return 2

    try:
        trash_root, moved, failed = execute_cleanup(assessment)
    except SafetyAbort as exc:
        # The assessment rows may have been updated before a metadata-sync
        # failure, so refresh the audit report before exiting.
        try:
            write_tsv(audit_path, assessment.rows)
            write_run_log(log_path, assessment, mode)
        except Exception:
            pass
        print(f"\nSAFETY ABORT: {exc}", file=sys.stderr)
        return 2

    # Refresh audit with execute-mode outcomes.
    write_tsv(audit_path, assessment.rows)
    write_run_log(log_path, assessment, mode)

    print("\nCleanup execution complete")
    print("--------------------------")
    print(f"Moved to trash: {len(moved)}")
    print(f"Move failures:  {len(failed)}")
    if moved:
        print(f"Trash archive:  {trash_root}")
        print(f"Metadata backup: {trash_root / 'metadata_backup'}")
    print(f"Audit report:   {audit_path}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
