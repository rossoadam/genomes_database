#!/usr/bin/env python3
"""Manage remote and manually downloaded genome assemblies.

This script preserves the project-manifest and metadata-management behavior of
``genomes.py`` while introducing a source-adapter architecture. Repository-
specific code is isolated in ``NCBISourceAdapter`` and ``GWHSourceAdapter``;
``LocalSourceAdapter`` imports files obtained through any repository browser.
Downstream metadata, manifest, retry, and rebuild operations consume normalized
``AssemblyRecord`` objects.

Version 3.1.0 file-recovery behavior
------------------------------------
Existing assembly directories are reconciled file-by-file. A valid canonical
genome FASTA or GFF is preserved, while a repository-provided missing file is
recovered independently. Missing optional GFF annotations do not make an
otherwise valid genome assembly a failed download.

Python requirements
-------------------
Python 3.9 or newer. GWH support uses only the Python standard library. NCBI
operations require the NCBI Datasets CLI. Taxonomic lineage enrichment is
optional and uses TaxonKit when it is available in PATH.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import gzip
import json
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any, Callable, Iterable, Optional


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SCRIPT_VERSION = "3.1.0"
SOURCE_METADATA_SCHEMA_VERSION = 1
USER_AGENT = f"GenomeManagerHybrid/{SCRIPT_VERSION} (+genome-repository-manager)"
DEFAULT_CONNECT_TIMEOUT = 30.0
DEFAULT_READ_TIMEOUT = 300.0
SUPPORTED_SOURCES = ("ncbi", "gwh", "local")
SOURCE_PREFIXES: dict[str, tuple[str, ...]] = {
    "ncbi": ("GCA_", "GCF_"),
    "gwh": ("GWH",),
}


# ---------------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------------


def now_iso() -> str:
    """Return a local ISO-8601 timestamp with second precision."""

    return datetime.datetime.now().isoformat(timespec="seconds")


def require_genomes_dir(path: Path) -> Path:
    """Validate and resolve the project genome directory."""

    resolved = path.expanduser().resolve()
    if resolved.name != "genomes":
        raise ValueError(
            f"Expected a directory ending with 'genomes', got: {resolved}"
        )
    return resolved


def safe_tag(text: str) -> str:
    """Convert a search term to a conservative filename component."""

    cleaned = text.strip().replace("/", "_").replace("\\", "_")
    cleaned = "_".join(cleaned.split())
    return cleaned or "query"


def is_nonempty_file(path: Path) -> bool:
    """Return True only when *path* is a regular file containing data."""

    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def root_acc(accession: str) -> str:
    """Remove an accession version suffix without assuming a repository."""

    return accession.split(".", 1)[0] if accession else accession


def build_species_key(organism_name: str) -> str:
    """Build the existing underscore-delimited species key."""

    return "_".join(organism_name.strip().split())


def detect_source(term: str) -> str:
    """Detect the repository associated with an accession or search term.

    GCA_/GCF_ accessions use NCBI, GWH accessions use Genome Warehouse, and
    taxonomic names or other non-accession search terms default to NCBI.
    """

    normalized = term.strip().upper()
    for source, prefixes in SOURCE_PREFIXES.items():
        if normalized.startswith(prefixes):
            return source
    return "ncbi"


def is_accession_term(term: str) -> bool:
    """Return True for a supported NCBI or GWH assembly accession."""

    normalized = term.strip().upper()
    return any(
        normalized.startswith(prefix)
        for prefixes in SOURCE_PREFIXES.values()
        for prefix in prefixes
    )


def resolve_source(term: str, requested_source: str = "auto") -> str:
    """Resolve an explicit or automatically detected repository choice."""

    if requested_source == "auto":
        return detect_source(term)
    if requested_source not in SUPPORTED_SOURCES:
        raise ValueError(
            f"Unsupported source '{requested_source}'. "
            f"Choose from: auto, {', '.join(SUPPORTED_SOURCES)}"
        )
    return requested_source


def normalize_repository_name(value: str) -> str:
    """Normalize a provenance repository label used by local imports."""

    normalized = value.strip().lower().replace("_", "-")
    aliases = {
        "cncb-gwh": "gwh",
        "cncb": "gwh",
        "genome-warehouse": "gwh",
        "genome warehouse": "gwh",
        "refseq": "ncbi",
        "genbank": "ncbi",
    }
    return aliases.get(normalized, normalized)


def source_database_for_repository(repository: str, accession: str = "") -> str:
    """Return a stable source-database label for normalized metadata."""

    normalized = normalize_repository_name(repository)
    if normalized == "gwh":
        return "CNCB-GWH"
    if normalized == "ncbi":
        return "RefSeq" if accession.upper().startswith("GCF_") else "GenBank"
    return repository.strip() or "Local import"


def default_source_record_url(repository: str, accession: str) -> str:
    """Build a provenance URL when a local import does not supply one."""

    normalized = normalize_repository_name(repository)
    quoted = urllib.parse.quote(accession, safe="")
    if normalized == "gwh":
        return f"https://ngdc.cncb.ac.cn/gwh/api/public/assembly/{quoted}"
    if normalized == "ncbi":
        return f"https://www.ncbi.nlm.nih.gov/datasets/genome/{quoted}/"
    return ""


def normalize_assembly_level(value: Any) -> str:
    """Normalize repository-specific assembly-level descriptions."""

    text = str(value or "").strip().lower()
    if "chromosome" in text:
        return "chromosome"
    if "complete" in text:
        return "complete"
    if "scaffold" in text:
        return "scaffold"
    if "contig" in text:
        return "contig"
    return text


def normalize_date(value: Any) -> str:
    """Normalize a repository date/time value to YYYY-MM-DD when possible."""

    text = str(value or "").strip()
    if not text:
        return ""
    if len(text) >= 10 and text[4:5] == "-" and text[7:8] == "-":
        return text[:10]
    return text


def parse_bool(value: Any, default: bool = False) -> bool:
    """Interpret common CSV and JSON boolean representations."""

    if isinstance(value, bool):
        return value
    if value is None:
        return default
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "y", "preferred", "current"}:
        return True
    if text in {"false", "0", "no", "n", "available", "superseded"}:
        return False
    return default


def stringify_value(value: Any) -> str:
    """Convert scalar/list API values to a stable metadata string."""

    if value is None:
        return ""
    if isinstance(value, (list, tuple, set)):
        return ";".join(str(item).strip() for item in value if str(item).strip())
    if isinstance(value, dict):
        for key in ("accession", "name", "value", "title"):
            if key in value and value[key] not in (None, ""):
                return str(value[key]).strip()
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value).strip()


def nested_get(mapping: dict[str, Any], dotted_path: str, default: Any = "") -> Any:
    """Safely read a dotted path from dictionaries and indexed lists."""

    current: Any = mapping
    for component in dotted_path.split("."):
        if isinstance(current, dict):
            if component not in current:
                return default
            current = current[component]
            continue
        if isinstance(current, list) and component.isdigit():
            index = int(component)
            if index >= len(current):
                return default
            current = current[index]
            continue
        return default
    return current


def first_nonempty(mapping: dict[str, Any], *paths: str, default: str = "") -> str:
    """Return the first nonempty value among several nested field paths."""

    for path in paths:
        value = nested_get(mapping, path, None)
        text = stringify_value(value)
        if text:
            return text
    return default


def extract_gwh_assembly_payload(payload: dict[str, Any]) -> Optional[dict[str, Any]]:
    """Locate a GWH assembly object inside a possibly wrapped API response.

    The public GWH documentation shows assembly fields at the JSON top level,
    but some deployments return envelopes such as ``data``, ``result``, or a
    one-element list. This helper searches those structures without tying the
    adapter to one response wrapper.
    """

    queue: list[tuple[Any, int]] = [(payload, 0)]
    visited: set[int] = set()
    preferred_keys = (
        "data",
        "result",
        "results",
        "assembly",
        "assemblies",
        "content",
        "payload",
        "response",
        "record",
        "records",
    )

    while queue:
        current, depth = queue.pop(0)
        if depth > 8:
            continue
        identity = id(current)
        if identity in visited:
            continue
        visited.add(identity)

        if isinstance(current, dict):
            accession = first_nonempty(
                current,
                "assemblyAccession",
                "assembly_accession",
                "assembly accession",
            )
            if accession.upper().startswith("GWH"):
                return current

            generic_accession = first_nonempty(current, "accession")
            if generic_accession.upper().startswith("GWH") and any(
                key in current
                for key in (
                    "ftpPathDna",
                    "ftp_path_dna",
                    "assemblyLevel",
                    "assembly_level",
                    "organism",
                )
            ):
                copied = dict(current)
                copied.setdefault("assemblyAccession", generic_accession)
                return copied

            for key in preferred_keys:
                if key in current:
                    queue.append((current[key], depth + 1))
            for key, value in current.items():
                if key not in preferred_keys and isinstance(value, (dict, list)):
                    queue.append((value, depth + 1))

        elif isinstance(current, list):
            for value in current:
                if isinstance(value, (dict, list)):
                    queue.append((value, depth + 1))

    return None


def summarize_gwh_response(payload: dict[str, Any], limit: int = 800) -> str:
    """Return a compact diagnostic summary of an unusable GWH response."""

    details: list[str] = []
    for key in (
        "message",
        "msg",
        "error",
        "errorMessage",
        "detail",
        "status",
        "code",
        "success",
    ):
        value = payload.get(key)
        if value not in (None, "", [], {}):
            details.append(f"{key}={stringify_value(value)!r}")

    keys = ", ".join(sorted(str(key) for key in payload.keys())) or "<none>"
    try:
        preview = json.dumps(payload, ensure_ascii=False, sort_keys=True)
    except (TypeError, ValueError):
        preview = repr(payload)
    if len(preview) > limit:
        preview = preview[:limit] + "..."

    prefix = "; ".join(details)
    if prefix:
        prefix += "; "
    return f"{prefix}top-level keys=[{keys}]; response={preview}"


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    """Read a CSV file or return an empty list when it does not exist."""

    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv_rows(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    """Atomically write a CSV table while ignoring unexpected legacy keys."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
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
    """Read a JSON object from disk and validate its top-level type."""

    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return payload


def write_json(path: Path, payload: dict[str, Any]) -> None:
    """Atomically write an indented UTF-8 JSON object."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def iter_json_lines(path: Path) -> Iterable[dict[str, Any]]:
    """Yield JSON objects from a JSON Lines file."""

    if not path.exists():
        return
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip():
                continue
            try:
                row = json.loads(raw_line)
            except json.JSONDecodeError as exc:
                raise ValueError(
                    f"Invalid JSON on line {line_number} of {path}: {exc}"
                ) from exc
            if isinstance(row, dict):
                yield row


def append_json_line(path: Path, payload: dict[str, Any]) -> None:
    """Append one compact JSON object to a JSON Lines file."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, ensure_ascii=False, sort_keys=True))
        handle.write("\n")


def canonical_output_paths(genomes_dir: Path, accession: str) -> dict[str, Path]:
    """Return repository-independent local paths for an assembly."""

    accession_root = genomes_dir / accession
    return {
        "root": accession_root,
        "fna": accession_root / f"{accession}_genomic.fna",
        "gff": accession_root / f"{accession}_genomic.gff",
        "lift_gff": accession_root / f"{accession}_lift.gff",
        "source_metadata": accession_root / "source_metadata.json",
    }


def _set_response_read_timeout(response: Any, read_timeout: float) -> None:
    """Best-effort read-timeout configuration for urllib responses.

    ``urlopen(timeout=...)`` applies a socket timeout to connection operations.
    CPython does not expose a public separate read-timeout argument, so this
    helper locates the underlying socket when possible and changes its timeout
    before the response body is consumed.
    """

    candidate_paths = (
        ("fp", "raw", "_sock"),
        ("fp", "_sock"),
        ("raw", "_sock"),
    )
    for path in candidate_paths:
        current = response
        try:
            for attribute in path:
                current = getattr(current, attribute)
            current.settimeout(read_timeout)
            return
        except (AttributeError, OSError):
            continue


def fetch_json(
    url: str,
    *,
    connect_timeout: float = DEFAULT_CONNECT_TIMEOUT,
    read_timeout: float = DEFAULT_READ_TIMEOUT,
    user_agent: str = USER_AGENT,
) -> dict[str, Any]:
    """Retrieve and validate a JSON object over HTTPS using the standard library."""

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": user_agent,
            "Accept": "application/json",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=connect_timeout) as response:
            _set_response_read_timeout(response, read_timeout)
            status = getattr(response, "status", 200)
            if status < 200 or status >= 300:
                raise RuntimeError(f"HTTP {status} while requesting {url}")
            raw = response.read()
    except urllib.error.HTTPError as exc:
        detail = ""
        try:
            detail = exc.read(500).decode("utf-8", errors="replace").strip()
        except Exception:
            detail = ""
        suffix = f": {detail}" if detail else ""
        raise RuntimeError(
            f"HTTP {exc.code} while requesting {url}{suffix}"
        ) from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Could not retrieve {url}: {exc.reason}") from exc
    except (TimeoutError, socket.timeout) as exc:
        raise RuntimeError(f"Timed out while requesting {url}") from exc

    try:
        payload = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Response from {url} was not valid UTF-8 JSON") from exc
    if not isinstance(payload, dict):
        raise RuntimeError(f"Expected a JSON object from {url}")
    return payload


def fetch_gwh_assembly(
    accession: str,
    *,
    fetcher: Callable[..., dict[str, Any]] = fetch_json,
    connect_timeout: float = DEFAULT_CONNECT_TIMEOUT,
    read_timeout: float = DEFAULT_READ_TIMEOUT,
) -> tuple[dict[str, Any], str]:
    """Fetch a GWH assembly, falling back to the unversioned root when needed.

    The exact user-requested accession is preserved later during normalization,
    even when the API request succeeds only for the unversioned root.
    """

    base = "https://ngdc.cncb.ac.cn/gwh/api/public/assembly"
    attempts = [accession]
    unversioned = root_acc(accession)
    if unversioned and unversioned != accession:
        attempts.append(unversioned)

    errors: list[str] = []
    for candidate in attempts:
        url = f"{base}/{urllib.parse.quote(candidate, safe='')}"
        try:
            payload = fetcher(
                url,
                connect_timeout=connect_timeout,
                read_timeout=read_timeout,
                user_agent=USER_AGENT,
            )
        except RuntimeError as exc:
            errors.append(f"{candidate}: {exc}")
            continue

        assembly_payload = extract_gwh_assembly_payload(payload)
        if assembly_payload is None:
            errors.append(
                f"{candidate}: response did not contain a usable GWH assembly "
                f"record; {summarize_gwh_response(payload)}"
            )
            continue

        # Preserve a nonstandard API envelope for troubleshooting and cache
        # rebuilding while keeping normalized assembly fields at the top level.
        if assembly_payload is not payload:
            assembly_payload = dict(assembly_payload)
            assembly_payload["_gwh_api_envelope"] = payload
        return assembly_payload, url

    joined = "; ".join(errors) if errors else "no response"
    raise RuntimeError(
        f"GWH metadata retrieval failed for {accession}. Attempts: {joined}"
    )


def ensure_absolute_url(url: str, base_url: str = "https://ngdc.cncb.ac.cn") -> str:
    """Normalize absolute, protocol-relative, and repository-relative URLs."""

    value = url.strip()
    if not value:
        return ""
    if value.startswith("//"):
        return "https:" + value
    return urllib.parse.urljoin(base_url, value)


def download_and_decompress(
    url: str,
    destination: Path,
    *,
    connect_timeout: float = DEFAULT_CONNECT_TIMEOUT,
    read_timeout: float = DEFAULT_READ_TIMEOUT,
    user_agent: str = USER_AGENT,
) -> None:
    """Stream a remote file to a temporary path and atomically install it.

    Files whose URL path ends in ``.gz`` are transparently decompressed. The
    entire assembly is never loaded into memory. Partial files are removed on
    failure, and zero-byte downloads are rejected.
    """

    if not url:
        raise ValueError("Cannot download an empty URL")

    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    partial.unlink(missing_ok=True)
    request = urllib.request.Request(
        url,
        headers={"User-Agent": user_agent, "Accept": "*/*"},
    )

    try:
        with urllib.request.urlopen(request, timeout=connect_timeout) as response:
            _set_response_read_timeout(response, read_timeout)
            status = getattr(response, "status", 200)
            if status < 200 or status >= 300:
                raise RuntimeError(f"HTTP {status} while downloading {url}")

            is_gzip = urllib.parse.urlparse(url).path.lower().endswith(".gz")
            source_stream: Any = response
            gzip_stream: Optional[gzip.GzipFile] = None
            if is_gzip:
                gzip_stream = gzip.GzipFile(fileobj=response, mode="rb")
                source_stream = gzip_stream

            try:
                with partial.open("wb") as output:
                    shutil.copyfileobj(source_stream, output, length=1024 * 1024)
            finally:
                if gzip_stream is not None:
                    gzip_stream.close()

        if not partial.exists() or partial.stat().st_size == 0:
            raise RuntimeError(f"Downloaded file was empty: {url}")
        os.replace(partial, destination)
    except urllib.error.HTTPError as exc:
        partial.unlink(missing_ok=True)
        raise RuntimeError(f"HTTP {exc.code} while downloading {url}") from exc
    except urllib.error.URLError as exc:
        partial.unlink(missing_ok=True)
        raise RuntimeError(f"Could not download {url}: {exc.reason}") from exc
    except (TimeoutError, socket.timeout) as exc:
        partial.unlink(missing_ok=True)
        raise RuntimeError(f"Timed out while downloading {url}") from exc
    except Exception:
        partial.unlink(missing_ok=True)
        raise


def copy_or_decompress_local(source: Path, destination: Path) -> None:
    """Copy a local file to a canonical destination, decompressing gzip files."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    partial.unlink(missing_ok=True)
    try:
        if source.suffix.lower() == ".gz":
            with gzip.open(source, "rb") as input_handle, partial.open("wb") as output:
                shutil.copyfileobj(input_handle, output, length=1024 * 1024)
        else:
            with source.open("rb") as input_handle, partial.open("wb") as output:
                shutil.copyfileobj(input_handle, output, length=1024 * 1024)
        if partial.stat().st_size == 0:
            raise RuntimeError(f"Source file was empty: {source}")
        os.replace(partial, destination)
    except Exception:
        partial.unlink(missing_ok=True)
        raise


# ---------------------------------------------------------------------------
# Normalized data structures
# ---------------------------------------------------------------------------


@dataclass
class AssemblyRecord:
    """Repository-independent assembly metadata written to genomes_metadata.csv."""

    accession: str = ""
    accession_root: str = ""

    repository: str = ""
    source_database: str = ""
    source_record_url: str = ""
    source_query: str = ""

    organism_name: str = ""
    taxon_id: str = ""
    species_key: str = ""

    assembly_name: str = ""
    assembly_level: str = ""
    assembly_method: str = ""
    sequencing_and_coverage: str = ""
    release_date: str = ""

    bioproject_accession: str = ""
    biosample_accession: str = ""
    raw_reads_accessions: str = ""

    refseq_category: str = ""
    is_refseq: bool = False

    is_preferred_for_species: bool = True
    preferred_over_accession: str = ""
    preference_reason: str = ""

    phylum: str = ""
    superorder: str = ""
    order: str = ""
    family: str = ""
    genus: str = ""
    genus_species: str = ""

    path_to_fna: str = ""
    path_to_gff: str = ""
    path_to_lift_gff: str = ""

    drive_to_fna: str = ""
    drive_to_gff: str = ""
    drive_to_lift_gff: str = ""

    assembly_folder: str = ""
    source_metadata_path: str = ""

    downloaded_at: str = ""
    metadata_last_updated: str = ""
    notes: str = ""

    # Backward-compatible columns retained from the original script.
    # is_current_for_species now mirrors is_preferred_for_species.
    is_current_for_species: bool = True
    # This field is populated only when replacement is repository-confirmed.
    replaced_by_accession: str = ""
    # Retained for older downstream scripts; source_query is the replacement.
    ncbi_search_term: str = ""

    @classmethod
    def from_mapping(cls, mapping: dict[str, Any]) -> "AssemblyRecord":
        """Construct a record from CSV/cache data while ignoring extra keys."""

        allowed = {field.name for field in fields(cls)}
        values = {key: value for key, value in mapping.items() if key in allowed}
        for boolean_field in (
            "is_refseq",
            "is_preferred_for_species",
            "is_current_for_species",
        ):
            if boolean_field in values:
                values[boolean_field] = parse_bool(values[boolean_field])
        return cls(**values)


@dataclass
class FetchedAssembly:
    """A normalized record paired with source metadata and transfer details."""

    record: AssemblyRecord
    raw_metadata: dict[str, Any]
    download_urls: dict[str, str]
    source_adapter: str = ""


@dataclass
class LocalImportConfig:
    """User-supplied metadata and file paths for a local assembly import."""

    repository: str
    local_fasta: Path
    local_gff: Optional[Path] = None
    organism_name: str = ""
    taxon_id: str = ""
    assembly_level: str = ""
    bioproject_accession: str = ""
    biosample_accession: str = ""
    raw_reads_accessions: str = ""
    release_date: str = ""
    assembly_method: str = ""
    assembly_name: str = ""
    sequencing_and_coverage: str = ""
    source_record_url: str = ""
    notes: str = ""


# ---------------------------------------------------------------------------
# Repository normalization
# ---------------------------------------------------------------------------


def _base_record(
    *,
    accession: str,
    repository: str,
    source_database: str,
    source_record_url: str,
    source_query: str,
    organism_name: str,
    taxon_id: str,
    genomes_dir: Path,
) -> AssemblyRecord:
    paths = canonical_output_paths(genomes_dir, accession)
    timestamp = now_iso()
    species_key = build_species_key(organism_name)
    return AssemblyRecord(
        accession=accession,
        accession_root=str(paths["root"]),
        repository=repository,
        source_database=source_database,
        source_record_url=source_record_url,
        source_query=source_query,
        organism_name=organism_name,
        taxon_id=taxon_id,
        species_key=species_key,
        assembly_folder=str(paths["root"]),
        source_metadata_path=str(paths["source_metadata"]),
        downloaded_at=timestamp,
        metadata_last_updated=timestamp,
        is_preferred_for_species=True,
        is_current_for_species=True,
        ncbi_search_term=source_query if repository == "ncbi" else "",
    )


def normalize_ncbi_entry(
    payload: dict[str, Any],
    requested_accession: str,
    *,
    genomes_dir: Path,
    source_query: str,
) -> AssemblyRecord:
    """Normalize one NCBI Datasets summary object to AssemblyRecord."""

    payload_accession = first_nonempty(payload, "accession")
    accession = payload_accession or requested_accession
    organism_name = first_nonempty(payload, "organism.organism_name")
    taxon_id = first_nonempty(payload, "organism.tax_id")
    source_database = "RefSeq" if accession.startswith("GCF_") else "GenBank"
    record = _base_record(
        accession=accession,
        repository="ncbi",
        source_database=source_database,
        source_record_url=f"https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/",
        source_query=source_query,
        organism_name=organism_name,
        taxon_id=taxon_id,
        genomes_dir=genomes_dir,
    )

    record.assembly_name = first_nonempty(
        payload,
        "assembly_info.assembly_name",
        "assembly_name",
    )
    record.assembly_level = normalize_assembly_level(
        first_nonempty(payload, "assembly_info.assembly_level", "assembly_level")
    )
    record.assembly_method = first_nonempty(
        payload,
        "assembly_info.assembly_method",
        "assembly_info.assembly_method_name",
        "assembly_method",
    )
    record.sequencing_and_coverage = first_nonempty(
        payload,
        "assembly_info.sequencing_tech",
        "assembly_info.sequencing_technology",
        "sequencing_tech",
    )
    record.release_date = normalize_date(
        first_nonempty(
            payload,
            "assembly_info.release_date",
            "assembly_info.submission_date",
            "release_date",
        )
    )
    record.bioproject_accession = first_nonempty(
        payload,
        "assembly_info.bioproject_accession",
        "bioproject_accession",
        "bioprojects.0.accession",
    )
    record.biosample_accession = first_nonempty(
        payload,
        "assembly_info.biosample.accession",
        "assembly_info.biosample_accession",
        "biosample.accession",
        "biosample_accession",
    )
    record.raw_reads_accessions = first_nonempty(
        payload,
        "assembly_info.sra_accessions",
        "sra_accessions",
        "raw_reads_accessions",
    )
    record.refseq_category = first_nonempty(
        payload,
        "assembly_info.refseq_category",
        "refseq_category",
    )
    record.is_refseq = accession.startswith("GCF_")

    # Only accept fields that explicitly state replacement/supersession.
    record.replaced_by_accession = first_nonempty(
        payload,
        "assembly_info.replaced_by_accession",
        "assembly_info.current_accession",
        "replaced_by_accession",
    )
    return record


def normalize_gwh_entry(
    payload: dict[str, Any],
    requested_accession: str,
    *,
    genomes_dir: Path,
    source_query: str,
    source_record_url: str = "",
) -> AssemblyRecord:
    """Normalize one GWH Assembly API response to AssemblyRecord."""

    extracted = extract_gwh_assembly_payload(payload)
    if extracted is None:
        raise ValueError(
            "GWH payload does not contain a usable assembly record: "
            + summarize_gwh_response(payload)
        )
    payload = extracted
    returned_accession = first_nonempty(
        payload,
        "assemblyAccession",
        "assembly_accession",
        "accession",
    )

    # Preserve the exact versioned accession requested by the user. For a taxon
    # or blank request (not currently supported by GWH), use the API accession.
    accession = (
        requested_accession
        if requested_accession.upper().startswith("GWH")
        else returned_accession
    )
    organism_name = first_nonempty(payload, "organism")
    taxon_id = first_nonempty(payload, "taxId")
    record_url = source_record_url or (
        "https://ngdc.cncb.ac.cn/gwh/api/public/assembly/"
        + urllib.parse.quote(returned_accession, safe="")
    )
    record = _base_record(
        accession=accession,
        repository="gwh",
        source_database="CNCB-GWH",
        source_record_url=record_url,
        source_query=source_query,
        organism_name=organism_name,
        taxon_id=taxon_id,
        genomes_dir=genomes_dir,
    )

    record.assembly_name = first_nonempty(payload, "assemblyName")
    record.assembly_level = normalize_assembly_level(
        first_nonempty(payload, "assemblyLevel")
    )
    record.assembly_method = first_nonempty(
        payload,
        "assemblyMethod",
        "assemblyMethods",
        "assembly_method",
    )
    record.sequencing_and_coverage = first_nonempty(
        payload,
        "sequencingAndCoverage",
        "sequencingCoverage",
        "sequencing_and_coverage",
        "sequencingTechnology",
    )
    record.release_date = normalize_date(first_nonempty(payload, "releaseTime"))
    record.bioproject_accession = first_nonempty(payload, "bioprojectAccession")
    record.biosample_accession = first_nonempty(payload, "biosampleAccession")
    record.raw_reads_accessions = first_nonempty(
        payload,
        "gsaRawReads",
        "gsaRawRead",
        "rawReadsAccessions",
        "rawReads",
        "gsaAccession",
    )
    record.refseq_category = ""
    record.is_refseq = False
    record.replaced_by_accession = first_nonempty(
        payload,
        "replacedByAccession",
        "replaced_by_accession",
    )
    return record


def gwh_download_urls(payload: dict[str, Any]) -> dict[str, str]:
    """Extract and normalize file URLs from a GWH Assembly API response."""

    extracted = extract_gwh_assembly_payload(payload) or payload
    return {
        "genome": ensure_absolute_url(
            first_nonempty(extracted, "ftpPathDna", "ftp_path_dna", "dnaUrl")
        ),
        "gff3": ensure_absolute_url(
            first_nonempty(extracted, "ftpPathGff", "ftp_path_gff", "gffUrl")
        ),
        "rna": ensure_absolute_url(
            first_nonempty(extracted, "ftpPathRna", "ftp_path_rna", "rnaUrl")
        ),
        "cds": ensure_absolute_url(
            first_nonempty(extracted, "ftpPathCds", "ftp_path_cds", "cdsUrl")
        ),
        "protein": ensure_absolute_url(
            first_nonempty(
                extracted, "ftpPathProtein", "ftp_path_protein", "proteinUrl"
            )
        ),
        "feature": ensure_absolute_url(
            first_nonempty(
                extracted, "ftpPathFeature", "ftp_path_feature", "featureUrl"
            )
        ),
    }


# ---------------------------------------------------------------------------
# Source adapters
# ---------------------------------------------------------------------------


class GenomeSourceAdapter(ABC):
    """Interface for repository metadata retrieval, normalization, and download."""

    repository: str
    source_database: str

    def __init__(
        self,
        *,
        genomes_dir: Path,
        records_dir: Path,
        logger: Callable[[str], None],
        connect_timeout: float = DEFAULT_CONNECT_TIMEOUT,
        read_timeout: float = DEFAULT_READ_TIMEOUT,
    ) -> None:
        self.genomes_dir = genomes_dir
        self.records_dir = records_dir
        self.logger = logger
        self.connect_timeout = connect_timeout
        self.read_timeout = read_timeout

    @abstractmethod
    def fetch(self, search_term: str) -> list[FetchedAssembly]:
        """Retrieve and normalize assemblies matching a repository query."""

    @abstractmethod
    def fetch_accession(self, accession: str) -> FetchedAssembly:
        """Retrieve and normalize exactly one accession."""

    @abstractmethod
    def normalize_cached(
        self,
        raw_metadata: dict[str, Any],
        *,
        accession: str,
        source_query: str,
        source_record_url: str = "",
    ) -> FetchedAssembly:
        """Normalize repository metadata read from a local cache."""

    @abstractmethod
    def download(self, assembly: FetchedAssembly, destination: Path) -> None:
        """Download repository files and create canonical local filenames."""


class NCBISourceAdapter(GenomeSourceAdapter):
    """NCBI Datasets CLI implementation of GenomeSourceAdapter."""

    repository = "ncbi"
    source_database = "NCBI"

    def require_dependencies(self) -> None:
        """Require the Datasets CLI only for an NCBI operation."""

        if shutil.which("datasets") is None:
            raise RuntimeError(
                "NCBI Datasets CLI was not found in PATH. "
                "Install it before running an NCBI query or download."
            )

    def _run_summary(
        self,
        search_term: str,
        *,
        is_accession: bool,
        apply_level_filters: bool,
    ) -> list[dict[str, Any]]:
        self.require_dependencies()
        query_kind = "accession" if is_accession else "taxon"
        command = ["datasets", "summary", "genome", query_kind, search_term]
        if apply_level_filters:
            command.extend(
                [
                    "--assembly-level",
                    "chromosome",
                    "--assembly-level",
                    "complete",
                ]
            )
        command.append("--as-json-lines")

        output_path = self.records_dir / f"{safe_tag(search_term)}_updates.jsonl"
        self.logger(
            f"Checking NCBI for {query_kind} '{search_term}' with command:\n"
            f"{' '.join(command)}\n"
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", encoding="utf-8") as output_handle:
            completed = subprocess.run(
                command,
                stdout=output_handle,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
        if completed.returncode != 0:
            stderr = completed.stderr.strip()
            raise RuntimeError(
                f"NCBI summary failed for {search_term} "
                f"(exit {completed.returncode}): {stderr}"
            )
        rows = list(iter_json_lines(output_path))
        if not rows:
            raise RuntimeError(f"NCBI returned no assemblies for {search_term}")
        return rows

    def _normalize(
        self,
        payload: dict[str, Any],
        *,
        requested_accession: str,
        source_query: str,
    ) -> FetchedAssembly:
        record = normalize_ncbi_entry(
            payload,
            requested_accession,
            genomes_dir=self.genomes_dir,
            source_query=source_query,
        )
        return FetchedAssembly(
            record=record,
            raw_metadata=payload,
            download_urls={},
            source_adapter="ncbi",
        )

    def fetch(self, search_term: str) -> list[FetchedAssembly]:
        if search_term.upper().startswith("GWH"):
            raise ValueError(
                f"'{search_term}' is a GWH accession and cannot be queried with "
                "the NCBI adapter. Use --source gwh or --source auto."
            )
        accession_query = search_term.upper().startswith(("GCA_", "GCF_"))
        payloads = self._run_summary(
            search_term,
            is_accession=accession_query,
            apply_level_filters=True,
        )
        assemblies: list[FetchedAssembly] = []
        for payload in payloads:
            payload_accession = first_nonempty(payload, "accession")
            requested = search_term if accession_query else payload_accession
            assemblies.append(
                self._normalize(
                    payload,
                    requested_accession=requested,
                    source_query=search_term,
                )
            )
        return assemblies

    def fetch_accession(self, accession: str) -> FetchedAssembly:
        if not accession.upper().startswith(("GCA_", "GCF_")):
            raise ValueError(f"Not an NCBI assembly accession: {accession}")
        payloads = self._run_summary(
            accession,
            is_accession=True,
            apply_level_filters=False,
        )
        for payload in payloads:
            if first_nonempty(payload, "accession") == accession:
                return self._normalize(
                    payload,
                    requested_accession=accession,
                    source_query=accession,
                )
        # NCBI occasionally normalizes an accession version. Preserve the exact
        # returned record rather than manufacturing metadata.
        if len(payloads) == 1:
            return self._normalize(
                payloads[0],
                requested_accession=accession,
                source_query=accession,
            )
        raise RuntimeError(f"NCBI did not return accession {accession}")

    def normalize_cached(
        self,
        raw_metadata: dict[str, Any],
        *,
        accession: str,
        source_query: str,
        source_record_url: str = "",
    ) -> FetchedAssembly:
        del source_record_url
        return self._normalize(
            raw_metadata,
            requested_accession=accession,
            source_query=source_query,
        )

    @staticmethod
    def _find_downloaded_file(
        root: Path,
        *,
        suffixes: tuple[str, ...],
        excluded_terms: tuple[str, ...] = (),
    ) -> Optional[Path]:
        candidates: list[Path] = []
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            lower_name = path.name.lower()
            if any(term in lower_name for term in excluded_terms):
                continue
            if lower_name.endswith(suffixes):
                candidates.append(path)
        if not candidates:
            return None
        candidates.sort(
            key=lambda path: (
                0 if "genomic" in path.name.lower() else 1,
                len(path.parts),
                str(path),
            )
        )
        return candidates[0]

    def download(self, assembly: FetchedAssembly, destination: Path) -> None:
        """Recover only NCBI files that are missing locally.

        A valid canonical genome FASTA or GFF is never overwritten.  When only
        the GFF is missing, the Datasets request asks only for ``gff3``.  A
        repository that returns a valid package without a GFF is treated as an
        assembly with no repository annotation rather than as a failed genome
        download.
        """

        self.require_dependencies()
        accession = assembly.record.accession
        destination.mkdir(parents=True, exist_ok=True)
        canonical = canonical_output_paths(self.genomes_dir, accession)

        need_genome = not is_nonempty_file(canonical["fna"])
        need_gff = not is_nonempty_file(canonical["gff"])
        if not need_genome and not need_gff:
            self.logger(
                f"NCBI reconciliation for {accession}: canonical genome FASTA "
                "and GFF are already present and nonempty; no files were "
                "redownloaded.\n"
            )
            return

        include_items: list[str] = []
        if need_genome:
            include_items.append("genome")
        if need_gff:
            include_items.append("gff3")
        include_spec = ",".join(include_items)

        zip_path = self.genomes_dir / f".{accession}.download.zip"
        zip_path.unlink(missing_ok=True)

        def command_variants(include_value: str) -> list[list[str]]:
            base = [
                "datasets",
                "download",
                "genome",
                "accession",
                accession,
            ]
            variants = [
                base
                + [
                    "--reference",
                    "--include",
                    include_value,
                    "--filename",
                    str(zip_path),
                ]
            ]
            if accession.startswith("GCA_"):
                variants.append(
                    base
                    + [
                        "--include",
                        include_value,
                        "--filename",
                        str(zip_path),
                    ]
                )
            return variants

        def request_package(include_value: str) -> tuple[bool, list[str]]:
            errors: list[str] = []
            for command in command_variants(include_value):
                zip_path.unlink(missing_ok=True)
                self.logger(
                    f"Downloading missing NCBI file(s) for {accession}: "
                    f"{' '.join(command)}\n"
                )
                completed = subprocess.run(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=False,
                )
                if (
                    completed.returncode == 0
                    and is_nonempty_file(zip_path)
                ):
                    return True, errors
                errors.append(
                    f"exit {completed.returncode}: {completed.stderr.strip()}"
                )
            zip_path.unlink(missing_ok=True)
            return False, errors

        success, errors = request_package(include_spec)
        optional_gff_request_failed = False

        # If genome+GFF retrieval fails, make one genome-only attempt so an
        # optional annotation problem cannot turn a recoverable genome into a
        # failed assembly.
        if not success and need_genome and need_gff:
            self.logger(
                f"WARNING: NCBI combined genome/GFF retrieval failed for "
                f"{accession}; retrying the required genome FASTA alone.\n"
            )
            success, genome_errors = request_package("genome")
            errors.extend(genome_errors)
            optional_gff_request_failed = success
            if success:
                include_spec = "genome"

        if not success:
            if need_genome:
                raise RuntimeError(
                    f"NCBI download failed for required genome {accession}: "
                    + "; ".join(errors)
                )
            # The genome already exists locally, so inability to retrieve the
            # optional GFF does not make the assembly a failed download.
            self.logger(
                f"WARNING: optional NCBI GFF retrieval failed for {accession}: "
                f"{'; '.join(errors)}. The existing genome FASTA was retained "
                "and the assembly remains valid.\n"
            )
            return

        with tempfile.TemporaryDirectory(
            prefix=f".{accession}.unpack_",
            dir=str(self.genomes_dir),
        ) as temporary_directory:
            unpack_root = Path(temporary_directory)
            try:
                shutil.unpack_archive(str(zip_path), str(unpack_root))
            finally:
                zip_path.unlink(missing_ok=True)

            source_fna = self._find_downloaded_file(
                unpack_root,
                suffixes=(
                    "genomic.fna",
                    "genomic.fna.gz",
                    ".genome.fasta",
                    ".genome.fasta.gz",
                ),
            )
            source_gff = self._find_downloaded_file(
                unpack_root,
                suffixes=(".gff", ".gff3", ".gff.gz", ".gff3.gz"),
                excluded_terms=("lift",),
            )

            if need_genome:
                if source_fna is None:
                    raise RuntimeError(
                        f"NCBI archive for {accession} did not contain the "
                        "required genome FASTA"
                    )
                if not is_nonempty_file(canonical["fna"]):
                    copy_or_decompress_local(source_fna, canonical["fna"])
                    self.logger(
                        f"Recovered missing NCBI genome FASTA for {accession}: "
                        f"{canonical['fna']}\n"
                    )

            if need_gff:
                if source_gff is not None:
                    if not is_nonempty_file(canonical["gff"]):
                        copy_or_decompress_local(source_gff, canonical["gff"])
                        self.logger(
                            f"Recovered missing NCBI GFF for {accession}: "
                            f"{canonical['gff']}\n"
                        )
                elif optional_gff_request_failed:
                    self.logger(
                        f"WARNING: the required genome FASTA was recovered for "
                        f"{accession}, but the optional NCBI GFF could not be "
                        "retrieved. path_to_gff will remain blank.\n"
                    )
                else:
                    self.logger(
                        f"NCBI provides no repository GFF for {accession}; "
                        "path_to_gff will remain blank.\n"
                    )

            # Preserve an existing raw NCBI package rather than replacing it.
            # Canonical recovery files above are the authoritative local files.
            raw_root = destination / "ncbi_dataset_package"
            if not raw_root.exists():
                shutil.copytree(unpack_root, raw_root)
            else:
                self.logger(
                    f"Existing NCBI package directory for {accession} was left "
                    "unchanged.\n"
                )


class GWHSourceAdapter(GenomeSourceAdapter):
    """CNCB Genome Warehouse Assembly API implementation."""

    repository = "gwh"
    source_database = "CNCB-GWH"

    def _normalize(
        self,
        payload: dict[str, Any],
        *,
        requested_accession: str,
        source_query: str,
        source_record_url: str,
    ) -> FetchedAssembly:
        record = normalize_gwh_entry(
            payload,
            requested_accession,
            genomes_dir=self.genomes_dir,
            source_query=source_query,
            source_record_url=source_record_url,
        )
        return FetchedAssembly(
            record=record,
            raw_metadata=payload,
            download_urls=gwh_download_urls(payload),
            source_adapter="gwh",
        )

    def fetch(self, search_term: str) -> list[FetchedAssembly]:
        if not search_term.upper().startswith("GWH"):
            raise ValueError(
                "GWH currently supports accession searches only. "
                f"The supplied term was: {search_term!r}"
            )
        return [self.fetch_accession(search_term)]

    def fetch_accession(self, accession: str) -> FetchedAssembly:
        if not accession.upper().startswith("GWH"):
            raise ValueError(f"Not a GWH assembly accession: {accession}")
        self.logger(f"Querying CNCB-GWH assembly API for {accession}\n")
        payload, record_url = fetch_gwh_assembly(
            accession,
            connect_timeout=self.connect_timeout,
            read_timeout=self.read_timeout,
        )
        return self._normalize(
            payload,
            requested_accession=accession,
            source_query=accession,
            source_record_url=record_url,
        )

    def normalize_cached(
        self,
        raw_metadata: dict[str, Any],
        *,
        accession: str,
        source_query: str,
        source_record_url: str = "",
    ) -> FetchedAssembly:
        return self._normalize(
            raw_metadata,
            requested_accession=accession,
            source_query=source_query,
            source_record_url=source_record_url,
        )

    def download(self, assembly: FetchedAssembly, destination: Path) -> None:
        """Recover only GWH genome/GFF files that are missing locally."""

        accession = assembly.record.accession
        paths = canonical_output_paths(self.genomes_dir, accession)
        destination.mkdir(parents=True, exist_ok=True)

        need_genome = not is_nonempty_file(paths["fna"])
        need_gff = not is_nonempty_file(paths["gff"])

        genome_url = assembly.download_urls.get("genome", "").strip()
        if need_genome:
            if not genome_url:
                raise RuntimeError(
                    f"GWH assembly {accession} does not provide the required "
                    "ftpPathDna/genome download URL"
                )
            self.logger(
                f"Downloading missing GWH genome FASTA for {accession}: "
                f"{genome_url}\n"
            )
            download_and_decompress(
                genome_url,
                paths["fna"],
                connect_timeout=self.connect_timeout,
                read_timeout=self.read_timeout,
            )
        else:
            self.logger(
                f"GWH reconciliation for {accession}: existing nonempty genome "
                "FASTA retained; it was not redownloaded.\n"
            )

        gff_url = assembly.download_urls.get("gff3", "").strip()
        if need_gff:
            if gff_url:
                self.logger(
                    f"Downloading missing GWH GFF for {accession}: {gff_url}\n"
                )
                try:
                    download_and_decompress(
                        gff_url,
                        paths["gff"],
                        connect_timeout=self.connect_timeout,
                        read_timeout=self.read_timeout,
                    )
                except Exception as exc:
                    self.logger(
                        f"WARNING: optional GWH GFF download failed for "
                        f"{accession}: {type(exc).__name__}: {exc}. The genome "
                        "FASTA was retained and the assembly remains valid; "
                        "path_to_gff will remain blank.\n"
                    )
            else:
                self.logger(
                    f"GWH provides no repository annotation/GFF for {accession}; "
                    "path_to_gff will remain blank.\n"
                )
        else:
            self.logger(
                f"GWH reconciliation for {accession}: existing nonempty GFF "
                "retained; it was not redownloaded.\n"
            )

        other_available = [
            key
            for key in ("rna", "cds", "protein", "feature")
            if assembly.download_urls.get(key)
        ]
        if other_available:
            self.logger(
                f"GWH provides optional {', '.join(other_available)} files for "
                f"{accession}. Their URLs are preserved in source_metadata.json; "
                "the workflow reconciles genome and GFF files only.\n"
            )


class LocalSourceAdapter(GenomeSourceAdapter):
    """Import user-downloaded files while preserving repository provenance.

    This adapter performs no repository network request. The assembly accession
    remains the positional search term, while ``LocalImportConfig.repository``
    records the repository from which the user manually obtained the files.
    """

    repository = "local"
    source_database = "Local import"

    def __init__(
        self,
        *,
        genomes_dir: Path,
        records_dir: Path,
        logger: Callable[[str], None],
        config: LocalImportConfig,
        connect_timeout: float = DEFAULT_CONNECT_TIMEOUT,
        read_timeout: float = DEFAULT_READ_TIMEOUT,
    ) -> None:
        super().__init__(
            genomes_dir=genomes_dir,
            records_dir=records_dir,
            logger=logger,
            connect_timeout=connect_timeout,
            read_timeout=read_timeout,
        )
        self.config = config

    @staticmethod
    def _validate_local_file(path: Path, label: str) -> Path:
        resolved = path.expanduser().resolve()
        if not resolved.exists():
            raise FileNotFoundError(f"{label} does not exist: {resolved}")
        if not resolved.is_file():
            raise ValueError(f"{label} is not a regular file: {resolved}")
        if resolved.stat().st_size == 0:
            raise ValueError(f"{label} is empty: {resolved}")
        return resolved

    def _provenance_note(self, repository: str) -> str:
        if repository == "gwh":
            reason = (
                "Files were manually downloaded through the CNCB Genome "
                "Warehouse browser because the public GWH Assembly API returned "
                "an all-null/unusable record for this accession."
            )
        else:
            reason = (
                "Files were manually downloaded outside this script and imported "
                "through the repository-independent local adapter."
            )
        if self.config.notes.strip():
            reason = f"{reason} {self.config.notes.strip()}"
        return reason

    def _build_assembly(self, accession: str) -> FetchedAssembly:
        fasta = self._validate_local_file(self.config.local_fasta, "Local FASTA")
        gff: Optional[Path] = None
        if self.config.local_gff is not None:
            gff = self._validate_local_file(self.config.local_gff, "Local GFF")

        repository = normalize_repository_name(self.config.repository)
        if not repository:
            raise ValueError("--repository must not be blank for --source local")
        source_database = source_database_for_repository(repository, accession)
        source_record_url = (
            self.config.source_record_url.strip()
            or default_source_record_url(repository, accession)
        )
        record = _base_record(
            accession=accession,
            repository=repository,
            source_database=source_database,
            source_record_url=source_record_url,
            source_query=accession,
            organism_name=self.config.organism_name.strip(),
            taxon_id=self.config.taxon_id.strip(),
            genomes_dir=self.genomes_dir,
        )
        record.assembly_name = self.config.assembly_name.strip()
        record.assembly_level = normalize_assembly_level(
            self.config.assembly_level
        )
        record.assembly_method = self.config.assembly_method.strip()
        record.sequencing_and_coverage = (
            self.config.sequencing_and_coverage.strip()
        )
        record.release_date = normalize_date(self.config.release_date)
        record.bioproject_accession = self.config.bioproject_accession.strip()
        record.biosample_accession = self.config.biosample_accession.strip()
        record.raw_reads_accessions = self.config.raw_reads_accessions.strip()
        record.refseq_category = ""
        record.is_refseq = (
            repository == "ncbi" and accession.upper().startswith("GCF_")
        )
        record.notes = self._provenance_note(repository)

        raw_metadata: dict[str, Any] = {
            "assemblyAccession": accession,
            "organism": record.organism_name,
            "taxId": record.taxon_id,
            "assemblyName": record.assembly_name,
            "assemblyLevel": record.assembly_level,
            "assemblyMethod": record.assembly_method,
            "sequencingAndCoverage": record.sequencing_and_coverage,
            "bioprojectAccession": record.bioproject_accession,
            "biosampleAccession": record.biosample_accession,
            "gsaRawReads": record.raw_reads_accessions,
            "releaseTime": record.release_date,
            "_local_import": {
                "mode": "manual-local-import",
                "repository": repository,
                "source_database": source_database,
                "source_fasta": str(fasta),
                "source_gff": str(gff) if gff else "",
                "source_record_url": source_record_url,
                "original_search_term": accession,
                "import_reason": record.notes,
                "metadata_supplied_by_user": True,
                "imported_at": now_iso(),
            },
        }
        return FetchedAssembly(
            record=record,
            raw_metadata=raw_metadata,
            download_urls={},
            source_adapter="local",
        )

    def fetch(self, search_term: str) -> list[FetchedAssembly]:
        if not search_term.strip():
            raise ValueError("A nonempty accession/identifier is required")
        return [self.fetch_accession(search_term)]

    def fetch_accession(self, accession: str) -> FetchedAssembly:
        self.logger(
            f"Preparing repository-independent local import for {accession}; "
            f"provenance repository={self.config.repository}.\n"
        )
        return self._build_assembly(accession)

    def normalize_cached(
        self,
        raw_metadata: dict[str, Any],
        *,
        accession: str,
        source_query: str,
        source_record_url: str = "",
    ) -> FetchedAssembly:
        del source_query, source_record_url
        local = raw_metadata.get("_local_import")
        if not isinstance(local, dict):
            raise ValueError("Cached local-import metadata lacks _local_import")
        repository = normalize_repository_name(
            str(local.get("repository", "") or self.config.repository)
        )
        record = _base_record(
            accession=accession,
            repository=repository,
            source_database=str(
                local.get("source_database", "")
                or source_database_for_repository(repository, accession)
            ),
            source_record_url=str(local.get("source_record_url", "") or ""),
            source_query=accession,
            organism_name=first_nonempty(raw_metadata, "organism"),
            taxon_id=first_nonempty(raw_metadata, "taxId"),
            genomes_dir=self.genomes_dir,
        )
        record.assembly_name = first_nonempty(raw_metadata, "assemblyName")
        record.assembly_level = normalize_assembly_level(
            first_nonempty(raw_metadata, "assemblyLevel")
        )
        record.assembly_method = first_nonempty(raw_metadata, "assemblyMethod")
        record.sequencing_and_coverage = first_nonempty(
            raw_metadata, "sequencingAndCoverage"
        )
        record.release_date = normalize_date(
            first_nonempty(raw_metadata, "releaseTime")
        )
        record.bioproject_accession = first_nonempty(
            raw_metadata, "bioprojectAccession"
        )
        record.biosample_accession = first_nonempty(
            raw_metadata, "biosampleAccession"
        )
        record.raw_reads_accessions = first_nonempty(
            raw_metadata, "gsaRawReads"
        )
        record.notes = str(local.get("import_reason", "") or "")
        return FetchedAssembly(
            record=record,
            raw_metadata=raw_metadata,
            download_urls={},
            source_adapter="local",
        )

    def download(self, assembly: FetchedAssembly, destination: Path) -> None:
        accession = assembly.record.accession
        paths = canonical_output_paths(self.genomes_dir, accession)
        destination.mkdir(parents=True, exist_ok=True)
        local = assembly.raw_metadata.get("_local_import")
        if not isinstance(local, dict):
            raise ValueError(
                f"Local-import transfer metadata is missing for {accession}"
            )
        fasta_text = str(local.get("source_fasta", "") or "").strip()
        if not fasta_text:
            raise ValueError(f"Local FASTA path is missing for {accession}")
        fasta = self._validate_local_file(Path(fasta_text), "Local FASTA")
        if not is_nonempty_file(paths["fna"]):
            self.logger(
                f"Importing local genome FASTA for {accession}: {fasta} -> "
                f"{paths['fna']}\n"
            )
            copy_or_decompress_local(fasta, paths["fna"])
        else:
            self.logger(
                f"Existing nonempty genome FASTA for {accession} was retained; "
                "the local import did not overwrite it.\n"
            )

        gff_text = str(local.get("source_gff", "") or "").strip()
        if gff_text:
            gff = self._validate_local_file(Path(gff_text), "Local GFF")
            if not is_nonempty_file(paths["gff"]):
                self.logger(
                    f"Importing optional local GFF for {accession}: {gff} -> "
                    f"{paths['gff']}\n"
                )
                copy_or_decompress_local(gff, paths["gff"])
            else:
                self.logger(
                    f"Existing nonempty GFF for {accession} was retained; the "
                    "local import did not overwrite it.\n"
                )
        else:
            self.logger(
                f"No local GFF was supplied for {accession}; path_to_gff will "
                "remain blank unless an annotation already exists locally.\n"
            )


# ---------------------------------------------------------------------------
# Preference semantics helpers
# ---------------------------------------------------------------------------


def _metadata_order(rows: list[dict[str, Any]]) -> dict[str, int]:
    return {
        str(row.get("accession", "")): index
        for index, row in enumerate(rows)
        if row.get("accession")
    }


def apply_species_preference_semantics(
    records: list[AssemblyRecord],
    *,
    existing_rows: Optional[list[dict[str, Any]]] = None,
    logger: Optional[Callable[[str], None]] = None,
) -> None:
    """Preserve one preferred assembly without inventing replacement events.

    Existing preference fields take priority. For a species represented by a
    single assembly, that assembly is preferred. When multiple assemblies have
    no previously selected preferred record, none is silently promoted; a
    warning is logged so the user can make an explicit choice later.
    """

    existing_rows = existing_rows or []
    existing_by_accession = {
        str(row.get("accession", "")): row
        for row in existing_rows
        if row.get("accession")
    }
    order = _metadata_order(existing_rows)

    for record in records:
        prior = existing_by_accession.get(record.accession)
        if prior:
            preferred_value = prior.get("is_preferred_for_species")
            if preferred_value in (None, ""):
                preferred_value = prior.get("is_current_for_species")
            record.is_preferred_for_species = parse_bool(
                preferred_value,
                default=record.is_preferred_for_species,
            )
            record.preferred_over_accession = str(
                prior.get("preferred_over_accession", "") or ""
            )
            record.preference_reason = str(
                prior.get("preference_reason", "") or ""
            )
            # The old script generated replacement links whenever two
            # assemblies shared a species. Those links have no provenance.
            # Preserve only version-to-version legacy links; newly normalized
            # repository-confirmed links remain untouched on new records.
            legacy_replacement = str(prior.get("replaced_by_accession", "") or "")
            if legacy_replacement and root_acc(legacy_replacement) == root_acc(record.accession):
                record.replaced_by_accession = legacy_replacement
            elif legacy_replacement and logger:
                record.replaced_by_accession = ""
                logger(
                    f"Cleared unverified legacy replacement link "
                    f"{record.accession} -> {legacy_replacement}; the accessions "
                    "do not represent versions of the same assembly.\n"
                )

    grouped: dict[str, list[AssemblyRecord]] = {}
    for record in records:
        key = record.species_key or root_acc(record.accession)
        grouped.setdefault(key, []).append(record)

    for species_key, group in grouped.items():
        if len(group) == 1:
            group[0].is_preferred_for_species = True
            group[0].is_current_for_species = True
            if not group[0].preference_reason:
                group[0].preference_reason = "Only assembly currently stored for species"
            continue

        # Defaults on freshly normalized records are not an explicit user
        # preference. Only metadata rows that existed before this operation can
        # designate a preferred assembly in a multi-assembly species group.
        preferred = [
            record
            for record in group
            if record.accession in existing_by_accession
            and record.is_preferred_for_species
        ]

        if len(preferred) > 1:
            preferred.sort(key=lambda record: order.get(record.accession, 10**9))
            selected = preferred[0]
            for record in group:
                record.is_preferred_for_species = record is selected
                record.is_current_for_species = record is selected
                if record is not selected and not record.preference_reason:
                    record.preference_reason = (
                        "Not automatically preferred; an existing preferred "
                        "assembly was preserved"
                    )
            if logger:
                logger(
                    f"Multiple preferred assemblies were found for {species_key}; "
                    f"preserved {selected.accession} based on existing metadata order.\n"
                )
        elif len(preferred) == 1:
            selected = preferred[0]
            for record in group:
                record.is_preferred_for_species = record is selected
                record.is_current_for_species = record is selected
                if record is not selected and not record.preference_reason:
                    record.preference_reason = (
                        f"Existing preferred assembly preserved: {selected.accession}"
                    )
        else:
            for record in group:
                record.is_preferred_for_species = False
                record.is_current_for_species = False
                if not record.preference_reason:
                    record.preference_reason = (
                        "No assembly was automatically preferred because multiple "
                        "independent assemblies are stored for this species"
                    )
            if logger:
                logger(
                    f"No preferred assembly is designated for {species_key}; "
                    "multiple independent assemblies were preserved.\n"
                )


# ---------------------------------------------------------------------------
# Genome manager
# ---------------------------------------------------------------------------


class GenomeManagerHybrid:
    """Coordinate repository adapters and repository-independent workflows."""

    def __init__(
        self,
        genomes_dir: str,
        search_term: str,
        *,
        source: str = "auto",
        project_name: Optional[str] = None,
        project_policy: str = "exploratory",
        project_root: Optional[str] = None,
        manifest_notes: str = "",
        retry_failures: bool = False,
        rebuild_metadata: bool = False,
        dry_run_rebuild: bool = False,
        project_accessions_csv: Optional[str] = None,
        project_accessions_txt: Optional[str] = None,
        project_accessions: Optional[list[str]] = None,
        connect_timeout: float = DEFAULT_CONNECT_TIMEOUT,
        read_timeout: float = DEFAULT_READ_TIMEOUT,
        local_import_config: Optional[LocalImportConfig] = None,
    ) -> None:
        self.genomes_dir = require_genomes_dir(Path(genomes_dir))
        self.records_dir = self.genomes_dir / "records"
        self.project_manifests_dir = self.records_dir / "project_manifests"

        self.search_term = search_term.strip()
        self.requested_source = source
        self.selected_source = resolve_source(self.search_term, source)
        self.source_was_auto = source == "auto"
        self.safe_tag = safe_tag(self.search_term)

        self.project_name = project_name
        self.project_policy = project_policy
        self.project_root = str(Path(project_root).expanduser().resolve()) if project_root else ""
        self.manifest_notes = manifest_notes
        self.retry_failures = retry_failures
        self.rebuild_metadata = rebuild_metadata
        self.dry_run_rebuild = dry_run_rebuild
        self.project_accessions_csv = project_accessions_csv
        self.project_accessions_txt = project_accessions_txt
        self.project_accessions = project_accessions or []
        self.root_local = self.genomes_dir.parent
        self.root_drive = Path("/content/drive/Othercomputers/macbook/projects")
        self.connect_timeout = connect_timeout
        self.read_timeout = read_timeout
        self.local_import_config = local_import_config

        self.genomes_catalog_jsonl = self.records_dir / "genomes_catalog.jsonl"
        self.genomes_metadata_csv = self.records_dir / "genomes_metadata.csv"
        self.species_assembly_index = self.records_dir / "species_assembly_index.csv"
        self.update_history_jsonl = self.records_dir / "assembly_update_history.jsonl"
        self.project_registry_csv = self.records_dir / "project_registry.csv"
        self.records_log_txt = self.records_dir / "genome_auto_update_records_log.txt"
        self.download_failures_csv = self.records_dir / "download_failures.csv"

        self.data: list[FetchedAssembly] = []
        self.successfully_downloaded_data: list[str] = []
        self.failed_downloads: list[dict[str, str]] = []
        self.project_selected_accessions: list[str] = []
        self.rebuild_preview: dict[str, Any] = {
            "local_accessions": [],
            "rebuilt_accessions": [],
            "broken_accessions": [],
            "would_trash": [],
            "missing_from_rebuild": [],
            "new_rows": [],
        }
        self._taxonkit_warning_emitted = False
        self.adapters: dict[str, GenomeSourceAdapter] = {
            "ncbi": NCBISourceAdapter(
                genomes_dir=self.genomes_dir,
                records_dir=self.records_dir,
                logger=self.log,
                connect_timeout=connect_timeout,
                read_timeout=read_timeout,
            ),
            "gwh": GWHSourceAdapter(
                genomes_dir=self.genomes_dir,
                records_dir=self.records_dir,
                logger=self.log,
                connect_timeout=connect_timeout,
                read_timeout=read_timeout,
            ),
        }
        if local_import_config is not None:
            self.adapters["local"] = LocalSourceAdapter(
                genomes_dir=self.genomes_dir,
                records_dir=self.records_dir,
                logger=self.log,
                config=local_import_config,
                connect_timeout=connect_timeout,
                read_timeout=read_timeout,
            )

    def initialize_records(self) -> None:
        self.genomes_dir.mkdir(parents=True, exist_ok=True)
        self.records_dir.mkdir(parents=True, exist_ok=True)
        self.project_manifests_dir.mkdir(parents=True, exist_ok=True)
        self.genomes_catalog_jsonl.touch(exist_ok=True)
        self.records_log_txt.touch(exist_ok=True)

    def log(self, message: str) -> None:
        """Write a timestamped message to the persistent records log."""

        self.records_dir.mkdir(parents=True, exist_ok=True)
        with self.records_log_txt.open("a", encoding="utf-8") as handle:
            handle.write(f"[{now_iso()}] {message}")
            if not message.endswith("\n"):
                handle.write("\n")

    def get_adapter(self, source: str) -> GenomeSourceAdapter:
        try:
            return self.adapters[source]
        except KeyError as exc:
            raise ValueError(f"No adapter is registered for source '{source}'") from exc

    def adapter_for_accession(self, accession: str) -> GenomeSourceAdapter:
        return self.get_adapter(detect_source(accession))

    # Retained method name for compatibility with prior code/users.
    def require_datasets_cli(self) -> None:
        adapter = self.get_adapter("ncbi")
        assert isinstance(adapter, NCBISourceAdapter)
        adapter.require_dependencies()
        self.log("NCBI Datasets CLI detected.\n")

    def accession_root(self, accession: str) -> Path:
        return self.genomes_dir / accession

    def source_metadata_file(self, accession: str) -> Path:
        return canonical_output_paths(self.genomes_dir, accession)["source_metadata"]

    def _existing_catalog_accessions(self) -> set[str]:
        accessions: set[str] = set()
        for row in iter_json_lines(self.genomes_catalog_jsonl):
            accession = first_nonempty(
                row,
                "accession",
                "normalized_record.accession",
            )
            if accession:
                accessions.add(accession)
        return accessions

    def query_source(self) -> None:
        """Query the selected adapter and reconcile every returned accession.

        Catalog membership is deliberately not a download-skip condition.  An
        accession may already be known to the project while still be missing a
        repository-provided GFF or even its canonical genome FASTA.
        """

        mode = "automatically detected" if self.source_was_auto else "explicitly selected"
        self.log(
            f"Repository {self.selected_source} was {mode} for search term "
            f"'{self.search_term}'.\n"
        )
        adapter = self.get_adapter(self.selected_source)
        self.data = adapter.fetch(self.search_term)
        existing = self._existing_catalog_accessions()
        reconcile_count = sum(
            1 for item in self.data if item.record.accession in existing
        )
        if reconcile_count:
            self.log(
                f"Rechecking {reconcile_count} catalogued accession(s) file-by-file; "
                "catalog membership does not imply that local files are complete.\n"
            )

    # Compatibility wrappers for the original NCBI-specific method names.
    def query_ncbi(self) -> None:
        self.selected_source = "ncbi"
        self.data = self.get_adapter("ncbi").fetch(self.search_term)

    def build_ncbi_summary_command(self) -> str:
        kind = "accession" if self.search_term.upper().startswith(("GCA_", "GCF_")) else "taxon"
        return (
            f"datasets summary genome {kind} {self.search_term!r} "
            "--assembly-level chromosome --assembly-level complete --as-json-lines"
        )

    def log_species_with_multiple_accessions(self) -> None:
        seen: dict[str, list[str]] = {}
        for assembly in self.data:
            seen.setdefault(assembly.record.organism_name, []).append(
                assembly.record.accession
            )
        for organism, accessions in seen.items():
            if organism and len(accessions) > 1:
                self.log(
                    f"Multiple accessions returned for {organism}: "
                    f"{', '.join(accessions)}\n"
                )

    def save_source_metadata(self, assembly: FetchedAssembly) -> None:
        """Persist raw and normalized repository metadata for offline rebuilds."""

        record = assembly.record
        path = self.source_metadata_file(record.accession)
        timestamp = now_iso()
        payload = {
            "schema_version": SOURCE_METADATA_SCHEMA_VERSION,
            "script_version": SCRIPT_VERSION,
            "repository": record.repository,
            "ingest_source": assembly.source_adapter or record.repository,
            "source_database": record.source_database,
            "source_record_url": record.source_record_url,
            "source_query": record.source_query,
            "requested_accession": record.accession,
            "downloaded_at": record.downloaded_at or timestamp,
            "metadata_last_updated": timestamp,
            "raw_metadata": assembly.raw_metadata,
            "download_urls": assembly.download_urls,
            "normalized_record": asdict(record),
        }
        write_json(path, payload)
        record.source_metadata_path = str(path)
        record.metadata_last_updated = timestamp

    def download_new_data(self) -> None:
        """Reconcile fetched assemblies against canonical local files.

        Genome FASTA and GFF are evaluated independently.  Existing nonempty
        canonical files are preserved, while missing files are requested from
        the selected repository adapter.  A missing or unavailable optional GFF
        never makes an otherwise valid genome assembly a failed download.
        """

        for assembly in self.data:
            record = assembly.record
            accession = record.accession
            adapter = self.get_adapter(assembly.source_adapter or record.repository)
            destination = self.accession_root(accession)
            created_new_directory = not destination.exists()
            destination.mkdir(parents=True, exist_ok=True)
            canonical = canonical_output_paths(self.genomes_dir, accession)
            before_fna = is_nonempty_file(canonical["fna"])
            before_gff = is_nonempty_file(canonical["gff"])

            try:
                adapter.download(assembly, destination)

                after_fna = is_nonempty_file(canonical["fna"])
                after_gff = is_nonempty_file(canonical["gff"])
                if not after_fna:
                    raise RuntimeError(
                        f"Required canonical genome FASTA is still missing or "
                        f"empty after reconciliation: {canonical['fna']}"
                    )

                added: list[str] = []
                if not before_fna and after_fna:
                    added.append("genome FASTA")
                if not before_gff and after_gff:
                    added.append("GFF")

                self.enrich_record_from_local_files(record)
                self.enrich_record_lineage(record)
                self.save_source_metadata(assembly)
                self.successfully_downloaded_data.append(accession)

                if added:
                    self.log(
                        f"Successfully reconciled {accession}; added missing "
                        f"{', '.join(added)} without replacing valid existing "
                        "files. Metadata was refreshed.\n"
                    )
                else:
                    self.log(
                        f"Successfully reconciled {accession}; the existing genome "
                        "FASTA remains valid and no repository file was added. "
                        "Metadata was refreshed.\n"
                    )
            except Exception as exc:
                message = f"{type(exc).__name__}: {exc}"
                self.log(
                    f"ERROR reconciling {accession} from {record.repository}: "
                    f"{message}\n"
                )
                self.append_download_failure(
                    accession=accession,
                    repository=record.repository,
                    ingest_source=assembly.source_adapter or record.repository,
                    stage="download",
                    command=f"{assembly.source_adapter or record.repository} adapter transfer",
                    error_note=message,
                )
                self.failed_downloads.append(
                    {
                        "accession": accession,
                        "repository": record.repository,
                        "stage": "download",
                        "error": message,
                    }
                )
                # Remove an empty directory created for a failed first attempt,
                # but never delete an existing user directory containing data.
                if created_new_directory and destination.exists():
                    try:
                        if not any(destination.iterdir()):
                            destination.rmdir()
                    except OSError:
                        pass

    def append_download_failure(
        self,
        accession: str,
        stage: str,
        command: str,
        error_note: str,
        repository: str = "",
        ingest_source: str = "",
    ) -> None:
        row = {
            "timestamp": now_iso(),
            "accession": accession,
            "repository": repository or detect_source(accession),
            "ingest_source": ingest_source or repository or detect_source(accession),
            "stage": stage,
            "command": command,
            "error_note": error_note,
        }
        existing = read_csv_rows(self.download_failures_csv)
        existing.append(row)
        write_csv_rows(
            self.download_failures_csv,
            [
                "timestamp",
                "accession",
                "repository",
                "ingest_source",
                "stage",
                "command",
                "error_note",
            ],
            existing,
        )

    def fetch_summary_for_accession(self, accession: str) -> FetchedAssembly:
        """Fetch a normalized accession using the centrally detected adapter."""

        adapter = self.adapter_for_accession(accession)
        self.log(
            f"Remote metadata fallback: querying {adapter.repository} for "
            f"{accession}.\n"
        )
        return adapter.fetch_accession(accession)

    def load_failed_accessions_for_retry(self) -> None:
        self.data = []
        if not self.download_failures_csv.exists():
            self.log("No download_failures.csv found; nothing to retry.\n")
            return
        rows = read_csv_rows(self.download_failures_csv)
        seen: set[str] = set()
        for row in rows:
            accession = row.get("accession", "").strip()
            if not accession or accession in seen:
                continue
            seen.add(accession)
            ingest_source = str(
                row.get("ingest_source", "")
                or row.get("repository", "")
                or detect_source(accession)
            )
            try:
                if ingest_source == "local":
                    if "local" not in self.adapters:
                        self.log(
                            f"Skipping local-import retry for {accession}; rerun "
                            "with --source local and the original local file arguments.\n"
                        )
                        continue
                    self.data.append(
                        self.get_adapter("local").fetch_accession(accession)
                    )
                else:
                    self.data.append(self.fetch_summary_for_accession(accession))
            except Exception as exc:
                self.log(
                    f"Could not prepare retry accession {accession}: {exc}\n"
                )
        self.log(f"Loaded {len(self.data)} failed accession(s) for retry.\n")

    def remove_resolved_failures(self, resolved_accessions: list[str]) -> None:
        if not resolved_accessions or not self.download_failures_csv.exists():
            return
        resolved = set(resolved_accessions)
        rows = read_csv_rows(self.download_failures_csv)
        kept = [
            row
            for row in rows
            if row.get("accession", "").strip() not in resolved
        ]
        write_csv_rows(
            self.download_failures_csv,
            [
                "timestamp",
                "accession",
                "repository",
                "ingest_source",
                "stage",
                "command",
                "error_note",
            ],
            kept,
        )
        self.log(
            f"Removed {len(rows) - len(kept)} resolved failure row(s).\n"
        )

    def find_accession_files(self, accession: str) -> tuple[str, str, str]:
        """Locate canonical files first, then legacy/repository package paths."""

        root = self.accession_root(accession)
        if not root.exists():
            return "", "", ""
        canonical = canonical_output_paths(self.genomes_dir, accession)
        found_fna = str(canonical["fna"]) if is_nonempty_file(canonical["fna"]) else ""
        found_gff = str(canonical["gff"]) if is_nonempty_file(canonical["gff"]) else ""
        found_lift = (
            str(canonical["lift_gff"]) if is_nonempty_file(canonical["lift_gff"]) else ""
        )

        fasta_candidates: list[Path] = []
        gff_candidates: list[Path] = []
        lift_candidates: list[Path] = []
        for path in root.rglob("*"):
            if (
                not path.is_file()
                or path.name.endswith(".part")
                or not is_nonempty_file(path)
            ):
                continue
            lower = path.name.lower()
            if not found_lift and "lift" in lower and lower.endswith((".gff", ".gff3")):
                lift_candidates.append(path)
                continue
            if not found_gff and "lift" not in lower and lower.endswith((".gff", ".gff3")):
                gff_candidates.append(path)
            if not found_fna and lower.endswith((".fna", ".fa", ".fasta")):
                if any(token in lower for token in ("genomic", "genome")):
                    fasta_candidates.append(path)

        if not found_fna and fasta_candidates:
            fasta_candidates.sort(key=lambda path: (len(path.parts), str(path)))
            found_fna = str(fasta_candidates[0])
        if not found_gff and gff_candidates:
            gff_candidates.sort(key=lambda path: (len(path.parts), str(path)))
            found_gff = str(gff_candidates[0])
        if not found_lift and lift_candidates:
            lift_candidates.sort(key=lambda path: (len(path.parts), str(path)))
            found_lift = str(lift_candidates[0])
        return found_fna, found_gff, found_lift

    def inspect_accession_root(self, accession: str) -> dict[str, Any]:
        found_fna, found_gff, found_lift = self.find_accession_files(accession)
        return {
            "accession": accession,
            "accession_root": self.accession_root(accession),
            "found_fna": found_fna,
            "found_gff": found_gff,
            "found_lift": found_lift,
            "is_broken": not bool(found_fna),
        }

    def to_local_and_drive(self, absolute_path: str) -> tuple[str, str]:
        if not absolute_path:
            return "", ""
        path = Path(absolute_path).expanduser().resolve()
        try:
            relative = path.relative_to(self.root_local)
        except ValueError:
            return str(path), ""
        return str(self.root_local / relative), str(self.root_drive / relative)

    def get_lineage(self, taxon_id: str) -> dict[str, str]:
        """Use TaxonKit when available; lineage enrichment is nonfatal."""

        empty = {
            "phylum": "",
            "superorder": "",
            "order": "",
            "family": "",
            "genus": "",
            "genus_species": "",
        }
        if not taxon_id:
            return empty
        if shutil.which("taxonkit") is None:
            if not self._taxonkit_warning_emitted:
                self.log(
                    "TaxonKit was not found in PATH; lineage columns will remain "
                    "blank where they are not already cached.\n"
                )
                self._taxonkit_warning_emitted = True
            return empty

        lineage = subprocess.run(
            ["taxonkit", "lineage"],
            input=f"{taxon_id}\n",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if lineage.returncode != 0:
            self.log(
                f"TaxonKit lineage failed for taxon {taxon_id}: "
                f"{lineage.stderr.strip()}\n"
            )
            return empty
        reformatted = subprocess.run(
            [
                "taxonkit",
                "reformat",
                "-f",
                "{p}\\t{c}\\t{o}\\t{f}\\t{g}\\t{s}",
            ],
            input=lineage.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if reformatted.returncode != 0:
            self.log(
                f"TaxonKit reformat failed for taxon {taxon_id}: "
                f"{reformatted.stderr.strip()}\n"
            )
            return empty
        values = reformatted.stdout.strip().split("\t")
        values += [""] * (6 - len(values))
        cleaned = [value.strip("{}") for value in values[:6]]
        return dict(zip(empty.keys(), cleaned))

    def enrich_record_lineage(self, record: AssemblyRecord) -> None:
        if all(
            (
                record.phylum,
                record.order,
                record.family,
                record.genus,
                record.genus_species,
            )
        ):
            return
        lineage = self.get_lineage(record.taxon_id)
        for key, value in lineage.items():
            if value and not getattr(record, key):
                setattr(record, key, value)

    def enrich_record_from_local_files(self, record: AssemblyRecord) -> None:
        found_fna, found_gff, found_lift = self.find_accession_files(record.accession)
        local_fna, drive_fna = self.to_local_and_drive(found_fna)
        local_gff, drive_gff = self.to_local_and_drive(found_gff)
        local_lift, drive_lift = self.to_local_and_drive(found_lift)
        paths = canonical_output_paths(self.genomes_dir, record.accession)
        record.accession_root = str(paths["root"])
        record.assembly_folder = str(paths["root"])
        record.source_metadata_path = str(paths["source_metadata"])
        record.path_to_fna = local_fna
        record.path_to_gff = local_gff
        record.path_to_lift_gff = local_lift
        record.drive_to_fna = drive_fna
        record.drive_to_gff = drive_gff
        record.drive_to_lift_gff = drive_lift

    # Retained name; accepts normalized objects rather than raw repository JSON.
    def build_assembly_record(self, entry: Any) -> AssemblyRecord:
        if isinstance(entry, FetchedAssembly):
            record = entry.record
        elif isinstance(entry, AssemblyRecord):
            record = entry
        elif isinstance(entry, dict) and "normalized_record" in entry:
            record = AssemblyRecord.from_mapping(entry["normalized_record"])
        elif isinstance(entry, dict):
            source = detect_source(first_nonempty(entry, "accession"))
            adapter = self.get_adapter(source)
            accession = first_nonempty(entry, "accession", "assemblyAccession")
            record = adapter.normalize_cached(
                entry,
                accession=accession,
                source_query=self.search_term,
            ).record
        else:
            raise TypeError(f"Unsupported assembly entry type: {type(entry).__name__}")
        self.enrich_record_from_local_files(record)
        self.enrich_record_lineage(record)
        return record

    def append_update_event(
        self,
        species_key: str,
        organism_name: str,
        old_accession: str,
        new_accession: str,
        reason: str,
    ) -> None:
        append_json_line(
            self.update_history_jsonl,
            {
                "timestamp": now_iso(),
                "species_key": species_key,
                "organism_name": organism_name,
                "old_accession": old_accession,
                "new_accession": new_accession,
                "reason": reason,
            },
        )

    def _upgrade_existing_row(self, row: dict[str, Any]) -> dict[str, Any]:
        upgraded = dict(row)
        accession = str(upgraded.get("accession", "") or "")
        repository = str(upgraded.get("repository", "") or "") or detect_source(accession)
        upgraded["repository"] = repository
        if not upgraded.get("source_database"):
            if repository == "gwh":
                upgraded["source_database"] = "CNCB-GWH"
            elif accession.startswith("GCF_"):
                upgraded["source_database"] = "RefSeq"
            else:
                upgraded["source_database"] = "GenBank"
        if not upgraded.get("source_query"):
            upgraded["source_query"] = upgraded.get("ncbi_search_term", "")
        if upgraded.get("is_preferred_for_species", "") == "":
            upgraded["is_preferred_for_species"] = upgraded.get(
                "is_current_for_species", "True"
            )
        if upgraded.get("is_current_for_species", "") == "":
            upgraded["is_current_for_species"] = upgraded.get(
                "is_preferred_for_species", "True"
            )
        paths = canonical_output_paths(self.genomes_dir, accession)
        upgraded.setdefault("source_metadata_path", str(paths["source_metadata"]))
        upgraded.setdefault("assembly_folder", str(paths["root"]))
        upgraded.setdefault("accession_root", str(paths["root"]))
        return upgraded

    def update_metadata_tables(self) -> None:
        existing_rows = [
            self._upgrade_existing_row(row)
            for row in read_csv_rows(self.genomes_metadata_csv)
        ]
        rows_by_accession = {
            str(row.get("accession", "")): row
            for row in existing_rows
            if row.get("accession")
        }

        records_to_add: list[AssemblyRecord] = []
        for assembly in self.data:
            record = assembly.record
            if record.accession not in self.successfully_downloaded_data:
                continue
            self.enrich_record_from_local_files(record)
            self.enrich_record_lineage(record)
            prior = rows_by_accession.get(record.accession)
            if prior and prior.get("downloaded_at"):
                record.downloaded_at = str(prior["downloaded_at"])
            records_to_add.append(record)

        combined_records: list[AssemblyRecord] = []
        incoming_by_accession = {record.accession: record for record in records_to_add}
        for row in existing_rows:
            accession = str(row.get("accession", ""))
            if accession in incoming_by_accession:
                combined_records.append(incoming_by_accession.pop(accession))
            else:
                combined_records.append(AssemblyRecord.from_mapping(row))
        combined_records.extend(incoming_by_accession.values())

        # New assemblies for an already represented species are not preferred
        # automatically. Existing metadata controls the selection.
        existing_species = {
            str(row.get("species_key", ""))
            for row in existing_rows
            if row.get("species_key")
        }
        for record in records_to_add:
            if record.species_key in existing_species and record.accession not in rows_by_accession:
                record.is_preferred_for_species = False
                record.is_current_for_species = False
                if not record.preference_reason:
                    record.preference_reason = (
                        "Not automatically preferred because another independent "
                        "assembly for this species was already stored"
                    )

        apply_species_preference_semantics(
            combined_records,
            existing_rows=existing_rows,
            logger=self.log,
        )
        rows = [asdict(record) for record in combined_records]
        rows.sort(key=lambda row: str(row.get("accession", "")))
        write_csv_rows(
            self.genomes_metadata_csv,
            [field.name for field in fields(AssemblyRecord)],
            rows,
        )
        self.update_species_index()

    def update_species_index(self) -> None:
        metadata_rows = read_csv_rows(self.genomes_metadata_csv)
        species_rows: list[dict[str, Any]] = []
        for row in metadata_rows:
            preferred = parse_bool(row.get("is_preferred_for_species"))
            replaced = str(row.get("replaced_by_accession", "") or "")
            if replaced:
                status = "repository-superseded"
            elif preferred:
                status = "preferred"
            else:
                status = "available"
            species_rows.append(
                {
                    "species_key": row.get("species_key", ""),
                    "organism_name": row.get("organism_name", ""),
                    "accession": row.get("accession", ""),
                    "accession_root": row.get("accession_root", ""),
                    "repository": row.get("repository", ""),
                    "assembly_name": row.get("assembly_name", ""),
                    "source_database": row.get("source_database", ""),
                    "assembly_level": row.get("assembly_level", ""),
                    "refseq_category": row.get("refseq_category", ""),
                    "downloaded_at": row.get("downloaded_at", ""),
                    "status": status,
                    "is_preferred_for_species": row.get(
                        "is_preferred_for_species", ""
                    ),
                    "preferred_over_accession": row.get(
                        "preferred_over_accession", ""
                    ),
                    "preference_reason": row.get("preference_reason", ""),
                    "replaced_by_accession": replaced,
                    "project_count": "0",
                    "notes": row.get("notes", ""),
                }
            )
        write_csv_rows(
            self.species_assembly_index,
            [
                "species_key",
                "organism_name",
                "accession",
                "accession_root",
                "repository",
                "assembly_name",
                "source_database",
                "assembly_level",
                "refseq_category",
                "downloaded_at",
                "status",
                "is_preferred_for_species",
                "preferred_over_accession",
                "preference_reason",
                "replaced_by_accession",
                "project_count",
                "notes",
            ],
            species_rows,
        )

    def append_successful_entries_to_catalog(self) -> None:
        if not self.successfully_downloaded_data:
            return
        successful = set(self.successfully_downloaded_data)
        existing = self._existing_catalog_accessions()
        for assembly in self.data:
            if assembly.record.accession not in successful:
                continue
            if assembly.record.accession in existing:
                continue
            append_json_line(
                self.genomes_catalog_jsonl,
                {
                    "schema_version": SOURCE_METADATA_SCHEMA_VERSION,
                    "repository": assembly.record.repository,
                    "ingest_source": assembly.source_adapter or assembly.record.repository,
                    "accession": assembly.record.accession,
                    "source_query": assembly.record.source_query,
                    "raw_metadata": assembly.raw_metadata,
                    "download_urls": assembly.download_urls,
                    "normalized_record": asdict(assembly.record),
                },
            )

    def project_manifest_path(self) -> Optional[Path]:
        if not self.project_name:
            return None
        return self.project_manifests_dir / f"{self.project_name}_manifest.csv"

    def update_project_registry(self, manifest_path: Path) -> None:
        rows = read_csv_rows(self.project_registry_csv)
        by_name = {
            row.get("project_name", ""): row
            for row in rows
            if row.get("project_name")
        }
        new_row = {
            "project_name": self.project_name or "",
            "project_root": self.project_root,
            "manifest_path": str(manifest_path),
            "policy": self.project_policy,
            "freeze_date": datetime.date.today().isoformat(),
            "notes": self.manifest_notes,
        }
        if self.project_name in by_name:
            by_name[self.project_name].update(new_row)
        else:
            rows.append(new_row)
        write_csv_rows(
            self.project_registry_csv,
            [
                "project_name",
                "project_root",
                "manifest_path",
                "policy",
                "freeze_date",
                "notes",
            ],
            rows,
        )

    def project_manifest_path_or_raise(self) -> Path:
        path = self.project_manifest_path()
        if path is None:
            raise RuntimeError("No project name was supplied")
        return path

    def load_accessions_from_csv(self, path: Path) -> list[str]:
        return [
            row.get("accession", "").strip()
            for row in read_csv_rows(path)
            if row.get("accession", "").strip()
        ]

    def load_accession_from_txt(self, path: Path) -> list[str]:
        accessions: list[str] = []
        with path.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                accession = raw_line.strip()
                if accession and not accession.startswith("#"):
                    accessions.append(accession)
        return accessions

    def normalize_accession_list(self, items: list[str]) -> list[str]:
        seen: set[str] = set()
        cleaned: list[str] = []
        for item in items:
            accession = item.strip()
            if accession and accession not in seen:
                seen.add(accession)
                cleaned.append(accession)
        return cleaned

    def resolve_project_accessions(self) -> list[str]:
        selected: list[str] = []
        if self.project_accessions_csv:
            selected.extend(
                self.load_accessions_from_csv(Path(self.project_accessions_csv))
            )
        if self.project_accessions_txt:
            selected.extend(
                self.load_accession_from_txt(Path(self.project_accessions_txt))
            )
        selected.extend(self.project_accessions)
        selected = self.normalize_accession_list(selected)
        if not selected:
            selected = list(self.successfully_downloaded_data)
        self.project_selected_accessions = selected
        return selected

    def write_project_manifest(self) -> None:
        manifest_path = self.project_manifest_path()
        if manifest_path is None:
            return
        selected_accessions = set(self.resolve_project_accessions())
        if not selected_accessions:
            self.log("No accessions selected for the project manifest.\n")
            return
        missing = [
            accession
            for accession in sorted(selected_accessions)
            if not self.accession_root(accession).exists()
        ]
        if missing:
            self.log(
                "WARNING: project manifest includes accessions not found on disk:\n"
                + "\n".join(missing)
                + "\n"
            )

        selected_rows: list[dict[str, Any]] = []
        for row in read_csv_rows(self.genomes_metadata_csv):
            accession = row.get("accession", "")
            if accession not in selected_accessions:
                continue
            selected_rows.append(
                {
                    "project_name": self.project_name,
                    "project_policy": self.project_policy,
                    "manifest_version": "2",
                    "manifest_date": datetime.date.today().isoformat(),
                    "species_key": row.get("species_key", ""),
                    "organism_name": row.get("organism_name", ""),
                    "accession": accession,
                    "accession_root": row.get("accession_root", ""),
                    "assembly_name": row.get("assembly_name", ""),
                    "project_role": "analysis",
                    "repository": row.get("repository", ""),
                    "source_database": row.get("source_database", ""),
                    "source_record_url": row.get("source_record_url", ""),
                    "is_preferred_for_species_at_freeze": row.get(
                        "is_preferred_for_species", ""
                    ),
                    # Backward-compatible manifest column.
                    "is_current_for_species_at_freeze": row.get(
                        "is_current_for_species", ""
                    ),
                    "analysis_notes": self.manifest_notes,
                }
            )
        write_csv_rows(
            manifest_path,
            [
                "project_name",
                "project_policy",
                "manifest_version",
                "manifest_date",
                "species_key",
                "organism_name",
                "accession",
                "accession_root",
                "assembly_name",
                "project_role",
                "repository",
                "source_database",
                "source_record_url",
                "is_preferred_for_species_at_freeze",
                "is_current_for_species_at_freeze",
                "analysis_notes",
            ],
            selected_rows,
        )
        self.update_project_registry(manifest_path)

    def looks_like_accession_dir(self, path: Path) -> bool:
        """Recognize known accession prefixes or local-import cache directories."""

        if not path.is_dir():
            return False
        if is_accession_term(path.name):
            return True
        return (path / "source_metadata.json").is_file()

    def discover_local_accessions(self) -> list[str]:
        accessions: list[str] = []
        if not self.genomes_dir.exists():
            return accessions
        for child in sorted(self.genomes_dir.iterdir()):
            if child.name in {"records", "trash"}:
                continue
            if self.looks_like_accession_dir(child):
                accessions.append(child.name)
        return accessions

    def load_catalog_entries_by_accession(self) -> dict[str, dict[str, Any]]:
        entries: dict[str, dict[str, Any]] = {}
        for row in iter_json_lines(self.genomes_catalog_jsonl):
            accession = first_nonempty(
                row,
                "accession",
                "normalized_record.accession",
            )
            if accession:
                entries[accession] = row
        return entries

    def load_cached_assembly(self, accession: str) -> Optional[FetchedAssembly]:
        """Normalize source_metadata.json without contacting a repository."""

        path = self.source_metadata_file(accession)
        if not path.exists():
            return None
        try:
            cache = read_json(path)
        except Exception as exc:
            self.log(f"Could not read cached metadata {path}: {exc}\n")
            return None

        repository = str(cache.get("repository", "") or "") or detect_source(accession)
        ingest_source = str(cache.get("ingest_source", "") or repository)
        cached_normalized = cache.get("normalized_record")
        if ingest_source == "local" and isinstance(cached_normalized, dict):
            record = AssemblyRecord.from_mapping(cached_normalized)
            record.accession = accession
            self.enrich_record_from_local_files(record)
            self.log(
                f"Using cached local-import source_metadata.json for {accession}.\n"
            )
            return FetchedAssembly(
                record=record,
                raw_metadata=(
                    cache.get("raw_metadata")
                    if isinstance(cache.get("raw_metadata"), dict)
                    else {}
                ),
                download_urls={
                    str(key): str(value or "")
                    for key, value in (cache.get("download_urls") or {}).items()
                },
                source_adapter="local",
            )
        if repository not in self.adapters:
            self.log(
                f"Cached metadata for {accession} names unsupported repository "
                f"'{repository}'.\n"
            )
            return None
        raw = cache.get("raw_metadata")
        if not isinstance(raw, dict):
            # Older/newer caches may still contain a complete normalized row.
            normalized = cache.get("normalized_record")
            if isinstance(normalized, dict):
                record = AssemblyRecord.from_mapping(normalized)
                record.accession = accession
                return FetchedAssembly(
                    record=record,
                    raw_metadata={},
                    download_urls={
                        str(key): str(value)
                        for key, value in (cache.get("download_urls") or {}).items()
                    },
                    source_adapter=ingest_source,
                )
            self.log(f"Cached metadata for {accession} contains no usable record.\n")
            return None

        adapter = self.get_adapter(repository)
        try:
            assembly = adapter.normalize_cached(
                raw,
                accession=accession,
                source_query=str(cache.get("source_query", "") or "rebuild_from_cache"),
                source_record_url=str(cache.get("source_record_url", "") or ""),
            )
            if not assembly.source_adapter:
                assembly.source_adapter = ingest_source
            if isinstance(cache.get("download_urls"), dict):
                assembly.download_urls.update(
                    {
                        str(key): str(value or "")
                        for key, value in cache["download_urls"].items()
                    }
                )
            cached_normalized = cache.get("normalized_record")
            if isinstance(cached_normalized, dict):
                cached_record = AssemblyRecord.from_mapping(cached_normalized)
                for field in (
                    "downloaded_at",
                    "notes",
                    "preference_reason",
                    "preferred_over_accession",
                    "replaced_by_accession",
                ):
                    value = getattr(cached_record, field)
                    if value:
                        setattr(assembly.record, field, value)
            self.log(f"Using cached source_metadata.json for {accession}.\n")
            return assembly
        except Exception as exc:
            self.log(f"Cached metadata normalization failed for {accession}: {exc}\n")
            return None

    def _assembly_from_catalog_entry(
        self,
        accession: str,
        entry: dict[str, Any],
    ) -> Optional[FetchedAssembly]:
        repository = str(entry.get("repository", "") or "") or detect_source(accession)
        ingest_source = str(entry.get("ingest_source", "") or repository)
        normalized = entry.get("normalized_record")
        if ingest_source == "local" and isinstance(normalized, dict):
            record = AssemblyRecord.from_mapping(normalized)
            record.accession = accession
            return FetchedAssembly(
                record=record,
                raw_metadata=(
                    entry.get("raw_metadata")
                    if isinstance(entry.get("raw_metadata"), dict)
                    else {}
                ),
                download_urls={
                    str(key): str(value or "")
                    for key, value in (entry.get("download_urls") or {}).items()
                },
                source_adapter="local",
            )
        adapter = self.get_adapter(repository)
        raw = entry.get("raw_metadata")
        if not isinstance(raw, dict):
            # Legacy genomes_catalog.jsonl contained raw NCBI objects directly.
            if "organism" in entry and "accession" in entry:
                raw = entry
                repository = "ncbi"
                adapter = self.get_adapter(repository)
            else:
                return None
        try:
            return adapter.normalize_cached(
                raw,
                accession=accession,
                source_query=str(entry.get("source_query", "") or "rebuild_from_catalog"),
                source_record_url=str(entry.get("source_record_url", "") or ""),
            )
        except Exception as exc:
            self.log(f"Catalog metadata failed for {accession}: {exc}\n")
            return None

    def trash_dir(self) -> Path:
        path = self.genomes_dir / "trash"
        path.mkdir(parents=True, exist_ok=True)
        return path

    def move_accession_root_to_trash(self, accession: str, reason: str) -> Path:
        source = self.accession_root(accession)
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        destination = self.trash_dir() / f"{root_acc(accession)}__{timestamp}"
        if not source.exists():
            self.log(
                f"Requested trash move for missing accession root {source}: {reason}\n"
            )
            return destination
        shutil.move(str(source), str(destination))
        self.log(f"Moved {source} -> {destination}: {reason}\n")
        return destination

    def handle_broken_accession_root(self, accession: str, reason: str) -> None:
        if self.dry_run_rebuild:
            self.rebuild_preview["broken_accessions"].append(accession)
            self.rebuild_preview["would_trash"].append(
                {"accession": accession, "reason": reason}
            )
            self.log(
                f"DRY-RUN would move {accession} to trash because {reason}.\n"
            )
            return
        self.move_accession_root_to_trash(accession, reason)

    def backup_metadata(self) -> None:
        if not self.genomes_metadata_csv.exists():
            return
        if self.dry_run_rebuild:
            self.log("DRY-RUN would back up genomes_metadata.csv.\n")
            return
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        destination = self.trash_dir() / f"genomes_metadata__{timestamp}.csv"
        shutil.copy2(self.genomes_metadata_csv, destination)
        self.log(f"Backed up genomes_metadata.csv to {destination}.\n")

    def current_metadata_accessions(self) -> list[str]:
        return [
            row.get("accession", "").strip()
            for row in read_csv_rows(self.genomes_metadata_csv)
            if row.get("accession", "").strip()
        ]

    def summarize_rebuild_diff(self, rebuilt_records: list[dict[str, Any]]) -> None:
        current = set(self.current_metadata_accessions())
        rebuilt = {
            str(row.get("accession", "")).strip()
            for row in rebuilt_records
            if row.get("accession")
        }
        self.rebuild_preview["rebuilt_accessions"] = sorted(rebuilt)
        self.rebuild_preview["missing_from_rebuild"] = sorted(current - rebuilt)
        self.rebuild_preview["new_rows"] = sorted(rebuilt - current)

    def report_rebuild_preview(self) -> None:
        preview = self.rebuild_preview
        self.log("===== DRY-RUN REBUILD PREVIEW =====\n")
        self.log(
            f"Local accession directories found: {len(preview['local_accessions'])}\n"
        )
        self.log(f"Rows that would be written: {len(preview['rebuilt_accessions'])}\n")
        self.log(f"Broken accession roots: {len(preview['broken_accessions'])}\n")
        self.log(
            "Current metadata rows that would disappear: "
            f"{len(preview['missing_from_rebuild'])}\n"
        )
        self.log(f"Accessions newly added: {len(preview['new_rows'])}\n")
        for label, key in (
            ("Broken roots", "broken_accessions"),
            ("Rows missing from rebuild", "missing_from_rebuild"),
            ("New rows", "new_rows"),
        ):
            values = preview[key]
            if values:
                self.log(f"{label}:\n")
                for value in values[:25]:
                    self.log(f"  - {value}\n")
                if len(values) > 25:
                    self.log(f"  ... and {len(values) - 25} more\n")
        self.log("===== END DRY-RUN PREVIEW =====\n")

    def build_local_only_record(
        self,
        accession: str,
        *,
        existing_row: Optional[dict[str, Any]] = None,
    ) -> Optional[AssemblyRecord]:
        """Build a minimal record when cached and remote metadata are unavailable."""

        inspection = self.inspect_accession_root(accession)
        if inspection["is_broken"]:
            self.handle_broken_accession_root(
                accession,
                "the primary genome FASTA is missing during metadata rebuild",
            )
            return None

        if existing_row:
            record = AssemblyRecord.from_mapping(self._upgrade_existing_row(existing_row))
        else:
            repository = detect_source(accession)
            source_database = (
                "CNCB-GWH"
                if repository == "gwh"
                else ("RefSeq" if accession.startswith("GCF_") else "GenBank")
            )
            record = _base_record(
                accession=accession,
                repository=repository,
                source_database=source_database,
                source_record_url="",
                source_query="rebuild_from_local",
                organism_name="",
                taxon_id="",
                genomes_dir=self.genomes_dir,
            )
            record.species_key = root_acc(accession)
            record.notes = (
                "Rebuilt from local files without usable cached, catalog, or "
                "remote repository metadata"
            )
        self.enrich_record_from_local_files(record)
        return record

    def rebuild_genomes_metadata(self) -> None:
        """Rebuild normalized metadata, preferring local source metadata caches."""

        local_accessions = self.discover_local_accessions()
        self.rebuild_preview["local_accessions"] = list(local_accessions)
        catalog = self.load_catalog_entries_by_accession()
        existing_rows = [
            self._upgrade_existing_row(row)
            for row in read_csv_rows(self.genomes_metadata_csv)
        ]
        existing_by_accession = {
            str(row.get("accession", "")): row
            for row in existing_rows
            if row.get("accession")
        }
        rebuilt: list[AssemblyRecord] = []

        for accession in local_accessions:
            inspection = self.inspect_accession_root(accession)
            if inspection["is_broken"]:
                self.handle_broken_accession_root(
                    accession,
                    "the primary genome FASTA is missing during metadata rebuild",
                )
                continue

            assembly = self.load_cached_assembly(accession)
            if assembly is None and accession in catalog:
                assembly = self._assembly_from_catalog_entry(accession, catalog[accession])
                if assembly is not None:
                    self.log(f"Using genomes_catalog.jsonl metadata for {accession}.\n")

            if assembly is None:
                try:
                    assembly = self.fetch_summary_for_accession(accession)
                    self.log(f"Remote metadata fallback succeeded for {accession}.\n")
                except Exception as exc:
                    self.log(
                        f"WARNING: remote metadata fallback failed for {accession}: "
                        f"{exc}. Preserving the assembly with local-only metadata.\n"
                    )

            if assembly is not None:
                record = assembly.record
                prior = existing_by_accession.get(accession)
                if prior and prior.get("downloaded_at"):
                    record.downloaded_at = str(prior["downloaded_at"])
                self.enrich_record_from_local_files(record)
                self.enrich_record_lineage(record)
            else:
                record = self.build_local_only_record(
                    accession,
                    existing_row=existing_by_accession.get(accession),
                )
            if record is not None:
                rebuilt.append(record)

        apply_species_preference_semantics(
            rebuilt,
            existing_rows=existing_rows,
            logger=self.log,
        )
        rebuilt.sort(key=lambda record: record.accession)
        rebuilt_rows = [asdict(record) for record in rebuilt]
        self.summarize_rebuild_diff(rebuilt_rows)

        if self.dry_run_rebuild:
            self.report_rebuild_preview()
            return

        self.backup_metadata()
        write_csv_rows(
            self.genomes_metadata_csv,
            [field.name for field in fields(AssemblyRecord)],
            rebuilt_rows,
        )
        self.update_species_index()
        self.log(
            f"Rebuilt metadata from {len(local_accessions)} local accession "
            f"directories; wrote {len(rebuilt_rows)} rows.\n"
        )

    def run(self) -> None:
        self.initialize_records()

        if self.dry_run_rebuild:
            self.rebuild_genomes_metadata()
            return
        if self.rebuild_metadata:
            self.rebuild_genomes_metadata()
            if self.project_name:
                self.write_project_manifest()
            return

        if self.retry_failures:
            self.load_failed_accessions_for_retry()
            if not self.data:
                self.log("No failed accessions could be prepared for retry.\n")
                return
        else:
            self.query_source()
            self.log_species_with_multiple_accessions()

        if not self.data:
            self.log("No new assemblies require downloading.\n")
            if self.project_name:
                self.write_project_manifest()
            return

        self.download_new_data()
        if self.successfully_downloaded_data:
            self.append_successful_entries_to_catalog()
            self.update_metadata_tables()
        else:
            self.log(
                "No downloads completed successfully; metadata tables were not modified.\n"
            )

        if self.retry_failures:
            self.remove_resolved_failures(self.successfully_downloaded_data)
        if self.project_name:
            self.write_project_manifest()


# ---------------------------------------------------------------------------
# Built-in focused tests
# ---------------------------------------------------------------------------


def run_self_tests() -> int:
    """Run focused offline tests with mocked repository/network behavior."""

    import unittest
    from unittest import mock

    class GenomeManagerTests(unittest.TestCase):
        def setUp(self) -> None:
            self.tempdir = tempfile.TemporaryDirectory()
            self.root = Path(self.tempdir.name) / "genomes"
            self.root.mkdir()

        def tearDown(self) -> None:
            self.tempdir.cleanup()

        def test_source_detection(self) -> None:
            self.assertEqual(detect_source("GCA_012345678.1"), "ncbi")
            self.assertEqual(detect_source("GCF_012345678.1"), "ncbi")
            self.assertEqual(detect_source("GWHFILM00000000.1"), "gwh")
            self.assertEqual(detect_source("Squamata"), "ncbi")

        def test_explicit_override(self) -> None:
            self.assertEqual(resolve_source("GCA_1.1", "gwh"), "gwh")
            self.assertEqual(resolve_source("GWHAAAA00000000", "ncbi"), "ncbi")

        def test_accession_and_root(self) -> None:
            self.assertTrue(is_accession_term("GWHFILM00000000.1"))
            self.assertEqual(root_acc("GWHFILM00000000.1"), "GWHFILM00000000")
            self.assertEqual(root_acc("GCA_012345678.2"), "GCA_012345678")

        def test_assembly_level(self) -> None:
            self.assertEqual(
                normalize_assembly_level("Draft genome in chromosome level"),
                "chromosome",
            )
            self.assertEqual(normalize_assembly_level("Scaffold"), "scaffold")

        def test_gwh_normalization_and_no_gff(self) -> None:
            payload = {
                "assemblyAccession": "GWHFILM00000000",
                "organism": "Indotyphlops braminus",
                "taxId": "42160",
                "assemblyName": "test",
                "assemblyLevel": "Draft genome in chromosome level",
                "bioprojectAccession": "PRJCA030601",
                "biosampleAccession": "SAMC4572868",
                "releaseTime": "2025-01-16 00:00:00",
                "ftpPathDna": "https://example.org/genome.fasta.gz",
                "ftpPathGff": "",
            }
            record = normalize_gwh_entry(
                payload,
                "GWHFILM00000000.1",
                genomes_dir=self.root,
                source_query="GWHFILM00000000.1",
            )
            self.assertEqual(record.accession, "GWHFILM00000000.1")
            self.assertEqual(record.repository, "gwh")
            self.assertEqual(record.assembly_level, "chromosome")
            self.assertEqual(record.bioproject_accession, "PRJCA030601")
            self.assertEqual(gwh_download_urls(payload)["gff3"], "")

        def test_wrapped_gwh_response(self) -> None:
            wrapped = {
                "code": 200,
                "data": {
                    "assemblyAccession": "GWHFILM00000000",
                    "organism": "Indotyphlops braminus",
                    "ftpPathDna": "https://example.org/genome.fasta.gz",
                },
            }
            extracted = extract_gwh_assembly_payload(wrapped)
            self.assertIsNotNone(extracted)
            assert extracted is not None
            self.assertEqual(extracted["assemblyAccession"], "GWHFILM00000000")
            record = normalize_gwh_entry(
                wrapped,
                "GWHFILM00000000.1",
                genomes_dir=self.root,
                source_query="GWHFILM00000000.1",
            )
            self.assertEqual(record.organism_name, "Indotyphlops braminus")
            self.assertTrue(gwh_download_urls(wrapped)["genome"].endswith(".gz"))

        def test_canonical_paths_and_directory_recognition(self) -> None:
            paths = canonical_output_paths(self.root, "GWHFILM00000000.1")
            self.assertEqual(paths["fna"].name, "GWHFILM00000000.1_genomic.fna")
            manager = GenomeManagerHybrid(
                str(self.root),
                "rebuild_local",
                source="auto",
            )
            accession_dir = self.root / "GWHFILM00000000.1"
            accession_dir.mkdir()
            self.assertTrue(manager.looks_like_accession_dir(accession_dir))

        def test_missing_required_gwh_url(self) -> None:
            adapter = GWHSourceAdapter(
                genomes_dir=self.root,
                records_dir=self.root / "records",
                logger=lambda message: None,
            )
            record = AssemblyRecord(
                accession="GWHAAAA00000000.1",
                repository="gwh",
            )
            assembly = FetchedAssembly(record, {}, {"genome": "", "gff3": ""})
            with self.assertRaisesRegex(RuntimeError, "required"):
                adapter.download(assembly, self.root / record.accession)

        def test_gwh_version_fallback_is_mocked(self) -> None:
            calls: list[str] = []

            def fake_fetcher(url: str, **kwargs: Any) -> dict[str, Any]:
                del kwargs
                calls.append(url)
                if url.endswith(".1"):
                    raise RuntimeError("not found")
                return {"assemblyAccession": "GWHFILM00000000"}

            payload, url = fetch_gwh_assembly(
                "GWHFILM00000000.1",
                fetcher=fake_fetcher,
            )
            self.assertEqual(payload["assemblyAccession"], "GWHFILM00000000")
            self.assertTrue(url.endswith("GWHFILM00000000"))
            self.assertEqual(len(calls), 2)

        def test_cached_metadata_used_without_remote(self) -> None:
            manager = GenomeManagerHybrid(str(self.root), "rebuild_local")
            accession = "GWHFILM00000000.1"
            paths = canonical_output_paths(self.root, accession)
            paths["root"].mkdir()
            paths["fna"].write_text(">x\nACGT\n", encoding="utf-8")
            raw = {
                "assemblyAccession": "GWHFILM00000000",
                "organism": "Indotyphlops braminus",
                "taxId": "42160",
                "ftpPathDna": "https://example.org/a.gz",
            }
            write_json(
                paths["source_metadata"],
                {
                    "repository": "gwh",
                    "source_query": accession,
                    "raw_metadata": raw,
                    "download_urls": {"genome": "https://example.org/a.gz"},
                },
            )
            with mock.patch.object(
                manager,
                "fetch_summary_for_accession",
                side_effect=AssertionError("remote fetch should not be called"),
            ):
                cached = manager.load_cached_assembly(accession)
            self.assertIsNotNone(cached)
            assert cached is not None
            self.assertEqual(cached.record.organism_name, "Indotyphlops braminus")

        def test_local_import_end_to_end(self) -> None:
            input_dir = Path(self.tempdir.name) / "downloads"
            input_dir.mkdir()
            fasta = input_dir / "GWHFILM00000000.1.genome.fasta.gz"
            gff = input_dir / "GWHFILM00000000.1.gff.gz"
            with gzip.open(fasta, "wt", encoding="utf-8") as handle:
                handle.write(">chr1\nACGTACGT\n")
            with gzip.open(gff, "wt", encoding="utf-8") as handle:
                handle.write("##gff-version 3\nchr1\ttest\tregion\t1\t8\t.\t+\t.\tID=chr1\n")

            config = LocalImportConfig(
                repository="gwh",
                local_fasta=fasta,
                local_gff=gff,
                organism_name="Indotyphlops braminus",
                taxon_id="42172",
                assembly_level="chromosome",
                bioproject_accession="PRJCA030601",
                biosample_accession="SAMC4572868",
                raw_reads_accessions="CRR1549896",
                release_date="2025-01-16",
                assembly_method="hifiasm 0.19.9-r616",
            )
            manager = GenomeManagerHybrid(
                str(self.root),
                "GWHFILM00000000.1",
                source="local",
                local_import_config=config,
            )
            manager.run()

            paths = canonical_output_paths(self.root, "GWHFILM00000000.1")
            self.assertEqual(paths["fna"].read_text(encoding="utf-8"), ">chr1\nACGTACGT\n")
            self.assertTrue(paths["gff"].read_text(encoding="utf-8").startswith("##gff-version 3"))
            cache = read_json(paths["source_metadata"] )
            self.assertEqual(cache["repository"], "gwh")
            self.assertEqual(cache["ingest_source"], "local")
            self.assertIn("all-null/unusable", cache["normalized_record"]["notes"] )

            rows = read_csv_rows(self.root / "records" / "genomes_metadata.csv")
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(row["source_database"], "CNCB-GWH")
            self.assertEqual(row["organism_name"], "Indotyphlops braminus")
            self.assertEqual(row["assembly_method"], "hifiasm 0.19.9-r616")
            self.assertEqual(row["raw_reads_accessions"], "CRR1549896")

        def test_local_import_cli_configuration(self) -> None:
            args = parse_args([
                str(self.root),
                "GWHFILM00000000.1",
                "--source", "local",
                "--repository", "gwh",
                "--local-fasta", "~/Downloads/example.fasta.gz",
                "--organism-name", "Indotyphlops braminus",
                "--taxon-id", "42172",
                "--assembly-level", "chromosome",
                "--bioproject", "PRJCA030601",
                "--biosample", "SAMC4572868",
                "--raw-reads", "CRR1549896",
                "--release-date", "2025-01-16",
                "--assembly-method", "hifiasm 0.19.9-r616",
            ])
            config = build_local_import_config(args)
            self.assertIsNotNone(config)
            assert config is not None
            self.assertEqual(config.repository, "gwh")
            self.assertEqual(config.raw_reads_accessions, "CRR1549896")
            self.assertEqual(config.local_fasta.name, "example.fasta.gz")

        def test_gwh_recovers_only_missing_gff(self) -> None:
            accession = "GWHFILM00000000.1"
            paths = canonical_output_paths(self.root, accession)
            paths["root"].mkdir()
            paths["fna"].write_text(">chr1\nACGT\n", encoding="utf-8")
            original_fasta = paths["fna"].read_bytes()

            adapter = GWHSourceAdapter(
                genomes_dir=self.root,
                records_dir=self.root / "records",
                logger=lambda message: None,
            )
            assembly = FetchedAssembly(
                AssemblyRecord(accession=accession, repository="gwh"),
                {},
                {
                    "genome": "https://example.org/genome.fasta.gz",
                    "gff3": "https://example.org/genes.gff.gz",
                },
                source_adapter="gwh",
            )
            calls: list[tuple[str, Path]] = []

            def fake_download(url: str, destination: Path, **kwargs: Any) -> None:
                del kwargs
                calls.append((url, destination))
                destination.write_text(
                    "##gff-version 3\nchr1\ttest\tgene\t1\t4\t.\t+\t.\tID=g1\n",
                    encoding="utf-8",
                )

            with mock.patch(
                f"{__name__}.download_and_decompress",
                side_effect=fake_download,
            ):
                adapter.download(assembly, paths["root"])

            self.assertEqual(paths["fna"].read_bytes(), original_fasta)
            self.assertTrue(is_nonempty_file(paths["gff"]))
            self.assertEqual(calls, [(assembly.download_urls["gff3"], paths["gff"])])

        def test_ncbi_recovers_only_missing_gff(self) -> None:
            import zipfile

            accession = "GCF_000000001.1"
            paths = canonical_output_paths(self.root, accession)
            paths["root"].mkdir()
            paths["fna"].write_text(">chr1\nACGT\n", encoding="utf-8")
            original_fasta = paths["fna"].read_bytes()
            commands: list[list[str]] = []

            adapter = NCBISourceAdapter(
                genomes_dir=self.root,
                records_dir=self.root / "records",
                logger=lambda message: None,
            )
            assembly = FetchedAssembly(
                AssemblyRecord(accession=accession, repository="ncbi"),
                {},
                {},
                source_adapter="ncbi",
            )

            def fake_run(command: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
                del kwargs
                commands.append(command)
                zip_path = Path(command[command.index("--filename") + 1])
                with zipfile.ZipFile(zip_path, "w") as archive:
                    archive.writestr(
                        f"ncbi_dataset/data/{accession}/genomic.gff",
                        "##gff-version 3\nchr1\ttest\tgene\t1\t4\t.\t+\t.\tID=g1\n",
                    )
                return subprocess.CompletedProcess(command, 0, "", "")

            with mock.patch.object(adapter, "require_dependencies"), mock.patch(
                f"{__name__}.subprocess.run", side_effect=fake_run
            ):
                adapter.download(assembly, paths["root"])

            self.assertEqual(paths["fna"].read_bytes(), original_fasta)
            self.assertTrue(is_nonempty_file(paths["gff"]))
            self.assertEqual(len(commands), 1)
            include_value = commands[0][commands[0].index("--include") + 1]
            self.assertEqual(include_value, "gff3")

        def test_ncbi_missing_repository_gff_is_not_failure(self) -> None:
            import zipfile

            accession = "GCF_000000002.1"
            paths = canonical_output_paths(self.root, accession)
            paths["root"].mkdir()
            paths["fna"].write_text(">chr1\nACGT\n", encoding="utf-8")
            logs: list[str] = []

            adapter = NCBISourceAdapter(
                genomes_dir=self.root,
                records_dir=self.root / "records",
                logger=logs.append,
            )
            assembly = FetchedAssembly(
                AssemblyRecord(accession=accession, repository="ncbi"),
                {},
                {},
                source_adapter="ncbi",
            )

            def fake_run(command: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
                del kwargs
                zip_path = Path(command[command.index("--filename") + 1])
                with zipfile.ZipFile(zip_path, "w") as archive:
                    archive.writestr("README.md", "No annotation available\n")
                return subprocess.CompletedProcess(command, 0, "", "")

            with mock.patch.object(adapter, "require_dependencies"), mock.patch(
                f"{__name__}.subprocess.run", side_effect=fake_run
            ):
                adapter.download(assembly, paths["root"])

            self.assertTrue(is_nonempty_file(paths["fna"]))
            self.assertFalse(paths["gff"].exists())
            self.assertTrue(
                any("provides no repository GFF" in message for message in logs)
            )

        def test_catalogued_accession_is_still_reconciled(self) -> None:
            manager = GenomeManagerHybrid(str(self.root), "GCF_000000003.1")
            manager.initialize_records()
            accession = "GCF_000000003.1"
            append_json_line(
                manager.genomes_catalog_jsonl,
                {"accession": accession},
            )
            fetched = FetchedAssembly(
                AssemblyRecord(accession=accession, repository="ncbi"),
                {},
                {},
                source_adapter="ncbi",
            )
            with mock.patch.object(
                manager.get_adapter("ncbi"), "fetch", return_value=[fetched]
            ):
                manager.query_source()
            self.assertEqual([item.record.accession for item in manager.data], [accession])

        def test_recovered_gff_refreshes_metadata_path(self) -> None:
            accession = "GWHMETA00000000.1"
            manager = GenomeManagerHybrid(str(self.root), accession, source="gwh")
            manager.initialize_records()
            paths = canonical_output_paths(self.root, accession)
            paths["root"].mkdir()
            paths["fna"].write_text(">chr1\nACGT\n", encoding="utf-8")

            record = AssemblyRecord(
                accession=accession,
                repository="gwh",
                source_database="CNCB-GWH",
                organism_name="Example species",
                species_key="Example_species",
            )
            assembly = FetchedAssembly(
                record,
                {"assemblyAccession": root_acc(accession)},
                {
                    "genome": "https://example.org/genome.fasta.gz",
                    "gff3": "https://example.org/genes.gff.gz",
                },
                source_adapter="gwh",
            )
            manager.data = [assembly]

            def fake_download(item: FetchedAssembly, destination: Path) -> None:
                del item, destination
                paths["gff"].write_text(
                    "##gff-version 3\nchr1\ttest\tgene\t1\t4\t.\t+\t.\tID=g1\n",
                    encoding="utf-8",
                )

            with mock.patch.object(
                manager.get_adapter("gwh"), "download", side_effect=fake_download
            ), mock.patch.object(manager, "enrich_record_lineage"):
                manager.download_new_data()
                manager.update_metadata_tables()

            rows = read_csv_rows(manager.genomes_metadata_csv)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["path_to_gff"], str(paths["gff"]))
            self.assertTrue(is_nonempty_file(Path(rows[0]["path_to_gff"])))

        def test_multiple_assemblies_do_not_replace_each_other(self) -> None:
            first = AssemblyRecord(
                accession="GCA_000000001.1",
                repository="ncbi",
                organism_name="Example species",
                species_key="Example_species",
                is_preferred_for_species=True,
            )
            second = AssemblyRecord(
                accession="GWHAAAA00000000.1",
                repository="gwh",
                organism_name="Example species",
                species_key="Example_species",
                is_preferred_for_species=False,
            )
            apply_species_preference_semantics(
                [first, second],
                existing_rows=[
                    {
                        "accession": first.accession,
                        "is_preferred_for_species": "True",
                    }
                ],
            )
            self.assertTrue(first.is_preferred_for_species)
            self.assertFalse(second.is_preferred_for_species)
            self.assertEqual(first.replaced_by_accession, "")
            self.assertEqual(second.replaced_by_accession, "")

    suite = unittest.defaultTestLoader.loadTestsFromTestCase(GenomeManagerTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Manage NCBI, CNCB-GWH, and manually downloaded genome assemblies "
            "using normalized metadata, canonical local files, and project manifests."
        )
    )
    parser.add_argument(
        "genomes_dir",
        nargs="?",
        help="Path to the genomes directory (its final component must be 'genomes')",
    )
    parser.add_argument(
        "search_term",
        nargs="?",
        default="rebuild_local",
        help="NCBI taxon name or an assembly accession/identifier",
    )
    parser.add_argument(
        "--source",
        choices=["auto", "ncbi", "gwh", "local"],
        default="auto",
        help=(
            "Ingest adapter. auto selects NCBI for GCA_/GCF_/taxon terms and "
            "GWH for GWH accessions; local imports user-downloaded files."
        ),
    )
    local_group = parser.add_argument_group(
        "repository-independent local import"
    )
    local_group.add_argument(
        "--repository",
        default=None,
        help=(
            "Provenance repository for --source local, for example gwh or ncbi. "
            "This is stored in metadata; it is not the transfer adapter."
        ),
    )
    local_group.add_argument(
        "--local-fasta",
        default=None,
        help="Required local genome FASTA for --source local; .gz is supported",
    )
    local_group.add_argument(
        "--local-gff",
        default=None,
        help="Optional local GFF/GFF3 for --source local; .gz is supported",
    )
    local_group.add_argument("--organism-name", default=None)
    local_group.add_argument("--taxon-id", default=None)
    local_group.add_argument("--assembly-level", default=None)
    local_group.add_argument("--bioproject", default=None)
    local_group.add_argument("--biosample", default=None)
    local_group.add_argument("--raw-reads", default=None)
    local_group.add_argument("--release-date", default=None)
    local_group.add_argument("--assembly-method", default=None)
    local_group.add_argument("--assembly-name", default=None)
    local_group.add_argument("--sequencing-and-coverage", default=None)
    local_group.add_argument("--source-record-url", default=None)
    local_group.add_argument(
        "--local-import-notes",
        default=None,
        help="Additional provenance notes for a local import",
    )
    parser.add_argument("--project-name", default=None)
    parser.add_argument(
        "--project-policy",
        choices=["frozen", "rolling", "exploratory"],
        default="exploratory",
    )
    parser.add_argument("--project-root", default=None)
    parser.add_argument("--manifest-notes", default="")
    parser.add_argument(
        "--retry-failures",
        action="store_true",
        help="Retry accessions in records/download_failures.csv",
    )
    parser.add_argument(
        "--rebuild-metadata",
        action="store_true",
        help="Rebuild metadata from local GCA_, GCF_, and GWH directories",
    )
    parser.add_argument(
        "--dry-run-rebuild",
        action="store_true",
        help="Preview metadata rebuilding without writing or moving directories",
    )
    parser.add_argument("--project-accessions-csv", default=None)
    parser.add_argument("--project-accessions-txt", default=None)
    parser.add_argument("--project-accessions", nargs="*", default=None)
    parser.add_argument(
        "--connect-timeout",
        type=float,
        default=DEFAULT_CONNECT_TIMEOUT,
        help="Connection timeout in seconds for HTTPS operations",
    )
    parser.add_argument(
        "--read-timeout",
        type=float,
        default=DEFAULT_READ_TIMEOUT,
        help="Read timeout in seconds for HTTPS operations",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run focused offline unit tests and exit",
    )
    parser.add_argument("--version", action="version", version=SCRIPT_VERSION)
    return parser.parse_args(argv)


def build_local_import_config(
    args: argparse.Namespace,
) -> Optional[LocalImportConfig]:
    """Validate CLI-local import fields and construct a typed configuration."""

    option_names = (
        "repository",
        "local_fasta",
        "local_gff",
        "organism_name",
        "taxon_id",
        "assembly_level",
        "bioproject",
        "biosample",
        "raw_reads",
        "release_date",
        "assembly_method",
        "assembly_name",
        "sequencing_and_coverage",
        "source_record_url",
        "local_import_notes",
    )
    supplied = [
        name for name in option_names
        if getattr(args, name, None) not in (None, "")
    ]
    if args.source != "local":
        if supplied:
            rendered = ", ".join(f"--{name.replace('_', '-')}" for name in supplied)
            raise ValueError(
                f"Local-import arguments ({rendered}) require --source local"
            )
        return None

    # Rebuild modes use saved source_metadata.json and do not need the original
    # local paths or user-supplied metadata again.
    if args.rebuild_metadata or args.dry_run_rebuild:
        return None

    required = {
        "repository": args.repository,
        "local-fasta": args.local_fasta,
        "organism-name": args.organism_name,
        "assembly-level": args.assembly_level,
    }
    missing = [name for name, value in required.items() if not str(value or "").strip()]
    if missing:
        raise ValueError(
            "--source local requires: "
            + ", ".join(f"--{name}" for name in missing)
        )
    if not args.search_term or args.search_term == "rebuild_local":
        raise ValueError(
            "--source local requires the positional assembly accession/identifier"
        )

    return LocalImportConfig(
        repository=str(args.repository),
        local_fasta=Path(str(args.local_fasta)).expanduser(),
        local_gff=(
            Path(str(args.local_gff)).expanduser()
            if args.local_gff
            else None
        ),
        organism_name=str(args.organism_name or ""),
        taxon_id=str(args.taxon_id or ""),
        assembly_level=str(args.assembly_level or ""),
        bioproject_accession=str(args.bioproject or ""),
        biosample_accession=str(args.biosample or ""),
        raw_reads_accessions=str(args.raw_reads or ""),
        release_date=str(args.release_date or ""),
        assembly_method=str(args.assembly_method or ""),
        assembly_name=str(args.assembly_name or ""),
        sequencing_and_coverage=str(args.sequencing_and_coverage or ""),
        source_record_url=str(args.source_record_url or ""),
        notes=str(args.local_import_notes or ""),
    )


def main(argv: Optional[list[str]] = None) -> None:
    args = parse_args(argv)
    if args.self_test:
        raise SystemExit(run_self_tests())
    if not args.genomes_dir:
        raise ValueError("genomes_dir is required unless --self-test is used")
    if args.rebuild_metadata and args.dry_run_rebuild:
        raise ValueError("Use either --rebuild-metadata or --dry-run-rebuild, not both")
    if args.connect_timeout <= 0 or args.read_timeout <= 0:
        raise ValueError("Timeout values must be greater than zero")

    local_import_config = build_local_import_config(args)

    manager = GenomeManagerHybrid(
        genomes_dir=args.genomes_dir,
        search_term=args.search_term,
        source=args.source,
        project_name=args.project_name,
        project_policy=args.project_policy,
        project_root=args.project_root,
        manifest_notes=args.manifest_notes,
        retry_failures=args.retry_failures,
        rebuild_metadata=args.rebuild_metadata,
        dry_run_rebuild=args.dry_run_rebuild,
        project_accessions_csv=args.project_accessions_csv,
        project_accessions_txt=args.project_accessions_txt,
        project_accessions=args.project_accessions,
        connect_timeout=args.connect_timeout,
        read_timeout=args.read_timeout,
        local_import_config=local_import_config,
    )
    manager.run()


if __name__ == "__main__":
    main()
