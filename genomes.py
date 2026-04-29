#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional

import jsonlines

def require_genomes_dir(path: Path) -> Path:
    path = path.resolve()
    if path.name != "genomes":
        raise ValueError(f"Expected a directory ending with 'genomes', got: {path}")
    return path

def is_accession_term(term: str) -> bool:
    return term.startswith(("GCA_", "GCF_"))
    
def safe_tag(text: str) -> str:
    return text.replace("/", "_").replace("\\", "_").replace(" ", "_")

def root_acc(accession: str) -> str:
    return accession.split(".")[0] if accession else accession

def build_species_key(organism_name: str) -> str:
    return "_".join(organism_name.strip().split())
    
def now_iso() -> str:
    return datetime.datetime.now().isoformat(timespec="seconds")

def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with open(path, "r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))

def write_csv_rows(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

@dataclass
class AssemblyRecord:
    accession: str
    accession_root: str
    organism_name: str
    taxon_id: str
    assembly_name: str
    assembly_level: str
    refseq_category: str
    source_database: str
    is_refseq: bool
    is_current_for_species: bool
    replaced_by_accession: str
    species_key: str
    phylum: str
    superorder: str
    order: str
    family: str
    genus: str
    genus_species: str
    path_to_fna: str
    path_to_gff: str
    path_to_lift_gff: str
    drive_to_fna: str
    drive_to_gff: str
    drive_to_lift_gff: str
    assembly_folder: str
    downloaded_at: str
    metadata_last_updated: str
    ncbi_search_term: str
    notes: str
    
class GenomeManagerHybrid:
    def __init__(
        self,
        genomes_dir: str,
        search_term: str,
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
    ):
        self.genomes_dir = require_genomes_dir(Path(genomes_dir))
        self.records_dir = self.genomes_dir / "records"
        self.project_manifests_dir = self.records_dir / "project_manifests"
        
        self.search_term = search_term.strip()
        self.is_accession = is_accession_term(self.search_term)
        # what is safe_tag used for?
        self.safe_tag = safe_tag(self.search_term) # is this correct? if the search term is squamate all species have the same safe_tag?
        
        self.project_name = project_name
        self.project_policy = project_policy
        self.project_root = str(Path(project_root).resolve()) if project_root else ""
        self.manifest_notes = manifest_notes
        self.retry_failures = retry_failures
        self.rebuild_metadata = rebuild_metadata
        self.dry_run_rebuild = dry_run_rebuild
        self.project_accessions_csv = project_accessions_csv
        self.project_accessions_txt = project_accessions_txt
        self.project_accessions = project_accessions or []
        self.on_colab = ("google.colab" in sys.modules) # I don't understand what this is doing
        self.root_local = self.genomes_dir.parent
        self.root_drive = Path("/content/drive/Othercomputers/macbook/projects")
        
        self.genomes_catalog_jsonl = self.records_dir / "genomes_catalog.jsonl"
        self.genomes_metadata_csv = self.records_dir / "genomes_metadata.csv"
        self.species_assembly_index = self.records_dir /  "species_assembly_index.csv"
        self.update_history_jsonl = self.records_dir / "assembly_update_history.jsonl"
        self.project_registry_csv = self.records_dir / "project_registry.csv"
        self.records_log_txt = self.records_dir / "genome_auto_update_records_log.txt"
        self.updates_jsonl = self.records_dir / f"{self.safe_tag}_updates.jsonl"
        self.download_failures_csv = self.records_dir / "download_failures.csv"
        
        self.data: list[dict] = []
        self.successfully_downloaded_data: list[str] = []
        self.failed_downloads: list[dict[str,str]]=[]
        self.project_selected_accessions: list[str] = []
        self.rebuild_preview: dict = {
            "local_accessions": [],
            "rebuilt_accessions": [],
            "broken_accessions": [],
            "would_trash": [],
            "missing_from_rebuild": [],
            "new_rows": [],
            }
        
        
    def initialize_records(self) -> None:
        self.records_dir.mkdir(parents=True, exist_ok=True)
        self.project_manifests_dir.mkdir(parents=True, exist_ok=True)
        self.genomes_catalog_jsonl.touch(exist_ok=True)
        self.records_log_txt.touch(exist_ok=True)

    def log(self, message: str) -> None:
        with open(self.records_log_txt, "a", encoding="utf-8") as handle:
            handle.write(message)

    def require_datasets_cli(self) -> None:
        exit_val = os.system("datasets --help > /dev/null 2>&1")
        if exit_val != 0:
            self.log(f"\n\n\n{now_iso()} ***ERROR*** datasets CLI not found in PATH\n")
            raise RuntimeError("NCBI datasets CLI not found in PATH") # this is better than sys.exit - which is considered unnecessary/extreme
        self.log(f"\n\n\n{now_iso()} datasets CLI detected\n")
    
    def build_ncbi_summary_command(self) -> str:
        if self.is_accession:
            query_bit = f"accession {shlex.quote(self.search_term)}"
        else:
            query_bit = f"taxon {shlex.quote(self.search_term)}"
       
        return (
            f"datasets summary genome {query_bit} "
            f"--assembly-level chromosome --assembly-level complete "
            f"--as-json-lines > {shlex.quote(str(self.updates_jsonl))}" # store the summary info in the search_term_updates_jsonl
            )
    
    def query_ncbi(self) -> None:
        cmd = self.build_ncbi_summary_command()
        self.log(f"\nChecking NCBI with command:\n{cmd}\n")
        
        os.chdir(self.records_dir)
        exit_val = os.system(cmd)
        if exit_val != 0:
            self.log(f"\n***ERROR*** failed to retrieve NCBI updates: {self.updates_jsonl}\n")
            raise RuntimeError("Failed to retrieve updates from NCBI")
            
        self.log("\n--Successfully retrieved updated information from NCBI--\n")
        self._compact_jsonl(self.updates_jsonl)
        self._load_new_entries_from_updates()

    def _compact_jsonl(self, path: Path) -> None:
        tmp = path.with_suffix(".compact.jsonl")
        with jsonlines.open(path, mode = "r") as reader, jsonlines.open(tmp, mode="w", compact=True) as writer:
            for line in reader:
                writer.write(line)
        os.replace(tmp, path) # used to rename a file or dir from src (tmp) -> dst (path), replaces dst file if it exists
    
    def _load_new_entries_from_updates(self) -> None:
        try:
            with open(self.genomes_catalog_jsonl, "r", encoding="utf-8") as handle:
                local_text = handle.read()
        except FileNotFoundError:
            local_text = ""
        
        self.data = [] # this will hold the accessions returned by the search term summary command
        with open(self.updates_jsonl, "r", encoding="utf-8") as handle:
            for raw_line in handle:
                if not raw_line.strip(): # if this is a blank line
                    continue # continue onto the next
                row = json.loads(raw_line)
                acc = row.get("accession", "")
                if acc and acc in local_text: # prevent local genomes from being in self.data
                    continue
                self.data.append(row)
                
    def log_species_with_multiple_accessions(self) -> None:
        seen: dict[str, list[str]] = {} # this will be a dict that is able to keep track of all accessions in the given search term. key[value[key[value]]] organism[value[organism_name[""]]]
        for genome in self.data:
            organism = genome.get("organism", {}).get("organism_name", "")
            acc = genome.get("accession", "")
            seen.setdefault(organism, []).append(acc) # dict.setdefault() - give me a value for this key (organism), if the key isn't there add it first so I don't get an error
        #self.log(
        #    f"\n\n\nTROUBLESHOOTING --dictionary seen:\n"
        #    f"{seen}"
        #)
        for organism, accessions in seen.items():
            if len(accessions) > 1:
                self.log(
                    f"Multiple accessions returned for organism {organism}: "
                    f"{', '.join(accessions)}\n"
                )
    
    def accession_root(self, accession: str) -> Path:
        return self.genomes_dir / accession
    
    def download_new_data(self) -> None:
        for entry in self.data:
            accession = entry.get("accession", "")
            if not accession:
                continue
            
            acc_root = self.accession_root(accession)
            acc_root.mkdir(parents=True, exist_ok=True)
            
            zip_path = self.genomes_dir / f"{accession}.zip"
            cmd = (
                f"datasets download genome accession {accession} "
                f"--reference --include genome,gff3 --filename {shlex.quote(str(zip_path))}" # shlex.quote(str()) allows this to be used in the command line argument
            )

            self.log(f"\nDownloading {accession}\n{cmd}\n")
            exit_val = os.system(cmd)
            if exit_val != 0:
                self.log(f"***ERROR*** failed to download {accession}\n")
                self.append_download_failure(
                    accession=accession,
                    stage="download",
                    command=cmd,
                    error_note=f"datasets download failed with exit code {exit_val}",
                )
                self.failed_downloads.append({
                    "accession": accession,
                    "stage": "download",
                    "error": f"exit code {exit_val}"
                })                
                continue
            try:
                shutil.unpack_archive(str(zip_path), str(acc_root))
                self.successfully_downloaded_data.append(accession)
                self.log(f"Successfully downloaded and unpacked {accession} -> {acc_root}\n")
            except Exception as exc:
                self.log(f"***ERROR*** failed to unpack {accession}: {exc}\n")
                if acc_root.exists(): # clean up the folder that was made during the attempt.
                    self.log(f"WARNING I'm removing the folder {acc_root} I initially created.")
                    shutil.rmtree(acc_root)
                self.append_download_failure(
                    accession=accession,
                    stage="unpack",
                    command=cmd,
                    error_note=str(exc),
                    )
                self.failed_downloads.append({
                    "accession": accession,
                    "stage": "unpack",
                    "error": str(exc)  
                    })
            finally:
                if zip_path.exists():
                    zip_path.unlink()

    def inspect_accession_root(self, accession: str) -> dict:
        accession_root = self.accession_root(accession)
        found_fna, found_gff, found_lift = self.find_accession_files(accession)

        return {
            "accession": accession,
            "accession_root": accession_root,
            "found_fna": found_fna,
            "found_gff": found_gff,
            "found_lift": found_lift,
            "is_broken": not bool(found_fna),
        }

    def handle_broken_accession_root(self, accession: str, reason: str) -> None:
        if self.dry_run_rebuild:
            self.rebuild_preview["broken_accessions"].append(accession)
            self.rebuild_preview["would_trash"].append({
                "accession": accession,
                "reason": reason,
            })
            self.log(
                f"DRY-RUN would move accession root to trash: {accession} | reason: {reason}\n"
            )
            return

        self.move_accession_root_to_trash(accession, reason=reason)

    def append_download_failure(
        self,
        accession: str,
        stage: str,
        command: str,
        error_note: str,
        ) -> None:
        row = {
            "timestamp": now_iso(),
            "accession": accession,
            "stage": stage,
            "command": command,
            "error_note": error_note,
            }
        file_exists = self.download_failures_csv.exists()
        with open(self.download_failures_csv, "a", newline="", encoding="utf-8") as handle:
            fieldnames = ["timestamp", "accession", "stage", "command", "error_note"]
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            if not file_exists:
                writer.writeheader()
            writer.writerow(row)
            
    def load_failed_accessions_for_retry(self) -> None:
        self.data = []
        if not self.download_failures_csv.exists():
            self.log("No download_failures.csv found; nothing to retry.\n")
            return

        rows = read_csv_rows(self.download_failures_csv)
        seen = set()

        for row in rows:
            accession = row.get("accession", "").strip()
            if not accession or accession in seen:
                continue
            seen.add(accession)
            enriched = self.fetch_summary_for_accession(accession)
            self.data.append(enriched)

        self.log(f"Loaded {len(self.data)} failed accessions for retry.\n")

    
    def remove_resolved_failures(self, resolved_accessions: list[str]) -> None:
        if not resolved_accessions:
            return
        if not self.download_failures_csv.exists():
            return

        resolved = set(resolved_accessions)
        rows = read_csv_rows(self.download_failures_csv)
        kept_rows = [
            row for row in rows
            if row.get("accession", "").strip() not in resolved
        ]

        fieldnames = ["timestamp", "accession", "stage", "command", "error_note"]
        write_csv_rows(self.download_failures_csv, fieldnames, kept_rows)
        
        self.log(
            f"Removed {len(rows) - len(kept_rows)} resolved rows from "
            f"download_failures.csv\n"
        )

    def find_accession_files(self, accession: str) -> tuple[str, str, str]:
        acc_root = self.accession_root(accession)
        if not acc_root.exists(): # remember ... if false fasle, if block triggers - if false true, if block does not trigger
            return "","",""
        
        fna_path = ""
        gff_path = ""
        lift_gff_path = ""
        
        for path in acc_root.rglob("*"):
            if not path.is_file(): # if false false i.e. the path is not a file skip
                continue
            
            name = path.name
            if name.endswith(".fna") and accession in name:
                fna_path = str(path)
            elif name.endswith(".gff") and "lift" not in name:
                gff_path = str(path)
            elif name.endswith("lift.gff"):
                lift_gff_path = str(path)
        return fna_path, gff_path, lift_gff_path
        
    def to_local_and_drive(self, abs_path: str) -> tuple[str, str]:
        if not abs_path:
            return "", ""
        abs_path_obj = Path(abs_path).resolve()
        try:
            rel = abs_path_obj.relative_to(self.root_local)
        except Exception:
            return str(abs_path_obj), ""
        
        local_ver = self.root_local / rel
        drive_ver = self.root_drive / rel
        return str(local_ver), str(drive_ver)
        
    def get_lineage(self, taxon_id: str) -> dict[str, str]:
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

        lineage_cmd = (
            f'echo {taxon_id} | taxonkit lineage | '
            f'taxonkit reformat -f "{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}"'
        )
        try:
            output = subprocess.check_output(lineage_cmd, shell=True, text=True).strip()
            fields = output.split("\t")
            cleaned = [fields[i].strip("{}") if i < len(fields) and fields[i] else "" for i in range(2, 8)]
            return dict(zip(empty.keys(), cleaned))
        except subprocess.CalledProcessError:
            return empty

    def build_assembly_record(self, entry: dict) -> AssemblyRecord:
        accession = entry.get("accession", "")
        organism_name = entry.get("organism", {}).get("organism_name", "")
        taxon_id = str(entry.get("organism", {}).get("tax_id", ""))
        assembly_info = entry.get("assembly_info", {})

        assembly_name = assembly_info.get("assembly_name", "")
        assembly_level = assembly_info.get("assembly_level", "") or "chromosome"
        refseq_category = assembly_info.get("refseq_category", "")
        source_database = "RefSeq" if accession.startswith("GCF_") else "GenBank"
        species_key = build_species_key(organism_name)
        lineage = self.get_lineage(taxon_id)
        
        found_fna, found_gff, found_lift = self.find_accession_files(accession)
        local_fna, drive_fna = self.to_local_and_drive(found_fna)
        local_gff, drive_gff = self.to_local_and_drive(found_gff)
        local_lift, drive_lift = self.to_local_and_drive(found_lift)
        
        timestamp = now_iso()
        
        return AssemblyRecord(
            accession=accession,
            accession_root=str(self.accession_root(accession)),
            organism_name=organism_name,
            taxon_id=taxon_id,
            assembly_name=assembly_name,
            assembly_level=assembly_level,
            refseq_category=refseq_category,
            source_database=source_database,
            is_refseq=accession.startswith("GCF_"),
            is_current_for_species=True,
            replaced_by_accession="",
            species_key=species_key,
            phylum=lineage["phylum"],
            superorder=lineage["superorder"],
            order=lineage["order"],
            family=lineage["family"],
            genus=lineage["genus"],
            genus_species=lineage["genus_species"],
            path_to_fna=local_fna,
            path_to_gff=local_gff,
            path_to_lift_gff=local_lift,
            drive_to_fna=drive_fna,
            drive_to_gff=drive_gff,
            drive_to_lift_gff=drive_lift,
            assembly_folder=str(self.accession_root(accession)),
            downloaded_at=timestamp,
            metadata_last_updated=timestamp,
            ncbi_search_term=self.search_term,
            notes="",
        )

    def append_update_event(
        self,
        species_key: str,
        organism_name: str,
        old_accession: str,
        new_accession: str,
        reason: str,
        ) -> None:
        payload = {
            "timestamp": now_iso(),
            "species_key": species_key,
            "organism_name": organism_name,
            "old_accession": old_accession,
            "new_accession": new_accession,
            "reason": reason,
            }
        with jsonlines.open(self.update_history_jsonl, mode="a", compact=True) as writer:
            writer.write(payload)
            
    def update_species_index(self) -> None:
        metadata_rows = read_csv_rows(self.genomes_metadata_csv)
        species_rows = []
        
        for row in metadata_rows:
            status = ( "current" if str(row.get("is_current_for_species", "")) == "True" else "superseded" )
            species_rows.append({
                "species_key": row.get("species_key", ""),
                "organism_name": row.get("organism_name", ""),
                "accession": row.get("accession", ""),
                "accession_root": row.get("accession_root", ""),
                "assembly_name": row.get("assembly_name", ""),
                "source_database": row.get("source_database", ""),
                "assembly_level": row.get("assembly_level", ""),
                "refseq_category": row.get("refseq_category", ""),
                "downloaded_at": row.get("downloaded_at", ""),
                "status": status,
                "replaced_by_accession": row.get("replaced_by_accession", ""),
                "project_count": "0",
                "notes": "",
                })
        fieldnames = [
            "species_key",
            "organism_name",
            "accession",
            "accession_root",
            "assembly_name",
            "source_database",
            "assembly_level",
            "refseq_category",
            "downloaded_at",
            "status",
            "replaced_by_accession",
            "project_count",
            "notes",
        ]                
        write_csv_rows(self.species_assembly_index, fieldnames, species_rows)        
                
    def update_metadata_tables(self) -> None:
        existing_rows = read_csv_rows(self.genomes_metadata_csv)
        rows_by_accession = {row["accession"]: row for row in existing_rows}
        
        for entry in self.data:
            accession = entry.get("accession", "")
            if accession not in self.successfully_downloaded_data:
                continue

            record = self.build_assembly_record(entry)
            new_row = asdict(record)
            species_key = record.species_key
            
            for row in existing_rows:
                if row.get("species_key") == species_key and row.get("accession") != accession:
                    row["is_current_for_species"] = "False"
                    row["replaced_by_accession"] = accession

                    self.append_update_event(
                        species_key=species_key,
                        organism_name=record.organism_name,
                        old_accession=row.get("accession", ""),
                        new_accession=accession,
                        reason="New assembly for same species discovered"
                    )

            if accession in rows_by_accession:
                rows_by_accession[accession].update(new_row)
            else:
                existing_rows.append(new_row)
                rows_by_accession[accession] = new_row
        fieldnames = list(AssemblyRecord.__dataclass_fields__.keys())
        write_csv_rows(self.genomes_metadata_csv, fieldnames, existing_rows)
        self.update_species_index()

    def append_successful_entries_to_catalog(self) -> None:
        if not self.successfully_downloaded_data:
            return

        with jsonlines.open(self.genomes_catalog_jsonl, mode="a", compact=True) as writer:
            for entry in self.data:
                accession = entry.get("accession", "")
                if accession in self.successfully_downloaded_data:
                    writer.write(entry)

    def project_manifest_path(self) -> Optional[Path]:
        if not self.project_name:
            return None
        return self.project_manifests_dir / f"{self.project_name}_manifest.csv"

    def update_project_registry(self, manifest_path: Path) -> None:
        rows = read_csv_rows(self.project_registry_csv)
        by_name = {row["project_name"]: row for row in rows if row.get("project_name")}

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

        fieldnames = [
            "project_name",
            "project_root",
            "manifest_path",
            "policy",
            "freeze_date",
            "notes",
        ]
        write_csv_rows(self.project_registry_csv, fieldnames, rows)

    def write_project_manifest(self) -> None:
        manifest_path = self.project_manifest_path()
        if manifest_path is None:
            return

        selected_accessions = set(self.resolve_project_accessions())
        missing = [acc for acc in sorted(selected_accessions) if not self.accession_root(acc).exists()]
        if missing:
            self.log(
                "***WARNING*** project manifest incluedes accessions not found on disk:\n"
                + "\n".join(missing)
                + "\n"
                )
        if not selected_accessions:
            self.log("No accessions selected for project manifest. \n")
            return
        
        metadata_rows = read_csv_rows(self.genomes_metadata_csv)
        selected_rows = []

        for row in metadata_rows:
            accession = row.get("accession", "")
            if accession not in selected_accessions:
                continue
            selected_rows.append({
                "project_name": self.project_name,
                "project_policy": self.project_policy,
                "manifest_version": "1",
                "manifest_date": datetime.date.today().isoformat(),
                "species_key": row.get("species_key", ""),
                "organism_name": row.get("organism_name", ""),
                "accession": accession,
                "accession_root": row.get("accession_root", ""),
                "assembly_name": row.get("assembly_name", ""),
                "project_role": "analysis",
                "source_database": row.get("source_database", ""),
                "is_current_for_species_at_freeze": row.get("is_current_for_species", ""),
                "analysis_notes": self.manifest_notes,
            })

        fieldnames = [
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
            "source_database",
            "is_current_for_species_at_freeze",
            "analysis_notes",
        ]
        write_csv_rows(manifest_path, fieldnames, selected_rows)
        self.update_project_registry(manifest_path)

    def looks_like_accession_dir(self, path: Path) -> bool:
        if not path.is_dir():
            return False
        return path.name.startswith(("GCA_", "GCF_"))

    def discover_local_accessions(self) -> list[str]:
        accessions = []
        for child in sorted(self.genomes_dir.iterdir()):
            if child.name == "records" or child.name == "trash":
                continue
            if self.looks_like_accession_dir(child):
                accessions.append(child.name)
        return accessions

    def load_accessions_from_csv(self, path: Path) -> list[str]:
        rows = read_csv_rows(path)
        accessions = []
        for row in rows:
            accession = row.get("accession", "").strip()
            if accession:
                accessions.append(accession)
        return accessions
        
    def load_accession_from_txt(self, path: Path) -> list[str]:
        accessions = []
        with open(path, 'r', encoding="utf-8") as handle:
            for raw_line in handle:
                accession = raw_line.strip()
                if not accession or accession.startswith('#'):
                    continue
                accessions.append(accession)
        return accessions
    
    def normalize_accession_list(self, items: list[str]) -> list[str]:
        seen = set()
        cleaned = []
        for item in items:
            acc = item.strip()
            if not acc:
                continue
            if acc in seen:
                continue
            seen.add(acc)
            cleaned.append(acc)
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
        if self.project_accessions:
            selected.extend(self.project_accessions)
        
        selected = self.normalize_accession_list(selected)
        
        if selected:
            self.project_selected_accessions = selected
            return selected
        # fallback behavior: if no explicit list is given, use newly downloaded accessions
        self.project_selected_accessions = list(self.successfully_downloaded_data)
        return self.project_selected_accessions
        
    def load_catalog_entries_by_accession(self) -> dict[str, dict]:
        by_accession: dict[str, dict] = {}
        if not self.genomes_catalog_jsonl.exists():
            return by_accession
        
        with open(self.genomes_catalog_jsonl, "r", encoding="utf-8") as handle:
            for raw_line in handle:
                if not raw_line.strip():
                    continue
                row = json.loads(raw_line)
                accession = row.get("accession", "").strip()
                if accession:
                    by_accession[accession] = row
        return by_accession
        
    def fetch_summary_for_accession(self, accession: str) -> dict:
        tmp_jsonl = self.records_dir / f"{accession}_retry_summary.jsonl"
        cmd = (
            f"datasets summary genome accession {shlex.quote(accession)} "
            f"--reference --as-json-lines > {shlex.quote(str(tmp_jsonl))}"
        )
        exit_val = os.system(cmd)
        if exit_val != 0 or not tmp_jsonl.exists():
            return {"accession": accession}

        try:
            with open(tmp_jsonl, "r", encoding="utf-8") as handle:
                for raw_line in handle:
                    if not raw_line.strip():
                        continue
                    row = json.loads(raw_line)
                    if row.get("accession", "") == accession:
                        return row
        finally:
            if tmp_jsonl.exists():
                tmp_jsonl.unlink()

        return {"accession": accession}

    def trash_dir(self) -> Path:
        path = self.genomes_dir / "trash"
        path.mkdir(parents=True, exist_ok=True)
        return path

    def move_accession_root_to_trash(self, accession: str, reason: str) -> Path:
        src = self.accession_root(accession)
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        dst = self.trash_dir() / f"{root_acc(accession)}__{timestamp}"
        
        if not src.exists():
            self.log(
                f"Requested trash move for missing accession root: {src} | reason: {reason}"
                )
            return dst
        shutil.move(str(src), str(dst))
        self.log(
            f"Moved accession root to trash: {src} -> {dst} | reason: {reason}")
        return dst
        
    def backup_metadata(self) -> None:
        if not self.genomes_metadata_csv.exists():
            return

        if self.dry_run_rebuild:
            self.log(
                "DRY-RUN would back up genomes_metadata.csv to genomes/trash/ before rewrite\n"
            )
            return

        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_path = self.trash_dir() / f"genomes_metadata__{timestamp}.csv"
        shutil.copy2(self.genomes_metadata_csv, backup_path)
        self.log(f"Backed up genomes_metadata.csv -> {backup_path}\n")

    def current_metadata_accessions(self) -> list[str]:
        if not self.genomes_metadata_csv.exists():
            return []

        rows = read_csv_rows(self.genomes_metadata_csv)
        accessions = []
        seen = set()

        for row in rows:
            accession = row.get("accession", "").strip()
            if not accession or accession in seen:
                continue
            seen.add(accession)
            accessions.append(accession)

        return accessions

    def summarize_rebuild_diff(self, rebuilt_records: list[dict]) -> None:
        current_accessions = set(self.current_metadata_accessions())
        rebuilt_accessions = set(
            row.get("accession", "").strip()
            for row in rebuilt_records
            if row.get("accession", "").strip()
        )

        self.rebuild_preview["rebuilt_accessions"] = sorted(rebuilt_accessions)
        self.rebuild_preview["missing_from_rebuild"] = sorted(
            current_accessions - rebuilt_accessions
        )
        self.rebuild_preview["new_rows"] = sorted(
            rebuilt_accessions - current_accessions
        )

    def report_rebuild_preview(self) -> None:
        local_accessions = self.rebuild_preview.get("local_accessions", [])
        rebuilt_accessions = self.rebuild_preview.get("rebuilt_accessions", [])
        broken_accessions = self.rebuild_preview.get("broken_accessions", [])
        missing_from_rebuild = self.rebuild_preview.get("missing_from_rebuild", [])
        new_rows = self.rebuild_preview.get("new_rows", [])

        self.log("\n===== DRY-RUN REBUILD PREVIEW =====\n")
        self.log(f"Local accession directories found: {len(local_accessions)}\n")
        self.log(f"Rows that would be written: {len(rebuilt_accessions)}\n")
        self.log(f"Broken accession roots: {len(broken_accessions)}\n")
        self.log(f"Current metadata rows that would disappear: {len(missing_from_rebuild)}\n")
        self.log(f"Accessions newly added to metadata: {len(new_rows)}\n")

        if broken_accessions:
            self.log("\nBroken accession roots that would be moved to trash:\n")
            for accession in broken_accessions[:25]:
                self.log(f"  - {accession}\n")
            if len(broken_accessions) > 25:
                self.log(f"  ... and {len(broken_accessions) - 25} more\n")

        if missing_from_rebuild:
            self.log("\nCurrent metadata rows that would be removed by rebuild:\n")
            for accession in missing_from_rebuild[:25]:
                self.log(f"  - {accession}\n")
            if len(missing_from_rebuild) > 25:
                self.log(f"  ... and {len(missing_from_rebuild) - 25} more\n")

        if new_rows:
            self.log("\nAccessions that would newly appear in rebuilt metadata:\n")
            for accession in new_rows[:25]:
                self.log(f"  - {accession}\n")
            if len(new_rows) > 25:
                self.log(f"  ... and {len(new_rows) - 25} more\n")

        self.log("===== END DRY-RUN PREVIEW =====\n\n")
 
    def build_local_only_record(self, accession: str) -> Optional[AssemblyRecord]:
        inspection = self.inspect_accession_root(accession)
        accession_root = inspection["accession_root"]
        found_fna = inspection["found_fna"]
        found_gff = inspection["found_gff"]
        found_lift = inspection["found_lift"]

        if inspection["is_broken"]:
            self.handle_broken_accession_root(
                accession,
                reason="missing primary genome FASTA during metadata rebuild",
            )
            return None

        local_fna, drive_fna = self.to_local_and_drive(found_fna)
        local_gff, drive_gff = self.to_local_and_drive(found_gff)
        local_lift, drive_lift = self.to_local_and_drive(found_lift)
        timestamp = now_iso()

        # if datasets is available, try to recover the summary metadata suing
        # the accession before falling all the way back to minimal local-only record.
        summary = None
        try:
            self.require_datasets_cli()
            summary = self.fetch_summary_for_accession(accession)
        except Exception:
            summary = None
        
        if summary:
            record = self.build_assembly_record(summary)
            
            # Preserve the file paths actually found on disk during rebuild
            record.path_to_fna = local_fna
            record.path_to_gff = local_gff
            record.path_to_lift_gff = local_lift
            record.drive_to_fna = drive_fna
            record.drive_to_gff = drive_gff
            record.drive_to_lift_gff = drive_lift
            record.assembly_folder = str(accession_root)
            record.ncbi_search_term = "rebuild_from_local"
            record.notes = "rebuild_from_local"
            return record

        # final fallback - keep the record, but make it explicit that this row was built
        # locally and could not be hydrated from datasets summary / catalog
            
        return AssemblyRecord(
            accession = accession,
            accession_root = str(accession_root),
            organism_name="",
            taxon_id="",
            assembly_name="",
            assembly_level="",
            refseq_category="",
            source_database="Refseq" if accession.startswith("GCF_") else "GenBank",
            is_refseq=accession.startswith("GCF_"),
            is_current_for_species=True,
            replaced_by_accession="",
            species_key=root_acc(accession),
            phylum="",
            superorder="",
            order="",
            family="",
            genus="",
            genus_species="",
            path_to_fna=local_fna,
            path_to_gff=local_gff,
            path_to_lift_gff=local_lift,
            drive_to_fna=drive_fna,
            drive_to_gff=drive_gff,
            drive_to_lift_gff=drive_lift,
            assembly_folder=str(accession_root),
            downloaded_at=timestamp,
            metadata_last_updated=timestamp,
            ncbi_search_term="rebuild_from_local",
            notes="rebuilt from local files without matching catalog entry",
            )

    def rebuild_genomes_metadata(self) -> None:
        local_accessions = self.discover_local_accessions()
        self.rebuild_preview["local_accessions"] = list(local_accessions)

        catalog_by_accession = self.load_catalog_entries_by_accession()
        rebuilt_records: list[dict] = []

        for accession in local_accessions:
            entry = catalog_by_accession.get(accession)

            if entry:
                record = self.build_assembly_record(entry)
            else:
                record = self.build_local_only_record(accession)

            if record is None:
                continue

            rebuilt_records.append(asdict(record))

        grouped: dict[str, list[dict]] = {}
        for row in rebuilt_records:
            species_key = row.get("species_key", "").strip()
            if not species_key:
                species_key = root_acc(row.get("accession", ""))
                row["species_key"] = species_key
            grouped.setdefault(species_key, []).append(row)

        for species_key, rows in grouped.items():
            if len(rows) == 1:
                rows[0]["is_current_for_species"] = True
                rows[0]["replaced_by_accession"] = ""
                continue

            rows_sorted = sorted(
                rows,
                key=lambda r: (
                    1 if str(r.get("is_refseq", "")) in ("True", "true", "1") or r.get("accession", "").startswith("GCF_") else 0,
                    r.get("accession", "")
                ),
                reverse=True,
            )

            current_acc = rows_sorted[0]["accession"]
            for row in rows:
                if row["accession"] == current_acc:
                    row["is_current_for_species"] = True
                    row["replaced_by_accession"] = ""
                else:
                    row["is_current_for_species"] = False
                    row["replaced_by_accession"] = current_acc

        rebuilt_records = sorted(rebuilt_records, key=lambda r: r.get("accession", ""))
        self.summarize_rebuild_diff(rebuilt_records)

        if self.dry_run_rebuild:
            self.report_rebuild_preview()
            return

        self.backup_metadata()
        fieldnames = list(AssemblyRecord.__dataclass_fields__.keys())
        write_csv_rows(self.genomes_metadata_csv, fieldnames, rebuilt_records)
        self.update_species_index()

        self.log(
            f"Rebuilt genomes metadata from local accession folders. "
            f"Found {len(local_accessions)} accession directories.\n"
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

        self.require_datasets_cli()

        if self.retry_failures:
            self.load_failed_accessions_for_retry()
            if not self.data:
                self.log("No failed accessions available for retry.\n")
                return
        else:
            self.query_ncbi()
            self.log_species_with_multiple_accessions()

        self.download_new_data()
        self.append_successful_entries_to_catalog()
        self.update_metadata_tables()

        if self.retry_failures:
            self.remove_resolved_failures(self.successfully_downloaded_data)

        if self.project_name:
            self.write_project_manifest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Manage genome assemblies in a mass-project-friendly hybrid model"
            "compatible with genomes/records/compleasm/."
        )
    )
    
    parser.add_argument(
        "genomes_dir",
        help="Path to genomes directory (must end with 'genomes')"
    )
    parser.add_argument(
        "search_term",
        help="NCBI taxon search term or accession"
    )
    parser.add_argument(
        "--project-name",
        default=None,
        help="Optional project name for writing a project manifest"
    )
    parser.add_argument(
        "--project-policy",
        choices=["frozen", "rolling", "exploratory"],
        default="exploratory",
        help="Manifest policy if --project-name is provided"
    )
    parser.add_argument(
        "--project-root",
        default=None,
        help="Optional project root path to record in project_registry.csv"
    )
    parser.add_argument(
        "--manifest-notes",
        default="",
        help="Optional notes written to the manifest rows"
    )
    parser.add_argument(
        "--retry-failures",
        action="store_true",
        help="Retry accessions listed in records/download_failuers.csv instead of running a new NCBI summary query"
    )
    parser.add_argument(
        "--rebuild-metadata",
        action="store_true",
        help=(
            "Rebuild genomes_metadata.csv and species_assembly_index.csv from "
            "accession folders already present under genomes/"
        )
    )
    parser.add_argument(
        "--project-accessions-csv",
        default=None,
        help=(
            "Optional CSV file containing accessions to include in the project "
            "manifest. Supports a column named accession."
        )
    )
    parser.add_argument(
        "--project-accessions",
        nargs="*",
        default=None,
        help=(
            "Optional explicit accession list to include in the project manifest, "
            "for example: --project-accessions GCF_123 GCA_456"
        )
    )
    parser.add_argument(
        "--project-accessions-txt",
        default=None,
        help=(
            "Optional plain text file containing one accession per line for the "
            "project manifest."
        )
    )

    parser.add_argument(
        "--dry-run-rebuild",
        action="store_true",
        help=(
            "Preview what rebuild-metadata would do without rewriting metadata, "
            "rewriting species index, or moving accession directories to trash"
        )
    )

    return parser.parse_args()
    
    
        

def main() -> None:
    args = parse_args()
    if args.rebuild_metadata and args.dry_run_rebuild:
        raise ValueError(
            "Use either --rebuild-metadata or --dry-run-rebuild, not both."
        )
    print("**Warning this version supports chromosome and complete assembly filters only. \n \
    It assumes you want only one genome per species and passes --reference to datasets***")

    mgr = GenomeManagerHybrid(
        genomes_dir=args.genomes_dir,
        search_term=args.search_term,
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
    )
    mgr.run()

if __name__ == "__main__":
    main()