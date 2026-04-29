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
            f"--reference --assembly-level chromosome --assembly-level complete "
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
                if not raw_line.strip():
                    continue
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
                f"--include genome,gff3 --filename {shlex.quote(str(zip_path))}" # shlex.quote(str()) allows this to be used in the command line argument
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
                self.data.append({"accession": accession})
            
            self.log(f"loaded {len(self.data)} failed accession for retry.\n")
    
    

class GenomeManagerHPC:
    def __init__(self, genomes_folder_arg: str, search_term: str):
        self.csv_path = ''

        # --- Validate genomes folder ends with "genomes" ---
        genomes_folder_arg = os.path.abspath(genomes_folder_arg)
        if os.path.basename(os.path.normpath(genomes_folder_arg)) != "genomes":
            print('current directory needs to labeled "genomes"')
            sys.exit(1)

        # --- Detect Colab for convenience (no behavior change required here) ---
        self.on_colab = ('google.colab' in sys.modules)

        # --- Canonical roots (unchanged) ---
        # self.root_local = "/Users/rossoaa/Desktop/dummy"
        # self.root_drive = "/content/drive/Othercomputers/macbook/projects"
        self.root_local = Path(genomes_folder_arg).parent
        self.root_drive = "/content/drive/Othercomputers/macbook/projects"

        # --- Active folders based on user-provided genomes folder ---
        # genomes folder is argv[1], records live inside this folder
        self.local_genomes_folder = genomes_folder_arg
        self.active_root = os.path.dirname(os.path.normpath(self.local_genomes_folder))
        self.local_records_folder = os.path.join(self.local_genomes_folder, 'records')

        # --- Search term + mode ---
        # If the term starts with GCA_ or GCF_ -> accession mode, else taxon mode
        self.search_term = search_term.strip()
        self.is_accession = self.search_term.startswith(("GCA_", "GCF_"))

        # A safe tag for filenames
        self.safe_tag = (
            self.search_term.replace("/", "_")
                            .replace("\\", "_")
                            .replace(" ", "_")
        )

        # Local JSONL catalog and update paths
        self.local_genomes_jsonl_pathway = os.path.join(self.local_records_folder, 'genomes_catalog.jsonl')

        # Updates file (summary output)
        self.updates_file_path = [
            os.path.join(self.local_records_folder, f'{self.safe_tag}_updates.jsonl')
        ]

        # NCBI datasets command (taxon vs accession)
        # Keep your reference and assembly-level filters
        if self.is_accession:
            # accession mode
            query_bit = f"accession {shlex.quote(self.search_term)}"
        else:
            # taxon mode
            query_bit = f"taxon {shlex.quote(self.search_term)}"

        self.ncbi_tags = (
            f'datasets summary genome {query_bit} '
            f'--reference --assembly-level chromosome --assembly-level complete '
            f'--as-json-lines > {shlex.quote(self.updates_file_path[0])}'
        )

        # Ensure directories exist
        Path(self.local_genomes_folder).mkdir(parents=True, exist_ok=True)
        Path(self.local_records_folder).mkdir(parents=True, exist_ok=True)

        self.data = []
        self.successfully_downloaded_data = []

    # --- Small helper: build local/drive twins from a path found under active_root ---
    def _to_local_and_drive(self, abs_active_path: str):
        """
        Given an absolute path under the active root (found via os.walk),
        return (local_version, drive_version) by reconstructing both roots.
        If abs_active_path is empty, returns ('','').
        """
        if not abs_active_path:
            return '', ''
        try:
            rel = os.path.relpath(abs_active_path, start=self.active_root)
        except Exception:
            # Fallback: try prefix swaps
            if abs_active_path.startswith(self.root_local):
                rel = os.path.relpath(abs_active_path, start=self.root_local)
            elif abs_active_path.startswith(self.root_drive):
                rel = os.path.relpath(abs_active_path, start=self.root_drive)
            else:
                # Unknown base; cannot map reliably
                return '', ''

        local_ver = os.path.normpath(os.path.join(self.root_local, rel))
        drive_ver = os.path.normpath(os.path.join(self.root_drive, rel))
        return local_ver, drive_ver

    def root_acc(self, x: str) -> str:
        return x.split('.')[0] if isinstance(x, str) else x

    def records_writing(self, r_data):
        os.makedirs(self.local_records_folder, exist_ok=True)
        with open(os.path.join(self.local_records_folder, 'genome_auto_updates_records.txt'), 'a') as file:
            file.write(r_data)

    def check_local_jsonl(self):
        Path(self.local_genomes_jsonl_pathway).touch(exist_ok=True)

    def check_records(self):
        records_path = os.path.join(self.local_records_folder, 'genome_auto_updates_records.txt')
        if not os.path.isfile(records_path):
            with open(records_path, 'w') as file:
                file.write('')

    def check_datasets_cli(self):
        if os.system('datasets --help > /dev/null 2>&1') != 0:
            r_data = f"\n\n\n{datetime.datetime.now()} ***ERROR*** PROGRAM TERMINATED. NCBI CLI 'datasets' not in $SHELL path.\n"
            self.records_writing(r_data)
            sys.exit()
        else:
            r_data = f"\n\n\n{datetime.datetime.now()}"
            self.records_writing(r_data)

    def check_for_new_data(self):
        r_data = f'\n\nChecking the NCBI Database for new updates using these tags:\n{self.ncbi_tags}\n'
        self.records_writing(r_data)

        self.check_local_jsonl()

        os.chdir(self.local_records_folder)
        exit_val = os.system(self.ncbi_tags)
        if exit_val != 0:
            r_data = f'\n***ERROR*** PROGRAM TERMINATED. FAILED TO RETRIEVE UPDATES FILE FROM NCBI, see: {self.updates_file_path[0]}\n'
            self.records_writing(r_data)
            sys.exit()
        else:
            r_data = '\n--Successfully retrieved updated information from NCBI--'
            self.records_writing(r_data)

        updates_file_path = self.updates_file_path[0]
        tmp_updates = updates_file_path.replace('.jsonl', '.compact.jsonl')
        with jsonlines.open(updates_file_path, mode='r') as reader, \
                jsonlines.open(tmp_updates, mode='w', compact=True) as writer:
            for line in reader:
                writer.write(line)
        os.replace(tmp_updates, updates_file_path)

        tmp_local = os.path.join(self.local_records_folder, 'genomes_catalog.compact.jsonl')
        wrote_any = False
        with jsonlines.open(self.local_genomes_jsonl_pathway, mode='r') as reader, \
              jsonlines.open(tmp_local, mode='w', compact=True) as writer:
            for line in reader:
                acc = line.get('accession', '')
                # ONLY proceed if this accession was just downloaded in this run
                if acc not in self.successfully_downloaded_data:
                    continue
                writer.write(line)
                wrote_any = True
        if wrote_any:
            os.replace(tmp_local, self.local_genomes_jsonl_pathway)
        elif os.path.exists(tmp_local):
            os.remove(tmp_local)

        try:
            with open(self.local_genomes_jsonl_pathway, 'r') as f:
                local_text = f.read()
        except FileNotFoundError:
            local_text = ''

        with open(updates_file_path, 'r') as k:
            for line in k:
                if not line.strip():
                    continue
                lk = json.loads(line)
                if lk.get('accession') and lk['accession'] in local_text:
                    continue
                self.data.append(lk)

    def find_duplicates(self):
        seen = {}
        duplicates = {}
        for genome in self.data:
            organism = genome.get('organism', {}).get('organism_name', '')
            acc = genome.get('accession', '')
            if organism in seen:
                if organism not in duplicates:
                    duplicates[organism] = [seen[organism]]
                duplicates[organism].append(acc)
            else:
                seen[organism] = acc
        return duplicates

    def fix_dupes(self):
        """
        Deduplicates self.data by organism name, prioritizing GCF (RefSeq) over GCA (GenBank).
        """
        # 1. Sort data so GCF comes before GCA for the same organism
        # GCF starts with 'GCF', GCA starts with 'GCA'. Alphabetically 'GCA' < 'GCF',
        # so we sort descending to put 'GCF' first.
        self.data.sort(key=lambda x: x.get('accession', ''), reverse=True)

        fixed = []
        seen_organisms = set()
        
        for genome in self.data:
            acc = genome.get('accession', '')
            organism = genome.get('organism', {}).get('organism_name', '')
            
            if organism not in seen_organisms:
                fixed.append(genome)
                seen_organisms.add(organism)
            else:
                self.records_writing(f"Removing duplicate/inferior record: {acc} ({organism})\n")
        
        self.data = fixed
        self.records_writing(f"Deduplication complete. {len(self.data)} unique organisms remaining.\n")

    def download_new_data(self):
        for entry in self.data:
            acc = entry.get('accession')
            zip_path = os.path.join(self.local_genomes_folder, f'{acc}.zip')
            dl_cmd = f'datasets download genome accession {acc} --include genome,gff3 --filename {shlex.quote(zip_path)}'

            if os.system(dl_cmd) != 0:
                self.records_writing(f'\nDownload failed for {acc}\n')
                continue

            self.records_writing(f'\nDownload succeeded for {acc}\n')
            self.successfully_downloaded_data.append(acc)

            extract_path = os.path.join(self.local_genomes_folder, acc)
            os.makedirs(extract_path, exist_ok=True)
            os.system(f'unzip -o {shlex.quote(zip_path)} -d {shlex.quote(extract_path)}')
            os.remove(zip_path)

            with jsonlines.open(self.local_genomes_jsonl_pathway, 'a', compact=True) as writer:
                writer.write(entry)

    def make_csv(self):
        csv_path = os.path.join(self.local_records_folder, 'genomes_metadata.csv')
        # Check if the file already exists to decide if we need a header
        file_exists = os.path.isfile(csv_path)
        with open(csv_path, 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            header = [
                'accession', 'organism_name', 'assembly_level',
                'path_to_fna', 'path_to_gff', 'path_to_lift_gff',
                'drive_to_fna', 'drive_to_gff', 'drive_to_lift_gff',
                'phylum', 'superorder', 'order', 'family', 'genus', 'genus_species'
            ]
            # Only write the header if the file is being created for the first time
            if not file_exists:
                writer.writerow(header)

            with jsonlines.open(self.local_genomes_jsonl_pathway, 'r') as reader:
                for line in reader:
                    acc = line.get('accession', '')
                    organism = line.get('organism', {}).get('organism_name', '')
                    taxon_id = line.get('organism', {}).get('tax_id', '')

                    # Paths discovered under the active root (via os.walk of argv[1])
                    found_fna_active = ''
                    found_gff_active = ''
                    found_lift_active = ''

                    for root, dirs, files in os.walk(self.local_genomes_folder):
                        dirs[:]=dirs
                        for f in files:
                            if f.endswith('.fna') and acc in f:
                                found_fna_active = os.path.join(root, f)
                            if f.endswith('.gff') and 'lift' not in f and acc in root and "records" not in root:
                                found_gff_active = os.path.join(root, f)
                                print(found_gff_active)
                            if f.endswith('lift.gff') and acc in root:
                                found_lift_active = os.path.join(root, f)
                                print(found_lift_active)

                    # Map the active-found paths to BOTH roots
                    local_fna, drive_fna = self._to_local_and_drive(found_fna_active)
                    local_gff, drive_gff = self._to_local_and_drive(found_gff_active)
                    local_lift, drive_lift = self._to_local_and_drive(found_lift_active)

                    # Lineage (unchanged)
                    lineage_dict = {k: '' for k in ['phylum', 'superorder', 'order', 'family', 'genus', 'genus_species']}
                    if taxon_id:
                        lineage_cmd = f'echo {taxon_id} | taxonkit lineage | taxonkit reformat -f "{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}" -F'
                        try:
                            lineage_output = subprocess.check_output(lineage_cmd, shell=True, text=True).strip()
                            fields = lineage_output.split('\t')
                            cleaned_fields = [(fields[i].strip('{}') if i < len(fields) and fields[i] else '') for i in range(2, 8)]
                            lineage_dict = dict(zip(['phylum', 'superorder', 'order', 'family', 'genus', 'genus_species'], cleaned_fields))
                        except subprocess.CalledProcessError:
                            lineage_dict = {k: '' for k in ['phylum', 'superorder', 'order', 'family', 'genus', 'genus_species']}
                            pass

                    writer.writerow([
                        acc, organism, "chromosome",
                        local_fna, local_gff, local_lift,
                        drive_fna, drive_gff, drive_lift,
                        lineage_dict['phylum'], lineage_dict['superorder'], lineage_dict['order'],
                        lineage_dict['family'], lineage_dict['genus'], lineage_dict['genus_species']
                    ])

"""
def _usage_and_exit():
    print("Usage: python3 script.py /abs/path/to/genomes 'search_term_or_accession'")
    print("  - First arg: genomes folder (must end with 'genomes')")
    print("  - Second arg: NCBI search term. If it starts with GCA_ or GCF_, accession mode is used; otherwise taxon mode.")
    print("**Warning this version supports chromosome assembly_level only**")
    sys.exit(1)
"""    
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
        "--retry-failures",
        action="store_true",
        help="Retry accessions listed in records/download_failuers.csv instead of running a new NCBI summary query"
    )
    return parser.parse_args()
    
    
        

if __name__ == '__main__':
    #if len(sys.argv) < 3:
    #    _usage_and_exit()
    # print("**Warning this version supports chromosome assembly_level only**")
    # genomes_folder = sys.argv[1]
    # search_term = sys.argv[2]
    args = parse_args()
    gmh = GenomeManagerHybrid(args.genomes_dir, args.search_term)
    gmh.initialize_records()
    gmh.log("check transfer to practice")
    ncbi_command = gmh.build_ncbi_summary_command()
    print(ncbi_command)
    gmh.query_ncbi()
    gmh.log_species_with_multiple_accessions()
    print(gmh.download_failures_csv)
    print(gmh.retry_failures)
    # genomes_folder = Path(args.genomes_dir).resolve()
    # search_term = args.search_term
    # print("PRINTING GENOMES FOLDER, SEARCH TERM:", genomes_folder, search_term)
    #mgr = GenomeManagerHPC(genomes_folder, search_term)
    #mgr.check_records()
    #mgr.check_datasets_cli()
    #mgr.check_for_new_data()
    #mgr.find_duplicates()
    #mgr.fix_dupes()
    #mgr.download_new_data()
    #mgr.make_csv()
