#!/usr/bin/env python3

import os
import sys
import csv
import json
import shutil
import datetime
from datetime import date
from pathlib import Path
import jsonlines
import subprocess
import shlex


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
        self.root_local = "/Users/rossoaa/projects"
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
        fixed = []
        seen = set()
        for genome in self.data:
            acc = genome.get('accession', '')
            organism = genome.get('organism', {}).get('organism_name', '')
            if organism in seen:
                continue
            seen.add(organism)
            fixed.append(genome)
        self.data = fixed

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
        with open(csv_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            header = [
                'accession', 'organism_name',
                'path_to_fna', 'path_to_gff', 'path_to_lift_gff',
                'drive_to_fna', 'drive_to_gff', 'drive_to_lift_gff',
                'phylum', 'superorder', 'order', 'family', 'genus', 'genus_species'
            ]
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
                        for f in files:
                            if f.endswith('.fna') and acc in f:
                                found_fna_active = os.path.join(root, f)
                            if f.endswith('.gff') and 'lift' not in f and acc in root:
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
                        acc, organism,
                        local_fna, local_gff, local_lift,
                        drive_fna, drive_gff, drive_lift,
                        lineage_dict['phylum'], lineage_dict['superorder'], lineage_dict['order'],
                        lineage_dict['family'], lineage_dict['genus'], lineage_dict['genus_species']
                    ])


def _usage_and_exit():
    print("Usage: python3 script.py /abs/path/to/genomes 'search_term_or_accession'")
    print("  - First arg: genomes folder (must end with 'genomes')")
    print("  - Second arg: NCBI search term. If it starts with GCA_ or GCF_, accession mode is used; otherwise taxon mode.")
    sys.exit(1)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        _usage_and_exit()

    genomes_folder = sys.argv[1]
    search_term = sys.argv[2]

    mgr = GenomeManagerHPC(genomes_folder, search_term)
    mgr.check_records()
    mgr.check_datasets_cli()
    mgr.check_for_new_data()
    mgr.find_duplicates()
    mgr.fix_dupes()
    mgr.download_new_data()
    mgr.make_csv()
