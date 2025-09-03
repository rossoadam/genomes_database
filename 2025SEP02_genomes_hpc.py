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


class GenomeManagerHPC:
    def __init__(self):
        self.csv_path = ''

        # Root path to genomes and records
        self.local_genomes_folder = r'/Users/rossoaa/projects/genomes'
        self.local_records_folder = os.path.join(self.local_genomes_folder, 'records')

        # Taxon name for current run
        self.taxon = 'lacertidae'

        # Local JSONL catalog and update paths
        self.local_genomes_jsonl_pathway = os.path.join(self.local_records_folder, 'genomes_catalog.jsonl')

        # Updates file (summary output)
        self.updates_file_path = [
            os.path.join(self.local_records_folder, f'{self.taxon}_updates.jsonl')
        ]

        # NCBI datasets command
        self.ncbi_tags = (
            f'datasets summary genome taxon {self.taxon} '
            f'--reference --assembly-level chromosome --assembly-level complete '
            f'--as-json-lines > {self.updates_file_path[0]}'
        )

        Path(self.local_genomes_folder).mkdir(parents=True, exist_ok=True)
        Path(self.local_records_folder).mkdir(parents=True, exist_ok=True)

        self.data = []
        self.successfully_downloaded_data = []

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
            dl_cmd = f'datasets download genome accession {acc} --include genome,gff3 --filename {zip_path}'

            if os.system(dl_cmd) != 0:
                self.records_writing(f'\nDownload failed for {acc}\n')
                continue

            self.records_writing(f'\nDownload succeeded for {acc}\n')
            self.successfully_downloaded_data.append(acc)

            extract_path = os.path.join(self.local_genomes_folder, acc)
            os.makedirs(extract_path, exist_ok=True)
            os.system(f'unzip -o {zip_path} -d {extract_path}')
            os.remove(zip_path)

            with jsonlines.open(self.local_genomes_jsonl_pathway, 'a', compact=True) as writer:
                writer.write(entry)

    def make_csv(self):
        csv_path = os.path.join(self.local_records_folder, 'genomes_metadata.csv')
        with open(csv_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            header = [
                'accession', 'organism_name', 'path_to_fna', 'path_to_gff', 'path_to_lift_gff',
                'phylum', 'superorder', 'order', 'family', 'genus', 'genus_species'
            ]
            writer.writerow(header)

            with jsonlines.open(self.local_genomes_jsonl_pathway, 'r') as reader:
                for line in reader:
                    acc = line.get('accession', '')
                    organism = line.get('organism', {}).get('organism_name', '')
                    taxon_id = line.get('organism', {}).get('tax_id', '')
                    path_to_fna = ''
                    path_to_gff = ''
                    path_to_lift = ''
                    for root, dirs, files in os.walk(self.local_genomes_folder):
                        for f in files:
                            if f.endswith('.fna') and acc in f:
                                path_to_fna = os.path.join(root, f)
                            if f.endswith('.gff') and 'lift' not in f and acc in root:
                                path_to_gff = os.path.join(root, f)
                                print(path_to_gff)
                            if f.endswith('lift.gff') and acc in root:
                                path_to_lift = os.path.join(root, f)
                                print(path_to_lift)

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
                        acc, organism, path_to_fna, path_to_gff, path_to_lift,
                        lineage_dict['phylum'], lineage_dict['superorder'], lineage_dict['order'],
                        lineage_dict['family'], lineage_dict['genus'], lineage_dict['genus_species']
                    ])


if __name__ == '__main__':
    today = GenomeManagerHPC()
    today.check_records()
    today.check_datasets_cli()
    today.check_for_new_data()
    today.find_duplicates()
    today.fix_dupes()
    today.download_new_data()
    today.make_csv()
