#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import math
import re
from collections import Counter
from pathlib import Path

import pandas as pd

################
# EXAMPLE USAGE:
#
# 1. Standard Run (80% occupancy, only Chromosome and Scaffold assemblies):
#    python 03_single_orthologs_csv.py ../../../genomes /path/to/project-manifest.csv "Complete Chromosome" "Scaffold"
#
# 2. Strict Run (100% occupancy for a high-quality phylogenomic matrix):
#    python 03_single_orthologs_csv.py ../../../genomes /path/to/project-manifest.csv "Complete Chromosome" --threshold 1.0
#
# 3. Relaxed Run (60% occupancy to maximize gene count):
#    python 03_single_orthologs_csv.py ../../../genomes /path/to/project-manifest.csv "Complete Chromosome" "Scaffold" --threshold 0.6
#
# ARGUMENTS:
#   genomes_dir    Path to the root project directory containing 'records/'.
#   accessions     CSV manifest with an accession column OR txt file with one accession per line.
#   targets        One or more assembly levels to include. Common values:
#                  "Complete Chromosome", "Scaffold", "Contig".
#                  (Separate multiple targets with spaces).
#
# OPTIONAL FLAGS:
#   --threshold    The fraction of genomes (0.0 - 1.0) that must have a gene
#                  marked as 'Single' to include it. Default is 0.80 (80%).
################


def root_acc(accession: str) -> str:
    accession = str(accession).strip()
    return accession.split('.')[0] if accession else accession


class SingleOrthologFinder:
    def __init__(self, genomes_dir, accessions_file, targets, threshold_pct=0.80):
        self.genomes_dir = Path(genomes_dir).resolve()
        self.accessions_file = Path(accessions_file).resolve()
        self.targets = targets or []
        self.threshold_pct = threshold_pct

        self.gene_id_dir = self.genomes_dir / 'records/compleasm/mb_downloads/sauropsida_odb12/hmms'
        self.genomes_meta = self.genomes_dir / 'records/genomes_metadata.csv'
        self.compleasm_meta = self.genomes_dir / 'records/compleasm/records/metadata.csv'
        self.records_dir = self.genomes_dir / 'records/compleasm/records'

        if self.gene_id_dir.exists():
            self.gene_id_list = sorted(f.stem for f in self.gene_id_dir.glob('*.hmm'))
        else:
            print(f"Error: HMM directory not found at {self.gene_id_dir}")
            self.gene_id_list = []

        self.allowed_accessions, self.allowed_roots = self._load_allowed_accessions()
        self.cohort_label = self._sanitize_label(self.accessions_file.stem)
        self.output_path = self.records_dir / f'shared_single_genes_{self.cohort_label}.csv'
        self.single_gene_ids = []

    def _sanitize_label(self, text):
        return re.sub(r'[^A-Za-z0-9._-]+', '_', text).strip('_')

    def _read_accessions_txt(self, txt_path):
        accessions = []
        with open(txt_path, 'r') as handle:
            for line in handle:
                value = line.strip()
                if not value or value.startswith('#'):
                    continue
                accessions.append(value.split()[0])
        return sorted(set(accessions))

    def _read_manifest_accessions(self, manifest_path):
        df = pd.read_csv(manifest_path)
        if 'accession' not in df.columns:
            raise ValueError(f"Manifest is missing required 'accession' column: {manifest_path}")
        accessions = df['accession'].dropna().astype(str).str.strip()
        accessions = accessions[accessions != '']
        return sorted(set(accessions))

    def _load_allowed_accessions(self):
        suffix = self.accessions_file.suffix.lower()
        if suffix == '.csv':
            accessions = self._read_manifest_accessions(self.accessions_file)
            print(f"Loaded {len(accessions)} accessions from manifest: {self.accessions_file}")
        else:
            accessions = self._read_accessions_txt(self.accessions_file)
            print(f"Loaded {len(accessions)} accessions from txt file: {self.accessions_file}")
        roots = {root_acc(a) for a in accessions}
        return set(accessions), roots

    def _repath(self, path_to_fix):
        path_str = str(path_to_fix)
        root_name = self.genomes_dir.name
        if root_name in path_str:
            suffix = path_str.split(root_name, 1)[-1].lstrip('/\\')
            return self.genomes_dir / suffix
        return Path(path_str)

    def _accession_allowed(self, accession_value: str) -> bool:
        accession_value = str(accession_value).strip()
        return accession_value in self.allowed_accessions or root_acc(accession_value) in self.allowed_roots

    def mask(self):
        if not self.gene_id_list:
            print('No BUSCO IDs found in HMM directory, aborting.')
            return []

        genomes_df = pd.read_csv(self.genomes_meta)
        compleasm_df = pd.read_csv(self.compleasm_meta)

        genomes_df['accession'] = genomes_df['accession'].astype(str).str.strip()
        compleasm_df['accession'] = compleasm_df['accession'].astype(str).str.strip()
        genomes_df['accession_root'] = genomes_df['accession'].map(root_acc)
        compleasm_df['accession_root'] = compleasm_df['accession'].map(root_acc)

        level_lookup = {}
        for _, row in genomes_df.iterrows():
            level_lookup[row['accession']] = row['assembly_level']
            level_lookup[row['accession_root']] = row['assembly_level']

        gene_occurrence_counter = Counter()
        included_count = 0
        skipped_for_cohort = 0
        skipped_for_target = 0
        missing_tsv = 0
        cohort_hits = []

        for _, row in compleasm_df.iterrows():
            acc = row['accession']
            acc_root = row['accession_root']
            level = level_lookup.get(acc, level_lookup.get(acc_root))

            if not self._accession_allowed(acc):
                skipped_for_cohort += 1
                continue

            cohort_hits.append(acc)

            if self.targets and level not in self.targets:
                skipped_for_target += 1
                continue

            tsv_path = self._repath(row['full_table'])
            if not tsv_path.exists():
                print(f"Warning: TSV file missing for {acc}: {tsv_path}")
                missing_tsv += 1
                continue

            try:
                df_tsv = pd.read_csv(
                    tsv_path,
                    sep='\t',
                    comment='#',
                    header=0,
                    usecols=[0, 1],
                    names=['Gene', 'Status'],
                )
                current_single = df_tsv[df_tsv['Status'] == 'Single']['Gene'].tolist()
                gene_occurrence_counter.update(current_single)
                print(f"Processed {acc}: {level}")
                included_count += 1
            except Exception as e:
                print(f"Error processing {acc}: {e}")

        if included_count == 0:
            print('No genomes passed the filters. Nothing to write.')
            print(f"Manifest/TXT accessions matched to compleasm rows: {len(set(cohort_hits))}")
            print(f"Rows skipped because not in cohort: {skipped_for_cohort}")
            print(f"Rows skipped because of assembly level: {skipped_for_target}")
            return []

        required_count = max(1, math.ceil(included_count * self.threshold_pct))
        self.single_gene_ids = sorted(
            gene
            for gene, count in gene_occurrence_counter.items()
            if count >= required_count and gene in self.gene_id_list
        )

        self.records_dir.mkdir(parents=True, exist_ok=True)
        with open(self.output_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(self.single_gene_ids)

        print('\nFinished.')
        print(f"Cohort filter file: {self.accessions_file}")
        print(f"Accessions requested: {len(self.allowed_accessions)}")
        print(f"Compleasm rows matching cohort: {len(set(cohort_hits))}")
        print(f"Assembly levels retained: {self.targets if self.targets else 'all'}")
        print(f"Rows skipped because not in cohort: {skipped_for_cohort}")
        print(f"Rows skipped because of assembly level: {skipped_for_target}")
        print(f"Rows skipped because full_table was missing: {missing_tsv}")
        print(f"Occupancy Threshold: {self.threshold_pct * 100:.1f}% ({required_count}/{included_count} genomes)")
        print(f"Genes Meeting Threshold: {len(self.single_gene_ids)}")
        print(f"Shared IDs written to: {self.output_path}")
        return self.single_gene_ids


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Find BUSCO/Compleasm single-copy orthologs across a selected genome cohort.'
    )
    parser.add_argument('genomes_dir', help='Path to genomes directory')
    parser.add_argument('accessions', help='CSV manifest with accession column or TXT with one accession per line')
    parser.add_argument(
        'targets',
        nargs='*',
        help='Optional assembly levels to keep, e.g. "Complete Chromosome" "Scaffold"',
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=0.80,
        help='Occupancy threshold (e.g. 0.80 for 80 percent)',
    )
    args = parser.parse_args()

    finder = SingleOrthologFinder(
        args.genomes_dir,
        args.accessions,
        args.targets,
        threshold_pct=args.threshold,
    )
    finder.mask()
