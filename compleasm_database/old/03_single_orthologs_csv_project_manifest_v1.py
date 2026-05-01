#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import math
import re
from collections import Counter
from pathlib import Path
import argparse
import pandas as pd


class SingleOrthologFinder:
    def __init__(
        self,
        genomes_dir,
        targets,
        threshold_pct=0.80,
        project_manifest=None,
        accessions_txt=None,
    ):
        self.genomes_dir = Path(genomes_dir).resolve()
        self.targets = targets or []
        self.threshold_pct = threshold_pct
        self.project_manifest = Path(project_manifest).resolve() if project_manifest else None
        self.accessions_txt = Path(accessions_txt).resolve() if accessions_txt else None

        self.gene_id_dir = self.genomes_dir / 'records/compleasm/mb_downloads/sauropsida_odb12/hmms'
        self.genomes_meta = self.genomes_dir / 'records/genomes_metadata.csv'
        self.compleasm_meta = self.genomes_dir / 'records/compleasm/records/metadata.csv'
        self.records_dir = self.genomes_dir / 'records/compleasm/records'

        if self.gene_id_dir.exists():
            self.gene_id_list = sorted(f.stem for f in self.gene_id_dir.glob('*.hmm'))
        else:
            print(f"Error: HMM directory not found at {self.gene_id_dir}")
            self.gene_id_list = []

        self.allowed_accessions = self._load_allowed_accessions()
        self.cohort_label = self._build_cohort_label()
        self.output_path = self._build_output_path()
        self.single_gene_ids = []

    def _sanitize_label(self, text):
        return re.sub(r'[^A-Za-z0-9._-]+', '_', text).strip('_')

    def _build_cohort_label(self):
        if self.project_manifest:
            return self._sanitize_label(self.project_manifest.stem)
        if self.accessions_txt:
            return self._sanitize_label(self.accessions_txt.stem)
        return None

    def _build_output_path(self):
        if self.cohort_label:
            return self.records_dir / f'shared_single_genes_{self.cohort_label}.csv'
        return self.records_dir / 'shared_single_genes.csv'

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
        if self.project_manifest and self.accessions_txt:
            raise ValueError('Use either --project-manifest or --accessions-txt, not both.')
        if self.project_manifest:
            accessions = self._read_manifest_accessions(self.project_manifest)
            print(f"Loaded {len(accessions)} accessions from manifest: {self.project_manifest}")
            return set(accessions)
        if self.accessions_txt:
            accessions = self._read_accessions_txt(self.accessions_txt)
            print(f"Loaded {len(accessions)} accessions from txt file: {self.accessions_txt}")
            return set(accessions)
        return None

    def _repath(self, path_to_fix):
        path_str = str(path_to_fix)
        root_name = self.genomes_dir.name
        if root_name in path_str:
            suffix = path_str.split(root_name, 1)[-1].lstrip('/\\')
            return self.genomes_dir / suffix
        return Path(path_str)

    def mask(self):
        if not self.gene_id_list:
            print('No BUSCO IDs found in HMM directory, aborting.')
            return []

        genomes_df = pd.read_csv(self.genomes_meta)
        compleasm_df = pd.read_csv(self.compleasm_meta)
        acc_to_level_dict = dict(zip(genomes_df['accession'], genomes_df['assembly_level']))

        gene_occurrence_counter = Counter()
        included_count = 0
        skipped_for_cohort = 0
        skipped_for_target = 0

        for _, row in compleasm_df.iterrows():
            acc = str(row['accession']).strip()
            level = acc_to_level_dict.get(acc)

            if self.allowed_accessions is not None and acc not in self.allowed_accessions:
                skipped_for_cohort += 1
                continue

            if self.targets and level not in self.targets:
                skipped_for_target += 1
                continue

            tsv_path = self._repath(row['full_table'])
            if not tsv_path.exists():
                print(f"Warning: TSV file missing for {acc}: {tsv_path}")
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
        if self.allowed_accessions is not None:
            print(f"Cohort filter file: {self.project_manifest or self.accessions_txt}")
            print(f"Accessions requested: {len(self.allowed_accessions)}")
            print(f"Rows skipped because not in cohort: {skipped_for_cohort}")
        if self.targets:
            print(f"Assembly levels retained: {self.targets}")
            print(f"Rows skipped because of assembly level: {skipped_for_target}")
        else:
            print('Assembly level filter: disabled')
        print(f"Occupancy Threshold: {self.threshold_pct * 100:.1f}% ({required_count}/{included_count} genomes)")
        print(f"Genes Meeting Threshold: {len(self.single_gene_ids)}")
        print(f"Shared IDs written to: {self.output_path}")
        return self.single_gene_ids


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Find BUSCO/Compleasm single-copy orthologs across a selected genome cohort.'
    )
    parser.add_argument('genomes_dir', help='Path to genomes directory')
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
    parser.add_argument(
        '--project-manifest',
        help='CSV manifest with an accession column defining the cohort to use',
    )
    parser.add_argument(
        '--accessions-txt',
        help='TXT file with one accession per line defining the cohort to use',
    )
    args = parser.parse_args()

    finder = SingleOrthologFinder(
        args.genomes_dir,
        args.targets,
        threshold_pct=args.threshold,
        project_manifest=args.project_manifest,
        accessions_txt=args.accessions_txt,
    )
    finder.mask()
