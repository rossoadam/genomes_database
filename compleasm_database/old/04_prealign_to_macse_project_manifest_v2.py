#!/usr/bin/env python3
import argparse
import logging
import random
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser


class MacseAligner:
    def __init__(
        self,
        genomes_dir,
        threads=1,
        test_mode=False,
        project_manifest=None,
        accessions_txt=None,
        shared_genes=None,
    ):
        self.genomes_dir = Path(genomes_dir).resolve()
        self.threads = threads
        self.test_mode = test_mode
        self.project_manifest = Path(project_manifest).resolve() if project_manifest else None
        self.accessions_txt = Path(accessions_txt).resolve() if accessions_txt else None
        self.shared_genes_arg = Path(shared_genes).resolve() if shared_genes else None

        self.records_dir = self.genomes_dir / 'records/compleasm/records'
        self.align_base = self.genomes_dir / 'records/compleasm/alignments'
        self.meta_path = self.records_dir / 'metadata.csv'
        self.pre_align_dir = self.align_base / '01_pre_alignments'
        self.macse_dir = self.align_base / '02_macse_alignments'
        self.log_file = self.records_dir / 'records_macse.log'

        self.pre_align_dir.mkdir(parents=True, exist_ok=True)
        self.macse_dir.mkdir(parents=True, exist_ok=True)

        logging.basicConfig(
            filename=self.log_file,
            level=logging.ERROR,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )

        self.allowed_accessions = self._load_allowed_accessions()
        self.cohort_label = self._build_cohort_label()
        self.shared_genes_path = self._resolve_shared_genes_path()
        self.orthologs = self._load_orthologs()

    def _sanitize_label(self, text):
        return re.sub(r'[^A-Za-z0-9._-]+', '_', text).strip('_')

    def _build_cohort_label(self):
        if self.project_manifest:
            return self._sanitize_label(self.project_manifest.stem)
        if self.accessions_txt:
            return self._sanitize_label(self.accessions_txt.stem)
        return None

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

    def _resolve_shared_genes_path(self):
        if self.shared_genes_arg:
            return self.shared_genes_arg
        if self.cohort_label:
            candidate = self.records_dir / f'shared_single_genes_{self.cohort_label}.csv'
            if candidate.exists():
                return candidate
        fallback = self.records_dir / 'shared_single_genes.csv'
        return fallback

    def _load_orthologs(self):
        if not self.shared_genes_path.exists():
            raise FileNotFoundError(f'Shared genes file not found: {self.shared_genes_path}')

        with open(self.shared_genes_path, 'r') as f:
            all_orthologs = [x for x in f.read().strip().split(',') if x]

        if self.test_mode:
            sample_size = min(300, len(all_orthologs))
            orthologs = set(random.sample(all_orthologs, sample_size))
            print(f"!!! TEST MODE ACTIVE: Only processing {len(orthologs)} genes !!!")
            return orthologs
        return set(all_orthologs)

    def _repath(self, original_path):
        path_obj = Path(original_path)
        if path_obj.exists():
            return path_obj
        parts = path_obj.parts
        if 'records' in parts:
            return self.genomes_dir / '/'.join(parts[parts.index('records'):])
        return path_obj

    def _run_single_macse(self, gene_id):
        fasta = self.pre_align_dir / f"{gene_id}_unaligned.fasta"
        out_nt = self.macse_dir / f"{gene_id}_NT.fasta"
        out_aa = self.macse_dir / f"{gene_id}_AA.fasta"

        if not fasta.exists():
            return

        if not out_nt.exists():
            try:
                subprocess.run([
                    'macse', '-prog', 'alignSequences',
                    '-seq', str(fasta),
                    '-out_NT', str(out_nt),
                    '-out_AA', str(out_aa)
                ], check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                error_msg = f"Gene: {gene_id} | Error: {e.stderr.strip()}"
                logging.error(error_msg)

    def run(self):
        if self.test_mode:
            print('Cleaning up old test files...')
            for f in self.pre_align_dir.glob('*.fasta'):
                f.unlink()

        print(f"--- Gathering Sequences for {len(self.orthologs)} genes ---")
        print(f"--- Shared genes source: {self.shared_genes_path} ---")
        df_meta = pd.read_csv(self.meta_path)

        if self.allowed_accessions is not None:
            before = len(df_meta)
            df_meta = df_meta[df_meta['accession'].astype(str).isin(self.allowed_accessions)].copy()
            print(f"--- Cohort filter retained {len(df_meta)} of {before} compleasm rows ---")

        for g in self.orthologs:
            with open(self.pre_align_dir / f"{g}_unaligned.fasta", 'w'):
                pass

        for _, row in df_meta.iterrows():
            cds_path = self._repath(row['cds_fasta'])
            if not cds_path.exists():
                continue
            with open(cds_path, 'r') as h:
                for head, seq in SimpleFastaParser(h):
                    gene_id = head.split()[0]
                    if gene_id in self.orthologs:
                        with open(self.pre_align_dir / f"{gene_id}_unaligned.fasta", 'a') as f:
                            f.write(f">{row['accession']}\n{seq}\n")

        print(f"--- Running MACSE in parallel (Threads: {self.threads}) ---")
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            executor.map(self._run_single_macse, list(self.orthologs))

        print(f"--- Finished. Check {self.log_file} for any failures. ---")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('genomes_dir')
    parser.add_argument('--project-manifest', help='CSV manifest with an accession column defining the cohort to use')
    parser.add_argument('--accessions-txt', help='TXT file with one accession per line defining the cohort to use')
    parser.add_argument('--shared-genes', help='Optional explicit path to the shared_single_genes CSV to use')
    parser.add_argument('--test', action='store_true')
    parser.add_argument('-j', '--threads', type=int, default=1)
    args = parser.parse_args()

    MacseAligner(
        args.genomes_dir,
        threads=args.threads,
        test_mode=args.test,
        project_manifest=args.project_manifest,
        accessions_txt=args.accessions_txt,
        shared_genes=args.shared_genes,
    ).run()
