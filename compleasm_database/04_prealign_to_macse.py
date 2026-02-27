#!/usr/bin/env python3
import subprocess
import pandas as pd
import argparse
import random
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser

class MacseAligner:
    def __init__(self, genomes_dir, test_mode=False):
        self.genomes_dir = Path(genomes_dir).resolve()
        self.test_mode = test_mode
        self.meta_path = self.genomes_dir / 'records/compleasm/records/metadata.csv'
        self.shared_genes_path = self.genomes_dir / 'records/compleasm/records/shared_single_genes.csv'
        self.pre_align_dir = self.genomes_dir / 'records/compleasm/alignments/01_pre_alignments'
        self.macse_dir = self.genomes_dir / 'records/compleasm/alignments/02_macse_alignments'

        self.pre_align_dir.mkdir(parents=True, exist_ok=True)
        self.macse_dir.mkdir(parents=True, exist_ok=True)

        with open(self.shared_genes_path, 'r') as f:
            all_orthologs = f.read().strip().split(',')
        self.orthologs = set(random.sample(all_orthologs, 30)) if test_mode else set(all_orthologs)

    def _repath(self, original_path):
        path_obj = Path(original_path)
        if path_obj.exists(): return path_obj
        parts = path_obj.parts
        return self.genomes_dir / "/".join(parts[parts.index('records'):]) if 'records' in parts else path_obj

    def run(self):
        print(f"--- Gathering Sequences ---")
        df_meta = pd.read_csv(self.meta_path)
        for g in self.orthologs:
            with open(self.pre_align_dir / f"{g}_unaligned.fasta", 'w') as f: pass

        for _, row in df_meta.iterrows():
            cds_path = self._repath(row['cds_fasta'])
            if not cds_path.exists(): continue
            with open(cds_path, 'r') as h:
                for head, seq in SimpleFastaParser(h):
                    gene_id = head.split()[0]
                    if gene_id in self.orthologs:
                        with open(self.pre_align_dir / f"{gene_id}_unaligned.fasta", 'a') as f:
                            f.write(f">{row['accession']}\n{seq}\n")

        print(f"--- Running MACSE ---")
        for fasta in self.pre_align_dir.glob("*_unaligned.fasta"):
            gene_id = fasta.name.replace("_unaligned.fasta", "")
            out_nt = self.macse_dir / f"{gene_id}_NT.fasta"
            if not out_nt.exists():
                subprocess.run(["macse", "-prog", "alignSequences", "-seq", str(fasta), 
                                "-out_NT", str(out_nt), "-out_AA", str(self.macse_dir / f"{gene_id}_AA.fasta")], check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes_dir")
    parser.add_argument("--test", action="store_true")
    args = parser.parse_args()
    MacseAligner(args.genomes_dir, args.test).run()