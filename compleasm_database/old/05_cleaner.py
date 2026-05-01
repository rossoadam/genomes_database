#!/usr/bin/env python3
import pandas as pd
import argparse
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser

class AlignmentCleaner:
    def __init__(self, genomes_dir, threshold, errors):
        self.genomes_dir = Path(genomes_dir).resolve()
        self.threshold = threshold
        self.errors = errors
        self.macse_dir = self.genomes_dir / 'records/compleasm/alignments/02_macse_alignments'
        self.fs_csv = self.genomes_dir / 'records/compleasm/records/records_frameshift.csv'
        self.sum_csv = self.genomes_dir / 'records/compleasm/records/records_filter_summary.csv'
        
        df_meta = pd.read_csv(self.genomes_dir / 'records/compleasm/records/metadata.csv')
        self.total_species = len(df_meta)

    def clean(self):
        print(f"--- Cleaning with Threshold: {self.threshold} ---")
        print(f"--- Cleaning with Errors: {self.errors} or less ---")
        fs_data, dropped_ids = [], []
        total_processed = 0

        for nt_file in self.macse_dir.glob("*_NT.fasta"):
            gene_id = nt_file.name.replace("_NT.fasta", "")
            total_processed += 1
            clean_seqs, clean_count = [], 0

            with open(nt_file, 'r') as h:
                for acc, seq in SimpleFastaParser(h):
                    fs_count = seq.count('!')
                    is_filtered = fs_count > int(self.errors)
                    fs_data.append({'threshold': self.threshold, 'gene_id': gene_id, 
                                    'accession': acc, 'frameshifts': fs_count, 'sequence_filtered': is_filtered})
                    
                    if not is_filtered:
                        clean_count += 1
                        clean_seqs.append(f">{acc}\n{seq.replace('!', '-')}\n")

            occupancy = clean_count / self.total_species
            if occupancy >= self.threshold:
                with open(self.macse_dir / f"{gene_id}_clean_thresh_{self.threshold}_errors_{self.errors}.fasta", 'w') as f:
                    f.writelines(clean_seqs)
            else:
                dropped_ids.append(gene_id)

        # Update CSVs (Append mode)
        pd.DataFrame(fs_data).to_csv(self.fs_csv, mode='a', index=False, header=not self.fs_csv.exists())
        
        summary = pd.DataFrame([{
            'threshold': self.threshold,
            'total_genes_analyzed': total_processed,
            'genes_dropped': len(dropped_ids),
            'genes_retained': total_processed - len(dropped_ids),
            'dropped_gene_ids': ",".join(dropped_ids)
        }])
        summary.to_csv(self.sum_csv, mode='a', index=False, header=not self.sum_csv.exists())
        print(f"Done. Retained {total_processed - len(dropped_ids)} genes.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes_dir")
    parser.add_argument("--threshold", type=float, default=.9)
    parser.add_argument("--errors", type=int, default=1)
    args = parser.parse_args()
    AlignmentCleaner(args.genomes_dir, args.threshold, args.errors).clean()