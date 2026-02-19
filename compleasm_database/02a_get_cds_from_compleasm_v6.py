# python script to get coding sequence from a genome from a compleasm output table
# python3 get_cds_from_compleasm_v5.py /Volumes/T9/Lepidodactylus_binning/triocanu/JUICER/out_JBAT_final.FINAL_moestus.fa /Volumes/T9/Lepidodactylus_binning/triocanu/lugubris_moestus_compleasm_after_juicer/sauropsida_odb12/full_table.tsv lugubris_moestus
# python3 get_cds_from_compleasm_v5.py /Volumes/T9/Lepidodactylus_binning/triocanu/JUICER/out_JBAT_final.FINAL_PANTAI.fa /Volumes/T9/Lepidodactylus_binning/triocanu/lugubris_pantai_compleasm_after_juicer/sauropsida_odb12/full_table.tsv lugubris_pantai

import sys
import pyfaidx
from pathlib import Path
from datetime import datetime

def load_genome(input_fasta: Path) -> pyfaidx.Fasta:
    # 1. Index the FASTA file with pyfaidx for fast access
    try:
        genome = pyfaidx.Fasta(
			str(input_fasta),
			as_raw=True,		# slices return plaint strings
			build_index=True 	# makes .fai if missing
        )
        print(f"Status: Genome FASTA indexed successfully: {input_fasta}")
        return genome

    except pyfaidx.FaidxException as e:
        print(f"CRITICAL ERROR: Failed to load genome FASTA with pyfaidx. Ensure the .fai file exists for {input_fasta}. Error: {e}")
        sys.exit(1)

    except Exception as e:
        print(f"CRITICAL ERROR: Failed to load genome FASTA. Error: {e}")
        sys.exit(1)

def reverse_complement(seq):
    complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
        'D': 'H', 'H': 'D', 'N': 'N',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
        'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
        'k': 'm', 'm': 'k', 'b': 'v', 'v': 'b',
        'd': 'h', 'h': 'd', 'n': 'n'
    }

    return ''.join(complement.get(base, base) for base in reversed(seq))

def complement(sequence):
    """Return the complement of a DNA sequence, including IUPAC ambiguity codes."""
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
        'H': 'D', 'V': 'B', 'N': 'N',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
        'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
        'k': 'm', 'm': 'k', 'b': 'v', 'd': 'h',
        'h': 'd', 'v': 'b', 'n': 'n'
    }
    return ''.join(complement.get(base, base) for base in sequence)

genome = load_genome(Path(sys.argv[1]))
species = sys.argv[3]
compleasm_path = Path(sys.argv[2])
outdir = compleasm_path.parent
cdsfilename = outdir / f"{species}_cds_compleasm.tsv"
cdsfasta = outdir / f"{species}_cds_compleasm.fasta"
record = outdir / f"{species}_cds_compleasm_records.txt"

#
with open(cdsfilename,"w") as outhandle, open(cdsfasta,"w") as fastahandle, open(compleasm_path, "r") as compleasmhandle, open(record,"w") as to_record:
	outhandle.write("Busco ID\tspecies\tlocation\tstrand\tsequence\n")
#
	for line in compleasmhandle:
		line = line.rstrip()
	
		sequence = ''
	
		if line.startswith("Gene"):
			now = datetime.now()
			print("running...")
			print(f"Time is {now}", file=to_record)
			# header, next
			continue
		else:
			data = line.split("\t")
			if len(data) < 13:
				print(f"WARNING: Skipping malformed line: {line}", file=to_record)
				continue
			status = data[1]
			buscoid = data[0]
			if status == "Single":
				chromosome = data[2]
				strand = data[5]
				coordinate_data = data[12].split("|")
				coordinate_data.sort(key=lambda x: int(x.split('_')[0])) # the - lines are in high to low order this will allow me to reverse complement at the end
				location = data[2]
				good_busco = True
				for coordinates in coordinate_data:
					parts = coordinates.split("_")
					start_1_based = int(parts[0])
					end_1_based = int(parts[1])
					start_0_based = start_1_based - 1
					end_0_based = end_1_based
					
					if len(parts) != 3:
						good_busco = False
						print(f"WARNING: {buscoid}: unexpected exon token '{coordinates}'", file=to_record)
						break 
					exon_strand = parts[2] 
					if exon_strand != strand:
						print(f"Skipping {buscoid} because exon strand != gene strand", file=to_record)
						good_busco = False
						break
					try:
						exon = genome[chromosome][start_0_based:end_0_based]
					except KeyError:
						print(f"WARNING: {buscoid}: contig '{chromosome}' not found in FASTA; skipping", file=to_record)
						good_busco = False
						break
					sequence += exon # The exons are ordered low to high see above
				if not good_busco:
					continue
				if strand == "-":
					sequence = reverse_complement(sequence)
				
				outhandle.write(buscoid + "\t" + species + "\t" + location + "\t" + strand + "\t" + sequence + "\n")
				fastahandle.write(">" + buscoid + "\n" + sequence + "\n")	
			else:
				print(f"{buscoid} is {status}.", file=to_record)


