#!/usr/bin/env python3
"""Rename accession-based FASTA headers using an explicit genome manifest.

Only the ``accession`` and ``species_key`` columns in --manifest are used
for taxon names. Species keys are normalized to ``Genus_species``: genus
capitalized, species lowercase, and any subspecies designation removed.

All FASTA files are validated before any output is written. Sequence lines
and their line endings are copied unchanged; only header lines are replaced.

Example
-------
python3 09_rename_aligned_fastas_v4.py \
    /Users/rossoaa/projects/genomes \
    --manifest /Users/rossoaa/projects/genomes/body-size-dna-dynamics_manifest.csv \
    --fastadir /path/to/input_fastas \
    --outdir /path/to/09_renamed_fastas_capitalized
"""

import argparse
import csv
import sys
import tempfile
from pathlib import Path


FASTA_SUFFIXES = {".fasta", ".fa", ".fas", ".fna"}


class ValidationError(Exception):
    """Raised when manifest or FASTA validation fails."""


class RenameAlignedFastas:
    def __init__(self, genomes_dir, manifest, fasta_dir, outdir):
        self.genomes_dir = Path(genomes_dir).resolve()
        self.manifest = Path(manifest).resolve()
        self.fasta_dir = Path(fasta_dir).resolve()
        self.outdir = Path(outdir).resolve()

        self.accession_to_organism = {}
        self.subspecies_accessions = set()
        self.fasta_list = []
        self.validated_files = 0
        self.accessions_with_subspecies_used = set()
        self.shortened_header_count = 0

    @staticmethod
    def normalize_species_key(species_key):
        """Convert Genus_species[_subspecies...] to Genus_species.

        Return (normalized_name, removed_subspecies). Invalid fields raise
        ValueError so malformed manifest entries cannot silently propagate.
        """
        parts = species_key.strip().split("_")
        if len(parts) < 2 or any(not part.strip() for part in parts):
            raise ValueError(
                f"expected Genus_species[_subspecies], got {species_key!r}"
            )
        genus, species = parts[0], parts[1]
        if any(character.isspace() for character in genus + species):
            raise ValueError(f"genus and species must not contain spaces: {species_key!r}")
        return f"{genus.capitalize()}_{species.lower()}", len(parts) > 2

    def load_manifest(self):
        """Build validated accession-to-species map from the supplied CSV."""
        errors = []
        accession_first_row = {}

        with self.manifest.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise ValidationError(f"Manifest is empty: {self.manifest}")

            # Strip only column-header whitespace (not values).
            reader.fieldnames = [name.strip() if name is not None else ""
                                 for name in reader.fieldnames]
            required = {"accession", "species_key"}
            missing = required.difference(reader.fieldnames)
            if missing:
                raise ValidationError(
                    f"Manifest missing required column(s): {', '.join(sorted(missing))}. "
                    f"Found: {', '.join(reader.fieldnames)}"
                )

            for row_number, row in enumerate(reader, start=2):
                accession = (row.get("accession") or "").strip()
                species_key = (row.get("species_key") or "").strip()

                if not accession:
                    errors.append(f"Manifest row {row_number}: empty accession")
                elif any(char.isspace() for char in accession):
                    errors.append(
                        f"Manifest row {row_number}: accession contains whitespace: "
                        f"{accession!r}"
                    )
                elif accession in accession_first_row:
                    errors.append(
                        f"Manifest row {row_number}: duplicate accession "
                        f"{accession!r} (first seen on row {accession_first_row[accession]})"
                    )
                else:
                    accession_first_row[accession] = row_number

                if not species_key:
                    errors.append(f"Manifest row {row_number}: empty species_key")
                    continue

                try:
                    normalized, removed_subspecies = self.normalize_species_key(species_key)
                except ValueError as exc:
                    errors.append(f"Manifest row {row_number}: {exc}")
                    continue

                if accession and accession_first_row.get(accession) == row_number:
                    self.accession_to_organism[accession] = normalized
                    if removed_subspecies:
                        self.subspecies_accessions.add(accession)

        if errors:
            raise ValidationError(
                f"Manifest validation failed ({len(errors)} issue(s)):\n"
                + "\n".join(f"  - {error}" for error in errors)
            )
        if not self.accession_to_organism:
            raise ValidationError("Manifest has no usable accession/species_key entries")

        print(f"Manifest: {self.manifest}")
        print(f"Unique manifest accessions loaded: {len(self.accession_to_organism)}")

    def load_fasta_list(self):
        """Collect supported FASTAs in deterministic order (not recursively)."""
        self.fasta_list = sorted(
            path for path in self.fasta_dir.iterdir()
            if path.is_file() and path.suffix in FASTA_SUFFIXES
        )
        if not self.fasta_list:
            raise ValidationError(
                f"No FASTA files (.fasta, .fa, .fas, .fna) found in: {self.fasta_dir}"
            )
        print(f"Input FASTA directory: {self.fasta_dir}")
        print(f"Output FASTA directory: {self.outdir}")
        print(f"FASTA files discovered: {len(self.fasta_list)}")

    @staticmethod
    def accession_from_header(line):
        """Extract accession as the first whitespace-delimited header token."""
        tokens = line[1:].strip().split()
        return tokens[0] if tokens else ""

    def validate_fastas(self):
        """Validate every input FASTA without creating any output files."""
        errors = []
        self.validated_files = 0
        self.shortened_header_count = 0
        self.accessions_with_subspecies_used.clear()

        for fasta_path in self.fasta_list:
            accessions_seen = {}  # accession -> first header line number
            names_seen = {}       # normalized species name -> (accession, line)
            encountered_header = False
            file_errors = 0

            with fasta_path.open("r", encoding="utf-8", newline="") as handle:
                for line_number, line in enumerate(handle, start=1):
                    if line.startswith(">"):
                        encountered_header = True
                        accession = self.accession_from_header(line)
                        if not accession:
                            errors.append(
                                f"{fasta_path.name}:{line_number}: empty FASTA header"
                            )
                            file_errors += 1
                            continue

                        if accession in accessions_seen:
                            errors.append(
                                f"{fasta_path.name}:{line_number}: duplicate accession "
                                f"{accession!r} (first at line {accessions_seen[accession]})"
                            )
                            file_errors += 1
                            continue
                        accessions_seen[accession] = line_number

                        species_name = self.accession_to_organism.get(accession)
                        if species_name is None:
                            errors.append(
                                f"{fasta_path.name}:{line_number}: accession "
                                f"{accession!r} not found in manifest {self.manifest.name}"
                            )
                            file_errors += 1
                            continue

                        if species_name in names_seen:
                            earlier_accession, earlier_line = names_seen[species_name]
                            errors.append(
                                f"{fasta_path.name}:{line_number}: duplicate normalized "
                                f"species name {species_name!r} from accessions "
                                f"{earlier_accession!r} (line {earlier_line}) and "
                                f"{accession!r} (line {line_number})"
                            )
                            file_errors += 1
                        else:
                            names_seen[species_name] = (accession, line_number)

                        if accession in self.subspecies_accessions:
                            self.shortened_header_count += 1
                            self.accessions_with_subspecies_used.add(accession)

                    elif line.strip() and not encountered_header:
                        errors.append(
                            f"{fasta_path.name}:{line_number}: sequence data "
                            "before first FASTA header"
                        )
                        file_errors += 1

            if not encountered_header:
                errors.append(f"{fasta_path.name}: no FASTA headers found")
                file_errors += 1
            if file_errors == 0:
                self.validated_files += 1

        print(f"FASTA files validated successfully: {self.validated_files}/{len(self.fasta_list)}")
        if errors:
            raise ValidationError(
                f"FASTA validation failed ({len(errors)} issue(s)); "
                "no output FASTA files were written:\n"
                + "\n".join(f"  - {error}" for error in errors)
            )

    @staticmethod
    def header_line_ending(original_line):
        """Keep input newline convention for renamed header lines."""
        if original_line.endswith("\r\n"):
            return "\r\n"
        if original_line.endswith("\n"):
            return "\n"
        if original_line.endswith("\r"):
            return "\r"
        return "\n"

    def write_new_fasta(self, fasta_path, target_dir):
        """Copy sequences verbatim, replacing only the accession header text."""
        output_path = target_dir / fasta_path.name
        with fasta_path.open("r", encoding="utf-8", newline="") as source, \
                output_path.open("w", encoding="utf-8", newline="") as dest:
            for line in source:
                if line.startswith(">"):
                    accession = self.accession_from_header(line)
                    new_header = self.accession_to_organism[accession]
                    dest.write(f">{new_header}{self.header_line_ending(line)}")
                else:
                    dest.write(line)
        return output_path

    def loop(self):
        """Stage conversions and publish them only after every FASTA passes QC."""
        self.load_manifest()
        self.load_fasta_list()
        self.validate_fastas()

        # Delay output-directory creation until after validation succeeds.
        self.outdir.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=".09_rename_aligned_fastas_v4_", dir=self.outdir
        ) as temporary_directory:
            staging_dir = Path(temporary_directory)
            for fasta_path in self.fasta_list:
                self.write_new_fasta(fasta_path, staging_dir)

            for fasta_path in self.fasta_list:
                target = self.outdir / fasta_path.name
                (staging_dir / fasta_path.name).replace(target)
                print(f"Wrote renamed FASTA: {target}")

        print("\nSummary")
        print(f"  Manifest: {self.manifest}")
        print(f"  Unique accessions loaded: {len(self.accession_to_organism)}")
        print(f"  Input directory: {self.fasta_dir}")
        print(f"  Output directory: {self.outdir}")
        print(f"  FASTA files discovered: {len(self.fasta_list)}")
        print(f"  FASTA files validated successfully: {self.validated_files}")
        print(
            "  Unique species labels shortened (accessions used): "
            f"{len(self.accessions_with_subspecies_used)}"
        )
        print(f"  FASTA headers shortened (all files): {self.shortened_header_count}")
        print(f"  FASTA files successfully written: {len(self.fasta_list)}")
        print("Done.")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Rename aligned FASTA headers using accession and species_key "
            "from a required project manifest. Species labels become Genus_species."
        )
    )
    parser.add_argument("genomes_dir", help="Path to genomes root directory")
    parser.add_argument(
        "--manifest", required=True,
        help="Required CSV manifest containing accession and species_key columns",
    )
    parser.add_argument(
        "--fastadir", required=True, help="Directory containing aligned input FASTA files",
    )
    parser.add_argument(
        "--outdir",
        help=(
            "Output directory (default: "
            "genomes/records/compleasm/alignments/07_renamed_fastas)"
        ),
    )
    return parser.parse_args()


def main():
    args = parse_args()
    genomes_dir = Path(args.genomes_dir).resolve()
    manifest = Path(args.manifest).resolve()
    fasta_dir = Path(args.fastadir).resolve()
    outdir = (
        Path(args.outdir).resolve() if args.outdir else
        genomes_dir / "records/compleasm/alignments/07_renamed_fastas"
    )

    if not genomes_dir.is_dir():
        raise ValidationError(f"Genomes directory not found: {genomes_dir}")
    if not manifest.is_file():
        raise ValidationError(f"Manifest file not found: {manifest}")
    if not fasta_dir.is_dir():
        raise ValidationError(f"FASTA input directory not found: {fasta_dir}")
    if fasta_dir == outdir:
        raise ValidationError(
            "Input and output directories must be different to protect source FASTAs"
        )
    if outdir.exists() and not outdir.is_dir():
        raise ValidationError(f"Output path exists but is not a directory: {outdir}")

    renamer = RenameAlignedFastas(
        genomes_dir=genomes_dir,
        manifest=manifest,
        fasta_dir=fasta_dir,
        outdir=outdir,
    )
    renamer.loop()


if __name__ == "__main__":
    try:
        main()
    except (ValidationError, OSError, csv.Error, UnicodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
