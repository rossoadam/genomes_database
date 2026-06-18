#!/usr/bin/env python3

"""
11_rename_treefile_keep_branch_lengths_updated.py

Rename accession labels in a Newick tree using Compleasm metadata.csv.

Default behavior is NHPhyML-friendly:
  - rename accession tip labels to genus_species
  - remove branch lengths
  - remove internal node bootstrap/support labels

For PGLS or other comparative methods that need branch lengths, add:
  --keep-branch-lengths

Example, topology-only tree for NHPhyML:

python3 11_rename_treefile_keep_branch_lengths_updated.py \
  /Users/rossoaa/projects/genomes \
  --tree /Users/rossoaa/projects/genomes/records/compleasm/alignments/06_rooted_tree_nonreversible/nonrev_dna.treefile

Example, renamed tree with branch lengths for PGLS:

python3 11_rename_treefile_keep_branch_lengths_updated.py \
  /Users/rossoaa/projects/genomes \
  --tree /Users/rossoaa/projects/genomes/records/compleasm/alignments/t1.0_e1_o7_with_s_punctatus/05_unrooted_tree_t1.0_e1_o7_with_s_punctatus/rev_dna.treefile \
  --keep-branch-lengths \
  --outdir /Users/rossoaa/projects/genomes/records/compleasm/alignments/renamed_trees_for_pgls
"""

import argparse
import csv
import re
from pathlib import Path


def load_metadata(metadata_csv: Path) -> dict[str, str]:
    """Return accession -> genus_species from Compleasm metadata.csv."""
    accession_to_species: dict[str, str] = {}

    with metadata_csv.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"accession", "genus_species"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Metadata file is missing required column(s): {sorted(missing)}. "
                f"Columns found: {reader.fieldnames}"
            )

        for row in reader:
            accession = (row.get("accession") or "").strip()
            genus_species = (row.get("genus_species") or "").strip().replace(" ", "_")
            if not accession or not genus_species:
                continue
            accession_to_species[accession] = genus_species

    if not accession_to_species:
        raise ValueError(f"No accession/genus_species mappings found in {metadata_csv}")

    return accession_to_species


def build_accession_pattern(accession_to_species: dict[str, str]) -> re.Pattern:
    """
    Match accessions as Newick labels.

    This captures labels such as:
      (GCA_123.1:0.01,GCF_456.1:0.02)
      (prefix|GCA_123.1:0.01)  [only accession portion is replaced if exact accession appears]
    """
    accessions = sorted(accession_to_species.keys(), key=len, reverse=True)
    escaped = [re.escape(acc) for acc in accessions]

    # Require a Newick label boundary before and after the accession so branch
    # length colons are preserved when --keep-branch-lengths is used.
    pattern = r"(?<=[(,\s])(" + "|".join(escaped) + r")(?=[:),;\s])"
    return re.compile(pattern)


def rename_tree_labels(tree_text: str, accession_to_species: dict[str, str], pattern: re.Pattern) -> tuple[str, int]:
    """Rename accessions and return renamed tree plus replacement count."""
    replacement_count = 0

    def replacer(match: re.Match) -> str:
        nonlocal replacement_count
        accession = match.group(1)
        replacement_count += 1
        return accession_to_species.get(accession, accession)

    renamed = pattern.sub(replacer, tree_text)
    return renamed, replacement_count


def remove_branch_lengths(tree_text: str) -> str:
    """
    Remove branch lengths like :0.123, :1e-05, :-0.02.
    """
    return re.sub(r":-?\d+(\.\d+)?([eE][-+]?\d+)?", "", tree_text)


def remove_bootstrap_values(tree_text: str) -> str:
    """
    Remove internal node support labels that appear immediately after ')'.

    Examples:
      )95:0.01      -> ):0.01
      )100          -> )
      )99.5:0.01    -> ):0.01
      )1e2:0.01     -> ):0.01

    Branch lengths are preserved because the regex removes only the support
    value between ')' and the following ':' / ',' / ')' / ';'.
    """
    return re.sub(r"\)(\d+(\.\d+)?([eE][-+]?\d+)?)(?=[:),;])", ")", tree_text)


def clean_tree(tree_text: str, keep_branch_lengths: bool) -> str:
    """
    Clean tree after accession renaming.

    Internal node support values are removed in both modes because they can be
    parsed as internal node labels by some tools. Branch lengths are removed
    only when keep_branch_lengths is False.
    """
    tree_text = remove_bootstrap_values(tree_text)
    if not keep_branch_lengths:
        tree_text = remove_branch_lengths(tree_text)
    return tree_text


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rename accession labels in a Newick tree using metadata.csv. "
            "By default, branch lengths are removed for NHPhyML. Use "
            "--keep-branch-lengths for PGLS/comparative analyses."
        )
    )
    parser.add_argument("genomes_dir", help="Path to genomes directory")
    parser.add_argument("--tree", required=True, help="Path to input tree file")
    parser.add_argument(
        "--metadata",
        default=None,
        help=(
            "Optional metadata CSV with accession and genus_species columns. "
            "Default: <genomes_dir>/records/compleasm/records/metadata.csv"
        ),
    )
    parser.add_argument(
        "--outdir",
        help=(
            "Optional output directory. Default: "
            "<genomes_dir>/records/compleasm/alignments/08_renamed_treefile"
        ),
    )
    parser.add_argument(
        "--keep-branch-lengths",
        action="store_true",
        help=(
            "Preserve branch lengths in the output Newick tree. "
            "Use this for PGLS. Internal support labels are still removed."
        ),
    )
    parser.add_argument(
        "--suffix",
        default=None,
        help=(
            "Optional suffix inserted before the output file extension. "
            "Default: '_renamed_keep_lengths' when --keep-branch-lengths is used, "
            "otherwise '_renamed_topology_only'."
        ),
    )
    return parser.parse_args()


def output_name(input_tree: Path, keep_branch_lengths: bool, suffix: str | None) -> str:
    if suffix is None:
        suffix = "_renamed_keep_lengths" if keep_branch_lengths else "_renamed_topology_only"

    # Many tree files use suffixes like .treefile or .nwk. Preserve the suffix
    # when present; otherwise append the suffix to the full filename.
    if input_tree.suffix:
        return f"{input_tree.stem}{suffix}{input_tree.suffix}"
    return f"{input_tree.name}{suffix}"


def main() -> None:
    args = parse_args()

    genomes_dir = Path(args.genomes_dir).expanduser().resolve()
    tree_file = Path(args.tree).expanduser().resolve()

    if args.metadata:
        metadata_csv = Path(args.metadata).expanduser().resolve()
    else:
        metadata_csv = genomes_dir / "records/compleasm/records/metadata.csv"

    if not metadata_csv.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_csv}")
    if not tree_file.exists():
        raise FileNotFoundError(f"Tree file not found: {tree_file}")

    if args.outdir:
        outdir = Path(args.outdir).expanduser().resolve()
    else:
        outdir = genomes_dir / "records/compleasm/alignments/08_renamed_treefile"

    outdir.mkdir(parents=True, exist_ok=True)
    output_tree = outdir / output_name(tree_file, args.keep_branch_lengths, args.suffix)

    accession_to_species = load_metadata(metadata_csv)
    pattern = build_accession_pattern(accession_to_species)

    tree_text = tree_file.read_text().strip()
    renamed_tree, replacement_count = rename_tree_labels(tree_text, accession_to_species, pattern)
    cleaned_tree = clean_tree(renamed_tree, keep_branch_lengths=args.keep_branch_lengths)

    output_tree.write_text(cleaned_tree + "\n")

    mode = "renamed tree with branch lengths preserved for PGLS" if args.keep_branch_lengths else "topology-only renamed tree for NHPhyML"
    print(f"Metadata: {metadata_csv}")
    print(f"Input tree: {tree_file}")
    print(f"Output tree: {output_tree}")
    print(f"Accessions available in metadata: {len(accession_to_species)}")
    print(f"Tip-label replacements made: {replacement_count}")
    print(f"Done. Wrote {mode}.")


if __name__ == "__main__":
    main()
