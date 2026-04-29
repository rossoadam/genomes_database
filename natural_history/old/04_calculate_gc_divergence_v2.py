#!/usr/bin/env python3
"""
02_calculate_gc_divergence.py

Calculate GC3 summaries and GC divergence values from nhPhyML outputs.

This is the second script in the refactored workflow. It uses:
  1. nhPhyML outputs (*.phylip, *.phylip_nhPhymlGC.tree, *.phylip_nhPhyml.lk)
  2. node_numbers_rebuilt.csv from 99_compout_tree_to_csv.py
  3. optional TSVs made by 01_make_master_input_tsvs.py:
       genome_sizes.tsv
       species_mass.tsv
  4. optional project manifest for final master_output.tsv ordering/filtering

Core calculations:
  d_ij   = sqrt(sum_k((GC3_i,k - GC3_j,k)^2))
  d_ianc = sqrt(sum_k((GC3_i,k - GC3_anc,k)^2))

Outputs:
  gc3_per_gene.tsv
  gc3_summary.tsv
  dij_results.json
  dij_results.tsv
  dianc_results.json
  dianc_results.tsv
  master_output.tsv

Example:
  python 02_calculate_gc_divergence.py /path/to/genomes \
    --nhphyml-dir /path/to/09_nhphyml_output \
    --node-numbers /path/to/node_numbers_rebuilt.csv \
    --manifest /path/to/project_manifest.csv
"""

import argparse
import csv
import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


def normalize_species_name(organism_name: str) -> str:
    """Normalize names using the same simple Genus_species rule planned for SQL."""
    organism_name = str(organism_name).strip()
    if "_" in organism_name:
        parts = [p for p in organism_name.split("_") if p]
    else:
        parts = organism_name.split()

    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    elif len(parts) == 1:
        return parts[0]
    return "unknown_species"


def read_table_auto(path: Path) -> pd.DataFrame:
    """Read CSV or TSV. Handles UTF-8 BOM."""
    suffixes = "".join(path.suffixes).lower()

    if suffixes.endswith(".tsv") or suffixes.endswith(".txt"):
        return pd.read_csv(path, sep="\t", encoding="utf-8-sig")

    if suffixes.endswith(".csv"):
        return pd.read_csv(path, encoding="utf-8-sig")

    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return pd.read_csv(handle, sep=dialect.delimiter)


def choose_first_column(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    cleaned_to_real = {str(c).strip().lower(): c for c in df.columns}
    for candidate in candidates:
        key = candidate.strip().lower()
        if key in cleaned_to_real:
            return cleaned_to_real[key]
    return None


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def write_json(obj: object, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=4, ensure_ascii=False), encoding="utf-8")


def load_manifest(manifest_path: Optional[Path]) -> Optional[pd.DataFrame]:
    if manifest_path is None:
        return None

    df = read_table_auto(manifest_path)

    accession_col = choose_first_column(
        df,
        [
            "accession",
            "assembly_accession",
            "assemblyAccession",
            "Assembly Accession",
            "genome_accession",
        ],
    )

    organism_col = choose_first_column(
        df,
        [
            "organism_name",
            "organismName",
            "Organism Name",
            "species",
            "species_name",
            "genus_species",
            "binomial",
            "binomial_2020",
        ],
    )

    if organism_col is None:
        raise ValueError(
            f"Could not find a species/organism column in manifest: {manifest_path}\n"
            f"Columns found: {list(df.columns)}"
        )

    out = pd.DataFrame()
    out["organism_name"] = df[organism_col].astype(str).str.strip()
    out["genus_species"] = out["organism_name"].map(normalize_species_name)

    if accession_col is not None:
        out["accession"] = df[accession_col].astype(str).str.strip()
    else:
        out["accession"] = pd.NA

    out = out[out["genus_species"] != "unknown_species"]
    out = out.drop_duplicates(subset=["genus_species"], keep="first")
    return out


def coerce_node_name(value: object) -> str:
    """
    Standardize node/species names from the rebuilt node table.

    Converts:
      44       -> node_44
      Node 44 -> node_44
      node_44 -> node_44
      ancestor/root -> ancestral
      Gekko kuhli -> Gekko_kuhli
    """
    if pd.isna(value):
        return ""

    text = str(value).strip()
    if not text:
        return ""

    low = text.lower().strip()

    if low in {"ancestor", "ancestral", "root"}:
        return "ancestral"

    m = re.fullmatch(r"(?:node[\s_\-]*)?(\d+)", low, flags=re.IGNORECASE)
    if m:
        return f"node_{int(m.group(1))}"

    m = re.fullmatch(r"node[\s_\-]*(\d+)", text, flags=re.IGNORECASE)
    if m:
        return f"node_{int(m.group(1))}"

    if " " in text and "_" not in text:
        return normalize_species_name(text)

    return text


def parse_node_relationships(path: Path) -> Dict[str, Dict[str, List[float]]]:
    """
    Parse node_numbers_rebuilt.csv from 99_compout_tree_to_csv.py.

    Preferred long format examples:
      node,parent
      daughter,parent
      child,parent
      node_id,parent_id
      child_node,parent_node

    Older wide format is also accepted:
      each column is a daughter node; non-empty values under that column are parent nodes.
    """
    df = read_table_auto(path)

    child_col = choose_first_column(
        df,
        [
            "daughter",
            "daughter_node",
            "child",
            "child_node",
            "node",
            "node_id",
            "descendant",
            "descendant_node",
            "son",
            "son_node",
        ],
    )
    parent_col = choose_first_column(
        df,
        [
            "parent",
            "parent_node",
            "ancestor",
            "ancestor_node",
            "mother",
            "parent_id",
        ],
    )

    relationships: Dict[str, Dict[str, List[float]]] = {}

    if child_col is not None and parent_col is not None and child_col != parent_col:
        for _, row in df.iterrows():
            child = coerce_node_name(row[child_col])
            parent = coerce_node_name(row[parent_col])

            if not child or not parent:
                continue

            relationships.setdefault(child, {})
            relationships[child].setdefault(parent, [])

    else:
        # Backward-compatible wide format from the older hand-built table.
        for daughter_raw in df.columns:
            daughter = coerce_node_name(daughter_raw)
            if not daughter:
                continue

            relationships.setdefault(daughter, {})
            for parent_raw in df[daughter_raw].dropna():
                parent = coerce_node_name(parent_raw)
                if parent:
                    relationships[daughter].setdefault(parent, [])

    if not relationships:
        raise ValueError(
            f"No node relationships could be parsed from {path}. "
            "Expected long format child/parent columns or older wide format."
        )

    return relationships


class NhPhymlGcDivergence:
    def __init__(
        self,
        nhphyml_dir: Path,
        node_numbers: Path,
        outdir: Path,
        input_tsv_dir: Optional[Path] = None,
        manifest: Optional[Path] = None,
        first_internal_node: int = 44,
        tree_suffix: str = ".phylip_nhPhymlGC.tree",
        lk_suffix: str = ".phylip_nhPhyml.lk",
        phylip_suffix: str = ".phylip",
        strict: bool = False,
    ) -> None:
        self.nhphyml_dir = nhphyml_dir
        self.node_numbers = node_numbers
        self.outdir = outdir
        self.input_tsv_dir = input_tsv_dir
        self.manifest = manifest
        self.first_internal_node = first_internal_node
        self.tree_suffix = tree_suffix
        self.lk_suffix = lk_suffix
        self.phylip_suffix = phylip_suffix
        self.strict = strict

        self.tree_files: Dict[str, Path] = {}
        self.lk_files: Dict[str, Path] = {}
        self.orthologs: List[str] = []

        self.node_relationship: Dict[str, Dict[str, List[float]]] = {}
        self.node_anc_relationship: Dict[str, List[float]] = defaultdict(list)
        self.node_gc3_total: Dict[str, List[float]] = defaultdict(list)
        self.gc3_per_gene_rows: List[Dict[str, object]] = []

        self.dij_dict: Dict[str, Dict[str, float]] = {}
        self.dianc_dict: Dict[str, float] = {}

    def index_nhphyml_files(self) -> None:
        if not self.nhphyml_dir.exists():
            raise FileNotFoundError(f"nhPhyML directory does not exist: {self.nhphyml_dir}")

        # nhPhyML output directories often contain only output files, not the
        # original .phylip inputs. A complete ortholog is therefore defined by
        # a matching pair of:
        #   <ortholog>.phylip_nhPhymlGC.tree
        #   <ortholog>.phylip_nhPhyml.lk
        # Any .phylip files are optional and used only for reporting.
        phylip_orthologs = {
            path.name[: -len(self.phylip_suffix)]
            for path in self.nhphyml_dir.iterdir()
            if path.is_file() and path.name.endswith(self.phylip_suffix)
        }

        self.tree_files = {
            path.name[: -len(self.tree_suffix)]: path
            for path in self.nhphyml_dir.iterdir()
            if path.is_file() and path.name.endswith(self.tree_suffix)
        }

        self.lk_files = {
            path.name[: -len(self.lk_suffix)]: path
            for path in self.nhphyml_dir.iterdir()
            if path.is_file() and path.name.endswith(self.lk_suffix)
        }

        complete = sorted(set(self.tree_files) & set(self.lk_files))

        if not complete:
            n_tree = len(self.tree_files)
            n_lk = len(self.lk_files)
            examples_tree = sorted(self.tree_files)[:3]
            examples_lk = sorted(self.lk_files)[:3]
            raise ValueError(
                "No complete nhPhyML tree/lk pair found. Expected matching files ending in:\n"
                f"  {self.tree_suffix}\n"
                f"  {self.lk_suffix}\n"
                f"Found {n_tree} tree files and {n_lk} lk files.\n"
                f"Example tree prefixes: {examples_tree}\n"
                f"Example lk prefixes: {examples_lk}"
            )

        missing_lk = sorted(set(self.tree_files) - set(self.lk_files))
        missing_tree = sorted(set(self.lk_files) - set(self.tree_files))

        if self.strict and (missing_tree or missing_lk):
            raise ValueError(
                f"Missing paired nhPhyML outputs. "
                f"tree_without_lk={len(missing_lk)}, lk_without_tree={len(missing_tree)}"
            )

        self.orthologs = complete
        print(
            f"[OK] Indexed nhPhyML outputs: trees={len(self.tree_files)}, "
            f"lk={len(self.lk_files)}, paired={len(self.orthologs)}, "
            f"optional_phylip_inputs_seen={len(phylip_orthologs)}"
        )

    def parse_gc_tree(self, tree_path: Path, ortholog: str) -> Dict[str, float]:
        """
        Parse GC values from an nhPhyML GC tree.

        This preserves the original traversal-based node numbering:
          first internal node -> node_<first_internal_node>
          next internal node  -> node_<first_internal_node + 1>
          etc.

        This must match node_numbers_rebuilt.csv from the same rooted tree/order.
        """
        tree = tree_path.read_text(encoding="utf-8", errors="replace").strip()
        node_gc: Dict[str, float] = {}
        node_counter = self.first_internal_node

        # Kept intentionally close to the original script's regex.
        pattern = re.compile(r"([^\(\):,\s]+).(\d+\.\d+):")

        for match in pattern.finditer(tree):
            match_string = match.group()
            raw_node = match.group(1)
            gc = float(match.group(2))

            if " " in match_string:
                node = normalize_species_name(raw_node)
            elif ")" in match_string:
                node = f"node_{node_counter}"
                node_counter += 1
            else:
                continue

            node_gc[node] = gc
            self.node_gc3_total[node].append(gc)
            self.gc3_per_gene_rows.append(
                {
                    "ortholog": ortholog,
                    "node_or_species": node,
                    "gc3": gc,
                }
            )

        if not node_gc:
            raise ValueError(f"No GC3 values were parsed from tree file: {tree_path}")

        return node_gc

    def parse_ancestral_gc(self, lk_path: Path, ortholog: str) -> float:
        lines = lk_path.read_text(encoding="utf-8", errors="replace").splitlines()

        for line in lines:
            if "Ancestral" in line:
                numbers = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)
                if numbers:
                    value = float(numbers[-1])
                    ancestral_gc = value * 100 if value <= 1 else value

                    self.node_gc3_total["ancestral"].append(ancestral_gc)
                    self.gc3_per_gene_rows.append(
                        {
                            "ortholog": ortholog,
                            "node_or_species": "ancestral",
                            "gc3": ancestral_gc,
                        }
                    )
                    return ancestral_gc

        raise ValueError(f"Could not find ancestral GC value in {lk_path}")

    def calculate_for_ortholog(self, ortholog: str) -> None:
        node_gc = self.parse_gc_tree(self.tree_files[ortholog], ortholog)
        node_gc["ancestral"] = self.parse_ancestral_gc(self.lk_files[ortholog], ortholog)

        missing_pairs = []

        # d_ij: node/species i versus parent j.
        for child, parent_dict in self.node_relationship.items():
            if child not in node_gc:
                missing_pairs.append((child, "CHILD_MISSING"))
                continue

            for parent in parent_dict:
                if parent not in node_gc:
                    missing_pairs.append((child, parent))
                    continue

                parent_dict[parent].append((node_gc[child] - node_gc[parent]) ** 2)

        if missing_pairs and self.strict:
            preview = ", ".join([f"{a}->{b}" for a, b in missing_pairs[:10]])
            raise KeyError(
                f"Missing node/species names for ortholog {ortholog}. "
                f"First missing pairs: {preview}"
            )

        # d_ianc: terminal taxa/species versus ancestral GC.
        for node in node_gc:
            if node == "ancestral" or node.startswith("node_"):
                continue
            self.node_anc_relationship[node].append(
                (node_gc[node] - node_gc["ancestral"]) ** 2
            )

    def run(self) -> None:
        self.index_nhphyml_files()
        self.node_relationship = parse_node_relationships(self.node_numbers)

        print(f"[OK] Found {len(self.orthologs)} complete orthologs")
        print(f"[OK] Parsed {len(self.node_relationship)} child nodes/species from {self.node_numbers}")

        for ortholog in self.orthologs:
            print(f"[INFO] Calculating {ortholog}")
            self.calculate_for_ortholog(ortholog)

        self.calculate_summaries()

    def calculate_summaries(self) -> None:
        for child, parent_dict in self.node_relationship.items():
            self.dij_dict[child] = {}
            for parent, values in parent_dict.items():
                self.dij_dict[child][parent] = math.sqrt(sum(values)) if values else math.nan

        for species, values in self.node_anc_relationship.items():
            self.dianc_dict[species] = math.sqrt(sum(values)) if values else math.nan

    def export_gc3(self) -> pd.DataFrame:
        per_gene = pd.DataFrame(self.gc3_per_gene_rows)
        write_tsv(per_gene, self.outdir / "gc3_per_gene.tsv")

        rows = []
        for node, values in sorted(self.node_gc3_total.items()):
            rows.append(
                {
                    "node_or_species": node,
                    "mean_gc3": sum(values) / len(values),
                    "sd_gc3": statistics.stdev(values) if len(values) > 1 else math.nan,
                    "n_orthologs": len(values),
                }
            )

        summary = pd.DataFrame(rows)
        write_tsv(summary, self.outdir / "gc3_summary.tsv")
        return summary

    def export_dij(self) -> pd.DataFrame:
        write_json(self.dij_dict, self.outdir / "dij_results.json")

        rows = []
        for child, parent_dict in self.dij_dict.items():
            for parent, dij in parent_dict.items():
                rows.append(
                    {
                        "node_or_species": child,
                        "parent": parent,
                        "dij": dij,
                    }
                )

        df = pd.DataFrame(rows)
        write_tsv(df, self.outdir / "dij_results.tsv")
        return df

    def export_dianc(self) -> pd.DataFrame:
        write_json(self.dianc_dict, self.outdir / "dianc_results.json")

        rows = [
            {
                "genus_species": normalize_species_name(species),
                "dianc": dianc,
            }
            for species, dianc in sorted(self.dianc_dict.items())
        ]
        df = pd.DataFrame(rows)
        write_tsv(df, self.outdir / "dianc_results.tsv")
        return df

    def read_input_tsv(self, filename: str) -> Optional[pd.DataFrame]:
        if self.input_tsv_dir is None:
            return None
        path = self.input_tsv_dir / filename
        if not path.exists():
            return None
        return pd.read_csv(path, sep="\t")

    def export_master(self, gc3_summary: pd.DataFrame, dianc_df: pd.DataFrame) -> None:
        master = gc3_summary.copy()
        master = master[
            (master["node_or_species"] != "ancestral")
            & (~master["node_or_species"].astype(str).str.startswith("node_"))
        ].copy()

        master = master.rename(columns={"node_or_species": "genus_species"})
        master["genus_species"] = master["genus_species"].map(normalize_species_name)

        master = master.merge(dianc_df, on="genus_species", how="left")

        genome_sizes = self.read_input_tsv("genome_sizes.tsv")
        if genome_sizes is not None:
            genome_sizes = genome_sizes.copy()
            genome_sizes["genus_species"] = genome_sizes["genus_species"].map(normalize_species_name)
            keep = [
                c
                for c in [
                    "genus_species",
                    "species_name",
                    "mean_c_value_pg",
                    "mean_genome_mb",
                    "mean_genome_gb",
                    "genome_size",
                ]
                if c in genome_sizes.columns
            ]
            master = master.merge(
                genome_sizes[keep].drop_duplicates("genus_species"),
                on="genus_species",
                how="left",
            )

        species_mass = self.read_input_tsv("species_mass.tsv")
        if species_mass is not None:
            species_mass = species_mass.copy()
            species_mass["genus_species"] = species_mass["genus_species"].map(normalize_species_name)
            keep = [
                c
                for c in [
                    "genus_species",
                    "mass_meiri",
                    "mass_title",
                    "mass_preferred",
                    "mass_source_preferred",
                ]
                if c in species_mass.columns
            ]
            master = master.merge(
                species_mass[keep].drop_duplicates("genus_species"),
                on="genus_species",
                how="left",
            )

        manifest_df = load_manifest(self.manifest)
        if manifest_df is not None:
            # Use manifest as final project species frame so missing GC values remain visible.
            master = manifest_df.merge(master, on="genus_species", how="left")
        else:
            if "organism_name" not in master.columns:
                master["organism_name"] = master["genus_species"].str.replace("_", " ", regex=False)
            if "accession" not in master.columns:
                master["accession"] = pd.NA

        preferred_order = [
            "genus_species",
            "organism_name",
            "accession",
            "species_name",
            "mean_gc3",
            "sd_gc3",
            "n_orthologs",
            "dianc",
            "mean_c_value_pg",
            "mean_genome_mb",
            "mean_genome_gb",
            "genome_size",
            "mass_meiri",
            "mass_title",
            "mass_preferred",
            "mass_source_preferred",
        ]
        ordered = [c for c in preferred_order if c in master.columns]
        remaining = [c for c in master.columns if c not in ordered]
        master = master[ordered + remaining]

        write_tsv(master, self.outdir / "master_output.tsv")

    def export_all(self) -> None:
        self.outdir.mkdir(parents=True, exist_ok=True)

        gc3_summary = self.export_gc3()
        self.export_dij()
        dianc_df = self.export_dianc()
        self.export_master(gc3_summary, dianc_df)

        print(f"[DONE] Wrote outputs to {self.outdir}")
        print(f"[INFO] Orthologs analyzed: {len(self.orthologs)}")
        print(f"[INFO] d_ij entries: {sum(len(v) for v in self.dij_dict.values())}")
        print(f"[INFO] d_ianc species: {len(self.dianc_dict)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate GC3 and GC divergence values from nhPhyML outputs."
    )
    parser.add_argument(
        "genomes_dir",
        type=Path,
        help="Path to the top-level genomes directory.",
    )
    parser.add_argument(
        "--nhphyml-dir",
        type=Path,
        required=True,
        help="Directory containing .phylip, .phylip_nhPhymlGC.tree, and .phylip_nhPhyml.lk files.",
    )
    parser.add_argument(
        "--node-numbers",
        type=Path,
        required=True,
        help="node_numbers_rebuilt.csv generated by 99_compout_tree_to_csv.py.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="Optional project manifest CSV/TSV for master_output.tsv ordering/filtering.",
    )
    parser.add_argument(
        "--input-tsv-dir",
        type=Path,
        default=None,
        help="Directory containing genome_sizes.tsv and species_mass.tsv. Default: <genomes_dir>/records/gc_metrics_inputs",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory. Default: <genomes_dir>/records/gc_metrics_outputs",
    )
    parser.add_argument(
        "--first-internal-node",
        type=int,
        default=44,
        help="First internal node number assigned while parsing GC trees. Default: 44.",
    )
    parser.add_argument(
        "--tree-suffix",
        default=".phylip_nhPhymlGC.tree",
        help="Suffix for nhPhyML GC tree files.",
    )
    parser.add_argument(
        "--lk-suffix",
        default=".phylip_nhPhyml.lk",
        help="Suffix for nhPhyML lk files.",
    )
    parser.add_argument(
        "--phylip-suffix",
        default=".phylip",
        help="Suffix for phylip files used to identify complete orthologs.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail if any node/species from node_numbers is missing in any parsed tree.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    genomes_dir = args.genomes_dir.expanduser().resolve()
    if not genomes_dir.exists():
        raise FileNotFoundError(f"genomes_dir does not exist: {genomes_dir}")

    input_tsv_dir = args.input_tsv_dir
    if input_tsv_dir is None:
        input_tsv_dir = genomes_dir / "records" / "gc_metrics_inputs"
    input_tsv_dir = input_tsv_dir.expanduser().resolve()

    outdir = args.outdir
    if outdir is None:
        outdir = genomes_dir / "records" / "gc_metrics_outputs"
    outdir = outdir.expanduser().resolve()

    manifest = args.manifest.expanduser().resolve() if args.manifest else None

    calculator = NhPhymlGcDivergence(
        nhphyml_dir=args.nhphyml_dir.expanduser().resolve(),
        node_numbers=args.node_numbers.expanduser().resolve(),
        outdir=outdir,
        input_tsv_dir=input_tsv_dir if input_tsv_dir.exists() else None,
        manifest=manifest,
        first_internal_node=args.first_internal_node,
        tree_suffix=args.tree_suffix,
        lk_suffix=args.lk_suffix,
        phylip_suffix=args.phylip_suffix,
        strict=args.strict,
    )

    calculator.run()
    calculator.export_all()


if __name__ == "__main__":
    main()
