#!/usr/bin/env python3
"""
15_calculate_gc_divergence_v2.py

Calculate GC3 summaries and GC divergence values from nhPhyML output files.

Important update:
  - --number-of-species is now a REQUIRED named argument.
  - first_internal_node is calculated as:
        first_internal_node = number_of_species + 1
    Example:
        43 species -> first internal node is node_44
        58 species -> first internal node is node_59

This matters when running different subsets of the same manifest.

Required inputs:
  1. genomes_dir
  2. --number-of-species
  3. --nhphyml-dir
  4. --node-numbers

Optional:
  --manifest
      Used for accession/organism metadata and optional project ordering.
      By default, the master output keeps only species actually found in the
      nhPhyML GC tree files. This avoids blank regression rows when the manifest
      contains species not included in the subset.
      Use --include-all-manifest-species if you want blank rows for missing
      species.

  --input-tsv-dir
      Default: <genomes_dir>/records/gc_metrics_inputs
      Reads genome_sizes.tsv, species_mass.tsv, and thermal_limits.tsv if present.

Outputs:
  <outdir>/gc3_per_gene.tsv
  <outdir>/gc3_summary.tsv
  <outdir>/dij_results.json
  <outdir>/dij_results.tsv
  <outdir>/dianc_results.json
  <outdir>/dianc_results.tsv
  <outdir>/master_output.tsv
  <outdir>/master_merge_diagnostics.tsv
  <outdir>/species_merge_report.tsv
  <outdir>/master_complete_cases.tsv

Example:
  python 04_calculate_gc_divergence.py /path/to/genomes \
    --number-of-species 58 \
    --nhphyml-dir /path/to/07_nhphyml \
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
from typing import Dict, List, Optional, Tuple

import pandas as pd


def normalize_species_name(organism_name: str) -> str:
    """
    Normalize species names to the SQL-style merge key used throughout this
    workflow.

    Important:
      - The returned key is lowercase.
      - Only the first two tokens are retained.
      - organism_name remains available separately for readable display.

    Examples:
      "Anolis carolinensis" -> "anolis_carolinensis"
      "Anolis_carolinensis" -> "anolis_carolinensis"
      "Aspidoscelis inornatus arizonae" -> "aspidoscelis_inornatus"
    """
    organism_name = str(organism_name).strip()
    if "_" in organism_name:
        parts = [p for p in organism_name.split("_") if p]
    else:
        parts = organism_name.split()

    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}".lower()
    elif len(parts) == 1:
        return parts[0].lower()
    return "unknown_species"


def read_table_auto(path: Path) -> pd.DataFrame:
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



def ensure_genus_species_key(df: pd.DataFrame, context: str) -> pd.DataFrame:
    """
    Ensure a dataframe has the legacy/internal genus_species merge key.

    Newer natural-history staging files use species_normalized as the canonical
    key. Older GC divergence code used genus_species internally. This helper
    accepts either spelling and normalizes the key for safe merges.
    """
    df = df.copy()

    if "genus_species" in df.columns:
        key_col = "genus_species"
    elif "species_normalized" in df.columns:
        key_col = "species_normalized"
    elif "sequence_id" in df.columns:
        key_col = "sequence_id"
    else:
        raise ValueError(
            f"Could not find a species merge key in {context}. Expected one of: "
            "genus_species, species_normalized, sequence_id. "
            f"Columns found: {list(df.columns)}"
        )

    df["genus_species"] = df[key_col].map(normalize_species_name)
    df = df[df["genus_species"] != "unknown_species"]
    return df


def summarize_gc4_species_file(path: Path) -> pd.DataFrame:
    """
    Summarize GC4 by species.

    Expected input contains:
      - species key as genus_species, species_normalized, or sequence_id
      - gc4 or mean_gc4

    If multiple rows per species are present, mean_gc4 and sd_gc4 are calculated
    across those rows. If one row per species is present, sd_gc4 is NA and
    n_gc4_orthologs is 1.
    """
    df = read_table_auto(path)
    df = ensure_genus_species_key(df, str(path))

    gc4_col = choose_first_column(df, ["gc4", "mean_gc4"])
    if gc4_col is None:
        raise ValueError(
            f"Could not find gc4 or mean_gc4 column in {path}. "
            f"Columns found: {list(df.columns)}"
        )

    df = df[["genus_species", gc4_col]].copy()
    df = df.rename(columns={gc4_col: "gc4"})
    df["gc4"] = pd.to_numeric(df["gc4"], errors="coerce")
    df = df.dropna(subset=["genus_species", "gc4"])

    summary = (
        df.groupby("genus_species", as_index=False)
        .agg(
            mean_gc4=("gc4", "mean"),
            sd_gc4=("gc4", lambda x: statistics.stdev(x) if len(x) > 1 else math.nan),
            n_gc4_orthologs=("gc4", "count"),
        )
    )

    return summary


def coerce_node_name(value: object) -> str:
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

    if " " in text and "_" not in text:
        return normalize_species_name(text)

    return text


def load_manifest(manifest_path: Optional[Path]) -> Optional[pd.DataFrame]:
    """
    Load project manifest and return accession/organism metadata.

    This is deliberately permissive because project manifests may evolve over
    time. The most important output columns are:

      genus_species
      organism_name
      accession

    For the merge key, species_key is preferred when present because it should
    already match the SQL-style Genus_species convention.
    """
    if manifest_path is None:
        return None

    df = read_table_auto(manifest_path)

    accession_col = choose_first_column(
        df,
        [
            "accession",
            "accession_root",
            "assembly_accession",
            "assemblyAccession",
            "Assembly Accession",
            "genome_accession",
            "genome_accession_id",
            "ncbi_accession",
        ],
    )

    species_key_col = choose_first_column(
        df,
        [
            "species_key",
            "genus_species",
            "genus_species_key",
            "normalized_species",
            "normalized_name",
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
            "binomial",
            "binomial_2020",
        ],
    )

    if species_key_col is None and organism_col is None:
        raise ValueError(
            f"Could not find species_key or species/organism column in manifest: {manifest_path}"
            f"Columns found: {list(df.columns)}"
        )

    out = pd.DataFrame()

    if species_key_col is not None:
        out["genus_species"] = df[species_key_col].astype(str).map(normalize_species_name)
    else:
        out["genus_species"] = df[organism_col].astype(str).map(normalize_species_name)

    if organism_col is not None:
        out["organism_name"] = df[organism_col].astype(str).str.strip()
    else:
        out["organism_name"] = out["genus_species"].str.replace("_", " ", regex=False)

    if accession_col is not None:
        out["accession"] = df[accession_col].astype(str).str.strip()
    else:
        out["accession"] = pd.NA

    # Keep a useful debug column so merge problems are easier to inspect.
    out["manifest_merge_key"] = out["genus_species"]

    out = out[out["genus_species"] != "unknown_species"]
    out = out[out["genus_species"].astype(str).str.len() > 0]
    out = out.drop_duplicates(subset=["genus_species"], keep="first")
    return out




def default_compleasm_metadata(genomes_dir: Path) -> Path:
    return genomes_dir / "records" / "compleasm" / "records" / "metadata.csv"


def load_species_metadata_from_manifest_or_compleasm(
    manifest_path: Optional[Path],
    genomes_dir: Path,
) -> Optional[pd.DataFrame]:
    """
    Load species metadata for accession/organism_name.

    Priority:
      1. --manifest if provided and readable
      2. <genomes_dir>/records/compleasm/records/metadata.csv if present

    This allows accession to remain populated even when the GC tree subset has
    species names but the manifest merge key differs slightly.
    """
    frames = []

    if manifest_path is not None and manifest_path.exists():
        frames.append(load_manifest(manifest_path))

    compleasm_path = default_compleasm_metadata(genomes_dir)
    if compleasm_path.exists():
        try:
            frames.append(load_manifest(compleasm_path))
        except Exception as exc:
            print(f"[WARN] Could not parse Compleasm metadata for accession fallback: {compleasm_path}: {exc}")

    frames = [f for f in frames if f is not None and not f.empty]
    if not frames:
        return None

    combined = pd.concat(frames, ignore_index=True)

    # Prefer rows that have accession populated.
    combined["_has_accession"] = combined["accession"].notna() & ~combined["accession"].astype(str).isin(["", "nan", "None"])
    combined = combined.sort_values("_has_accession", ascending=False)
    combined = combined.drop_duplicates(subset=["genus_species"], keep="first")
    combined = combined.drop(columns=["_has_accession"])
    return combined




def parse_node_relationships(path: Path) -> Dict[str, Dict[str, List[float]]]:
    """
    Parse node_numbers_rebuilt.csv from 99_compout_tree_to_csv.py.

    Supported long columns include:
      child,parent
      daughter,parent
      node,parent
      node_id,parent_id
      son,parent

    Older wide format is also supported:
      each column = daughter/child, values = parent(s)
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
        for child_raw in df.columns:
            child = coerce_node_name(child_raw)
            if not child:
                continue

            relationships.setdefault(child, {})
            for parent_raw in df[child_raw].dropna():
                parent = coerce_node_name(parent_raw)
                if parent:
                    relationships[child].setdefault(parent, [])

    if not relationships:
        raise ValueError(
            f"No node relationships could be parsed from {path}. "
            "Expected long child/parent columns or older wide format."
        )

    return relationships


class NhPhymlGcDivergence:
    def __init__(
        self,
        genomes_dir: Path,
        number_of_species: int,
        nhphyml_dir: Path,
        node_numbers: Path,
        outdir: Path,
        input_tsv_dir: Optional[Path] = None,
        manifest: Optional[Path] = None,
        include_all_manifest_species: bool = False,
        tree_suffix: str = ".phylip_nhPhymlGC.tree",
        lk_suffix: str = ".phylip_nhPhyml.lk",
        strict: bool = False,
        gc4_species_summary: Optional[Path] = None,
    ) -> None:
        self.genomes_dir = genomes_dir
        self.number_of_species = number_of_species
        self.first_internal_node = number_of_species + 1
        self.nhphyml_dir = nhphyml_dir
        self.node_numbers = node_numbers
        self.outdir = outdir
        self.input_tsv_dir = input_tsv_dir
        self.manifest = manifest
        self.include_all_manifest_species = include_all_manifest_species
        self.tree_suffix = tree_suffix
        self.lk_suffix = lk_suffix
        self.strict = strict
        self.gc4_species_summary = gc4_species_summary

        self.tree_files: Dict[str, Path] = {}
        self.lk_files: Dict[str, Path] = {}
        self.orthologs: List[str] = []

        self.node_relationship: Dict[str, Dict[str, List[float]]] = {}
        self.node_anc_relationship: Dict[str, List[float]] = defaultdict(list)
        self.node_gc3_total: Dict[str, List[float]] = defaultdict(list)
        self.gc3_per_gene_rows: List[Dict[str, object]] = []

        self.dij_dict: Dict[str, Dict[str, float]] = {}
        self.dianc_dict: Dict[str, float] = {}
        self.displacement_examples: List[Dict[str, object]] = []

    def index_nhphyml_files(self) -> None:
        if not self.nhphyml_dir.exists():
            raise FileNotFoundError(f"nhPhyML directory does not exist: {self.nhphyml_dir}")

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
            tree_preview = "\n".join(sorted([p.name for p in self.nhphyml_dir.iterdir() if p.is_file()])[:10])
            raise ValueError(
                "No complete ortholog set found. Expected paired files ending in:\n"
                f"  {self.tree_suffix}\n"
                f"  {self.lk_suffix}\n\n"
                f"First files seen in --nhphyml-dir:\n{tree_preview}"
            )

        missing_lk = sorted(set(self.tree_files) - set(self.lk_files))
        missing_tree = sorted(set(self.lk_files) - set(self.tree_files))

        if self.strict and (missing_tree or missing_lk):
            raise ValueError(
                f"Missing paired nhPhyML outputs. "
                f"tree_without_lk={len(missing_lk)}, lk_without_tree={len(missing_tree)}"
            )

        self.orthologs = complete

    def parse_gc_tree(self, tree_path: Path, ortholog: str) -> Dict[str, float]:
        """
        Parse GC values from an nhPhyML GC tree.

        Internal nodes are numbered by traversal order:
            node_<number_of_species + 1>, node_<number_of_species + 2>, ...
        """
        tree = tree_path.read_text(encoding="utf-8", errors="replace").strip()
        node_gc: Dict[str, float] = {}
        node_counter = self.first_internal_node

        # nhPhyML GC trees commonly annotate branch GC like:
        # taxon_name 42.01:0.123
        # )42.01:0.123
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

        for child, parent_dict in self.node_relationship.items():
            if child not in node_gc:
                missing_pairs.append((child, "CHILD_MISSING"))
                continue

            for parent in parent_dict:
                if parent not in node_gc:
                    missing_pairs.append((child, parent))
                    continue

                squared_displacement = (node_gc[child] - node_gc[parent]) ** 2
                parent_dict[parent].append(squared_displacement)

                if len(self.displacement_examples) < 200:
                    self.displacement_examples.append(
                        {
                            "ortholog": ortholog,
                            "child": child,
                            "parent": parent,
                            "child_gc3": node_gc[child],
                            "parent_gc3": node_gc[parent],
                            "difference_child_minus_parent": node_gc[child] - node_gc[parent],
                            "squared_displacement": squared_displacement,
                        }
                    )

        if missing_pairs and self.strict:
            preview = ", ".join([f"{a}->{b}" for a, b in missing_pairs[:10]])
            raise KeyError(
                f"Missing node/species names for ortholog {ortholog}. "
                f"First missing pairs: {preview}"
            )

        for node in node_gc:
            if node == "ancestral" or node.startswith("node_"):
                continue
            self.node_anc_relationship[node].append(
                (node_gc[node] - node_gc["ancestral"]) ** 2
            )

    def run(self) -> None:
        self.index_nhphyml_files()
        self.node_relationship = parse_node_relationships(self.node_numbers)

        print(f"[OK] number_of_species={self.number_of_species}")
        print(f"[OK] first_internal_node=node_{self.first_internal_node}")
        print(f"[OK] Found {len(self.orthologs)} paired nhPhyML tree/lk orthologs")
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

    def read_input_tsv(self, filename: str) -> Optional[pd.DataFrame]:
        if self.input_tsv_dir is None:
            return None
        path = self.input_tsv_dir / filename
        if not path.exists():
            return None
        return pd.read_csv(path, sep="\t")

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

    def export_master(self, gc3_summary: pd.DataFrame, dianc_df: pd.DataFrame) -> None:
        """
        Export a GC-only species summary.

        Final core columns:
          genus_species
          accession
          mean_gc3
          sd_gc3
          n_orthologs
          mean_gc4
          sd_gc4
          dianc
          manifest_merge_key

        Phenotypic, physiological, and genome-size traits are intentionally
        excluded here. They should be joined downstream using SQL-ready TSVs.
        """
        master = gc3_summary.copy()
        master = master[
            (master["node_or_species"] != "ancestral")
            & (~master["node_or_species"].astype(str).str.startswith("node_"))
        ].copy()

        master = master.rename(columns={"node_or_species": "genus_species"})
        master["genus_species"] = master["genus_species"].map(normalize_species_name)

        master = master.merge(dianc_df, on="genus_species", how="left")

        if self.gc4_species_summary is not None:
            if not self.gc4_species_summary.exists():
                raise FileNotFoundError(f"GC4 species summary does not exist: {self.gc4_species_summary}")

            gc4_summary = summarize_gc4_species_file(self.gc4_species_summary)
            write_tsv(gc4_summary, self.outdir / "gc4_summary.tsv")

            master = master.merge(gc4_summary, on="genus_species", how="left")
        else:
            master["mean_gc4"] = pd.NA
            master["sd_gc4"] = pd.NA
            master["n_gc4_orthologs"] = pd.NA

        manifest_df = load_species_metadata_from_manifest_or_compleasm(self.manifest, self.genomes_dir)
        if manifest_df is not None:
            manifest_keep = [
                c
                for c in ["genus_species", "accession", "manifest_merge_key"]
                if c in manifest_df.columns
            ]

            if self.include_all_manifest_species:
                master = manifest_df[manifest_keep].merge(master, on="genus_species", how="left")
            else:
                master = master.merge(
                    manifest_df[manifest_keep].drop_duplicates("genus_species"),
                    on="genus_species",
                    how="left",
                )
        else:
            master["accession"] = pd.NA
            master["manifest_merge_key"] = master["genus_species"]

        for col in ["accession", "manifest_merge_key"]:
            if col not in master.columns:
                master[col] = pd.NA
            master[col] = master[col].replace({"nan": pd.NA, "None": pd.NA, "": pd.NA})

        final_cols = [
            "genus_species",
            "accession",
            "mean_gc3",
            "sd_gc3",
            "n_orthologs",
            "mean_gc4",
            "sd_gc4",
            "n_gc4_orthologs",
            "dianc",
            "manifest_merge_key",
        ]

        for col in final_cols:
            if col not in master.columns:
                master[col] = pd.NA

        master = master[final_cols]

        write_tsv(master, self.outdir / "master_output.tsv")
        write_tsv(master, self.outdir / "gc_summary_output.tsv")

        self.export_merge_diagnostics(master, manifest_df)

    def export_merge_diagnostics(
        self,
        master: pd.DataFrame,
        manifest_df: Optional[pd.DataFrame],
    ) -> None:
        rows = []

        gc_species = set(
            normalize_species_name(x)
            for x in self.node_gc3_total
            if x != "ancestral" and not str(x).startswith("node_")
        )

        def add(category: str, value: object) -> None:
            rows.append({"diagnostic": category, "value": value})

        add("number_of_species_argument", self.number_of_species)
        add("first_internal_node_used", f"node_{self.first_internal_node}")
        add("orthologs_analyzed_for_gc3", len(self.orthologs))
        add("species_with_gc3_results", len(gc_species))
        add("rows_in_gc_summary_output", len(master))
        add("rows_with_accession", int(master["accession"].notna().sum()) if "accession" in master.columns else 0)
        add("rows_with_mean_gc3", int(master["mean_gc3"].notna().sum()) if "mean_gc3" in master.columns else 0)
        add("rows_with_sd_gc3", int(master["sd_gc3"].notna().sum()) if "sd_gc3" in master.columns else 0)
        add("rows_with_dianc", int(master["dianc"].notna().sum()) if "dianc" in master.columns else 0)
        add("rows_with_mean_gc4", int(master["mean_gc4"].notna().sum()) if "mean_gc4" in master.columns else 0)
        add("rows_with_sd_gc4", int(master["sd_gc4"].notna().sum()) if "sd_gc4" in master.columns else 0)
        add("rows_with_n_gc4_orthologs", int(master["n_gc4_orthologs"].notna().sum()) if "n_gc4_orthologs" in master.columns else 0)

        if manifest_df is not None:
            manifest_species = set(manifest_df["genus_species"].dropna().map(normalize_species_name))
            add("species_in_manifest_or_compleasm_metadata", len(manifest_species))
            add("gc_species_missing_from_manifest_or_compleasm_metadata", len(gc_species - manifest_species))
            add("manifest_or_compleasm_species_missing_gc_results", len(manifest_species - gc_species))

        write_tsv(pd.DataFrame(rows), self.outdir / "master_merge_diagnostics.tsv")


    def export_displacement_examples(self) -> None:
        """
        Write a readable text file with a few example per-gene displacement
        calculations.

        The examples are chosen from:
          - one internal node near the root/highest numbered internal nodes
          - one Anolis terminal node if present
          - otherwise the first few available relationships
        """
        if not self.displacement_examples:
            path = self.outdir / "displacement_examples.txt"
            path.write_text("No displacement examples were recorded.\n", encoding="utf-8")
            return

        df = pd.DataFrame(self.displacement_examples)

        chosen_frames = []

        # Near-root-ish internal node: with this numbering scheme, the root or
        # near-root internal nodes are usually among the highest node numbers.
        internal = df[df["child"].astype(str).str.startswith("node_")].copy()
        if not internal.empty:
            internal["node_number"] = (
                internal["child"]
                .astype(str)
                .str.extract(r"node_(\d+)", expand=False)
                .astype(float)
            )
            max_node = internal["node_number"].max()
            chosen_frames.append(internal[internal["node_number"] == max_node].head(3))

        # Anolis example if represented in this subset.
        anolis = df[df["child"].astype(str).str.lower().str.startswith("anolis_")].copy()
        if not anolis.empty:
            chosen_frames.append(anolis.head(3))

        if chosen_frames:
            chosen = pd.concat(chosen_frames, ignore_index=True).drop_duplicates(
                subset=["ortholog", "child", "parent"]
            )
        else:
            chosen = df.head(6).copy()

        lines = []
        lines.append("Example GC3 displacement calculations")
        lines.append("=" * 42)
        lines.append("")
        lines.append(f"number_of_species = {self.number_of_species}")
        lines.append(f"first_internal_node = node_{self.first_internal_node}")
        lines.append("")
        lines.append(
            "Formula per ortholog k: squared_displacement = (GC3_child,k - GC3_parent,k)^2"
        )
        lines.append("Final d_ij over all orthologs: sqrt(sum_k squared_displacement)")
        lines.append("")

        for i, row in chosen.iterrows():
            child = row["child"]
            parent = row["parent"]
            child_gc3 = float(row["child_gc3"])
            parent_gc3 = float(row["parent_gc3"])
            diff = float(row["difference_child_minus_parent"])
            sq = float(row["squared_displacement"])

            lines.append(f"Example {i + 1}")
            lines.append("-" * 20)
            lines.append(f"ortholog: {row['ortholog']}")
            lines.append(f"child i:  {child}")
            lines.append(f"parent j: {parent}")
            lines.append(f"GC3_i,k:  {child_gc3:.6f}")
            lines.append(f"GC3_j,k:  {parent_gc3:.6f}")
            lines.append(f"GC3_i,k - GC3_j,k = {child_gc3:.6f} - {parent_gc3:.6f} = {diff:.6f}")
            lines.append(f"(GC3_i,k - GC3_j,k)^2 = {sq:.6f}")

            if child in self.dij_dict and parent in self.dij_dict[child]:
                lines.append(f"Final d_ij for {child} -> {parent}: {self.dij_dict[child][parent]:.6f}")

            lines.append("")

        path = self.outdir / "displacement_examples.txt"
        path.write_text("\n".join(lines), encoding="utf-8")


    def export_all(self) -> None:
        self.outdir.mkdir(parents=True, exist_ok=True)

        gc3_summary = self.export_gc3()
        self.export_dij()
        dianc_df = self.export_dianc()
        self.export_master(gc3_summary, dianc_df)
        self.export_displacement_examples()

        print(f"[DONE] Wrote outputs to {self.outdir}")
        print(f"[INFO] Orthologs analyzed: {len(self.orthologs)}")
        print(f"[INFO] Species with GC results: {len(self.dianc_dict)}")


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
        "--number-of-species",
        dest="number_of_species",
        type=int,
        required=True,
        help=(
            "Number of terminal species/taxa in this nhPhyML subset. "
            "The first internal node will be node_(number_of_species + 1)."
        ),
    )
    parser.add_argument(
        "--nhphyml-dir",
        type=Path,
        required=True,
        help="Directory containing paired .phylip_nhPhymlGC.tree and .phylip_nhPhyml.lk files.",
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
        help="Optional project manifest CSV/TSV for accession/organism metadata.",
    )
    parser.add_argument(
        "--include-all-manifest-species",
        action="store_true",
        help=(
            "Include every species in the manifest in master_output.tsv. "
            "By default, only species found in the GC tree subset are included."
        ),
    )
    parser.add_argument(
        "--input-tsv-dir",
        type=Path,
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--gc4-species-summary",
        type=Path,
        default=None,
        help=(
            "Optional gc4_species_summary.tsv from 16_gc4_from_alignment.py. "
            "The script summarizes gc4 by species as mean_gc4, sd_gc4, and n_gc4_orthologs."
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory. Default: <genomes_dir>/records/gc_metrics_outputs",
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

    if args.number_of_species < 2:
        raise ValueError("number_of_species must be at least 2")

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
        genomes_dir=genomes_dir,
        number_of_species=args.number_of_species,
        nhphyml_dir=args.nhphyml_dir.expanduser().resolve(),
        node_numbers=args.node_numbers.expanduser().resolve(),
        outdir=outdir,
        input_tsv_dir=input_tsv_dir if input_tsv_dir.exists() else None,
        manifest=manifest,
        include_all_manifest_species=args.include_all_manifest_species,
        tree_suffix=args.tree_suffix,
        lk_suffix=args.lk_suffix,
        strict=args.strict,
        gc4_species_summary=args.gc4_species_summary.expanduser().resolve() if args.gc4_species_summary else None,
    )

    calculator.run()
    calculator.export_all()


if __name__ == "__main__":
    main()
