#!/usr/bin/env python3
"""15_calculate_gc_divergence_v6.py

Summarize nhPhyML GC3 and independently computed fourfold-degenerate GC4.

* nhPhyML GC3 exports use proportions (0..1); per-ortholog arithmetic means
  and sample standard deviations are unweighted.
* Historical d_ij and d_ianc = sqrt(sum of squared per-ortholog GC3
  displacements) retain the **percentage-point** scale of v5.
* Original nhPhyML terminal labels are preserved, never regenerated from a
  manifest or lowercased. The manifest provides independent validation and
  accession assignments.
* GC4 comes from both TSVs emitted by 14_gc4_from_alignment_v3.py; the
  per-ortholog file supplies unweighted sample SD and coverage diagnostics.
* The consolidated master_output.tsv has exactly 12 columns.

Usage:
    python 15_calculate_gc_divergence_v6.py GENOMES_DIR \\
       --number-of-species 80 --nhphyml-dir DIR --node-numbers FILE \\
       --manifest FILE --gc4-detailed detailed.tsv \\
       --gc4-species-summary species_summary.tsv --outdir DIR

Optional --alignment-dir enables direct per-ortholog FASTA/PHYLIP label
validation; it does not recalculate GC. --strict treats unexpected data
mismatches as errors. Valid no-fourfold-site GC4 exclusions are reported but
are not themselves considered unexpected.

Requires Python 3.9+ and pandas.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd

MASTER_COLS = [
    "nhphyml_run_pk", "run_name", "alignment_name", "species_key_short",
    "accession", "mean_gc3_nhphyml_unweighted", "sd_gc3_nhphyml_unweighted",
    "mean_gc4_unweighted", "sd_gc4_unweighted", "dianc", "n_orthologs",
    "created_at",
]
GC3_DETAIL_COLS = ["odb12_id", "node_or_species_key_short", "gc3_nhphyml"]
GC3_SUMMARY_COLS = ["node_or_species_key_short", "mean_gc3_nhphyml_unweighted",
                    "sd_gc3_nhphyml_unweighted", "n_orthologs"]
GC4_SUMMARY_COLS = ["species_key_short", "mean_gc4_unweighted", "sd_gc4_unweighted",
                    "mean_gc4_weighted", "n_gc4_values"]
DIAGNOSTIC_DETAIL_COLS = ["category", "species_key_short", "odb12_id",
                          "source_value", "reference_value", "details"]
NUMBER = r"(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
TERMINAL_GC = re.compile(r"(?<=[(,])\s*([^\s(),:;]+)\s+(" + NUMBER + r")\s*:")
INTERNAL_GC = re.compile(r"\s*(" + NUMBER + r")(?=\s*[:,;)])")
FLOAT_NUMBERS = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")
MISSING = {"", "nan", "none", "null", "<na>", "\\n", "\\N"}
GC4_REJECT_STATUSES = {"fail", "failed", "error", "invalid", "skipped", "no_sites",
                       "no_gc4", "no_gc4_sites", "no_retained_sites", "missing"}
GC4_GOOD_STATUSES = {"", "ok", "pass", "passed", "valid", "success", "partial"}
GC4_EXPECTED_NO_SITE_STATUSES = {"no_sites", "no_gc4", "no_gc4_sites", "no_retained_sites"}


def read_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Required input file does not exist: {path}")
    if path.suffix.lower() in {".tsv", ".txt"}:
        delimiter = "\t"
    elif path.suffix.lower() == ".csv":
        delimiter = ","
    else:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            sample = handle.read(4096)
        delimiter = csv.Sniffer().sniff(sample, delimiters=",\t").delimiter
    # Keep IDs and accessions exactly as authored, including leading zeros.
    return pd.read_csv(path, sep=delimiter, dtype=str, encoding="utf-8-sig",
                       keep_default_na=False)


def require_columns(df: pd.DataFrame, columns: List[str], source: Path) -> None:
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise ValueError(f"Missing column(s) {missing} in {source}; available: {list(df.columns)}")


def is_missing(value: Any) -> bool:
    return value is None or (isinstance(value, float) and math.isnan(value)) or str(value).strip().lower() in MISSING


def clean(value: Any) -> str:
    return "" if is_missing(value) else str(value).strip()


def numeric(value: Any, context: str, *, allow_missing: bool = False) -> float:
    if is_missing(value):
        if allow_missing:
            return math.nan
        raise ValueError(f"Missing numeric value: {context}")
    try:
        answer = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"Non-numeric value {value!r}: {context}") from exc
    if not math.isfinite(answer):
        raise ValueError(f"Non-finite value {value!r}: {context}")
    return answer


def proportion(value: Any, context: str, *, allow_missing: bool = False) -> float:
    result = numeric(value, context, allow_missing=allow_missing)
    if not math.isnan(result) and not (0.0 <= result <= 1.0):
        raise ValueError(f"Expected 0–1 GC proportion; found {result} at {context}")
    return result


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="\\N", float_format="%.12g")


def write_json(obj: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # Python's default NaN JSON is non-standard; export null instead.
    def make_safe(value: Any) -> Any:
        if isinstance(value, dict):
            return {str(k): make_safe(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [make_safe(v) for v in value]
        if isinstance(value, float) and not math.isfinite(value):
            return None
        return value
    path.write_text(json.dumps(make_safe(obj), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def is_internal(value: str) -> bool:
    return value == "ancestral" or re.fullmatch(r"node_\d+", value) is not None


def derive_species_key_short(value: str) -> str:
    parts = [x for x in re.split(r"[_\s]+", clean(value)) if x]
    return "_".join(parts[:2]) if len(parts) >= 2 else ""


def coerce_node_name(value: Any) -> str:
    """Preserve terminal names; only normalize known numeric/internal aliases."""
    item = clean(value)
    if not item:
        return ""
    if item.lower() in {"ancestral", "ancestor", "root"}:
        return "ancestral"
    match = re.fullmatch(r"(?:node[ _-]*)?(\d+)", item, flags=re.I)
    if match:
        return f"node_{int(match.group(1))}"
    return item


def parse_node_relationships(path: Path) -> Dict[str, Dict[str, List[float]]]:
    """Maintain v5's long- and wide-format node relationship interpretation."""
    df = read_table(path)
    columns = {str(c).lower().strip(): c for c in df.columns}
    def pick(options: List[str]) -> Optional[str]:
        return next((columns[c] for c in options if c in columns), None)
    child_col = pick(["daughter", "daughter_node", "child", "child_node", "node",
                      "node_id", "descendant", "descendant_node", "son", "son_node"])
    parent_col = pick(["parent", "parent_node", "ancestor", "ancestor_node",
                       "mother", "parent_id"])
    relationships: Dict[str, Dict[str, List[float]]] = {}
    if child_col and parent_col and child_col != parent_col:
        for _, row in df.iterrows():
            child, parent = coerce_node_name(row[child_col]), coerce_node_name(row[parent_col])
            if child and parent:
                relationships.setdefault(child, {}).setdefault(parent, [])
    else:
        for raw_child in df.columns:
            child = coerce_node_name(raw_child)
            if not child:
                continue
            relationships.setdefault(child, {})
            for raw_parent in df[raw_child]:
                parent = coerce_node_name(raw_parent)
                if parent:
                    relationships[child].setdefault(parent, [])
    if not relationships:
        raise ValueError(f"No parent/child relationships could be read from {path}")
    return relationships


def normalize_odb12_id(filename: str, suffix: str) -> str:
    prefix = filename[:-len(suffix)]
    # Handle 401at8457__401at8457.phylip_nhPhyml*, _nhPhyml*, etc.
    if prefix.endswith(".phylip"):
        prefix = prefix[:-len(".phylip")]
    if re.fullmatch(r"([^/]+)__\1", prefix):
        prefix = prefix.split("__", 1)[0]
    if not prefix or prefix in {".", ".."}:
        raise ValueError(f"Cannot identify ODB12 ID in {filename}")
    return prefix


def index_nhphyml_files(directory: Path, tree_suffix: Optional[str], lk_suffix: Optional[str],
                         strict: bool, problems: List[str]) -> Dict[str, Tuple[Path, Path]]:
    if not directory.is_dir():
        raise NotADirectoryError(f"nhPhyML directory not found: {directory}")
    tree_suffixes = ([tree_suffix] if tree_suffix else [".phylip_nhPhymlGC.tree", "_nhPhymlGC.tree"])
    lk_suffixes = ([lk_suffix] if lk_suffix else [".phylip_nhPhyml.lk", "_nhPhyml.lk"])
    trees: Dict[str, Path] = {}
    lks: Dict[str, Path] = {}
    for path in sorted(directory.iterdir()):
        if not path.is_file():
            continue
        for suffixes, index, label in [(tree_suffixes, trees, "GC tree"), (lk_suffixes, lks, "likelihood")]:
            match = next((suffix for suffix in suffixes if path.name.endswith(suffix)), None)
            if match is None:
                continue
            odb = normalize_odb12_id(path.name, match)
            if odb in index:
                raise ValueError(f"Ambiguous duplicate {label} for odb12_id={odb}: "
                                 f"{index[odb]} and {path}")
            index[odb] = path
            break
    paired = sorted(trees.keys() & lks.keys())
    if not paired:
        preview = ', '.join(p.name for p in sorted(directory.iterdir())[:8])
        raise ValueError(f"No paired nhPhyML GC tree/likelihood files in {directory}. "
                         f"First entries: {preview}")
    for label, items in [("GC tree without likelihood", sorted(trees.keys() - lks.keys())),
                         ("likelihood without GC tree", sorted(lks.keys() - trees.keys()))]:
        if items:
            problems.append(f"{label}: {len(items)} orthologs, e.g. {items[:8]}")
    if strict and problems:
        raise ValueError("Incomplete nhPhyML output pairs: " + "; ".join(problems))
    return {orth: (trees[orth], lks[orth]) for orth in paired}


class GCDivergenceV6:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.first_internal_node = args.number_of_species + 1
        self.created_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
        self.details: List[Dict[str, str]] = []
        self.strict_issues: List[str] = []
        self.notes: List[str] = []
        self.node_gc3_total: Dict[str, List[float]] = defaultdict(list)
        self.node_anc_relationship: Dict[str, List[float]] = defaultdict(list)
        self.gc3_per_gene_rows: List[Dict[str, Any]] = []
        self.gc3_orthologs: Dict[str, Set[str]] = defaultdict(set)
        self.gc4_orthologs: Dict[str, Set[str]] = defaultdict(set)
        self.gc4_recorded: Dict[str, Set[str]] = defaultdict(set)
        self.gc4_exclusions: Dict[Tuple[str, str], str] = {}
        self.terminal_sets: Dict[str, Set[str]] = {}
        self.displacement_examples: List[Dict[str, Any]] = []
        self.dianc_dict: Dict[str, float] = {}
        self.dij_dict: Dict[str, Dict[str, float]] = {}
        self.node_relationship = parse_node_relationships(args.node_numbers)
        self.paired = index_nhphyml_files(args.nhphyml_dir, args.tree_suffix, args.lk_suffix,
                                          args.strict, self.notes)
        self.run_meta: Dict[str, str] = {}
        self.manifest = pd.DataFrame()
        self.gc4_summary = pd.DataFrame()
        self.gc4_detailed = pd.DataFrame()
        self.gc4_summary_source = pd.DataFrame()
        self.gc4_mean_disagreements = 0
        self.gc3_species: Set[str] = set()
        self.gc4_species: Set[str] = set()
        self.manifest_species: Set[str] = set()

    def record(self, category: str, species: str = "", ortholog: str = "",
               source: Any = "", reference: Any = "", details: str = "",
               *, strict: bool = False) -> None:
        row = dict(zip(DIAGNOSTIC_DETAIL_COLS,
                       [category, species, ortholog, clean(source), clean(reference), details]))
        self.details.append(row)
        if strict:
            self.strict_issues.append(f"{category}: {species} {ortholog} {details}".strip())

    def load_manifest(self) -> None:
        df = read_table(self.args.manifest)
        require_columns(df, ["species_key", "accession"], self.args.manifest)
        df = df[["species_key", "accession"]].copy()
        df["species_key_short_manifest"] = df["species_key"].map(derive_species_key_short)
        if (df["species_key_short_manifest"] == "").any():
            bad = df[df["species_key_short_manifest"] == ""].head(8)
            raise ValueError(f"Manifest has incomplete species_key values: {bad.to_dict('records')}")
        if df["accession"].map(is_missing).any():
            bad = df[df["accession"].map(is_missing)].head(8)
            raise ValueError(f"Manifest contains blank accessions: {bad.to_dict('records')}")
        df["accession"] = df["accession"].map(clean)
        if df["accession"].duplicated().any():
            bad = df[df["accession"].duplicated(keep=False)]
            raise ValueError(f"Manifest accession is not unique: {bad.head(8).to_dict('records')}")
        self.manifest = df
        self.manifest_species = set(df["species_key_short_manifest"])
        duplicates = df[df["species_key_short_manifest"].duplicated(keep=False)]
        for label, rows in duplicates.groupby("species_key_short_manifest"):
            self.record("ambiguous_shortened_manifest_species", label, source=",".join(rows["accession"]),
                        details="Multiple distinct manifest records share this shortened species label",
                        strict=True)
        print(f"[OK] Loaded {len(df)} manifest entries; {len(self.manifest_species)} short species labels")

    def load_gc4(self) -> None:
        detailed = read_table(self.args.gc4_detailed)
        source = read_table(self.args.gc4_species_summary)
        require_columns(detailed, ["nhphyml_run_pk", "run_name", "alignment_name", "odb12_id",
                                    "species_key_short", "gc4", "valid_gc4_bases", "status"], self.args.gc4_detailed)
        require_columns(source, ["nhphyml_run_pk", "run_name", "alignment_name", "species_key_short",
                                  "gc4_unweighted", "gc4_weighted"], self.args.gc4_species_summary)
        for table, path in [(detailed, self.args.gc4_detailed), (source, self.args.gc4_species_summary)]:
            if table.empty:
                raise ValueError(f"Empty GC4 input: {path}")
            for field in ["nhphyml_run_pk", "run_name", "alignment_name"]:
                values = set(table[field].map(clean))
                if len(values) != 1 or "" in values:
                    raise ValueError(f"Expected exactly one nonempty {field} in {path}: {values}")
            for field in ["species_key_short"]:
                if any(table[field].map(is_missing)):
                    raise ValueError(f"Missing species labels in {path}")
        for field in ["nhphyml_run_pk", "run_name", "alignment_name"]:
            left = clean(source[field].iloc[0])
            right = clean(detailed[field].iloc[0])
            if left != right:
                raise ValueError(f"GC4 metadata conflict for {field}: species summary={left!r}, detailed={right!r}")
            self.run_meta[field] = left
        if not re.fullmatch(r"\d+", self.run_meta["nhphyml_run_pk"]) or int(self.run_meta["nhphyml_run_pk"]) < 1:
            raise ValueError(f"nhphyml_run_pk must be a positive integer, got {self.run_meta['nhphyml_run_pk']!r}")
        if source["species_key_short"].duplicated().any():
            raise ValueError("Duplicate species_key_short values in GC4 species_summary.tsv")
        if detailed.duplicated(subset=["odb12_id", "species_key_short"]).any():
            dup = detailed[detailed.duplicated(subset=["odb12_id", "species_key_short"], keep=False)]
            raise ValueError(f"Duplicate GC4 species/ortholog records: {dup[['odb12_id','species_key_short']].head(8).to_dict('records')}")
        self.gc4_species = set(source["species_key_short"])
        gc4_rows = defaultdict(list)
        for row in detailed.itertuples(index=False):
            species = clean(row.species_key_short)
            odb = clean(row.odb12_id)
            if not species or not odb:
                raise ValueError("GC4 detailed.tsv contains blank species_key_short/odb12_id")
            self.gc4_recorded[species].add(odb)
            raw_status = clean(row.status)
            status = raw_status.lower().replace("-", "_").replace(" ", "_")
            count = numeric(row.valid_gc4_bases, f"GC4 valid_gc4_bases {species}/{odb}", allow_missing=True)
            gc = proportion(row.gc4, f"GC4 detailed {species}/{odb}", allow_missing=True)
            if status not in GC4_GOOD_STATUSES and status not in GC4_REJECT_STATUSES:
                self.record("unrecognized_gc4_status", species, odb, source=raw_status,
                            details="Using numeric value and valid_gc4_bases to determine eligibility")
            reason = ""
            if status in GC4_REJECT_STATUSES:
                reason = f"GC4 status={raw_status}"
            elif math.isnan(count) or count <= 0:
                reason = "No valid fourfold-degenerate bases"
            elif math.isnan(gc):
                reason = "Missing/invalid GC4 despite usable bases"
            if reason:
                self.gc4_exclusions[(species, odb)] = reason
                self.record("gc4_no_usable_value", species, odb, source=row.gc4,
                            reference=row.valid_gc4_bases, details=reason)
                continue
            self.gc4_orthologs[species].add(odb)
            gc4_rows[species].append(gc)
        gc4_summaries: List[Dict[str, Any]] = []
        for record in source.to_dict("records"):
            species = clean(record["species_key_short"])
            unweighted = proportion(record["gc4_unweighted"], f"GC4 species unweighted {species}", allow_missing=True)
            weighted = proportion(record["gc4_weighted"], f"GC4 species weighted {species}", allow_missing=True)
            values = gc4_rows.get(species, [])
            if values:
                checked_mean = statistics.mean(values)
                if math.isnan(unweighted) or not math.isclose(checked_mean, unweighted, abs_tol=1e-8, rel_tol=1e-8):
                    self.gc4_mean_disagreements += 1
                    self.record("gc4_mean_disagreement", species,
                                source=checked_mean, reference=unweighted,
                                details=f"Detailed valid ortholog count={len(values)}; inspect filtering/status", strict=True)
            elif not math.isnan(unweighted):
                self.gc4_mean_disagreements += 1
                self.record("gc4_mean_without_detail", species, source=unweighted,
                            details="Summary reports a mean but detailed file has no valid values", strict=True)
            if math.isnan(unweighted) and values:
                self.record("gc4_summary_missing_mean", species, source=statistics.mean(values),
                            details="Detailed data have GC4 values, but species summary mean is missing", strict=True)
            gc4_summaries.append({
                "species_key_short": species,
                "mean_gc4_unweighted": unweighted,
                "sd_gc4_unweighted": statistics.stdev(values) if len(values) > 1 else math.nan,
                "mean_gc4_weighted": weighted,
                "n_gc4_values": len(values),
            })
        # Detect detail-only species without dropping their records silently.
        for species in sorted(set(self.gc4_recorded) - self.gc4_species):
            self.record("gc4_detail_species_missing_summary", species,
                        details="Found in GC4 detailed input but absent from species_summary", strict=True)
        for species in sorted(self.gc4_species - set(self.gc4_recorded)):
            self.record("gc4_summary_species_missing_detail", species,
                        details="Found in species_summary but absent from GC4 detailed input", strict=True)
        self.gc4_detailed = detailed
        self.gc4_summary_source = source
        self.gc4_summary = pd.DataFrame(gc4_summaries, columns=GC4_SUMMARY_COLS)
        print(f"[OK] GC4 run {self.run_meta['nhphyml_run_pk']}: {len(source)} species, "
              f"{len(detailed)} detailed records")

    def parse_gc_tree(self, path: Path, ortholog: str) -> Dict[str, float]:
        """Parse annotated GC tree into percent-scale values.

        Leaves: `(Genus_species 42.01:0.01,...)`.
        Internal: `)47.8:0.04` (or `)47.8;` at the root).
        Internal numbering follows the left-to-right closing-parenthesis order,
        starting at number_of_species + 1, as in v5's documented procedure.
        Unlike v5's ambiguous dot regex, the two cases are parsed separately.
        """
        tree = path.read_text(encoding="utf-8", errors="replace").strip()
        if not tree or tree.count("(") != tree.count(")"):
            raise ValueError(f"Malformed/unbalanced GC tree: {path}")
        gc: Dict[str, float] = {}
        leaves = list(TERMINAL_GC.finditer(tree))
        for match in leaves:
            label, raw_value = match.group(1), float(match.group(2))
            if is_internal(label) or label in gc:
                raise ValueError(f"Duplicate or reserved terminal species label {label!r}: {path}")
            gc[label] = raw_value
        if len(leaves) != self.args.number_of_species:
            raise ValueError(f"Expected {self.args.number_of_species} labeled terminals in {path}, "
                             f"parsed {len(leaves)}. Check nhPhyML GC tree formatting")
        n_annotated_internal = 0
        for i, closing in enumerate(re.finditer(r"\)", tree)):
            identifier = f"node_{self.first_internal_node + i}"
            m = INTERNAL_GC.match(tree, closing.end())
            if m:
                gc[identifier] = float(m.group(1))
                n_annotated_internal += 1
            else:
                self.record("unannotated_internal_node", "", ortholog, source=identifier,
                            details=f"No GC3 annotation after ')' in {path.name}", strict=True)
        raw_values = list(gc.values())
        if self.args.tree_gc_scale == "percent":
            if any(not 0 <= x <= 100 for x in raw_values):
                raise ValueError(f"Tree GC3 outside percentage range 0..100: {path}")
            if all(x <= 1 for x in raw_values):
                raise ValueError(f"All annotated GC3 values in {path.name} lie in 0..1; "
                                 "input may be proportions. Specify --tree-gc-scale proportion if confirmed")
        else:
            if any(not 0 <= x <= 1 for x in raw_values):
                raise ValueError(f"GC tree contains values >1, inconsistent with --tree-gc-scale proportion: {path}")
            gc = {key: value * 100 for key, value in gc.items()}
        self.terminal_sets[ortholog] = set(match.group(1) for match in leaves)
        if n_annotated_internal == 0:
            raise ValueError(f"No annotated internal GC3 nodes parsed from {path}")
        return gc

    def parse_ancestral_gc(self, path: Path, ortholog: str) -> float:
        """Preserve v5 last-number-on-Ancestral-line extraction, in percent units."""
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if "Ancestral" not in line:
                continue
            matches = FLOAT_NUMBERS.findall(line)
            if matches:
                raw = float(matches[-1])
                if raw < 0 or raw > 100:
                    raise ValueError(f"Ancestral GC outside 0..100 at {path}: {raw}")
                result = raw * 100 if raw <= 1 else raw
                if result > 100:
                    raise ValueError(f"Ancestral GC not within 0..100 at {path}")
                return result
        raise ValueError(f"Could not identify ancestral GC3 in likelihood file for {ortholog}: {path}")

    def run_gc3(self) -> None:
        first_species: Optional[Set[str]] = None
        for index, (ortholog, (tree, lk)) in enumerate(self.paired.items(), 1):
            if index == 1 or index % 100 == 0 or index == len(self.paired):
                print(f"[INFO] nhPhyML ortholog {index}/{len(self.paired)}: {ortholog}")
            node_gc = self.parse_gc_tree(tree, ortholog)
            if first_species is None:
                first_species = self.terminal_sets[ortholog]
            elif first_species != self.terminal_sets[ortholog]:
                for species in sorted(first_species - self.terminal_sets[ortholog]):
                    self.record("nhphyml_species_missing_from_ortholog", species, ortholog,
                                details="Species present in first ortholog but absent here", strict=True)
                for species in sorted(self.terminal_sets[ortholog] - first_species):
                    self.record("nhphyml_species_unexpected_in_ortholog", species, ortholog,
                                details="Not present in first ortholog", strict=True)
            node_gc["ancestral"] = self.parse_ancestral_gc(lk, ortholog)
            for node, percent in node_gc.items():
                frac = percent / 100.0
                if not (0 <= frac <= 1):
                    raise ValueError(f"GC3 proportion invalid for {node}, {ortholog}: {frac}")
                self.node_gc3_total[node].append(frac)
                self.gc3_per_gene_rows.append({
                    "odb12_id": ortholog, "node_or_species_key_short": node, "gc3_nhphyml": frac,
                })
                if not is_internal(node):
                    self.gc3_orthologs[node].add(ortholog)
            missing_pairs = []
            for child, parent_dict in self.node_relationship.items():
                if child not in node_gc:
                    missing_pairs.append((child, "CHILD_MISSING"))
                    continue
                for parent, values in parent_dict.items():
                    if parent not in node_gc:
                        missing_pairs.append((child, parent))
                        continue
                    # These remain percentage-point differences, exactly as in v5.
                    difference = node_gc[child] - node_gc[parent]
                    sq = difference ** 2
                    values.append(sq)
                    if len(self.displacement_examples) < 200:
                        self.displacement_examples.append({
                            "odb12_id": ortholog, "child": child, "parent": parent,
                            "child_gc3_percent": node_gc[child],
                            "parent_gc3_percent": node_gc[parent],
                            "difference_percent_points": difference,
                            "squared_displacement": sq,
                        })
            if missing_pairs:
                example = ", ".join(f"{a}->{b}" for a, b in missing_pairs[:8])
                self.record("nhphyml_node_relationship_missing", "", ortholog,
                            details=f"{len(missing_pairs)} missing pairs: {example}", strict=True)
            for species, percent in node_gc.items():
                if is_internal(species):
                    continue
                difference = percent - node_gc["ancestral"]
                self.node_anc_relationship[species].append(difference ** 2)
        self.gc3_species = set(self.gc3_orthologs)
        for child, parents in self.node_relationship.items():
            self.dij_dict[child] = {parent: math.sqrt(sum(values)) if values else math.nan
                                    for parent, values in parents.items()}
        self.dianc_dict = {species: math.sqrt(sum(values)) if values else math.nan
                           for species, values in self.node_anc_relationship.items()}
        print(f"[OK] nhPhyML: {len(self.paired)} orthologs, {len(self.gc3_species)} terminal labels; "
              f"first internal node node_{self.first_internal_node}")
        if len(self.gc3_species) != self.args.number_of_species:
            self.record("nhphyml_species_count_mismatch", source=len(self.gc3_species),
                        reference=self.args.number_of_species, strict=True)

    def validate_species(self) -> pd.DataFrame:
        for species in sorted(self.gc3_species - self.manifest_species):
            self.record("nhphyml_species_missing_manifest", species,
                        details="Original nhPhyML label does not exactly equal a manifest-derived label", strict=True)
        for species in sorted(self.manifest_species - self.gc3_species):
            self.record("manifest_species_not_in_subset", species,
                        details="May be expected when manifest includes other subset taxa")
        for species in sorted(self.gc3_species - self.gc4_species):
            self.record("nhphyml_species_missing_gc4_summary", species, strict=True)
        for species in sorted(self.gc4_species - self.gc3_species):
            self.record("gc4_species_missing_nhphyml", species, strict=True)
        for species in sorted(self.gc3_species - set(self.gc4_recorded)):
            self.record("nhphyml_species_missing_gc4_detailed", species, strict=True)
        for species in sorted(self.gc4_species - self.manifest_species):
            self.record("gc4_species_missing_manifest", species, strict=True)
        # Report a validation result for every terminal, not only mismatches.
        match_groups = self.manifest.groupby("species_key_short_manifest", sort=False)
        manifest_ids = {label: block["accession"].tolist() for label, block in match_groups}
        for species in sorted(self.gc3_species):
            accessions = manifest_ids.get(species, [])
            self.record(
                "species_label_validation", species,
                source=species, reference=species if accessions else "",
                details=(f"manifest={'match' if len(accessions) == 1 else 'ambiguous' if accessions else 'missing'}; "
                         f"gc4_summary={'match' if species in self.gc4_species else 'missing'}; "
                         f"gc4_detailed={'match' if species in self.gc4_recorded else 'missing'}; "
                         f"manifest_accession={','.join(accessions) if accessions else 'missing'}; "
                         "accession origin=manifest (not independently confirmed from sequence data)")
            )
        # A missing manifest match means an accession cannot be assigned
        # without guessing; fail even without --strict, preserving the labels.
        missing_manifest = sorted(self.gc3_species - self.manifest_species)
        if missing_manifest:
            case_matches = {
                label: sorted(x for x in self.manifest_species if x.casefold() == label.casefold())
                for label in missing_manifest
            }
            raise ValueError(
                "Cannot assign manifest accessions to nhPhyML labels (labels are preserved, "
                "not silently normalized): " + repr(missing_manifest[:20]) +
                f". Case-insensitive near-matches: {case_matches}"
            )
        # Resolve a unique manifest accession; never relabel a nhPhyML taxon.
        grouped = {name: block for name, block in self.manifest.groupby("species_key_short_manifest", sort=False)}
        selected: List[Dict[str, str]] = []
        for species in sorted(self.gc3_species):
            if species not in grouped:
                selected.append({"species_key_short": species, "accession": ""})
                continue
            block = grouped[species]
            if len(block) != 1:
                raise ValueError(f"Ambiguous accession for {species}: "
                                 f"{block[['species_key', 'accession']].to_dict('records')}")
            selected.append({"species_key_short": species, "accession": block.iloc[0]["accession"]})
        return pd.DataFrame(selected, columns=["species_key_short", "accession"])

    def validate_alignment_headers(self) -> None:
        if self.args.alignment_dir is None:
            self.notes.append("Direct alignment-header validation: not performed (--alignment-dir omitted)")
            return
        folder = self.args.alignment_dir
        if not folder.is_dir():
            raise NotADirectoryError(f"--alignment-dir does not exist: {folder}")
        extensions = (".fasta", ".fa", ".fas", ".fna", ".phy", ".phylip", ".aln")
        paths = [p for p in folder.iterdir() if p.is_file() and p.name.lower().endswith(extensions)]
        index: Dict[str, Path] = {}
        for path in sorted(paths):
            odb = path.name
            for ext in extensions:
                if odb.lower().endswith(ext):
                    odb = odb[:-len(ext)]
                    break
            if odb in index:
                raise ValueError(f"Duplicate alignments for {odb}: {index[odb]} and {path}")
            index[odb] = path
        for ortholog in self.paired:
            path = index.get(ortholog)
            if path is None:
                self.record("alignment_file_missing", "", ortholog,
                            details="Alignment not found for nhPhyML ortholog", strict=True)
                continue
            if path.suffix.lower() in {".phy", ".phylip"}:
                labels = self.read_phylip_labels(path)
            else:
                labels = self.read_fasta_labels(path)
            found = set(labels)
            duplicates = sorted([species for species, n in Counter(labels).items() if n > 1])
            if duplicates:
                self.record("alignment_duplicate_headers", "", ortholog,
                            source=",".join(duplicates), strict=True)
            for species in sorted(self.terminal_sets[ortholog] - found):
                self.record("nhphyml_species_absent_alignment", species, ortholog, strict=True)
            for species in sorted(found - self.terminal_sets[ortholog]):
                self.record("alignment_species_absent_nhphyml", species, ortholog, strict=True)
        print(f"[OK] Alignment validation examined {len(self.paired)} orthologs")

    @staticmethod
    def read_fasta_labels(path: Path) -> List[str]:
        labels = []
        with path.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(">"):
                    label = line[1:].strip().split()[0] if line[1:].strip() else ""
                    if label:
                        labels.append(label)
        if not labels:
            raise ValueError(f"No FASTA headers found: {path}")
        return labels

    @staticmethod
    def read_phylip_labels(path: Path) -> List[str]:
        with path.open(encoding="utf-8", errors="replace") as handle:
            header = next((line.strip() for line in handle if line.strip()), "")
            parts = header.split()
            if len(parts) < 2 or not all(x.isdigit() for x in parts[:2]):
                raise ValueError(f"Invalid PHYLIP heading in {path}: {header}")
            count = int(parts[0])
            labels: List[str] = []
            for line in handle:
                if not line.strip():
                    continue
                fields = line.split()
                if len(fields) >= 2 and re.fullmatch(r"[A-Za-z?*!.\-]+", fields[1]):
                    labels.append(fields[0])
                    if len(labels) == count:
                        break
            if len(labels) != count:
                raise ValueError(f"Could not read {count} PHYLIP labels from {path}; "
                                 "interleaved/strict PHYLIP formatting may require a dedicated parser")
            return labels

    def validate_ortholog_coverage(self) -> None:
        for species in sorted(self.gc3_species | self.gc4_species | set(self.gc4_recorded)):
            gc3 = self.gc3_orthologs.get(species, set())
            gc4 = self.gc4_orthologs.get(species, set())
            recorded = self.gc4_recorded.get(species, set())
            for ortholog in sorted(gc3 - gc4):
                if ortholog not in recorded:
                    kind = "gc3_ortholog_missing_gc4_record"
                    reason = "No record in GC4 detailed input"
                    is_unexpected = True
                else:
                    kind = "gc3_ortholog_without_usable_gc4"
                    reason = self.gc4_exclusions.get((species, ortholog), "GC4 record excluded")
                    is_unexpected = not (
                        "No valid fourfold-degenerate bases" in reason
                        or any(f"status={status}" in reason for status in GC4_EXPECTED_NO_SITE_STATUSES)
                    )
                self.record(kind, species, ortholog, details=reason, strict=is_unexpected)
            for ortholog in sorted(gc4 - gc3):
                self.record("gc4_ortholog_without_gc3", species, ortholog,
                            details="GC4 valid but nhPhyML GC3 unavailable", strict=True)
            for ortholog in sorted(recorded - gc3):
                if ortholog not in gc4:
                    self.record("gc4_record_without_gc3", species, ortholog,
                                details=self.gc4_exclusions.get((species, ortholog), "GC4 excluded"), strict=True)

    def export_gc3(self) -> pd.DataFrame:
        detail = pd.DataFrame(self.gc3_per_gene_rows, columns=GC3_DETAIL_COLS)
        if detail.duplicated(["odb12_id", "node_or_species_key_short"]).any():
            raise ValueError("Duplicate nhPhyML ortholog/node records detected")
        write_tsv(detail, self.args.outdir / "gc3_per_gene.tsv")
        rows = []
        for node, values in sorted(self.node_gc3_total.items()):
            rows.append({"node_or_species_key_short": node,
                         "mean_gc3_nhphyml_unweighted": statistics.mean(values),
                         "sd_gc3_nhphyml_unweighted": statistics.stdev(values) if len(values) > 1 else math.nan,
                         "n_orthologs": len(values)})
        summary = pd.DataFrame(rows, columns=GC3_SUMMARY_COLS)
        write_tsv(summary, self.args.outdir / "gc3_summary.tsv")
        return summary

    def export_divergence(self) -> None:
        write_json(self.dij_dict, self.args.outdir / "dij_results.json")
        dij_rows = [dict(node_or_species_key_short=child, parent=parent, dij=value)
                    for child, parent_vals in self.dij_dict.items() for parent, value in parent_vals.items()]
        write_tsv(pd.DataFrame(dij_rows, columns=["node_or_species_key_short", "parent", "dij"]),
                  self.args.outdir / "dij_results.tsv")
        write_json(self.dianc_dict, self.args.outdir / "dianc_results.json")
        write_tsv(pd.DataFrame([dict(species_key_short=s, dianc=d) for s, d in sorted(self.dianc_dict.items())],
                               columns=["species_key_short", "dianc"]),
                  self.args.outdir / "dianc_results.tsv")

    def export_displacement_examples(self) -> None:
        path = self.args.outdir / "displacement_examples.txt"
        header = ["Example nhPhyML GC3 displacement calculations", "=" * 48, "",
                  f"number_of_species = {self.args.number_of_species}",
                  f"first_internal_node = node_{self.first_internal_node}", "",
                  "GC3 in gc3_per_gene.tsv is a PROPORTION (0..1).",
                  "The unchanged v5 d_ij and d_ianc use GC3 in PERCENTAGE units (0..100).",
                  "d_ij = sqrt(sum_k (GC3_child,k - GC3_parent,k)^2)", ""]
        examples = self.displacement_examples
        df = pd.DataFrame(examples)
        if not df.empty:
            internal = df[df["child"].map(lambda label: str(label).startswith("node_"))].copy()
            selected = []
            if not internal.empty:
                internal["num"] = internal["child"].str.extract(r"node_(\d+)").astype(int)
                selected += internal[internal["num"] == internal["num"].max()].head(3).to_dict("records")
            anolis = df[df["child"].str.lower().str.startswith("anolis_")]
            selected += anolis.head(3).to_dict("records")
            if not selected:
                selected = df.head(6).to_dict("records")
            for i, row in enumerate(selected, start=1):
                child, parent = row["child"], row["parent"]
                a = row["child_gc3_percent"]
                b = row["parent_gc3_percent"]
                header += [f"Example {i}", "-" * 20,
                           f"odb12_id: {row['odb12_id']}", f"child: {child}", f"parent: {parent}",
                           f"GC3 child (proportion): {a / 100:.6f}",
                           f"GC3 parent (proportion): {b / 100:.6f}",
                           f"GC3 child (percent): {a:.6f}", f"GC3 parent (percent): {b:.6f}",
                           f"difference (percentage points): {a - b:.6f}",
                           f"squared displacement: {row['squared_displacement']:.6f}",
                           f"Final d_ij for {child} -> {parent}: "
                           f"{self.dij_dict.get(child, {}).get(parent, math.nan):.6f}", ""]
        else:
            header.append("No displacement examples were recorded.")
        path.write_text("\n".join(header) + "\n", encoding="utf-8")

    def export_master(self, gc3_summary: pd.DataFrame, manifest_matches: pd.DataFrame) -> pd.DataFrame:
        species = gc3_summary.loc[
            ~gc3_summary["node_or_species_key_short"].map(is_internal)
        ].rename(columns={"node_or_species_key_short": "species_key_short"}).copy()
        master = species.merge(manifest_matches, on="species_key_short", how="left", validate="one_to_one")
        master = master.merge(self.gc4_summary[["species_key_short", "mean_gc4_unweighted", "sd_gc4_unweighted"]],
                              on="species_key_short", how="left", validate="one_to_one")
        master["dianc"] = master["species_key_short"].map(self.dianc_dict)
        if self.args.include_all_manifest_species:
            absent = sorted(self.manifest_species - set(master["species_key_short"]))
            manifest_all = {key: rows for key, rows in self.manifest.groupby("species_key_short_manifest")}
            extra = []
            for key in absent:
                rows = manifest_all[key]
                if len(rows) != 1:
                    raise ValueError(f"Ambiguous manifest entry for --include-all-manifest-species: {key}")
                extra.append({"species_key_short": key, "accession": rows.iloc[0]["accession"]})
            if extra:
                master = pd.concat([master, pd.DataFrame(extra)], ignore_index=True)
        master["nhphyml_run_pk"] = int(self.run_meta["nhphyml_run_pk"])
        master["run_name"] = self.run_meta["run_name"]
        master["alignment_name"] = self.run_meta["alignment_name"]
        master["created_at"] = self.created_at
        master = master.reindex(columns=MASTER_COLS)
        master = master.sort_values("species_key_short").reset_index(drop=True)
        write_tsv(master, self.args.outdir / "master_output.tsv")
        return master

    def export_diagnostics(self, master: pd.DataFrame) -> None:
        gc3_only = sum(len(v - self.gc4_orthologs.get(k, set())) for k, v in self.gc3_orthologs.items())
        gc4_only = sum(len(v - self.gc3_orthologs.get(k, set())) for k, v in self.gc4_orthologs.items())
        shared = sum(len(v & self.gc4_orthologs.get(k, set())) for k, v in self.gc3_orthologs.items())
        diagnostics: List[Tuple[str, Any]] = [
            ("nhphyml_run_pk", self.run_meta["nhphyml_run_pk"]),
            ("run_name", self.run_meta["run_name"]),
            ("alignment_name", self.run_meta["alignment_name"]),
            ("number_of_species_argument", self.args.number_of_species),
            ("first_internal_node_used", f"node_{self.first_internal_node}"),
            ("paired_nhphyml_orthologs", len(self.paired)),
            ("species_with_gc3", len(self.gc3_species)),
            ("species_with_gc4_summary", len(self.gc4_species)),
            ("manifest_short_species", len(self.manifest_species)),
            ("species_matched_manifest", len(self.gc3_species & self.manifest_species)),
            ("species_matched_gc4", len(self.gc3_species & self.gc4_species)),
            ("gc3_ortholog_species_records", sum(map(len, self.gc3_orthologs.values()))),
            ("gc4_valid_ortholog_species_records", sum(map(len, self.gc4_orthologs.values()))),
            ("gc3_gc4_shared_ortholog_species_records", shared),
            ("gc3_only_ortholog_species_records", gc3_only),
            ("gc4_only_ortholog_species_records", gc4_only),
            ("gc4_mean_disagreements", self.gc4_mean_disagreements),
            ("rows_in_master_output", len(master)),
            ("rows_with_accession", master["accession"].map(lambda a: not is_missing(a)).sum()),
            ("rows_with_mean_gc3", master["mean_gc3_nhphyml_unweighted"].notna().sum()),
            ("rows_with_mean_gc4", master["mean_gc4_unweighted"].notna().sum()),
            ("rows_with_sd_gc3", master["sd_gc3_nhphyml_unweighted"].notna().sum()),
            ("rows_with_sd_gc4", master["sd_gc4_unweighted"].notna().sum()),
            ("rows_with_dianc", master["dianc"].notna().sum()),
            ("direct_alignment_validation", "performed" if self.args.alignment_dir else "not_performed"),
            ("detailed_validation_findings", len(self.details)),
            ("strict_validation_findings", len(self.strict_issues)),
        ]
        for note in self.notes:
            diagnostics.append(("note", note))
        write_tsv(pd.DataFrame(diagnostics, columns=["diagnostic", "value"]),
                  self.args.outdir / "master_merge_diagnostics.tsv")
        write_tsv(pd.DataFrame(self.details, columns=DIAGNOSTIC_DETAIL_COLS),
                  self.args.outdir / "validation_details.tsv")

    def run(self) -> None:
        self.args.outdir.mkdir(parents=True, exist_ok=True)
        legacy_output = self.args.outdir / "gc_summary_output.tsv"
        if legacy_output.exists():
            self.notes.append("Legacy gc_summary_output.tsv already exists from an earlier run; "
                              "v6 does not update or delete that file")
        print(f"[OK] first_internal_node=node_{self.first_internal_node}; "
              f"paired nhPhyML outputs={len(self.paired)}")
        self.load_manifest()
        self.load_gc4()
        self.run_gc3()
        matches = self.validate_species()
        self.validate_alignment_headers()
        self.validate_ortholog_coverage()
        gc3_summary = self.export_gc3()
        write_tsv(self.gc4_summary, self.args.outdir / "gc4_summary.tsv")
        self.export_divergence()
        self.export_displacement_examples()
        master = self.export_master(gc3_summary, matches)
        self.export_diagnostics(master)
        for note in self.notes[:15]:
            print(f"[WARN] {note}")
        counts = Counter(row["category"] for row in self.details)
        if counts.get("species_label_validation", 0):
            print(f"[OK] Documented species-label checks for {counts['species_label_validation']} terminal taxa")
        for category, count in counts.most_common(15):
            if category != "species_label_validation":
                print(f"[WARN] {category}: {count}")
        print(f"[DONE] Outputs written to {self.args.outdir}; master rows: {len(master)}")
        if self.args.strict and self.strict_issues:
            sample = "\n  - ".join(self.strict_issues[:12])
            raise ValueError(f"Strict validation failed ({len(self.strict_issues)} finding(s)); "
                             f"see validation_details.tsv and master_merge_diagnostics.tsv:\n  - {sample}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("genomes_dir", type=Path, help="Root directory of genomes project")
    parser.add_argument("--number-of-species", required=True, type=int, help="Number of terminal taxa")
    parser.add_argument("--nhphyml-dir", required=True, type=Path,
                        help="Directory of paired GC.tree and .lk nhPhyML outputs")
    parser.add_argument("--node-numbers", required=True, type=Path,
                        help="node_numbers_rebuilt.csv (long/wide parent-child relationships)")
    parser.add_argument("--manifest", required=True, type=Path,
                        help="Project manifest with species_key and accession columns")
    parser.add_argument("--gc4-detailed", required=True, type=Path,
                        help="detailed.tsv produced by 14_gc4_from_alignment_v3.py")
    parser.add_argument("--gc4-species-summary", required=True, type=Path,
                        help="species_summary.tsv produced by 14_gc4_from_alignment_v3.py")
    parser.add_argument("--alignment-dir", type=Path,
                        help="Optional source FASTA/PHYLIP alignments for label validation only")
    parser.add_argument("--outdir", type=Path,
                        help="Default: GENOMES_DIR/records/gc_metrics_outputs")
    parser.add_argument("--include-all-manifest-species", action="store_true",
                        help="Add blank-result rows for manifest taxa not in nhPhyML subset")
    parser.add_argument("--tree-suffix", default=None,
                        help="Override auto-recognized GC tree filename suffix")
    parser.add_argument("--lk-suffix", default=None,
                        help="Override auto-recognized likelihood filename suffix")
    parser.add_argument("--tree-gc-scale", choices=["percent", "proportion"], default="percent",
                        help="Scale of GC values annotated in nhPhyML GC trees (default: percent)")
    parser.add_argument("--strict", action="store_true", help="Fail on unexpected coverage/name mismatches")
    # Deprecated v5 hidden argument: unused natural-history metadata inputs.
    parser.add_argument("--input-tsv-dir", type=Path, default=None, help=argparse.SUPPRESS)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.genomes_dir = args.genomes_dir.expanduser().resolve()
    if not args.genomes_dir.is_dir():
        raise NotADirectoryError(f"genomes_dir not found: {args.genomes_dir}")
    if args.number_of_species < 2:
        raise ValueError("--number-of-species must be >= 2")
    for name in ["nhphyml_dir", "node_numbers", "manifest", "gc4_detailed", "gc4_species_summary", "alignment_dir"]:
        value = getattr(args, name)
        if value is not None:
            setattr(args, name, value.expanduser().resolve())
    args.outdir = (args.outdir or args.genomes_dir / "records" / "gc_metrics_outputs").expanduser().resolve()
    GCDivergenceV6(args).run()


if __name__ == "__main__":
    main()
