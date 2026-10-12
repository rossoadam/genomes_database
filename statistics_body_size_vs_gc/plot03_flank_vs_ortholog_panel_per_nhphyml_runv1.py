#!/usr/bin/env python3
"""Plot genome-mean GC3/GC4 versus paired genome-mean flank GC in four runs.

Dependencies: python -m pip install pandas numpy scipy matplotlib

Example:
  python plot_flank_vs_ortholog_panel_per_nhphyml_runv1.py \
    --flank-tsv /Users/rossoaa/projects/genomes/records/sql_tsvs/flank_compleasm.tsv \
    --nhphyml-tsv /Users/rossoaa/projects/genomes/records/sql_tsvs/nhphyml_ortholog_results.tsv

Each point represents ONE genome. Join on genome-specific ortholog_pk, then
within each (run, genome) average the two coordinates over the SAME orthologs.
GC3 and GC4 use separate paired subsets if their missing values differ.
Exclude flank GC that is zero, NULL, NaN or infinite. Missing coding GC is also
excluded, but coding GC=0 is retained. GC values retain their input units.

Run expectations are genome counts and maximum ortholog counts PER GENOME,
not counts of distinct genome-specific ortholog_pk values across a whole run.
These two tables cannot independently verify BUSCO/family identities; doing
that would require an ortholog-to-family mapping table.

Multiple flank sets require --flank-set-pk. Duplicate rows are never averaged
away. By default validation requires exactly the expected genome count before
filtering AND in both plots. --allow-incomplete-genomes permits fewer genomes,
explicitly reporting the actual counts; it never permits extra genomes or
excess orthologs. Missing orthologs below the per-genome limit are permitted.

Outputs: two 2x2 figures (PNG and PDF), genome_means.tsv,
run_validation.tsv, correlation_statistics.tsv, validation_issues.tsv.
--validate-only writes the audit tables and skips figures.
Pearson/regression statistics describe correlations ACROSS GENOME MEANS;
they do not account for phylogenetic non-independence.
"""

from pathlib import Path
import argparse
import sys

try:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LinearLocator, ScalarFormatter
    from scipy.stats import linregress
except ImportError as exc:
    raise SystemExit(f"Missing dependency: {exc}. Install pandas numpy scipy matplotlib.")


# ---------------------- USER-EDITABLE SETTINGS ----------------------
axes_width = 3
figure_size = (16.5, 16.5)  # cm, for EACH complete four-panel figure
regression_line_width = 5
point_diameter = 10        # points; plot(markersize=...) gives a diameter
point_color = "black"
axes_color = "black"
regression_line_color = "black"
regression_line_style = "--"
axes_font_size = 10
key_font_size = 10         # retained as requested; no key/legend is drawn
descriptive_statistics_font_size = 8
number_of_ticks = 4        # exactly this many major ticks on BOTH axes
output_dpi = 300
panel_labels = ("A", "B", "C", "D")

# nhphyml_run_pk: (expected genomes, maximum retained orthologs per genome)
RUN_EXPECTATIONS = {
    1: (80, 645),
    2: (70, 1189),
    3: (60, 1794),
    4: (50, 2609),
}
DEFAULT_OUTPUT_DIR = Path(
    "/Users/rossoaa/projects/genomes/records/statistics_body_size_vs_gc/"
    "correlation_flank_vs_ortholog_per_nhphyml_run"
)
# -------------------------------------------------------------------

NULL_VALUES = ["", "\\N", "NULL", "null", "None", "none", "NA", "N/A", "NaN", "nan"]
METRICS = ("gc3_nhphyml", "gc4")


def load_tsv(path, required):
    """Read plain or compressed TSVs, preserving large integer keys as text."""
    df = pd.read_csv(path, sep="\t", dtype="string", keep_default_na=False,
                     na_values=NULL_VALUES, encoding="utf-8-sig")
    df.columns = df.columns.str.strip()
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise ValueError(f"{path}: missing required columns: {', '.join(missing)}")
    for col in required:
        df[col] = df[col].str.strip().replace(NULL_VALUES, pd.NA)
    return df


def integer_key(df, column, source):
    s = df[column]
    bad = s.isna() | ~s.str.fullmatch(r"[0-9]+", na=False)
    if bad.any():
        raise ValueError(f"{source}: invalid/missing {column} in {int(bad.sum())} rows")
    # Canonicalize without floating-point conversions or rounding large keys.
    df[column] = s.map(lambda value: str(int(value)))
    if df[column].eq("0").any():
        raise ValueError(f"{source}: {column} must be a positive integer")


def numeric_gc(df, column, source):
    values = pd.to_numeric(df[column], errors="coerce").astype(float)
    malformed = df[column].notna() & values.isna()
    if malformed.any():
        examples = df.loc[malformed, column].head(3).tolist()
        raise ValueError(f"{source}: nonnumeric {column}: {examples}")
    finite = np.isfinite(values)
    if (finite & ((values < 0) | (values > 100))).any():
        raise ValueError(f"{source}: {column} outside [0, 100]; check GC units")
    df[column] = values.where(finite, np.nan)


def reject_duplicates(df, columns, label):
    mask = df.duplicated(columns, keep=False)
    if mask.any():
        sample = df.loc[mask, columns].head(5).to_dict("records")
        raise ValueError(f"Duplicate {label}: {sample}. Resolve duplicates before plotting.")


def prepare_inputs(args):
    n = load_tsv(args.nhphyml_tsv,
                 ["nhphyml_run_pk", "genome_pk", "ortholog_pk", *METRICS])
    for col in ["nhphyml_run_pk", "genome_pk", "ortholog_pk"]:
        integer_key(n, col, args.nhphyml_tsv)
    n["nhphyml_run_pk"] = n["nhphyml_run_pk"].map(int)
    ignored = sorted(set(n["nhphyml_run_pk"]) - set(RUN_EXPECTATIONS))
    if ignored:
        print(f"[INFO] Ignoring runs outside RUN_EXPECTATIONS: {ignored}")
    n = n[n["nhphyml_run_pk"].isin(RUN_EXPECTATIONS)].copy()
    if n.empty:
        raise ValueError("No rows match RUN_EXPECTATIONS")
    for col in METRICS:
        numeric_gc(n, col, args.nhphyml_tsv)
    reject_duplicates(n, ["nhphyml_run_pk", "genome_pk", "ortholog_pk"],
                      "(run, genome, ortholog) results")
    owners = n.groupby("ortholog_pk")["genome_pk"].nunique()
    if (owners > 1).any():
        raise ValueError("ortholog_pk belongs to multiple genomes. Joining on it alone "
                         "is unsafe; supply genome-specific ortholog IDs.")

    f = load_tsv(args.flank_tsv, ["ortholog_pk", "gc"])
    integer_key(f, "ortholog_pk", args.flank_tsv)
    # Restrict to relevant IDs before checking flank-set ambiguity.
    f = f[f["ortholog_pk"].isin(n["ortholog_pk"])].copy()
    if "flank_set_pk" in f:
        f["flank_set_pk"] = f["flank_set_pk"].str.strip().replace(NULL_VALUES, pd.NA)
        integer_key(f, "flank_set_pk", args.flank_tsv)
        sets = sorted(f["flank_set_pk"].unique(), key=int)
        if args.flank_set_pk is not None:
            f = f[f["flank_set_pk"].eq(str(args.flank_set_pk))].copy()
            if f.empty:
                raise ValueError("Selected flank_set_pk has no relevant flank rows")
        elif len(sets) > 1:
            raise ValueError(f"Multiple relevant flank sets {sets}; select --flank-set-pk")
    elif args.flank_set_pk is not None:
        raise ValueError("--flank-set-pk requires the flank_set_pk column")
    reject_duplicates(f, ["ortholog_pk"], "flank ortholog_pk rows")
    numeric_gc(f, "gc", args.flank_tsv)
    joined = n.merge(f[["ortholog_pk", "gc"]].rename(columns={"gc": "flank_gc"}),
                     on="ortholog_pk", how="left", validate="many_to_one", indicator=True)
    return joined


def summarize(joined, allow_incomplete):
    means, audits, issues = [], [], []

    def issue(run, metric, severity, message):
        issues.append(dict(nhphyml_run_pk=run, metric=metric,
                           severity=severity, message=message))

    for run, (expected_genomes, max_orthologs) in RUN_EXPECTATIONS.items():
        rows = joined[joined["nhphyml_run_pk"].eq(run)]
        raw_counts = rows.groupby("genome_pk")["ortholog_pk"].nunique()
        observed = len(raw_counts)
        if observed != expected_genomes:
            severity = "WARNING" if allow_incomplete and 0 < observed < expected_genomes else "ERROR"
            issue(run, "both", severity,
                  f"Input has {observed} genomes; expected {expected_genomes}")
        excessive = raw_counts[raw_counts > max_orthologs]
        if not excessive.empty:
            issue(run, "both", "ERROR", f"Ortholog limit {max_orthologs} exceeded: {excessive.to_dict()}")
        short = int((raw_counts < max_orthologs).sum())
        if short:
            issue(run, "both", "INFO", f"{short} genomes have fewer than {max_orthologs} input orthologs")
        valid_flank = rows["flank_gc"].notna() & rows["flank_gc"].gt(0)
        base = rows[["genome_pk"]].drop_duplicates().set_index("genome_pk")
        base["input_orthologs"] = raw_counts
        base["orthologs_without_flank_row"] = rows[rows["_merge"].eq("left_only")].groupby("genome_pk").size()
        base["orthologs_zero_or_null_flank"] = rows[rows["_merge"].eq("both") & ~valid_flank].groupby("genome_pk").size()
        base = base.fillna(0)
        for metric in METRICS:
            paired = rows[valid_flank & rows[metric].notna()]
            stats = paired.groupby("genome_pk").agg(
                mean_coding_gc=(metric, "mean"), mean_flank_gc=("flank_gc", "mean"),
                n_paired_orthologs=("ortholog_pk", "nunique"))
            out = base.join(stats)
            out["n_paired_orthologs"] = out["n_paired_orthologs"].fillna(0).astype(int)
            for col in ["input_orthologs", "orthologs_without_flank_row", "orthologs_zero_or_null_flank"]:
                out[col] = out[col].astype(int)
            out["orthologs_missing_coding_gc_with_valid_flank"] = (
                rows[valid_flank & rows[metric].isna()].groupby("genome_pk").size()
                .reindex(out.index, fill_value=0).astype(int))
            out["nhphyml_run_pk"] = run
            out["metric"] = metric
            out["expected_genomes"] = expected_genomes
            out["max_orthologs_per_genome"] = max_orthologs
            out["included_in_plot"] = out["n_paired_orthologs"].gt(0)
            means.append(out.reset_index())
            n_points = int(out["included_in_plot"].sum())
            if n_points != expected_genomes:
                severity = "WARNING" if allow_incomplete and 3 <= n_points < expected_genomes else "ERROR"
                issue(run, metric, severity,
                      f"After paired filtering, {n_points} genome points; expected {expected_genomes}")
            if n_points < 3:
                issue(run, metric, "ERROR", "At least three genome means are required")
            if (out["n_paired_orthologs"] > max_orthologs).any():
                issue(run, metric, "ERROR", "Paired ortholog count exceeds expected maximum")
            audits.append(dict(
                nhphyml_run_pk=run, metric=metric, expected_genomes=expected_genomes,
                input_genomes=observed, plotted_genomes=n_points,
                maximum_orthologs_per_genome=max_orthologs,
                input_ortholog_count_min=int(raw_counts.min()) if observed else 0,
                input_ortholog_count_max=int(raw_counts.max()) if observed else 0,
                paired_ortholog_count_min=int(stats.n_paired_orthologs.min()) if n_points else 0,
                paired_ortholog_count_max=int(stats.n_paired_orthologs.max()) if n_points else 0,
                missing_flank_rows=int(rows["_merge"].eq("left_only").sum()),
                zero_or_null_flank_rows=int((rows["_merge"].eq("both") & ~valid_flank).sum()),
                missing_coding_rows_with_valid_flank=int((valid_flank & rows[metric].isna()).sum()),
            ))
    return (pd.concat(means, ignore_index=True), pd.DataFrame(audits),
            pd.DataFrame(issues, columns=["nhphyml_run_pk", "metric", "severity", "message"]))


def correlation_table(means):
    records = []
    for run in RUN_EXPECTATIONS:
        for metric in METRICS:
            points = means[(means.nhphyml_run_pk == run) & (means.metric == metric)
                           & means.included_in_plot]
            x, y = points.mean_coding_gc.to_numpy(), points.mean_flank_gc.to_numpy()
            result = dict(nhphyml_run_pk=run, metric=metric, n_genomes=len(points),
                          slope=np.nan, intercept=np.nan, pearson_r=np.nan,
                          r_squared=np.nan, p_value=np.nan, slope_standard_error=np.nan,
                          status="undefined: fewer than 3 points or constant coordinate")
            if len(points) >= 3 and np.ptp(x) > 0 and np.ptp(y) > 0:
                fit = linregress(x, y)
                result.update(slope=fit.slope, intercept=fit.intercept,
                              pearson_r=fit.rvalue, r_squared=fit.rvalue ** 2,
                              p_value=fit.pvalue, slope_standard_error=fit.stderr, status="ok")
            records.append(result)
    return pd.DataFrame(records)


def padded_limits(values):
    lo, hi = float(np.min(values)), float(np.max(values))
    pad = (hi - lo) * 0.08 if hi > lo else max(abs(lo) * 0.05, 0.01)
    return lo - pad, hi + pad


def plot_figure(means, correlations, metric, output_dir):
    if number_of_ticks < 2:
        raise ValueError("number_of_ticks must be at least 2")
    fig, axes = plt.subplots(2, 2, figsize=tuple(v / 2.54 for v in figure_size))
    # Leave room between panels for labels and reserve an annotation strip.
    fig.subplots_adjust(left=0.14, right=0.94, bottom=0.11, top=0.98,
                        wspace=0.57, hspace=0.62)
    coding_label = "GC3" if metric == "gc3_nhphyml" else "GC4"
    for index, (ax, run) in enumerate(zip(axes.flat, RUN_EXPECTATIONS)):
        points = means[(means.nhphyml_run_pk == run) & (means.metric == metric)
                       & means.included_in_plot]
        stat = correlations[(correlations.nhphyml_run_pk == run)
                            & (correlations.metric == metric)].iloc[0]
        x, y = points.mean_coding_gc.to_numpy(), points.mean_flank_gc.to_numpy()
        ax.plot(x, y, linestyle="none", marker="o", markersize=point_diameter,
                color=point_color, markeredgewidth=0, zorder=3)
        ax.set_xlim(*padded_limits(x))
        ax.set_ylim(*padded_limits(y))
        if stat["status"] == "ok":
            line_x = np.array([x.min(), x.max()])
            ax.plot(line_x, stat.slope * line_x + stat.intercept,
                    color=regression_line_color, linestyle=regression_line_style,
                    linewidth=regression_line_width, zorder=2)
            description = f"r = {stat.pearson_r:.3f}; R² = {stat.r_squared:.3f}\np = {stat.p_value:.2g}"
        else:
            description = "Correlation undefined"
        expected, cap = RUN_EXPECTATIONS[run]
        pair_min = int(points.n_paired_orthologs.min())
        pair_max = int(points.n_paired_orthologs.max())
        annotation = (f"Run {run}; n = {len(points)}/{expected} genomes\n"
                      f"{pair_min}–{pair_max} paired orthologs (≤{cap})\n{description}")
        ax.text(0, 1.04, annotation, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=descriptive_statistics_font_size, color="black")
        ax.text(-0.27, 1.03, panel_labels[index], transform=ax.transAxes,
                ha="left", va="bottom", fontsize=axes_font_size, fontweight="bold")
        ax.set_xlabel(f"Mean ortholog {coding_label}", fontsize=axes_font_size)
        ax.set_ylabel("Mean flank GC", fontsize=axes_font_size)
        for name in ["bottom", "left"]:
            ax.spines[name].set_linewidth(axes_width)
            ax.spines[name].set_color(axes_color)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.minorticks_off()
        for axis in [ax.xaxis, ax.yaxis]:
            axis.set_major_locator(LinearLocator(number_of_ticks))
            formatter = ScalarFormatter(useOffset=False)
            formatter.set_scientific(False)
            axis.set_major_formatter(formatter)
        ax.tick_params(axis="both", which="major", width=axes_width, length=5,
                       colors=axes_color, labelsize=axes_font_size)
        ax.grid(False)
    # Annotation strips extend above axes; lower the top row to keep them in page.
    fig.subplots_adjust(top=0.85)
    stem = "correlation_flank_vs_" + coding_label.lower() + "_per_nhphyml_run"
    for ext in ["png", "pdf"]:
        path = output_dir / f"{stem}.{ext}"
        fig.savefig(path, dpi=output_dpi, facecolor="white")
        print(f"[OK] {path}")
    plt.close(fig)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--flank-tsv", type=Path, required=True)
    parser.add_argument("--nhphyml-tsv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--flank-set-pk", type=int, help="Select one flank set if multiple are present")
    parser.add_argument("--allow-incomplete-genomes", action="store_true",
                        help="Permit fewer than expected genomes; report actual n")
    parser.add_argument("--validate-only", action="store_true", help="Write audit tables without figures")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        if len(RUN_EXPECTATIONS) != 4:
            raise ValueError("RUN_EXPECTATIONS must define exactly four panels")
        if args.flank_set_pk is not None and args.flank_set_pk < 1:
            raise ValueError("--flank-set-pk must be positive")
        joined = prepare_inputs(args)
        means, audits, issues = summarize(joined, args.allow_incomplete_genomes)
        correlations = correlation_table(means)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for table, name in [(means, "genome_means.tsv"), (audits, "run_validation.tsv"),
                            (issues, "validation_issues.tsv"),
                            (correlations, "correlation_statistics.tsv")]:
            table.to_csv(args.output_dir / name, sep="\t", index=False, na_rep="\\N")
        for row in issues.itertuples(index=False):
            print(f"[{row.severity}] Run {row.nhphyml_run_pk} ({row.metric}): {row.message}")
        print(audits[["nhphyml_run_pk", "metric", "input_genomes", "plotted_genomes",
                      "paired_ortholog_count_min", "paired_ortholog_count_max"]].to_string(index=False))
        if issues.severity.eq("ERROR").any():
            print(f"[ERROR] Validation failed; inspect audit TSVs in {args.output_dir}. "
                  "No figures generated by this invocation.", file=sys.stderr)
            return 1
        if not args.validate_only:
            for metric in METRICS:
                plot_figure(means, correlations, metric, args.output_dir)
        print(f"[OK] Validation passed. Audit tables: {args.output_dir}")
        return 0
    except (ValueError, OSError, pd.errors.ParserError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
