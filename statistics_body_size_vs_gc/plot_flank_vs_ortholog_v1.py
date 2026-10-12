#!/usr/bin/env python3
"""Plot GC3/GC4 versus flank GC, retaining individual genome observations.

Dependencies: pandas, numpy, scipy, matplotlib
    python -m pip install pandas numpy scipy matplotlib

Example (the default output directory is already set below):
    python plot_flank_vs_ortholog_v1.py \
        --flank-tsv /Users/rossoaa/projects/genomes/records/sql_tsvs/flank_compleasm.tsv \
        --nhphyml-tsv /Users/rossoaa/projects/genomes/records/sql_tsvs/nhphyml_ortholog_results.tsv

Use --run-pk 1 --flank-set-pk 1 to select one run and flank set. Otherwise,
produce one two-panel figure per matched run/flank-set combination. The whole
figure is 16.5 x 16.5 cm, with GC3 above GC4. GC3/GC4 are x; flank GC is y.

Join on ortholog_pk within each flank set, retaining every nhPhyML row and
its genome_pk. No genome averaging, ortholog aggregation, or subsampling.
Each flank set must contain at most one flank row per ortholog_pk. A key
repeated across genomes in nhPhyML is retained and receives the matching
flank GC on every row. If those genomes actually have different flanks, the
supplied flank table must identify their distinct genome-specific ortholog
keys: it has no genome_pk column with which to resolve such a distinction.
In the project schema ortholog_pk identifies a genome-specific copy, while
odb12_id identifies an ortholog family shared across genomes.

SQL NULL (\\N), empty cells, NA/NaN, and NULL are accepted. Missing,
nonnumeric, infinite, or out-of-range GC values are excluded independently
in each panel and counted in the summary. Numeric GC units are preserved
as supplied (fractions or percentages); axis labels do not assert units.
Status values are retained without assumptions; optional status filters
are available. Pearson r/p and ordinary least-squares regression describe
pooled observations and do not adjust for phylogeny or repeated families.

Outputs: PNG/PDF figures, joined_observations.tsv.gz (including unmatched
nhPhyML rows), correlation_summary.tsv, and join_summary.tsv.
"""

# ------------------------- Editable settings -------------------------
axes_width = 3
figure_size = (16.5, 16.5)  # cm; entire two-panel figure
regression_line_width = 5
circular_point_diameter = 10  # points; plot markersize, not scatter area
circular_point_color = "black"
axes_color = "black"
axes_font_size = 12  # labels and tick labels
key_font_size = 10  # retained as requested; unused because there is no key
descriptive_statistics_font_size = 9
tick_count = 5  # exactly this many major ticks on both axes in every panel
output_dpi = 300
output_directory = (
    "/Users/rossoaa/projects/genomes/records/statistics_body_size_vs_gc/"
    "correlation_flank_vs_ortholog"
)
# ---------------------------------------------------------------------

import argparse
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullLocator, ScalarFormatter
    import numpy as np
    import pandas as pd
    from scipy.stats import linregress, pearsonr
except ImportError as exc:
    raise SystemExit(
        f"Missing dependency: {exc}. Install with: "
        "python -m pip install pandas numpy scipy matplotlib"
    ) from exc


NULL_VALUES = [r"\N", "", "NA", "N/A", "NaN", "nan", "NULL", "null", "None"]


def read_table(path, required):
    """Keep integer IDs as strings to avoid float precision loss."""
    frame = pd.read_csv(
        path, sep="\t", dtype="string", na_values=NULL_VALUES,
        keep_default_na=False, encoding="utf-8-sig", compression="infer",
    )
    frame.columns = frame.columns.str.strip()
    if frame.columns.duplicated().any():
        raise ValueError(f"{path}: duplicate column names")
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f"{path}: missing required columns: {', '.join(missing)}")
    for column in frame.columns:
        frame[column] = frame[column].str.strip().replace(NULL_VALUES, pd.NA)
    if frame.empty:
        raise ValueError(f"{path}: no data rows")
    return frame


def validate_ids(frame, columns, label):
    for column in columns:
        valid = frame[column].str.fullmatch(r"[0-9]+", na=False)
        if not valid.all():
            examples = frame.loc[~valid, column].head(5).tolist()
            raise ValueError(f"{label}: missing/invalid {column}: {examples}")
        # Canonicalize leading zeros without coercing BIGINTs to floats.
        frame[column] = frame[column].str.lstrip("0").replace("", "0")


def require_unique(frame, columns, label):
    duplicates = frame.duplicated(columns, keep=False)
    if duplicates.any():
        examples = frame.loc[duplicates, columns].head(5).to_dict("records")
        raise ValueError(
            f"{label}: duplicate {columns}; joining these could inflate point "
            f"counts. Example keys: {examples}"
        )


def filter_table(frame, column, selections, label):
    if not selections:
        return frame
    if column not in frame.columns:
        raise ValueError(f"{label}: cannot filter; column {column} is absent")
    missing = set(selections) - set(frame[column].dropna())
    if missing:
        raise ValueError(f"{label}: {column} value(s) not found: {sorted(missing)}")
    return frame.loc[frame[column].isin(selections)].copy()


def numeric_gc(series):
    values = pd.to_numeric(series, errors="coerce").to_numpy(
        dtype=float, na_value=np.nan
    )
    valid = np.isfinite(values) & (values >= 0) & (values <= 100)
    return values, valid


def padded_limits(values):
    if len(values) == 0:
        return (0.0, 1.0)
    low, high = float(np.min(values)), float(np.max(values))
    padding = (high - low) * 0.07 if high > low else max(abs(low) * 0.07, 0.05)
    return (low - padding, high + padding)


def format_axes(ax, x, y, x_label):
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(axes_width)
        ax.spines[side].set_color(axes_color)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlim(*padded_limits(x))
    ax.set_ylim(*padded_limits(y))
    ax.set_xticks(np.linspace(*ax.get_xlim(), tick_count))
    ax.set_yticks(np.linspace(*ax.get_ylim(), tick_count))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())
    for axis in (ax.xaxis, ax.yaxis):
        formatter = ScalarFormatter(useOffset=False)
        formatter.set_scientific(False)
        axis.set_major_formatter(formatter)
    ax.tick_params(
        axis="both", which="major", width=axes_width, length=6,
        colors=axes_color, labelsize=axes_font_size, direction="out",
    )
    ax.set_xlabel(x_label, fontsize=axes_font_size, color=axes_color)
    ax.set_ylabel("Flank GC", fontsize=axes_font_size, color=axes_color)
    ax.grid(False)


def plot_panel(ax, frame, metric, x_label):
    x_all, x_valid = numeric_gc(frame[metric])
    y_all, y_valid = numeric_gc(frame["flank_gc"])
    valid = x_valid & y_valid
    x, y = x_all[valid], y_all[valid]
    selected = frame.loc[valid]
    stats = {
        "metric": metric, "n_joined_rows": len(frame), "n_points": len(x),
        "n_excluded_rows": int((~valid).sum()),
        "n_invalid_or_missing_x": int((~x_valid).sum()),
        "n_invalid_or_missing_flank_gc": int((~y_valid).sum()),
        "n_genomes": selected["genome_pk"].nunique(),
        "n_ortholog_pks": selected["ortholog_pk"].nunique(),
        "pearson_r": np.nan, "pearson_p": np.nan, "r_squared": np.nan,
        "slope": np.nan, "intercept": np.nan, "slope_stderr": np.nan,
        "correlation_status": "insufficient_points",
    }
    ax.plot(
        x, y, linestyle="none", marker="o",
        markersize=circular_point_diameter,
        markerfacecolor=circular_point_color,
        markeredgecolor=circular_point_color, markeredgewidth=0,
        rasterized=True, zorder=2,
    )
    annotation = f"n = {len(x):,}\nCorrelation unavailable"
    if len(x) >= 2:
        if np.ptp(x) == 0:
            stats["correlation_status"] = "constant_x"
            annotation = f"n = {len(x):,}\nConstant x; correlation undefined"
        else:
            fit = linregress(x, y)
            stats.update(slope=fit.slope, intercept=fit.intercept,
                         slope_stderr=fit.stderr)
            xx = np.linspace(float(x.min()), float(x.max()), 200)
            ax.plot(xx, fit.intercept + fit.slope * xx, "--", color="black",
                    linewidth=regression_line_width, zorder=3)
            if np.ptp(y) == 0:
                stats["correlation_status"] = "constant_y"
                annotation = f"n = {len(x):,}\nConstant y; correlation undefined"
            else:
                r, p = pearsonr(x, y)
                stats.update(pearson_r=r, pearson_p=p, r_squared=r ** 2,
                             correlation_status="ok")
                p_text = "< 1e-300" if p == 0 else f"= {p:.3g}"
                sign = "+" if fit.intercept >= 0 else "−"
                annotation = (
                    f"y = {fit.slope:.3g}x {sign} {abs(fit.intercept):.3g}\n"
                    f"r = {r:.3f}; $R^2$ = {r ** 2:.3f}\n"
                    f"p {p_text}\nn = {len(x):,}"
                )
    format_axes(ax, x, y, x_label)
    ax.text(
        0.98, 0.98, annotation, transform=ax.transAxes,
        ha="right", va="top", fontsize=descriptive_statistics_font_size,
        color="black", zorder=4,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2},
    )
    # No title, legend/key, or minor ticks.
    return stats


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--flank-tsv", required=True, type=Path,
                        help="flank_compleasm.tsv (or flanks_compleasm.tsv); .gz accepted")
    parser.add_argument("--nhphyml-tsv", required=True, type=Path,
                        help="nhphyml_ortholog_results.tsv; .gz accepted")
    parser.add_argument("--output-dir", type=Path, default=Path(output_directory))
    parser.add_argument("--run-pk", type=int, nargs="+",
                        help="Select one or more nhphyml_run_pk values; default all")
    parser.add_argument("--flank-set-pk", type=int, nargs="+",
                        help="Select one or more flank_set_pk values; default all")
    parser.add_argument("--nhphyml-status", nargs="+",
                        help="Optional exact status values to include; default all")
    parser.add_argument("--flank-status", nargs="+",
                        help="Optional exact status values to include; default all")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    flank = read_table(args.flank_tsv, ["flank_pk", "ortholog_pk", "flank_set_pk", "gc"])
    nh = read_table(args.nhphyml_tsv, [
        "nhphyml_ortholog_result_pk", "nhphyml_run_pk", "ortholog_pk",
        "genome_pk", "gc3_nhphyml", "gc4",
    ])
    validate_ids(flank, ["flank_pk", "ortholog_pk", "flank_set_pk"], "Flanks")
    validate_ids(nh, ["nhphyml_ortholog_result_pk", "nhphyml_run_pk", "ortholog_pk",
                      "genome_pk"], "nhPhyML")
    require_unique(flank, ["flank_pk"], "Flanks")
    require_unique(nh, ["nhphyml_ortholog_result_pk"], "nhPhyML")
    require_unique(flank, ["ortholog_pk", "flank_set_pk"], "Flanks")
    require_unique(nh, ["nhphyml_run_pk", "ortholog_pk", "genome_pk"], "nhPhyML")
    nh = filter_table(nh, "nhphyml_run_pk",
                      [str(i) for i in args.run_pk] if args.run_pk else None, "nhPhyML")
    flank = filter_table(flank, "flank_set_pk",
                         [str(i) for i in args.flank_set_pk] if args.flank_set_pk else None,
                         "Flanks")
    nh = filter_table(nh, "status", args.nhphyml_status, "nhPhyML")
    flank = filter_table(flank, "status", args.flank_status, "Flanks")
    if nh.empty or flank.empty:
        raise ValueError("No input rows remain after filtering")
    if tick_count < 2:
        raise ValueError("tick_count must be at least 2")
    out = args.output_dir.expanduser()
    out.mkdir(parents=True, exist_ok=True)
    summary, join_summary = [], []
    joined_path = out / "joined_observations.tsv.gz"
    wrote_joined = False
    n_figures = 0
    for run_pk, run in nh.groupby("nhphyml_run_pk", sort=True):
        for set_pk, flank_set in flank.groupby("flank_set_pk", sort=True):
            # Explicitly many-to-one: repeated genome rows are never collapsed.
            renamed = flank_set.rename(columns={"gc": "flank_gc", "status": "flank_status"})
            joined = run.merge(renamed, on="ortholog_pk", how="left",
                               validate="many_to_one", indicator="flank_join_status")
            # Keep the selected set identifier even for unmatched observations.
            joined["flank_set_pk"] = set_pk
            joined.to_csv(joined_path, sep="\t", index=False, na_rep=r"\N",
                          compression="gzip", mode="a" if wrote_joined else "w",
                          header=not wrote_joined)
            wrote_joined = True
            matched = joined.loc[joined["flank_join_status"] == "both"].copy()
            join_summary.append({
                "nhphyml_run_pk": run_pk, "flank_set_pk": set_pk,
                "n_nhphyml_rows": len(run), "n_matched_rows": len(matched),
                "n_unmatched_nhphyml_rows": len(joined) - len(matched),
                "n_unmatched_flank_rows": int((~flank_set["ortholog_pk"].isin(
                    run["ortholog_pk"])).sum()),
            })
            print(f"[INFO] run={run_pk}, flank_set={set_pk}: "
                  f"{len(matched):,}/{len(run):,} nhPhyML rows matched")
            if matched.empty:
                print("[WARNING] No shared ortholog_pk values; no figure for this combination")
                continue
            fig, axes = plt.subplots(
                2, 1, figsize=tuple(cm / 2.54 for cm in figure_size),
                layout="constrained",
            )
            fig.set_constrained_layout_pads(w_pad=0.06, h_pad=0.06, hspace=0.10)
            stem = f"flank_vs_ortholog_run_{run_pk}_flank_set_{set_pk}"
            for ax, metric, label in zip(
                axes, ["gc3_nhphyml", "gc4"], ["GC3 (nhPhyML)", "GC4"]
            ):
                stats = plot_panel(ax, matched, metric, label)
                stats.update(nhphyml_run_pk=run_pk, flank_set_pk=set_pk, figure=stem)
                summary.append(stats)
                print(f"[INFO] {metric}: n={stats['n_points']:,}, "
                      f"excluded={stats['n_excluded_rows']:,}, "
                      f"status={stats['correlation_status']}")
            for extension in ("png", "pdf"):
                # Do not use bbox_inches='tight': retain the specified physical size.
                fig.savefig(out / f"{stem}.{extension}", dpi=output_dpi, facecolor="white")
            plt.close(fig)
            n_figures += 1
    pd.DataFrame(join_summary).to_csv(out / "join_summary.tsv", sep="\t", index=False)
    if summary:
        stats_frame = pd.DataFrame(summary)
        first = ["nhphyml_run_pk", "flank_set_pk", "metric"]
        stats_frame = stats_frame[first + [c for c in stats_frame.columns if c not in first]]
        stats_frame.to_csv(out / "correlation_summary.tsv", sep="\t", index=False, na_rep=r"\N")
    if n_figures == 0:
        raise ValueError(f"No matched run/flank-set combinations; inspect {out / 'join_summary.tsv'}")
    print(f"[OK] Saved {n_figures} two-panel figure(s), PNG/PDF, and audit TSVs to {out}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, pd.errors.ParserError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
