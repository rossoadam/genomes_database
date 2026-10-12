#!/usr/bin/env python3
"""Plot mean versus SD of GC3 and GC4 for nhPhyML runs 1–4.

Usage:
    python plot_nhphyml_mean_sd_correlations_v2.py nhphyml_species_summary.tsv
    python plot_nhphyml_mean_sd_correlations_v2.py INPUT.tsv --outdir figures

Dependencies: numpy, scipy, matplotlib
    python -m pip install numpy scipy matplotlib

Outputs: two four-panel figures (PNG and PDF) and correlation_statistics.tsv.
Default output directory: /genomes/records/statistics_body_size_vs_gc
Use --outdir to override this path. Panel annotations show equation, R², n,
and the two-sided Pearson p-value. n is the number of valid genome pairs
in each panel; it is not fixed at the number of panels (four).
Each point is one genome within a run. The x variable is mean GC; y is SD GC.
Values retain the input units. Sample size counts finite mean/SD pairs, not
orthologs or runs. Pearson tests are two-sided and p-values are unadjusted.
The dashed line is an ordinary least-squares fit with an intercept; its R²
is Pearson r squared. Missing/invalid pairs are excluded separately for GC3
and GC4 and their counts are reported. Duplicate genome/run pairs are errors.
"""

# ======================= EDIT PLOTTING SETTINGS HERE =======================
axes_width = 3                  # Matplotlib widths are in points.
figure_size = (16.5, 16.5)       # Full four-panel figure: width, height in cm.
regression_line_width = 3
regression_line_color = "black"
regression_line_style = "--"
point_diameter = 5             # Marker size in points; use plot(), not scatter s.
point_color = "black"
axes_color = "black"
axes_font_size = 9              # Axis labels and major tick labels.
key_font_size = 6.5
# Descriptive-statistics settings (the statistics are computed from the data).
descriptive_statistics_font_size = 6.5
descriptive_statistics_position = (0.97, 0.97)  # Axes-relative, top right.
descriptive_statistics_precision = 3          # Significant digits.
show_pearson_r = False          # Optional; r is always included in the results TSV.
show_p_value = True             # Always included in the results TSV.
statistics_background_alpha = 0.9
key_position = "lower left"
point_key_label = "genome"
regression_key_label = "regression line"
number_of_major_ticks = 4       # EXACT same count on x and y; no minor ticks.
tick_length = 5
tick_label_precision = 3
run_ids = (1, 2, 3, 4)          # Top left, top right, bottom left, bottom right.
output_dpi = 300
output_formats = ("png", "pdf")
output_directory = "/Users/rossoaa/projects/genomes/records/statistics_body_size_vs_gc/nhphyml_stdv_vs_mean"  # --outdir overrides.
# Extra y space leaves room for statistics and the key inside each panel.
y_padding_below = 0.35
y_padding_above = 0.65
x_padding = 0.10
# ========================================================================

import argparse
import csv
import sys
from pathlib import Path

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    from scipy.stats import linregress, pearsonr
except ImportError as exc:
    raise SystemExit("Missing dependency. Install with: python -m pip install numpy scipy matplotlib") from exc

REQUIRED_COLUMNS = {
    "nhphyml_species_summary_pk", "nhphyml_run_pk", "genome_pk",
    "mean_gc3", "mean_gc4", "sd_gc3", "sd_gc4", "n_orthologs", "created_at",
}


def numeric(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def read_input(path):
    groups = {run: [] for run in run_ids}
    seen = set()
    ignored = 0
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Missing columns: " + ", ".join(sorted(missing)))
        for line, row in enumerate(reader, 2):
            try:
                run = int(row["nhphyml_run_pk"])
                genome = int(row["genome_pk"])
            except (ValueError, TypeError):
                raise ValueError(f"Line {line}: run and genome IDs must be integers.")
            if run not in groups:
                ignored += 1
                continue
            pair = (run, genome)
            if pair in seen:
                raise ValueError(f"Line {line}: duplicate genome/run pair {pair}.")
            seen.add(pair)
            groups[run].append(row)
    if ignored:
        print(f"[INFO] Ignored {ignored} rows outside requested runs {run_ids}.")
    return groups


def limits(values, below, above):
    if not len(values):
        return (0.0, 1.0)
    low, high = float(np.min(values)), float(np.max(values))
    span = high - low
    if span == 0:
        span = max(abs(low) * 0.1, 0.01)
    return low - below * span, high + above * span


def analyze(rows, metric, run):
    pairs = np.array([(numeric(r[f"mean_{metric}"]), numeric(r[f"sd_{metric}"]))
                      for r in rows], dtype=float).reshape(-1, 2)
    valid = np.isfinite(pairs).all(axis=1)
    pairs = pairs[valid]
    if len(pairs) and (pairs[:, 1] < 0).any():
        raise ValueError(f"Run {run}, {metric}: standard deviation cannot be negative.")
    x, y = pairs[:, 0], pairs[:, 1]
    record = dict(metric=metric.upper(), nhphyml_run_pk=run, n_input=len(rows),
                  n=len(x), n_excluded=int((~valid).sum()), pearson_r=np.nan,
                  pearson_p_value=np.nan, r_squared=np.nan, slope=np.nan,
                  intercept=np.nan, slope_standard_error=np.nan, status="ok")
    if len(x) < 2:
        record["status"] = "fewer_than_two_valid_pairs"
    elif np.ptp(x) == 0:
        record["status"] = "constant_mean"
    else:
        fit = linregress(x, y)
        record.update(slope=fit.slope, intercept=fit.intercept,
                      slope_standard_error=fit.stderr)
        if np.ptp(y) == 0:
            record["status"] = "constant_sd_pearson_undefined"
        else:
            r, p = pearsonr(x, y)
            record.update(pearson_r=float(r), pearson_p_value=float(p),
                          r_squared=float(r * r))
    return x, y, record


def statistics_text(record):
    digits = descriptive_statistics_precision
    def fmt(value):
        return f"{value:.{digits}g}" if np.isfinite(value) else "NA"
    m, b = record["slope"], record["intercept"]
    equation = (f"y = {fmt(m)}x {'+' if b >= 0 else '−'} {fmt(abs(b))}"
                if np.isfinite(m) else "Regression unavailable")
    lines = [equation, f"R² = {fmt(record['r_squared'])}"]
    if show_pearson_r:
        lines.append(f"Pearson r = {fmt(record['pearson_r'])}")
    lines.append(f"n = {record['n']}")
    if show_p_value:
        lines.append(f"p = {fmt(record['pearson_p_value'])}")
    return "\n".join(lines)


def make_figure(groups, metric, outdir):
    fig, axes = plt.subplots(2, 2, figsize=tuple(v / 2.54 for v in figure_size))
    fig.subplots_adjust(left=0.12, right=0.97, bottom=0.12, top=0.97,
                        wspace=0.40, hspace=0.48)
    records = []
    for panel, (ax, run) in enumerate(zip(axes.flat, run_ids)):
        x, y, record = analyze(groups[run], metric, run)
        records.append(record)
        ax.plot(x, y, linestyle="none", marker="o", markersize=point_diameter,
                markerfacecolor=point_color, markeredgecolor=point_color,
                markeredgewidth=0, zorder=3)
        if np.isfinite(record["slope"]):
            line_x = np.array([np.min(x), np.max(x)])
            ax.plot(line_x, record["slope"] * line_x + record["intercept"],
                    color=regression_line_color, linestyle=regression_line_style,
                    linewidth=regression_line_width, zorder=2)
        ax.set_xlim(*limits(x, x_padding, x_padding))
        ax.set_ylim(*limits(y, y_padding_below, y_padding_above))
        for axis, bounds in [(ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())]:
            axis.set_major_locator(FixedLocator(np.linspace(*bounds, number_of_major_ticks)))
            axis.set_minor_locator(NullLocator())
            axis.set_major_formatter(FuncFormatter(
                lambda value, pos: f"{value:.{tick_label_precision}g}"))
        ax.tick_params(axis="both", which="major", color=axes_color,
                       labelcolor=axes_color, width=axes_width, length=tick_length,
                       labelsize=axes_font_size)
        for name, spine in ax.spines.items():
            spine.set_visible(name in ("left", "bottom"))
            spine.set_color(axes_color)
            spine.set_linewidth(axes_width)
        ax.set_xlabel(f"Mean {metric.upper()}\n({chr(65 + panel)}) nhphyml_run_pk = {run}",
                      fontsize=axes_font_size, color=axes_color)
        ax.set_ylabel(f"SD {metric.upper()}", fontsize=axes_font_size, color=axes_color)
        ax.text(*descriptive_statistics_position, statistics_text(record),
                transform=ax.transAxes, ha="right", va="top",
                fontsize=descriptive_statistics_font_size, color="black",
                bbox=dict(facecolor="white", edgecolor="none",
                          alpha=statistics_background_alpha, pad=1.5), zorder=5)
        handles = []
        if np.isfinite(record["slope"]):
            handles.append(Line2D([], [], color=regression_line_color,
                                  linestyle=regression_line_style,
                                  linewidth=regression_line_width,
                                  label=regression_key_label))
        handles.append(Line2D([], [], linestyle="none", marker="o",
                              markersize=point_diameter, markerfacecolor=point_color,
                              markeredgecolor=point_color, markeredgewidth=0,
                              label=point_key_label))
        ax.legend(handles=handles, loc=key_position, fontsize=key_font_size,
                  frameon=False, handlelength=2.0, borderaxespad=0.4)
        print(f"[INFO] {metric.upper()} run {run}: n={record['n']}, "
              f"excluded={record['n_excluded']}, status={record['status']}")
    for extension in output_formats:
        path = outdir / f"nhphyml_{metric}_mean_vs_sd.{extension}"
        # No bbox_inches='tight': preserve the exact configured figure dimensions.
        fig.savefig(path, dpi=output_dpi, facecolor="white")
        print(f"[OK] {path}")
    plt.close(fig)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_tsv", nargs="?", type=Path,
                        default=Path("nhphyml_species_summary.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path(output_directory),
                        help=f"Output directory (default: {output_directory}).")
    args = parser.parse_args()
    if len(run_ids) != 4 or len(set(run_ids)) != 4:
        parser.error("run_ids must contain four distinct run IDs.")
    if number_of_major_ticks < 2:
        parser.error("number_of_major_ticks must be at least 2.")
    try:
        groups = read_input(args.input_tsv)
        args.outdir.mkdir(parents=True, exist_ok=True)
        records = []
        for metric in ("gc3", "gc4"):
            records.extend(make_figure(groups, metric, args.outdir))
        path = args.outdir / "correlation_statistics.tsv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(records)
        print(f"[OK] {path}")
    except (OSError, ValueError) as exc:
        parser.exit(1, f"[ERROR] {exc}\n")


if __name__ == "__main__":
    main()
