#!/usr/bin/env python3
"""Plot the descriptive tradeoff between genomes included and orthologs retained.

Requirements: Python >= 3.9, numpy, matplotlib.
Example:
    python plot_genome_ortholog_tradeoff_v3.py genome_ortholog_tradeoff_example.tsv

Outputs default to the input TSV's directory. Use --output-dir to override.

Figure styling: 16.5 x 16.5 cm; black axes of width 3 points; five major
ticks on each axis; no minor ticks; black dashed regression line of width
5 points; black circular markers of diameter 10 points; no title. The legend
is at the bottom left and descriptive statistics are at the top right.
Edit axes_font_size, key_font_size, and descriptive_statistics_font_size below
to control axis labels/tick labels, legend text, and the statistics box.

Each TSV row represents an evaluated species/genome subset. Required headers are
number_of_species and orthologs_retained, or their human-readable equivalents
"Number of species" and "Orthologs retained". Additional count columns are
optional and checked for consistency when present. One genome per species is
assumed; use --x-column number_of_genomes for a genome-count input instead.

This program summarizes supplied counts; it does not compute ortholog
intersections from compleasm files or choose species subsets. The regression is
descriptive: overlapping/nested subsets are not independent replicates. No
inferential P values or confidence intervals are calculated. Its slope describes
the retention criterion, not biological gene loss or gain in individual genomes.
"""

import argparse
import csv
import math
import re
import sys
from pathlib import Path


FIGURE_SIZE_CM = 16.5
AXIS_WIDTH_PT = 3
REGRESSION_WIDTH_PT = 5
MARKER_DIAMETER_PT = 10
MAJOR_TICK_COUNT = 5

# Font sizes are in points. Axes size controls both axis labels and tick labels.
axes_font_size = 10
key_font_size = 10
descriptive_statistics_font_size = 10

OPTIONAL_COUNTS = (
    'total_eligible_species', 'number_of_species_omitted',
    'orthologs_identified_by_compleasm', 'orthologs_filtered',
)


def normalize_header(value):
    return re.sub(r'[^a-z0-9]+', '_', value.strip().lower()).strip('_')


def count(value, column, line_number):
    # Accept integers and integer-valued decimal notation without rounding.
    from decimal import Decimal, InvalidOperation
    try:
        number = Decimal(value.strip())
    except (InvalidOperation, AttributeError):
        raise ValueError(f'Line {line_number}: {column} must contain an integer; got {value!r}.') from None
    if not number.is_finite() or number < 0 or number != number.to_integral_value():
        raise ValueError(f'Line {line_number}: {column} must be a nonnegative integer; got {value!r}.')
    return int(number)


def read_counts(path, x_column, y_column):
    rows = []
    with path.open('r', encoding='utf-8-sig', newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        try:
            headers = next(reader)
        except StopIteration:
            raise ValueError('Input TSV is empty.') from None
        normalized = [normalize_header(h) for h in headers]
        if not all(normalized) or len(set(normalized)) != len(normalized):
            raise ValueError('Input headers must be nonempty and unique after normalization.')
        missing = {x_column, y_column} - set(normalized)
        if missing:
            raise ValueError('Missing required TSV column(s): ' + ', '.join(sorted(missing)))
        for cells in reader:
            line_number = reader.line_num
            if not cells or not any(c.strip() for c in cells):
                continue
            if len(cells) != len(headers):
                raise ValueError(f'Line {line_number}: expected {len(headers)} tab-separated fields, got {len(cells)}.')
            raw = dict(zip(normalized, cells))
            row = {'source_line': line_number}
            for column in dict.fromkeys((x_column, y_column, *OPTIONAL_COUNTS)):
                if column in raw:
                    row[column] = count(raw[column], column, line_number)
            if row[x_column] == 0:
                raise ValueError(f'Line {line_number}: {x_column} must be greater than zero.')
            # The eligible/omitted columns describe species, so only compare to
            # the selected x-axis when x explicitly counts species.
            if x_column == 'number_of_species':
                eligible = row.get('total_eligible_species')
                omitted = row.get('number_of_species_omitted')
                if eligible is not None and row[x_column] > eligible:
                    raise ValueError(f'Line {line_number}: included species exceed total eligible species.')
                if eligible is not None and omitted is not None and row[x_column] + omitted != eligible:
                    raise ValueError(f'Line {line_number}: included + omitted species must equal total eligible species.')
            if y_column == 'orthologs_retained':
                identified = row.get('orthologs_identified_by_compleasm')
                filtered = row.get('orthologs_filtered')
                if identified is not None and row[y_column] > identified:
                    raise ValueError(f'Line {line_number}: retained orthologs exceed identified orthologs.')
                if identified is not None and filtered is not None and row[y_column] + filtered != identified:
                    raise ValueError(f'Line {line_number}: retained + filtered orthologs must equal identified orthologs.')
            rows.append(row)
    if len(rows) < 2 or len({row[x_column] for row in rows}) < 2:
        raise ValueError('At least two observations with distinct genome/species counts are required.')
    return sorted(rows, key=lambda row: (row[x_column], row['source_line']))


def write_tsv(path, headers, rows):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(headers)
        writer.writerows(rows)


def run(args):
    try:
        import numpy as np
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FixedLocator, NullLocator, StrMethodFormatter
    except ImportError as exc:
        raise ValueError('Missing plotting dependency. Install with: python -m pip install numpy matplotlib') from exc

    x_column = normalize_header(args.x_column)
    y_column = normalize_header(args.y_column)
    if x_column == y_column:
        raise ValueError('--x-column and --y-column must be different.')
    rows = read_counts(args.input_tsv, x_column, y_column)
    x = np.array([row[x_column] for row in rows], dtype=float)
    y = np.array([row[y_column] for row in rows], dtype=float)
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError('Input counts exceed the numeric range supported for regression.')
    x_centered = x - x.mean()
    y_centered = y - y.mean()
    slope = float(np.dot(x_centered, y_centered) / np.dot(x_centered, x_centered))
    intercept = float(y.mean() - slope * x.mean())
    fitted = intercept + slope * x
    residuals = y - fitted
    sse = float(np.dot(residuals, residuals))
    sst = float(np.dot(y_centered, y_centered))
    r_squared = 1.0 - sse / sst if sst > 0 else math.nan
    rmse = float(np.sqrt(sse / len(rows)))

    warnings = []
    if len(rows) == 2:
        warnings.append('Only two observations: the fitted line cannot assess departures from linearity.')
    if len(set(x)) < len(x):
        warnings.append('Repeated genome/species counts: every input row receives equal weight in the fit.')
    if len({row.get('total_eligible_species') for row in rows}) > 1:
        warnings.append('Eligible species totals differ between rows; the fit combines different species pools.')
    if len({row.get('orthologs_identified_by_compleasm') for row in rows}) > 1:
        warnings.append('Identified ortholog totals differ between rows; the fit combines different ortholog pools.')
    if any(x[i + 1] > x[i] and y[i + 1] > y[i] for i in range(len(x) - 1)):
        warnings.append('Retained orthologs increase at some larger subset sizes. Check the subsets and retention criterion; independent selections need not be monotonic.')
    if sst == 0:
        warnings.append('Retained ortholog counts are constant; R-squared is undefined.')

    count_label = 'genome'
    units = f'retained orthologs per additional {count_label} included'
    if abs(slope) < 1e-10:
        change_sentence = f'The fitted change in retained shared orthologs was approximately zero per additional {count_label} included.'
    else:
        direction = 'decrease' if slope < 0 else 'increase'
        change_sentence = (
            f'Across the evaluated subsets ({int(x.min())}–{int(x.max())} {count_label}s included), '
            f'the descriptive linear fit estimated a {direction} of approximately '
            f'{abs(slope):.2f} retained shared orthologs per additional {count_label} included.'
        )
    interpretation = (
        'The slope describes how the number of orthologs satisfying the chosen '
        'shared-across-genomes criterion changes with subset size; it does not '
        'measure biological gene loss or gain. Subsets may overlap or be nested, '
        'so this fit is descriptive and no P values or confidence intervals are reported.'
    )

    output_dir = (args.output_dir.expanduser() if args.output_dir is not None
                  else args.input_tsv.expanduser().resolve().parent)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_dir / args.prefix
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': axes_font_size,
        'axes.labelsize': axes_font_size,
        'xtick.labelsize': axes_font_size, 'ytick.labelsize': axes_font_size,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.edgecolor': 'black', 'axes.linewidth': AXIS_WIDTH_PT,
        'axes.labelcolor': 'black', 'xtick.color': 'black', 'ytick.color': 'black',
        'pdf.fonttype': 42, 'svg.fonttype': 'none',
    })
    figure_size_inches = FIGURE_SIZE_CM / 2.54
    fig, ax = plt.subplots(figsize=(figure_size_inches, figure_size_inches),
                           layout='constrained')
    grid = np.linspace(x.min(), x.max(), 200)
    ax.plot(grid, intercept + slope * grid, color='black', linestyle='--',
            linewidth=REGRESSION_WIDTH_PT,
            label='Descriptive linear fit', zorder=2)
    ax.plot(x, y, linestyle='none', marker='o', markersize=MARKER_DIAMETER_PT,
            color='black', markerfacecolor='black', markeredgecolor='black',
            markeredgewidth=0, label='Evaluated genome subsets', zorder=3)
    if args.label_points:
        for xv, yv in zip(x, y):
            ax.annotate(f'{int(yv):,}', (xv, yv), xytext=(7, 7),
                        textcoords='offset points', fontsize=9)
    ax.set_xlabel(args.x_label)
    ax.set_ylabel(args.y_label)
    # Fixed locators guarantee exactly the same number of major ticks. For x,
    # select five evenly spaced integer positions that cover the observed range.
    x_tick_step = max(1, math.ceil((x.max() - x.min()) / (MAJOR_TICK_COUNT - 1)))
    x_ticks = x.min() + np.arange(MAJOR_TICK_COUNT) * x_tick_step
    # Use a rounded upper y tick to keep labels legible without clipping data.
    y_max = max(float(y.max()), float(fitted.max()), 1.0)
    y_scale = 10 ** math.floor(math.log10(y_max))
    y_tick_max = math.ceil(y_max / y_scale) * y_scale
    y_ticks = np.linspace(0, y_tick_max, MAJOR_TICK_COUNT)
    ax.xaxis.set_major_locator(FixedLocator(x_ticks))
    ax.yaxis.set_major_locator(FixedLocator(y_ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:g}'))
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax.grid(axis='y', alpha=0.16, zorder=0)
    x_span = x_ticks[-1] - x_ticks[0]
    ax.set_xlim(x_ticks[0] - 0.08 * x_span, x_ticks[-1] + 0.08 * x_span)
    y_min = min(0.0, float(fitted.min()))
    y_span = y_tick_max - y_min
    ax.set_ylim(y_min - 0.03 * y_span, y_tick_max + 0.12 * y_span)
    ax.tick_params(axis='both', which='major', colors='black',
                   width=AXIS_WIDTH_PT, length=7)
    ax.minorticks_off()
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(AXIS_WIDTH_PT)
    sign = '+' if intercept >= 0 else '−'
    r2_text = f'{r_squared:.4f}' if math.isfinite(r_squared) else 'undefined'
    annotation = (f'y = {slope:.2f}x {sign} {abs(intercept):.2f}\n'
                  f'R² = {r2_text}; n = {len(rows)} subsets\n'
                  f'Slope: {slope:+.2f} orthologs per genome included')
    ax.text(0.97, 0.97, annotation, transform=ax.transAxes, va='top', ha='right',
            fontsize=descriptive_statistics_font_size,
            bbox={'boxstyle': 'round,pad=0.5', 'facecolor': 'white',
                             'edgecolor': '#DDDDDD', 'alpha': 0.95})
    ax.legend(loc='lower left', frameon=False, fontsize=key_font_size)
    paths = []
    for suffix in ('png', 'pdf', 'svg'):
        path = Path(f'{prefix}.{suffix}')
        # Cropping with bbox_inches='tight' would change the requested physical
        # size. Constrained layout keeps all labels inside this fixed canvas.
        fig.savefig(path, dpi=args.dpi, facecolor='white', bbox_inches=None)
        paths.append(path)
    plt.close(fig)

    summary_path = Path(f'{prefix}_regression_summary.tsv')
    summary = [
        ('input_tsv', str(args.input_tsv.resolve())),
        ('x_column', x_column), ('y_column', y_column),
        ('n_observations', len(rows)), ('n_distinct_subset_sizes', len(set(x))),
        ('minimum_genomes_included', int(x.min())), ('maximum_genomes_included', int(x.max())),
        ('slope', format(slope, '.12g')), ('slope_units', units),
        ('intercept', format(intercept, '.12g')),
        ('r_squared', format(r_squared, '.12g') if math.isfinite(r_squared) else 'NA'),
        ('rmse_orthologs', format(rmse, '.12g')),
        ('fit_type', 'descriptive ordinary least squares; equal weight per input row'),
        ('statement', change_sentence), ('interpretation', interpretation),
        ('warnings', ' | '.join(warnings) or 'None'),
    ]
    write_tsv(summary_path, ('metric', 'value'), summary)
    paths.append(summary_path)
    points_path = Path(f'{prefix}_fitted_values.tsv')
    optional_present = [name for name in OPTIONAL_COUNTS if name in rows[0] and name not in {x_column, y_column}]
    point_headers = ['source_line', x_column, y_column, *optional_present,
                     'fitted_orthologs', 'residual_orthologs']
    write_tsv(points_path, point_headers, (
        [row['source_line'], row[x_column], row[y_column],
         *[row[name] for name in optional_present],
         format(pred, '.12g'), format(resid, '.12g')]
        for row, pred, resid in zip(rows, fitted, residuals)
    ))
    paths.append(points_path)
    statement_path = Path(f'{prefix}_statement.txt')
    statement_path.write_text(change_sentence + '\n\n' + interpretation + '\n', encoding='utf-8')
    paths.append(statement_path)
    print(change_sentence)
    print(f'R² = {r2_text}; RMSE = {rmse:.2f} orthologs.')
    for warning in warnings:
        print('[WARNING] ' + warning, file=sys.stderr)
    for path in paths:
        print('[SAVED] ' + str(path.resolve()))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_tsv', type=Path, help='Tab-separated subset count table.')
    parser.add_argument('--output-dir', type=Path, default=None,
                        help="Output directory (default: input TSV's directory).")
    parser.add_argument('--prefix', default='genome_ortholog_tradeoff', help='Output filename stem (no directories).')
    parser.add_argument('--x-column', default='number_of_species')
    parser.add_argument('--y-column', default='orthologs_retained')
    parser.add_argument('--x-label', default='Number of genomes included (one per species)')
    parser.add_argument('--y-label', default='Shared orthologs retained')
    parser.add_argument('--label-points', action='store_true', help='Label each point with its retained count.')
    parser.add_argument('--dpi', type=int, default=300)
    args = parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', args.prefix):
        parser.error('--prefix must be a filename stem containing letters, numbers, underscores, dots, or hyphens.')
    if args.dpi < 1:
        parser.error('--dpi must be a positive integer.')
    try:
        run(args)
    except (ValueError, OSError, OverflowError) as exc:
        print(f'[ERROR] {exc}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
