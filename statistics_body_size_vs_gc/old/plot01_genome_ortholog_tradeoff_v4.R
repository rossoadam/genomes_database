#!/usr/bin/env Rscript
# plot01_genome_ortholog_tradeoff_v4.R
# R translation of plot01_genome_ortholog_tradeoff_v3.py.
# Requires R >= 4.0; uses only packages distributed with R.
#
# RUN FROM A TERMINAL (replace the input path with your own):
# Rscript plot01_genome_ortholog_tradeoff_v4.R genome_ortholog_tradeoff.tsv
# Rscript plot01_genome_ortholog_tradeoff_v4.R genome_ortholog_tradeoff.tsv \
#   --output-dir /path/to/genomes/records/figures --label-points --diagnostics
# Rscript plot01_genome_ortholog_tradeoff_v4.R --help
#
# Input and output conventions, validation, and main figure styling follow v3.
# Defaults: 16.5 x 16.5 cm; black 3-point axes; five major ticks per axis;
# black dashed 5-point fitted line; approximately 10-point black circles;
# no title; key at bottom left; descriptive statistics at top right.
# Outputs default to the input TSV's directory. For this project, place the
# input under genomes/records/ or pass --output-dir beneath genomes/records/.
#
# WHAT THIS ANALYSIS DOES
# Each input row is an evaluated genome/species subset. It fits an unweighted
# ordinary least-squares line: retained orthologs ~ genomes included.
# It summarizes supplied counts; it does not select subsets or calculate
# ortholog intersections from Compleasm files. One genome per species is
# assumed with the default x column. The slope is a change in the shared set,
# not a measurement of biological gene loss/gain in individual genomes.
#
# READING GUIDE: Zuur, Ieno, Walker, Saveliev & Smith (2009),
# Mixed Effects Models and Extensions in Ecology with R, Springer.
# All page references below are PRINTED book pages, not PDF viewer pages.
#   Section 1.7, pp. 9-10: importing a table with read.table().
#   Section 2.2, pp. 17-19: response, predictor, intercept, slope, assumptions.
#   Appendix A.3, pp. 536-539: lm() and summary() output, including R-squared.
#   Section 2.3.5, pp. 21-22: dependence and its consequences for inference.
#   Section 2.3.6.1, pp. 23-25, and Appendix A.3.2, pp. 542-543:
#     graphical model validation, residuals, normality, homogeneity.
#   Section 2.3.7, pp. 26-27: leverage and influential observations.
#   Appendix A.3.3, pp. 543-546: interpreting and plotting fitted predictions.
#
# IMPORTANT DISTINCTION: fitting a line is not itself a hypothesis test.
# The book also demonstrates coefficient t-tests and model-comparison F-tests
# (Appendix A.3). This script intentionally does NOT report those tests or
# confidence intervals. Overlapping/nested subsets are not independent
# replicates, and four subset summaries do not provide a reliable basis for
# ordinary regression inference. This application of the book's independence
# discussion is specific to our data; it is not a worked example in the book.
# No mixed effects model is added: these aggregate rows do not specify a
# defensible random-effects or residual covariance structure.
#
# --diagnostics adds exploratory plots, following the book's graphical
# workflow. With four rows they have very limited resolution and cannot
# establish independence, normality, or constant variance. A high R-squared
# also does not establish those properties or prove the relation is linear.
# Chapter 3 discusses nonlinear relationships; no nonlinear model selection
# is performed here with only four evaluated subsets.

FIGURE_SIZE_CM <- 16.5
AXIS_WIDTH_PT <- 3
REGRESSION_WIDTH_PT <- 5
MARKER_DIAMETER_PT <- 10
MAJOR_TICK_COUNT <- 5L

# Font sizes in points. The first controls axis labels AND tick labels.
axes_font_size <- 10
key_font_size <- 10
descriptive_statistics_font_size <- 10

OPTIONAL_COUNTS <- c(
  "total_eligible_species", "number_of_species_omitted",
  "orthologs_identified_by_compleasm", "orthologs_filtered"
)

fail <- function(...) stop(..., call. = FALSE)

normalize_header <- function(x) {
  x <- gsub("[^a-z0-9]+", "_", tolower(trimws(x)))
  gsub("^_+|_+$", "", x)
}

help_text <- function() {
  cat(paste0(
    "Usage: Rscript plot01_genome_ortholog_tradeoff_v4.R INPUT.tsv [options]\n\n",
    "Required TSV headers (spaces/capitalization are normalized):\n",
    "  number_of_species, orthologs_retained\n\n",
    "Options:\n",
    "  --output-dir DIR     Default: input TSV's directory\n",
    "  --prefix STEM        Default: genome_ortholog_tradeoff\n",
    "  --x-column NAME      Default: number_of_species\n",
    "  --y-column NAME      Default: orthologs_retained\n",
    "  --x-label TEXT       Default: Number of genomes included (one per species)\n",
    "  --y-label TEXT       Default: Shared orthologs retained\n",
    "  --dpi INTEGER        PNG resolution (default: 300)\n",
    "  --label-points       Annotate retained counts\n",
    "  --diagnostics        Also save a four-panel diagnostic PNG and PDF\n",
    "  -h, --help           Show this help\n\n",
    "Outputs: PNG/PDF/SVG main figure, *_regression_summary.tsv,\n",
    "  *_fitted_values.tsv, *_statement.txt; optional *_diagnostics.png/.pdf.\n",
    "Uses base R; no contributed packages required. No P values or CIs.\n"
  ))
}

parse_args <- function(argv) {
  args <- list(
    input_tsv = NULL, output_dir = NULL, prefix = "genome_ortholog_tradeoff",
    x_column = "number_of_species", y_column = "orthologs_retained",
    x_label = "Number of genomes included (one per species)",
    y_label = "Shared orthologs retained", label_points = FALSE,
    diagnostics = FALSE, dpi = 300L
  )
  if (any(argv %in% c("-h", "--help"))) {
    help_text()
    return(NULL)
  }
  flags <- c("--label-points" = "label_points", "--diagnostics" = "diagnostics")
  options <- c(
    "--output-dir" = "output_dir", "--prefix" = "prefix",
    "--x-column" = "x_column", "--y-column" = "y_column",
    "--x-label" = "x_label", "--y-label" = "y_label", "--dpi" = "dpi"
  )
  i <- 1L
  while (i <= length(argv)) {
    token <- argv[i]
    if (token %in% names(flags)) {
      args[[flags[[token]]]] <- TRUE
    } else if (token %in% names(options)) {
      if (i == length(argv) || startsWith(argv[i + 1L], "--")) {
        fail("Missing value for ", token)
      }
      i <- i + 1L
      args[[options[[token]]]] <- argv[i]
    } else if (startsWith(token, "-")) {
      fail("Unknown option: ", token, ". Use --help.")
    } else if (is.null(args$input_tsv)) {
      args$input_tsv <- token
    } else {
      fail("Unexpected argument: ", token)
    }
    i <- i + 1L
  }
  if (is.null(args$input_tsv)) fail("An input TSV is required. Use --help.")
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", args$prefix)) {
    fail("--prefix must be a filename stem containing letters, numbers, _, ., or -.")
  }
  dpi <- suppressWarnings(as.numeric(args$dpi))
  if (length(dpi) != 1L || !is.finite(dpi) || dpi < 1 ||
      dpi != floor(dpi) || dpi > .Machine$integer.max) {
    fail("--dpi must be a positive integer supported by R.")
  }
  args$dpi <- as.integer(dpi)
  args$x_column <- normalize_header(args$x_column)
  args$y_column <- normalize_header(args$y_column)
  if (args$x_column == args$y_column) fail("The x and y columns must differ.")
  if ("source_line" %in% c(args$x_column, args$y_column)) {
    fail("source_line is reserved for input provenance; choose a count column.")
  }
  args
}

parse_count <- function(value, column, source_line) {
  # Validate integer-valued decimal/scientific notation BEFORE conversion to
  # R's double type so a fractional count cannot silently round to an integer.
  bad <- function() fail("Line ", source_line, ": ", column,
                         " must be a nonnegative integer; got '", value, "'.")
  text <- trimws(value)
  if (!grepl("^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$", text)) bad()
  pieces <- strsplit(tolower(text), "e", fixed = TRUE)[[1L]]
  exponent <- if (length(pieces) == 2L) as.numeric(pieces[2L]) else 0
  mantissa <- sub("^[+-]", "", pieces[1L])
  decimals <- if (grepl(".", mantissa, fixed = TRUE)) {
    nchar(sub("^[^.]*[.]", "", mantissa))
  } else 0L
  digits <- sub("^0+", "", gsub(".", "", mantissa, fixed = TRUE))
  if (!nzchar(digits)) return(0)
  if (startsWith(text, "-") || !is.finite(exponent)) bad()
  shift <- exponent - decimals
  if (shift < 0) {
    remove <- -shift
    if (remove >= nchar(digits)) bad()
    tail_digits <- substr(digits, nchar(digits) - remove + 1, nchar(digits))
    if (grepl("[^0]", tail_digits)) bad()
    digits <- substr(digits, 1L, nchar(digits) - remove)
  } else {
    if (nchar(digits) + shift > 16) bad()
    digits <- paste0(digits, strrep("0", shift))
  }
  number <- as.numeric(digits)
  if (!is.finite(number) || number > 2^53 - 1) {
    fail("Line ", source_line, ": ", column,
         " exceeds R's exactly representable integer range (2^53 - 1).")
  }
  number
}

read_counts <- function(path, x_column, y_column) {
  if (!file.exists(path) || dir.exists(path)) fail("Input TSV not found: ", path)
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  if (!length(lines)) fail("Input TSV is empty.")
  lines[1L] <- sub("^\ufeff", "", lines[1L])
  # Preserve physical line numbers, including any intervening blank lines.
  nonblank <- which(nzchar(trimws(lines)))
  if (!length(nonblank) || nonblank[1L] != 1L) {
    fail("The first TSV line must contain nonempty headers.")
  }
  split_fields <- function(line) {
    fields <- strsplit(paste0(line, "\t__END__"), "\t", fixed = TRUE)[[1L]]
    fields[-length(fields)]
  }
  headers <- normalize_header(split_fields(lines[1L]))
  if (any(!nzchar(headers)) || anyDuplicated(headers)) {
    fail("Input headers must be nonempty and unique after normalization.")
  }
  missing <- setdiff(c(x_column, y_column), headers)
  if (length(missing)) fail("Missing required TSV column(s): ", paste(missing, collapse = ", "))
  source_lines <- nonblank[-1L]
  if (length(source_lines) < 2L) fail("At least two observations are required.")
  for (line_no in source_lines) {
    if (length(split_fields(lines[line_no])) != length(headers)) {
      fail("Line ", line_no, ": expected ", length(headers), " tab-separated fields.")
    }
  }
  # Zuur Section 1.7, pp. 9-10: read.table() and an explicit data frame.
  # Here every column is initially character so invalid counts are rejected
  # explicitly rather than converted to NA and silently omitted from the fit.
  # Plain TSV convention: tabs delimit fields; quotation marks are literal.
  con <- textConnection(lines[nonblank])
  on.exit(close(con), add = TRUE)
  raw <- read.table(con, header = TRUE, sep = "\t", quote = "", comment.char = "",
                    colClasses = "character", check.names = FALSE,
                    na.strings = character(), row.names = NULL, fill = FALSE)
  names(raw) <- headers
  selected <- unique(c(x_column, y_column, intersect(OPTIONAL_COUNTS, headers)))
  Data <- data.frame(source_line = source_lines)
  for (column in selected) {
    Data[[column]] <- vapply(seq_len(nrow(raw)), function(i) {
      parse_count(raw[[column]][i], column, source_lines[i])
    }, numeric(1))
  }
  check_rows <- function(bad, message) {
    if (any(bad)) fail("Line ", Data$source_line[which(bad)[1L]], ": ", message)
  }
  check_rows(Data[[x_column]] == 0, paste0(x_column, " must be greater than zero."))
  if (x_column == "number_of_species" && "total_eligible_species" %in% names(Data)) {
    check_rows(Data[[x_column]] > Data$total_eligible_species,
               "included species exceed total eligible species.")
    if ("number_of_species_omitted" %in% names(Data)) {
      check_rows(Data[[x_column]] + Data$number_of_species_omitted != Data$total_eligible_species,
                 "included + omitted species must equal total eligible species.")
    }
  }
  if (y_column == "orthologs_retained" && "orthologs_identified_by_compleasm" %in% names(Data)) {
    check_rows(Data[[y_column]] > Data$orthologs_identified_by_compleasm,
               "retained orthologs exceed identified orthologs.")
    if ("orthologs_filtered" %in% names(Data)) {
      check_rows(Data[[y_column]] + Data$orthologs_filtered != Data$orthologs_identified_by_compleasm,
                 "retained + filtered orthologs must equal identified orthologs.")
    }
  }
  if (length(unique(Data[[x_column]])) < 2L) {
    fail("At least two observations with distinct genome/species counts are required.")
  }
  Data <- Data[order(Data[[x_column]], Data$source_line), , drop = FALSE]
  rownames(Data) <- NULL
  Data
}

write_tsv <- function(data, path) {
  write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE,
              na = "NA", fileEncoding = "UTF-8", eol = "\n")
}

open_device <- function(path, format, size_cm, dpi) {
  inches <- size_cm / 2.54
  if (format == "png") {
    if (!capabilities("png")) fail("This R build lacks PNG graphics support.")
    png(path, width = inches, height = inches, units = "in", res = dpi,
        pointsize = axes_font_size, bg = "white")
  } else if (format == "pdf") {
    # Cairo embeds the selected fonts for consistent viewing on other systems.
    if (capabilities("cairo")) {
      cairo_pdf(path, width = inches, height = inches, pointsize = axes_font_size,
                family = "sans", bg = "white")
    } else {
      pdf(path, width = inches, height = inches, pointsize = axes_font_size,
          family = "Helvetica", bg = "white", useDingbats = FALSE)
    }
  } else if (format == "svg") {
    if (!capabilities("cairo")) fail("This R build needs Cairo support for base R SVG output.")
    svg(path, width = inches, height = inches, pointsize = axes_font_size,
        family = "sans", bg = "white")
  }
}

draw_tradeoff <- function(path, format, args, x, y, D1, P1, slope, intercept, r2) {
  open_device(path, format, FIGURE_SIZE_CM, args$dpi)
  on.exit(dev.off(), add = TRUE)
  # R's lwd unit is 1/96 inch; convert the v3 widths from physical points.
  pt_to_lwd <- 96 / 72
  axis_lwd <- AXIS_WIDTH_PT * pt_to_lwd
  marker_cex <- MARKER_DIAMETER_PT / (0.75 * axes_font_size)
  par(mar = c(4.8, 5.4, 1.0, 1.0), mgp = c(3.2, 0.9, 0), las = 1,
      family = "sans", fg = "black", col.axis = "black", col.lab = "black",
      cex.axis = 1, cex.lab = 1, xaxs = "i", yaxs = "i", lend = "butt")
  x_step <- max(1, ceiling(diff(range(x)) / (MAJOR_TICK_COUNT - 1L)))
  x_ticks <- min(x) + seq.int(0L, MAJOR_TICK_COUNT - 1L) * x_step
  y_max <- max(y, P1, 1)
  y_scale <- 10^floor(log10(y_max))
  y_tick_max <- ceiling(y_max / y_scale) * y_scale
  y_ticks <- seq(0, y_tick_max, length.out = MAJOR_TICK_COUNT)
  x_span <- diff(range(x_ticks))
  y_min <- min(0, P1)
  y_span <- y_tick_max - y_min
  # Zuur Appendix A.3.3, pp. 543-546: observations plus predicted line.
  # No extrapolation: the line covers only the observed subset-size range.
  plot(x, y, type = "n", axes = FALSE, ann = FALSE,
       xlim = range(x_ticks) + c(-0.08, 0.08) * x_span,
       ylim = c(y_min - 0.03 * y_span, y_tick_max + 0.12 * y_span))
  abline(h = y_ticks, col = "#D6D6D6", lwd = 0.6)
  lines(D1$genomes, P1, col = "black", lty = "dashed",
        lwd = REGRESSION_WIDTH_PT * pt_to_lwd)
  points(x, y, pch = 19, col = "black", cex = marker_cex)
  if (args$label_points) {
    text(x, y, labels = format(y, big.mark = ",", scientific = FALSE, trim = TRUE),
         pos = 4, offset = 0.6, cex = 0.9)
  }
  axis(1, at = x_ticks, labels = format(x_ticks, scientific = FALSE, trim = TRUE),
       lwd = 0, lwd.ticks = axis_lwd, tcl = -0.4)
  axis(2, at = y_ticks, labels = format(round(y_ticks), big.mark = ",",
                                     scientific = FALSE, trim = TRUE),
       lwd = 0, lwd.ticks = axis_lwd, tcl = -0.4)
  box(bty = "l", lwd = axis_lwd)
  title(xlab = args$x_label, ylab = args$y_label)
  sign <- if (intercept >= 0) "+" else "-"
  annotation <- c(
    sprintf("y = %.2fx %s %.2f", slope, sign, abs(intercept)),
    sprintf("R\u00b2 = %s; n = %d subsets", r2, length(x)),
    sprintf("Slope: %+.2f orthologs per genome included", slope)
  )
  legend("topright", legend = annotation, bty = "o", bg = "white",
         box.col = "#DDDDDD", box.lwd = 0.7,
         cex = descriptive_statistics_font_size / axes_font_size,
         x.intersp = 0, y.intersp = 1.15, inset = 0.015)
  legend("bottomleft", legend = c("Descriptive linear fit", "Evaluated genome subsets"),
         lty = c(2, NA), lwd = c(REGRESSION_WIDTH_PT * pt_to_lwd, NA),
         pch = c(NA, 19), pt.cex = marker_cex, col = "black", bty = "n",
         cex = key_font_size / axes_font_size, inset = 0.015, seg.len = 3.5)
}

draw_diagnostics <- function(path, format, args, M1, ModelData) {
  open_device(path, format, 20, args$dpi)
  on.exit(dev.off(), add = TRUE)
  op <- par(mfrow = c(2, 2), mar = c(4.4, 4.4, 2.2, 1),
            oma = c(3, 0, 1, 0), family = "sans", las = 1)
  on.exit(par(op), add = TRUE, after = FALSE)
  # Zuur Section 2.3.6.1, pp. 23-25, and Appendix A.3.2, pp. 542-543:
  # par(mfrow=...), plot(M1, ...), residual-vs-predictor plots, and qqnorm().
  # Graphs diagnose possible patterns; they do not certify assumptions.
  E <- resid(M1)
  Fit <- fitted(M1)
  numerically_exact <- sum(E^2) <= .Machine$double.eps *
    max(1, sum((ModelData$retained - mean(ModelData$retained))^2))
  if (df.residual(M1) <= 0L || numerically_exact) {
    # A saturated or effectively exact line cannot support standardized
    # residual diagnostics. Show ordinary residuals without dividing by zero.
    plot(Fit, E, pch = 19, xlab = "Fitted orthologs", ylab = "Residuals",
         main = "Residuals vs fitted")
    abline(h = 0, lty = 2)
    plot(ModelData$genomes, E, pch = 19, xlab = args$x_label,
         ylab = "Residuals", main = "Residuals vs subset size")
    abline(h = 0, lty = 2)
    for (i in seq_len(2L)) {
      plot.new()
      text(0.5, 0.5, "Standardized diagnostics unavailable:\nzero residual df or effectively exact fit.", cex = 0.9)
    }
  } else {
    # The book uses this plot.lm() call for residuals vs fitted values.
    plot(M1, which = 1, add.smooth = FALSE, id.n = 0, pch = 19, sub.caption = "")
    E <- rstandard(M1)
    plot(ModelData$genomes, E, pch = 19, xlab = args$x_label,
         ylab = "Standardized residuals", main = "Residuals vs subset size")
    abline(h = 0, lty = 2)
    # Zuur Section 2.3.2, pp. 19-20: inspect residuals rather than a
    # pooled histogram of the unadjusted response to assess normality.
    qqnorm(E, pch = 19, main = "Normal Q-Q: residuals")
    qqline(E, lty = 2)
    # Zuur Section 2.3.7, pp. 26-27, and Appendix A.3.2, p. 543:
    # Cook's distance describes how strongly each row influences the fit.
    # Use a direct plot; do not automatically remove influential subsets.
    Cook <- cooks.distance(M1)
    plot(seq_along(Cook), Cook, type = "h", lwd = 2,
         xlab = "Row in sorted fitted-values table", ylab = "Cook's distance",
         main = "Influence of each subset")
    points(seq_along(Cook), Cook, pch = 19)
  }
  mtext("Exploratory only: few, overlapping subsets cannot validate independence.",
        side = 1, outer = TRUE, line = 1, cex = 0.9)
}

run <- function(args) {
  input <- normalizePath(path.expand(args$input_tsv), mustWork = TRUE)
  Data <- read_counts(input, args$x_column, args$y_column)
  # Zuur Section 2.2, pp. 17-19: Y is retained orthologs; X is subset size.
  # Use explicit data= (as in Appendix A.3), without attach() or globals.
  ModelData <- data.frame(genomes = Data[[args$x_column]], retained = Data[[args$y_column]])
  x <- ModelData$genomes
  y <- ModelData$retained

  # Zuur Appendix A.3, pp. 537-539: lm() fits the ordinary least-squares
  # model. The intercept is included by default. No weights are supplied,
  # so every subset row receives equal weight, including repeated sizes.
  M1 <- lm(retained ~ genomes, data = ModelData, na.action = na.fail)
  if (M1$rank != 2L) fail("Numerically singular fit; check the subset-size range.")
  intercept <- unname(coef(M1)["(Intercept)"])
  slope <- unname(coef(M1)["genomes"])
  Fit <- fitted(M1)
  E <- resid(M1)
  if (any(!is.finite(c(intercept, slope, Fit, E)))) fail("Nonfinite regression result.")

  # summary(M1), demonstrated in Appendix A.3, also calculates inferential
  # statistics internally. Extract ONLY the descriptive R-squared component;
  # do not print summary(M1), its coefficient tests, or its F-statistic.
  # If Y is constant, define R-squared as NA because total variation is zero.
  # Suppress R's 'essentially perfect fit' warning for valid two-point fits.
  r_squared <- if (sum((y - mean(y))^2) > 0) {
    suppressWarnings(summary(M1))$r.squared
  } else NA_real_
  # RMSE uses n, as in v3. It is NOT summary(M1)$sigma, which uses n-2.
  rmse <- sqrt(mean(E^2))
  r2 <- if (is.finite(r_squared)) sprintf("%.4f", r_squared) else "undefined"

  # Zuur Appendix A.3.3, pp. 543-546: predict(M1, newdata=D1) and lines().
  # D1 is sorted and stays within the observed range. No interval= is used.
  D1 <- data.frame(genomes = seq(min(x), max(x), length.out = 200L))
  P1 <- predict(M1, newdata = D1)

  warnings <- character()
  if (nrow(Data) == 2L) warnings <- c(warnings, "Only two observations: the fitted line cannot assess departures from linearity.")
  if (anyDuplicated(x)) warnings <- c(warnings, "Repeated genome/species counts: every input row receives equal weight in the fit.")
  for (column in c("total_eligible_species", "orthologs_identified_by_compleasm")) {
    if (column %in% names(Data) && length(unique(Data[[column]])) > 1L) {
      warnings <- c(warnings, paste0(column, " differs between rows; the fit combines different pools."))
    }
  }
  if (any(diff(x) > 0 & diff(y) > 0)) {
    warnings <- c(warnings, "Retained orthologs increase at some larger subset sizes. Check the subsets and retention criterion; independent selections need not be monotonic.")
  }
  if (!is.finite(r_squared)) warnings <- c(warnings, "Retained ortholog counts are constant; R-squared is undefined.")

  if (abs(slope) < 1e-10) {
    statement <- "The fitted change in retained shared orthologs was approximately zero per additional genome included."
  } else {
    statement <- sprintf(paste0(
      "Across the evaluated subsets (%.0f-%.0f genomes included), the descriptive ",
      "linear fit estimated a %s of approximately %.2f retained shared orthologs ",
      "per additional genome included."
    ), min(x), max(x), if (slope < 0) "decrease" else "increase", abs(slope))
  }
  interpretation <- paste(
    "The slope describes how the number of orthologs satisfying the chosen",
    "shared-across-genomes criterion changes with subset size; it does not",
    "measure biological gene loss or gain. Subsets may overlap or be nested,",
    "so this fit is descriptive and no P values or confidence intervals are reported."
  )
  # Independence: Zuur Section 2.3.5, pp. 21-22. No t-test, F-test,
  # confint(), ANOVA, Shapiro-Wilk test, or automatic model selection here.

  output_dir <- if (is.null(args$output_dir)) dirname(input) else path.expand(args$output_dir)
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    fail("Could not create output directory: ", output_dir)
  }
  prefix <- file.path(normalizePath(output_dir, mustWork = TRUE), args$prefix)
  paths <- character()
  for (format in c("png", "pdf", "svg")) {
    path <- paste0(prefix, ".", format)
    draw_tradeoff(path, format, args, x, y, D1, P1, slope, intercept, r2)
    paths <- c(paths, path)
  }
  num <- function(z) if (is.finite(z)) sprintf("%.12g", z) else "NA"
  summary <- c(
    input_tsv = input, x_column = args$x_column, y_column = args$y_column,
    n_observations = nrow(Data), n_distinct_subset_sizes = length(unique(x)),
    minimum_genomes_included = min(x), maximum_genomes_included = max(x),
    slope = num(slope), slope_units = "retained orthologs per additional genome included",
    intercept = num(intercept), r_squared = num(r_squared), rmse_orthologs = num(rmse),
    fit_type = "descriptive ordinary least squares; equal weight per input row",
    statistical_test = "None; descriptive fit only",
    r_function = "stats::lm(retained ~ genomes, data = ModelData)",
    statement = statement, interpretation = interpretation,
    warnings = if (length(warnings)) paste(warnings, collapse = " | ") else "None",
    reading_fit = "Zuur et al. (2009), Section 2.2 pp. 17-19; Appendix A.3 pp. 536-539",
    reading_independence = "Zuur et al. (2009), Section 2.3.5 pp. 21-22",
    reading_diagnostics = "Zuur et al. (2009), Section 2.3.6.1 pp. 23-25; Appendix A.3.2 pp. 542-543"
  )
  path <- paste0(prefix, "_regression_summary.tsv")
  write_tsv(data.frame(metric = names(summary), value = unname(summary)), path)
  paths <- c(paths, path)
  Data$fitted_orthologs <- unname(Fit)
  Data$residual_orthologs <- unname(E)
  path <- paste0(prefix, "_fitted_values.tsv")
  write_tsv(Data, path)
  paths <- c(paths, path)
  path <- paste0(prefix, "_statement.txt")
  writeLines(c(statement, "", interpretation), path, useBytes = TRUE)
  paths <- c(paths, path)
  if (args$diagnostics) {
    for (format in c("png", "pdf")) {
      path <- paste0(prefix, "_diagnostics.", format)
      draw_diagnostics(path, format, args, M1, ModelData)
      paths <- c(paths, path)
    }
  }
  cat(statement, "\n", sep = "")
  cat(sprintf("R^2 = %s; RMSE = %.2f orthologs.\n", r2, rmse))
  cat("Descriptive fit only; no hypothesis-test P values or confidence intervals.\n")
  for (message in warnings) cat("[WARNING] ", message, "\n", sep = "", file = stderr())
  for (path in paths) cat("[SAVED] ", path, "\n", sep = "")
  invisible(M1)
}

main <- function() {
  tryCatch({
    args <- parse_args(commandArgs(trailingOnly = TRUE))
    if (!is.null(args)) run(args)
    0L
  }, error = function(e) {
    cat("[ERROR] ", conditionMessage(e), "\n", sep = "", file = stderr())
    1L
  })
}

# Running with Rscript executes main(). source() only loads the functions;
# from RStudio, call run(parse_args(c("/path/to/input.tsv", "--diagnostics"))).
if (sys.nframe() == 0L) quit(save = "no", status = main())
