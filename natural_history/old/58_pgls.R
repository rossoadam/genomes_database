#!/usr/bin/env Rscript

###############################################################################
# 58_pgls.R
#
# Purpose:
#   Practical PGLS workflow for testing individual and additive hypotheses for
#   mean GC3 using the 58-species master_output.tsv table.
#
# Main biological model requested:
#   mean_gc3 ~ genome_size + ctmax + mass_preferred + ctmin
#
# Workflow:
#   1. Load data and tree
#   2. Clean numeric variables
#   3. Match species names between data and tree
#   4. Explore missingness and predictor correlations
#   5. Fit individual PGLS models
#   6. Fit the full additive PGLS model when enough complete data exist
#   7. Compare candidate additive models using AICc on the same complete dataset
#   8. Save model summaries and diagnostic plots
#
# Example usage:
#   Rscript 58_pgls.R \
#     --data master_output.tsv \
#     --tree tree_rooted_with_geckos_gc_final_oct_29_treeonly.tre \
#     --outdir pgls_58_results
###############################################################################

suppressPackageStartupMessages({
  library(ape)
  library(nlme)
})

# -----------------------------
# Small argument parser
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) {
    return(args[[idx + 1]])
  }
  return(default)
}

DATA_PATH <- get_arg("--data", "master_output.tsv")
TREE_PATH <- get_arg("--tree", "~/projects/mass_predicts_dna_dynamics/tree_rooted_with_geckos_gc_final_oct_29_treeonly.tre")
OUTDIR <- get_arg("--outdir", "pgls_58_results")
MANIFEST_PATH <- get_arg("--manifest", NULL)
MIN_N_PER_PARAMETER <- as.numeric(get_arg("--min-n-per-parameter", "3"))

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

sink(file.path(OUTDIR, "58_pgls_log.txt"), split = TRUE)
cat("58-species PGLS workflow\n")
cat("========================\n\n")
cat("Data:", DATA_PATH, "\n")
cat("Tree:", TREE_PATH, "\n")
cat("Output directory:", OUTDIR, "\n")
cat("Manifest:", ifelse(is.null(MANIFEST_PATH), "none", MANIFEST_PATH), "\n\n")

# -----------------------------
# Helper functions
# -----------------------------
required_columns <- c(
  "mean_gc3",
  "genome_size",
  "mass_preferred",
  "ctmax",
  "ctmin"
)

species_name_candidates <- c(
  "genus_species",
  "species_normalized",
  "species_key",
  "genus_species_key",
  "normalized_species",
  "normalized_name",
  "organism_name",
  "organismName",
  "Organism Name",
  "species",
  "species_name",
  "binomial",
  "binomial_2020"
)

accession_candidates <- c(
  "accession",
  "accession_id",
  "assembly_accession",
  "assemblyAccession",
  "Assembly Accession"
)

make_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

normalize_species_name <- function(value) {
  if (length(value) == 0 || is.na(value)) return("unknown_species")
  text <- trimws(as.character(value))
  if (text == "") return("unknown_species")

  text <- gsub("_", " ", text)
  text <- gsub("\\([^)]*\\)", " ", text)
  text <- gsub("[^A-Za-z×x[:space:].-]", " ", text)
  text <- gsub("\\s+", " ", text)
  text <- trimws(text)

  tokens <- unlist(strsplit(text, "\\s+"))
  tokens <- tolower(gsub("^\\.+|\\.+$", "", tokens))
  bad_tokens <- c("sp", "spp", "cf", "aff", "nr", "sp.", "spp.", "cf.", "aff.", "nr.")
  hybrid_tokens <- c("x", "×")
  tokens <- tokens[tokens != "" & !(tokens %in% bad_tokens) & !(tokens %in% hybrid_tokens)]

  if (length(tokens) >= 2) return(paste(tokens[1], tokens[2], sep = "_"))
  if (length(tokens) == 1) return(tokens[1])
  return("unknown_species")
}

extract_accession <- function(value) {
  if (length(value) == 0 || is.na(value)) return(NA_character_)
  text <- as.character(value)
  hit <- regmatches(text, regexpr("GC[AF]_[0-9]+\\.[0-9]+", text, perl = TRUE))
  if (length(hit) == 0 || hit == "") return(NA_character_)
  return(hit)
}

choose_first_column <- function(df, candidates) {
  lowered <- tolower(trimws(colnames(df)))
  for (candidate in candidates) {
    idx <- match(tolower(trimws(candidate)), lowered)
    if (!is.na(idx)) return(colnames(df)[idx])
  }
  return(NULL)
}

load_manifest_accession_map <- function(path) {
  if (is.null(path) || is.na(path) || path == "") {
    return(data.frame(accession = character(), species_from_manifest = character(), stringsAsFactors = FALSE))
  }

  manifest_path <- path.expand(path)
  if (!file.exists(manifest_path)) stop("Manifest file not found: ", manifest_path)

  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  acc_col <- choose_first_column(manifest, accession_candidates)
  species_col <- choose_first_column(manifest, c("species_key", "species_normalized", "genus_species", "organism_name", "organismName", "species", "species_name"))

  if (is.null(acc_col)) stop("Could not find an accession column in manifest. Columns: ", paste(colnames(manifest), collapse = ", "))
  if (is.null(species_col)) stop("Could not find a species column in manifest. Columns: ", paste(colnames(manifest), collapse = ", "))

  out <- data.frame(
    accession = vapply(manifest[[acc_col]], extract_accession, character(1)),
    species_from_manifest = vapply(manifest[[species_col]], normalize_species_name, character(1)),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$accession) & out$accession != "" &
             !is.na(out$species_from_manifest) & out$species_from_manifest != "unknown_species", , drop = FALSE]
  out <- out[!duplicated(out$accession), , drop = FALSE]
  return(out)
}

map_label_to_species <- function(label, manifest_map = NULL) {
  acc <- extract_accession(label)
  if (!is.na(acc) && !is.null(manifest_map) && nrow(manifest_map) > 0 && acc %in% manifest_map$accession) {
    return(manifest_map$species_from_manifest[match(acc, manifest_map$accession)])
  }
  return(normalize_species_name(label))
}

safe_aicc <- function(model) {
  # AICc = AIC + [2k(k+1)] / (n-k-1)
  n <- length(resid(model, type = "normalized"))
  k <- length(coef(model)) + 1  # fixed effects + residual variance parameter approximation
  if (is.na(n) || is.na(k) || (n - k - 1) <= 0) return(NA_real_)
  AIC(model) + (2 * k * (k + 1)) / (n - k - 1)
}

write_model_summary <- function(model, path) {
  con <- file(path, open = "wt")
  sink(con)
  print(summary(model))
  cat("\nANOVA table:\n")
  print(anova(model))
  cat("\nAIC:", AIC(model), "\n")
  cat("AICc:", safe_aicc(model), "\n")
  sink()
  close(con)
}

plot_pgls_diagnostics <- function(model, prefix) {
  pdf(file.path(OUTDIR, paste0(prefix, "_diagnostics.pdf")), width = 7, height = 7)
  par(mfrow = c(2, 2))
  fitted_vals <- fitted(model)
  norm_resid <- resid(model, type = "normalized")
  raw_resid <- resid(model, type = "response")

  plot(fitted_vals, norm_resid,
       xlab = "Fitted values", ylab = "Normalized residuals",
       main = paste(prefix, "residuals vs fitted"))
  abline(h = 0, lty = 2)

  qqnorm(norm_resid, main = paste(prefix, "normal Q-Q"))
  qqline(norm_resid)

  hist(norm_resid, main = paste(prefix, "normalized residuals"),
       xlab = "Normalized residuals")

  plot(fitted_vals, raw_resid,
       xlab = "Fitted values", ylab = "Raw residuals",
       main = paste(prefix, "raw residuals vs fitted"))
  abline(h = 0, lty = 2)
  dev.off()
}

fit_pgls <- function(formula_obj, data_obj, tree_obj, model_name, method = "ML") {
  vars <- all.vars(formula_obj)
  complete_data <- data_obj[complete.cases(data_obj[, vars, drop = FALSE]), , drop = FALSE]
  n <- nrow(complete_data)
  p <- length(vars) - 1
  min_needed <- max(p + 2, p * MIN_N_PER_PARAMETER)

  cat("\nModel:", model_name, "\n")
  cat("Formula:", deparse(formula_obj), "\n")
  cat("Complete rows:", n, "\n")
  cat("Predictors:", p, "\n")
  cat("Suggested minimum rows:", min_needed, "\n")

  if (n < min_needed) {
    cat("SKIPPED:", model_name, "has too few complete rows for reliable PGLS.\n")
    return(NULL)
  }

  tree_sub <- drop.tip(tree_obj, setdiff(tree_obj$tip.label, complete_data$genus_species))
  complete_data <- complete_data[tree_sub$tip.label, , drop = FALSE]

  cor_struct <- corBrownian(phy = tree_sub, form = ~ genus_species)
  model <- gls(formula_obj, data = complete_data, correlation = cor_struct, method = method)

  write_model_summary(model, file.path(OUTDIR, paste0(model_name, "_summary.txt")))
  plot_pgls_diagnostics(model, model_name)

  return(list(model = model, data = complete_data, tree = tree_sub, n = n, p = p))
}

# -----------------------------
# Load data and tree
# -----------------------------
if (!file.exists(DATA_PATH)) stop("Data file not found: ", DATA_PATH)
if (!file.exists(path.expand(TREE_PATH))) stop("Tree file not found: ", TREE_PATH)

data <- read.delim(DATA_PATH, stringsAsFactors = FALSE, check.names = FALSE)
tree <- read.tree(path.expand(TREE_PATH))

missing_cols <- setdiff(required_columns, colnames(data))
if (length(missing_cols) > 0) {
  stop("Missing required numeric/model columns in data: ", paste(missing_cols, collapse = ", "))
}

manifest_map <- load_manifest_accession_map(MANIFEST_PATH)
cat("Manifest accession mappings loaded:", nrow(manifest_map), "\n")

species_col <- choose_first_column(data, species_name_candidates)
accession_col <- choose_first_column(data, accession_candidates)

if (is.null(species_col) && is.null(accession_col)) {
  stop("Could not find a species-name or accession column in data. Columns: ", paste(colnames(data), collapse = ", "))
}

data$source_species_label <- if (!is.null(species_col)) as.character(data[[species_col]]) else NA_character_
data$source_accession_label <- if (!is.null(accession_col)) as.character(data[[accession_col]]) else NA_character_

# analysis_species is the single matching key used by both data and tree.
# If a row has an accession and that accession exists in the manifest, the manifest species is used.
# Otherwise, the best available species-name column is normalized.
data$genus_species <- mapply(
  function(species_label, accession_label) {
    acc <- extract_accession(accession_label)
    if (!is.na(acc) && nrow(manifest_map) > 0 && acc %in% manifest_map$accession) {
      return(manifest_map$species_from_manifest[match(acc, manifest_map$accession)])
    }
    if (!is.na(species_label) && species_label != "") {
      return(normalize_species_name(species_label))
    }
    return(map_label_to_species(accession_label, manifest_map))
  },
  data$source_species_label,
  data$source_accession_label,
  USE.NAMES = FALSE
)

# Standardize and clean key variables
data$genus_species <- as.character(data$genus_species)
data$mean_gc3 <- make_numeric(data$mean_gc3)
data$genome_size <- make_numeric(data$genome_size)
data$mass_preferred <- make_numeric(data$mass_preferred)
data$ctmax <- make_numeric(data$ctmax)
data$ctmin <- make_numeric(data$ctmin)

# Optional transformed columns for sensitivity analyses.
# The requested primary model below uses raw mass_preferred and genome_size.
data$log10_mass_preferred <- ifelse(data$mass_preferred > 0, log10(data$mass_preferred), NA_real_)
data$log10_genome_size <- ifelse(data$genome_size > 0, log10(data$genome_size), NA_real_)

# Keep only rows with species names and GC3
data <- data[!is.na(data$genus_species) & data$genus_species != "" & !is.na(data$mean_gc3), ]
rownames(data) <- data$genus_species

# -----------------------------
# Match tree and data
# -----------------------------
cat("Original data rows:", nrow(data), "\n")
cat("Original tree tips:", length(tree$tip.label), "\n\n")

tree_mapping <- data.frame(
  original_tree_tip = tree$tip.label,
  tree_accession = vapply(tree$tip.label, extract_accession, character(1)),
  mapped_species = vapply(tree$tip.label, map_label_to_species, character(1), manifest_map = manifest_map),
  stringsAsFactors = FALSE
)
tree_mapping$mapping_source <- ifelse(
  !is.na(tree_mapping$tree_accession) & tree_mapping$tree_accession %in% manifest_map$accession,
  "manifest_accession",
  "normalized_tree_tip"
)

data_mapping <- data.frame(
  row_number = seq_len(nrow(data)),
  source_species_label = data$source_species_label,
  source_accession_label = data$source_accession_label,
  mapped_species = data$genus_species,
  stringsAsFactors = FALSE
)

write.csv(tree_mapping, file.path(OUTDIR, "tree_tip_to_analysis_species_mapping.csv"), row.names = FALSE)
write.csv(data_mapping, file.path(OUTDIR, "data_to_analysis_species_mapping.csv"), row.names = FALSE)

unmapped_tree <- tree_mapping$mapped_species == "unknown_species" | is.na(tree_mapping$mapped_species) | tree_mapping$mapped_species == ""
if (any(unmapped_tree)) {
  write.csv(tree_mapping[unmapped_tree, , drop = FALSE], file.path(OUTDIR, "unmapped_tree_tips.csv"), row.names = FALSE)
}

# Drop tree tips that collapse to duplicated mapped species. This avoids invalid duplicated
# tip labels after accession -> species conversion. The first occurrence is retained.
duplicated_tree_species <- duplicated(tree_mapping$mapped_species) & !unmapped_tree
if (any(duplicated_tree_species)) {
  write.csv(tree_mapping[duplicated_tree_species, , drop = FALSE], file.path(OUTDIR, "duplicate_tree_mapped_species_dropped.csv"), row.names = FALSE)
  tree <- drop.tip(tree, tree_mapping$original_tree_tip[duplicated_tree_species])
  tree_mapping <- tree_mapping[!duplicated_tree_species, , drop = FALSE]
}

# Rename tree tips from accessions/raw labels to the same species key used in data.
tree$tip.label <- tree_mapping$mapped_species[match(tree$tip.label, tree_mapping$original_tree_tip)]

# Drop duplicated data species, if present, keeping the first row and writing a report.
duplicated_data_species <- duplicated(data$genus_species)
if (any(duplicated_data_species)) {
  write.csv(data[duplicated_data_species, , drop = FALSE], file.path(OUTDIR, "duplicate_data_mapped_species_dropped.csv"), row.names = FALSE)
  data <- data[!duplicated_data_species, , drop = FALSE]
}
rownames(data) <- data$genus_species

species_in_both <- intersect(tree$tip.label, data$genus_species)
species_only_tree <- setdiff(tree$tip.label, data$genus_species)
species_only_data <- setdiff(data$genus_species, tree$tip.label)

cat("Species shared by tree and data:", length(species_in_both), "\n")
cat("Species only in tree:", length(species_only_tree), "\n")
cat("Species only in data:", length(species_only_data), "\n\n")

writeLines(species_only_tree, file.path(OUTDIR, "species_only_in_tree.txt"))
writeLines(species_only_data, file.path(OUTDIR, "species_only_in_data.txt"))

if (length(species_in_both) < 3) {
  stop("Fewer than 3 species overlap between data and tree. Check species names. See tree_tip_to_analysis_species_mapping.csv and data_to_analysis_species_mapping.csv.")
}

tree <- drop.tip(tree, setdiff(tree$tip.label, species_in_both))
data <- data[tree$tip.label, , drop = FALSE]

# -----------------------------
# Missingness and exploration
# -----------------------------
model_vars <- c("mean_gc3", "genome_size", "mass_preferred", "ctmax", "ctmin")
missingness <- data.frame(
  variable = model_vars,
  n_missing = sapply(data[, model_vars, drop = FALSE], function(x) sum(is.na(x))),
  n_present = sapply(data[, model_vars, drop = FALSE], function(x) sum(!is.na(x)))
)
write.csv(missingness, file.path(OUTDIR, "missingness_summary.csv"), row.names = FALSE)
print(missingness)

complete_full <- data[complete.cases(data[, model_vars, drop = FALSE]), model_vars, drop = FALSE]
cat("\nComplete rows for full requested additive model:", nrow(complete_full), "\n")

# Correlations among predictors using pairwise complete observations
predictor_vars <- c("genome_size", "mass_preferred", "ctmax", "ctmin")
cor_mat <- cor(data[, predictor_vars, drop = FALSE], use = "pairwise.complete.obs", method = "pearson")
write.csv(cor_mat, file.path(OUTDIR, "predictor_pairwise_correlations.csv"))
cat("\nPairwise predictor correlations:\n")
print(cor_mat)

pdf(file.path(OUTDIR, "exploratory_pairs_plots.pdf"), width = 9, height = 9)
pairs(data[, model_vars, drop = FALSE], main = "PGLS variables: pairwise exploration")
dev.off()

# -----------------------------
# Individual PGLS hypotheses
# -----------------------------
individual_formulas <- list(
  pgls_mass = mean_gc3 ~ mass_preferred,
  pgls_genome_size = mean_gc3 ~ genome_size,
  pgls_ctmax = mean_gc3 ~ ctmax,
  pgls_ctmin = mean_gc3 ~ ctmin
)

individual_results <- list()
for (nm in names(individual_formulas)) {
  individual_results[[nm]] <- fit_pgls(individual_formulas[[nm]], data, tree, nm, method = "ML")
}

# -----------------------------
# Full additive PGLS hypothesis
# -----------------------------
full_formula <- mean_gc3 ~ genome_size + ctmax + mass_preferred + ctmin
full_result <- fit_pgls(full_formula, data, tree, "pgls_full_additive", method = "ML")

# -----------------------------
# Candidate additive model comparison on same complete dataset
# -----------------------------
# This section is intentionally restricted to rows complete for all predictors so
# AICc values are comparable across models. If too few rows exist, the script
# records that model comparison was skipped.
cat("\nCandidate additive model comparison on full-complete dataset\n")
cat("----------------------------------------------------------\n")

candidate_formulas <- list(
  intercept_only = mean_gc3 ~ 1,
  mass = mean_gc3 ~ mass_preferred,
  genome_size = mean_gc3 ~ genome_size,
  ctmax = mean_gc3 ~ ctmax,
  ctmin = mean_gc3 ~ ctmin,
  mass_genome_size = mean_gc3 ~ mass_preferred + genome_size,
  mass_ctmax = mean_gc3 ~ mass_preferred + ctmax,
  mass_ctmin = mean_gc3 ~ mass_preferred + ctmin,
  genome_size_ctmax = mean_gc3 ~ genome_size + ctmax,
  genome_size_ctmin = mean_gc3 ~ genome_size + ctmin,
  ctmax_ctmin = mean_gc3 ~ ctmax + ctmin,
  full_additive = mean_gc3 ~ genome_size + ctmax + mass_preferred + ctmin
)

complete_species <- rownames(data[complete.cases(data[, model_vars, drop = FALSE]), , drop = FALSE])
comparison_rows <- list()

if (length(complete_species) >= 8) {
  comp_tree <- drop.tip(tree, setdiff(tree$tip.label, complete_species))
  comp_data <- data[comp_tree$tip.label, , drop = FALSE]
  comp_cor <- corBrownian(phy = comp_tree, form = ~ genus_species)

  for (nm in names(candidate_formulas)) {
    f <- candidate_formulas[[nm]]
    p <- length(all.vars(f)) - 1
    if (nrow(comp_data) < max(p + 2, p * MIN_N_PER_PARAMETER)) {
      comparison_rows[[nm]] <- data.frame(
        model = nm,
        formula = deparse(f),
        n = nrow(comp_data),
        predictors = p,
        AIC = NA_real_,
        AICc = NA_real_,
        delta_AICc = NA_real_,
        weight_AICc = NA_real_,
        note = "too_few_rows"
      )
      next
    }
    m <- gls(f, data = comp_data, correlation = comp_cor, method = "ML")
    comparison_rows[[nm]] <- data.frame(
      model = nm,
      formula = deparse(f),
      n = nrow(comp_data),
      predictors = p,
      AIC = AIC(m),
      AICc = safe_aicc(m),
      delta_AICc = NA_real_,
      weight_AICc = NA_real_,
      note = "fit"
    )
  }

  model_comparison <- do.call(rbind, comparison_rows)
  if (any(!is.na(model_comparison$AICc))) {
    min_aicc <- min(model_comparison$AICc, na.rm = TRUE)
    model_comparison$delta_AICc <- model_comparison$AICc - min_aicc
    weights <- exp(-0.5 * model_comparison$delta_AICc)
    model_comparison$weight_AICc <- weights / sum(weights, na.rm = TRUE)
    model_comparison <- model_comparison[order(model_comparison$AICc), ]
  }
} else {
  model_comparison <- data.frame(
    model = names(candidate_formulas),
    formula = sapply(candidate_formulas, deparse),
    n = length(complete_species),
    predictors = sapply(candidate_formulas, function(f) length(all.vars(f)) - 1),
    AIC = NA_real_,
    AICc = NA_real_,
    delta_AICc = NA_real_,
    weight_AICc = NA_real_,
    note = "skipped_too_few_complete_rows_for_full_predictor_set"
  )
  cat("SKIPPED candidate model comparison: only", length(complete_species),
      "species are complete for all predictors.\n")
}

write.csv(model_comparison, file.path(OUTDIR, "candidate_model_comparison_AICc.csv"), row.names = FALSE)
print(model_comparison)

# -----------------------------
# Optional sensitivity model using log10 body mass
# -----------------------------
# Body mass is often log-scaled in comparative analyses. This is not the primary
# requested model, but it is useful to examine as a sensitivity check.
log_mass_formula <- mean_gc3 ~ genome_size + ctmax + log10_mass_preferred + ctmin
log_mass_result <- fit_pgls(log_mass_formula, data, tree, "pgls_full_additive_log10_mass_sensitivity", method = "ML")

cat("\nDone. Results written to:", OUTDIR, "\n")
sink()
