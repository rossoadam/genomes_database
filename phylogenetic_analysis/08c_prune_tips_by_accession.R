#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(optparse)
  library(readr)
})

option_list <- list(
  make_option(c("-t", "--tree"), type = "character", help = "Path to input tree file in Newick format."),
  make_option(c("-m", "--manifest"), type = "character", help = "Path to manifest CSV."),
  make_option(c("-a", "--accessions"), type = "character", help = "Path to text file with one accession per line."),
  make_option(c("-o", "--output"), type = "character", help = "Path to output pruned tree file in Newick format."),
  make_option(c("--manifest_accession_col"), type = "character", default = "accession", help = "Manifest column containing accession IDs [default %default]."),
  make_option(c("--manifest_tip_col"), type = "character", default = "species_key", help = "Manifest column containing species names or expected tree tip labels [default %default]."),
  make_option(c("--drop_unmatched"), action = "store_true", default = FALSE,
              help = "If set, silently ignore accessions that do not map to any manifest row or mapped tip not present in the tree. By default the script reports them."),
  make_option(c("--report"), type = "character", default = NULL,
              help = "Optional path to write a pruning report TSV."),
  make_option(c("--allow_substring_match"), action = "store_true", default = TRUE,
              help = "Allow normalized species names to match tree tips by substring when exact normalized matching fails [default %default]."),
  make_option(c("--no_substring_match"), action = "store_false", dest = "allow_substring_match",
              help = "Disable substring matching and require exact normalized matches only.")
)

parser <- OptionParser(
  usage = "%prog --tree tree.nwk --manifest manifest.csv --accessions drop_accessions.txt --output pruned_tree.nwk [options]",
  option_list = option_list
)
opt <- parse_args(parser)

required_args <- c("tree", "manifest", "accessions", "output")
missing_args <- required_args[vapply(required_args, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing_args) > 0) {
  print_help(parser)
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "), call. = FALSE)
}

read_accession_list <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  lines <- lines[!grepl("^#", lines)]
  unique(lines)
}

normalize_label <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("['\"]", "", x)
  x <- gsub("\\[.*?\\]|\\(.*?\\)", "", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

species_token_key <- function(x) {
  x <- normalize_label(x)
  pieces <- unlist(strsplit(x, "_+"))
  pieces <- pieces[pieces != ""]
  if (length(pieces) >= 2) {
    paste(pieces[1:2], collapse = "_")
  } else if (length(pieces) == 1) {
    pieces[1]
  } else {
    ""
  }
}

if (!file.exists(opt$tree)) {
  stop("Tree file not found: ", opt$tree, call. = FALSE)
}
if (!file.exists(opt$manifest)) {
  stop("Manifest file not found: ", opt$manifest, call. = FALSE)
}
if (!file.exists(opt$accessions)) {
  stop("Accession file not found: ", opt$accessions, call. = FALSE)
}

tree <- read.tree(opt$tree)
manifest <- read_csv(opt$manifest, show_col_types = FALSE)
accessions_to_drop <- read_accession_list(opt$accessions)

if (!(opt$manifest_accession_col %in% names(manifest))) {
  stop(
    "Manifest is missing accession column '", opt$manifest_accession_col,
    "'. Available columns: ", paste(names(manifest), collapse = ", "),
    call. = FALSE
  )
}
if (!(opt$manifest_tip_col %in% names(manifest))) {
  stop(
    "Manifest is missing tip/species column '", opt$manifest_tip_col,
    "'. Available columns: ", paste(names(manifest), collapse = ", "),
    call. = FALSE
  )
}

manifest[[opt$manifest_accession_col]] <- trimws(as.character(manifest[[opt$manifest_accession_col]]))
manifest[[opt$manifest_tip_col]] <- trimws(as.character(manifest[[opt$manifest_tip_col]]))
manifest <- manifest[manifest[[opt$manifest_accession_col]] != "" & manifest[[opt$manifest_tip_col]] != "", , drop = FALSE]

matched_manifest <- manifest[manifest[[opt$manifest_accession_col]] %in% accessions_to_drop, , drop = FALSE]
missing_accessions <- setdiff(accessions_to_drop, unique(matched_manifest[[opt$manifest_accession_col]]))

mapping <- unique(matched_manifest[, c(opt$manifest_accession_col, opt$manifest_tip_col), drop = FALSE])
names(mapping) <- c("accession", "manifest_tip")
mapping$manifest_tip_normalized <- normalize_label(mapping$manifest_tip)
mapping$manifest_species_key <- vapply(mapping$manifest_tip, species_token_key, character(1))

tree_lookup <- data.frame(
  tree_tip = tree$tip.label,
  tree_tip_normalized = normalize_label(tree$tip.label),
  tree_species_key = vapply(tree$tip.label, species_token_key, character(1)),
  stringsAsFactors = FALSE
)

resolved_tip <- rep(NA_character_, nrow(mapping))
match_method <- rep(NA_character_, nrow(mapping))
match_note <- rep(NA_character_, nrow(mapping))

for (i in seq_len(nrow(mapping))) {
  target_norm <- mapping$manifest_tip_normalized[i]
  target_key <- mapping$manifest_species_key[i]

  exact_idx <- which(tree_lookup$tree_tip_normalized == target_norm)
  if (length(exact_idx) == 1) {
    resolved_tip[i] <- tree_lookup$tree_tip[exact_idx]
    match_method[i] <- "exact_normalized"
    next
  }
  if (length(exact_idx) > 1) {
    resolved_tip[i] <- tree_lookup$tree_tip[exact_idx[1]]
    match_method[i] <- "exact_normalized_ambiguous"
    match_note[i] <- paste(tree_lookup$tree_tip[exact_idx], collapse = ";")
    next
  }

  key_idx <- which(tree_lookup$tree_species_key == target_key & target_key != "")
  if (length(key_idx) == 1) {
    resolved_tip[i] <- tree_lookup$tree_tip[key_idx]
    match_method[i] <- "species_key"
    next
  }
  if (length(key_idx) > 1) {
    resolved_tip[i] <- tree_lookup$tree_tip[key_idx[1]]
    match_method[i] <- "species_key_ambiguous"
    match_note[i] <- paste(tree_lookup$tree_tip[key_idx], collapse = ";")
    next
  }

  if (isTRUE(opt$allow_substring_match) && target_norm != "") {
    sub_idx <- which(
      grepl(target_norm, tree_lookup$tree_tip_normalized, fixed = TRUE) |
      grepl(tree_lookup$tree_tip_normalized, target_norm, fixed = TRUE)
    )
    if (length(sub_idx) == 1) {
      resolved_tip[i] <- tree_lookup$tree_tip[sub_idx]
      match_method[i] <- "substring_normalized"
      next
    }
    if (length(sub_idx) > 1) {
      resolved_tip[i] <- tree_lookup$tree_tip[sub_idx[1]]
      match_method[i] <- "substring_normalized_ambiguous"
      match_note[i] <- paste(tree_lookup$tree_tip[sub_idx], collapse = ";")
      next
    }
  }

  match_method[i] <- "unmatched"
}

mapping$tree_tip <- resolved_tip
mapping$match_method <- match_method
mapping$match_note <- match_note
mapping$tip_in_tree <- !is.na(mapping$tree_tip) & mapping$tree_tip %in% tree$tip.label

resolved_unique_tips <- unique(na.omit(mapping$tree_tip))
missing_manifest_tips <- unique(mapping$manifest_tip[is.na(mapping$tree_tip)])

cat("Input tree tips:", length(tree$tip.label), "\n")
cat("Accessions requested for pruning:", length(accessions_to_drop), "\n")
cat("Accessions matched in manifest:", length(unique(mapping$accession)), "\n")
cat("Tips matched in tree:", length(resolved_unique_tips), "\n")

method_counts <- sort(table(mapping$match_method), decreasing = TRUE)
if (length(method_counts) > 0) {
  cat("\nMatch methods used:\n")
  for (nm in names(method_counts)) {
    cat("-", nm, ":", method_counts[[nm]], "\n")
  }
}

if (!opt$drop_unmatched) {
  if (length(missing_accessions) > 0) {
    cat("\nAccessions not found in manifest:\n")
    cat(paste(missing_accessions, collapse = "\n"), "\n")
  }
  if (length(missing_manifest_tips) > 0) {
    cat("\nManifest species/tips that could not be matched to any tree tip:\n")
    cat(paste(missing_manifest_tips, collapse = "\n"), "\n")
  }

  ambiguous_mask <- !is.na(mapping$match_method) & grepl("ambiguous", mapping$match_method, fixed = TRUE)
  ambiguous_rows <- mapping[ambiguous_mask, , drop = FALSE]
  if (nrow(ambiguous_rows) > 0) {
    cat("\nAmbiguous matches (first candidate used):\n")
    amb_out <- paste0(
      ambiguous_rows$manifest_tip,
      " -> ", ambiguous_rows$tree_tip,
      " | candidates=", ambiguous_rows$match_note
    )
    cat(paste(amb_out, collapse = "\n"), "\n")
  }
}

if (length(resolved_unique_tips) == 0) {
  cat("\nFirst 20 tree tips for debugging:\n")
  cat(paste(head(tree$tip.label, 20), collapse = "\n"), "\n")
  stop("No matching tips were found in the tree, so nothing was pruned.", call. = FALSE)
}

pruned_tree <- drop.tip(tree, resolved_unique_tips)
write.tree(pruned_tree, file = opt$output)

cat("\nPruned", length(resolved_unique_tips), "tips.\n")
cat("Output tree written to:", opt$output, "\n")
cat("Remaining tips:", length(pruned_tree$tip.label), "\n")

if (!is.null(opt$report)) {
  write_tsv(mapping, opt$report)
  cat("Pruning report written to:", opt$report, "\n")
}
