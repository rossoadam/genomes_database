#!/usr/bin/env Rscript
# 08a_root_a_tree_v6.R
#
# Required positional arguments:
#   Rscript 08a_root_a_tree_v6.R \
#     <genomes_dir> \
#     <outgroup_list> \
#     <newick_file> \
#     <project_manifest.csv>
#
# Optional arguments (prompted interactively when omitted):
#   --promoted-tree-dir PATH
#   --tree-files-tsv PATH
#   --nhphyml-run-pk INTEGER

suppressPackageStartupMessages({
  library(phytools)
  library(ape)
})

# ── 1. Helpers ────────────────────────────────────────────────────────────────

print_usage <- function() {
  cat(
    paste0(
      "Usage:\n",
      "  Rscript 08a_root_a_tree_v6.R <genomes_dir> <outgroup_list> <newick_file> <project_manifest.csv> \\\n",
      "      [--promoted-tree-dir PATH] \\\n",
      "      [--tree-files-tsv PATH] \\\n",
      "      [--nhphyml-run-pk INTEGER]\n\n",
      "If an optional value is omitted, the script prompts for it interactively.\n"
    )
  )
}

clean_path <- function(path) {
  path <- trimws(path)
  if (!nzchar(path)) return(path)
  path <- path.expand(path)
  path <- gsub("\\\\", "/", path)
  path <- gsub("/{2,}", "/", path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

prompt_if_missing <- function(value, prompt_text) {
  if (is.null(value) || !nzchar(trimws(value))) {
    value <- readline(prompt = prompt_text)
  }
  value <- trimws(value)
  if (!nzchar(value)) {
    stop(sprintf("No value supplied for: %s", sub(":\\s*$", "", prompt_text)))
  }
  value
}

format_genus_species <- function(label) {
  normalized <- gsub("\\s+", "_", trimws(label))
  parts <- strsplit(normalized, "_", fixed = TRUE)[[1]]
  parts <- parts[nchar(parts) > 0]

  if (length(parts) < 2) {
    stop(sprintf(
      "Cannot convert organism name '%s' to Genus_species: fewer than two name components.",
      label
    ))
  }

  genus <- parts[1]
  species <- parts[2]

  genus <- paste0(
    toupper(substr(genus, 1, 1)),
    tolower(substr(genus, 2, nchar(genus)))
  )
  species <- tolower(species)

  paste(genus, species, sep = "_")
}

rename_tips_full <- function(tree, label_map) {
  current_tips <- tree$tip.label
  tree$tip.label <- unname(label_map[current_tips])
  cat(sprintf("Renamed %d tips using manifest organism_name values\n",
              length(current_tips)))
  tree
}

rename_tips_genus_species <- function(tree, label_map) {
  current_tips <- tree$tip.label
  full_names <- unname(label_map[current_tips])
  formatted <- vapply(full_names, format_genus_species, character(1), USE.NAMES = FALSE)

  duplicate_labels <- unique(formatted[duplicated(formatted)])
  if (length(duplicate_labels) > 0) {
    cat("\nERROR: Genus_species formatting produced duplicate tip labels:\n")
    cat(paste(" -", duplicate_labels), sep = "\n")
    cat("\nThese likely represent multiple accessions/subspecies that collapse to the same species label.\n")
    stop("Refusing to write ambiguous species-renamed trees.")
  }

  tree$tip.label <- formatted
  cat(sprintf("Formatted %d tips as Genus_species\n", length(formatted)))
  tree
}

plot_tree_pdf <- function(tree, pdf_path, title_str) {
  n_tips <- length(tree$tip.label)
  fig_h <- max(8, n_tips * 0.18)
  fnt_size <- min(0.7, 25 / n_tips)

  pdf(pdf_path, width = 12, height = fig_h)
  plotTree(
    tree,
    fsize = fnt_size,
    lwd = 0.8,
    mar = c(1, 0, 2, 1),
    direction = "rightwards"
  )
  title(main = title_str, cex.main = 0.9)
  dev.off()

  cat(sprintf("PDF written: %s\n", pdf_path))
}

read_tree_text_one_line <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Tree file does not exist: %s", path))
  }
  lines <- readLines(path, warn = FALSE)
  # Concatenate Newick records without changing the spelling of any taxa.
  # In particular, do NOT use the regex "[\\r\\n\\t]": depending on
  # interpretation it can delete literal r, n, and t from taxon names.
  text <- paste(trimws(lines), collapse = "")
  if (!nzchar(text)) {
    stop(sprintf("Tree file is empty: %s", path))
  }
  # A TSV field must be one physical line. Never silently mutate tree labels.
  if (any(vapply(c("\r", "\n", "\t"), function(char) {
    grepl(char, text, fixed = TRUE)
  }, logical(1)))) {
    stop(sprintf("Tree contains an embedded tab or line break: %s", path))
  }
  text
}

sanitize_tsv_field <- function(x) {
  x <- as.character(x)
  # Match ACTUAL control characters as fixed strings, not regex character
  # classes that can accidentally remove ordinary letters r, n, and t.
  x <- gsub("\r", " ", x, fixed = TRUE)
  x <- gsub("\n", " ", x, fixed = TRUE)
  x <- gsub("\t", " ", x, fixed = TRUE)
  x
}

# Catch any future accidental loss of letters before writing user data.
stopifnot(
  identical(sanitize_tsv_field("unrooted"), "unrooted"),
  identical(sanitize_tsv_field("rooted"), "rooted"),
  identical(sanitize_tsv_field("substitution_per_site"), "substitution_per_site"),
  identical(sanitize_tsv_field("Crotalus_viridis"), "Crotalus_viridis"),
  identical(sanitize_tsv_field("Sphenodon_punctatus"), "Sphenodon_punctatus"),
  identical(sanitize_tsv_field("A\tB\r\nC"), "A B  C")
)

make_tsv_line <- function(record, header_fields) {
  values <- vapply(
    header_fields,
    function(field) {
      if (field %in% names(record)) sanitize_tsv_field(record[[field]]) else ""
    },
    character(1)
  )
  paste(values, collapse = "\t")
}

read_tree_pk_values <- function(tsv_path, header_fields) {
  lines <- readLines(tsv_path, warn = FALSE)
  if (length(lines) <= 1) return(integer(0))

  tree_pk_idx <- match("tree_pk", header_fields)
  data_lines <- lines[-1]
  data_lines <- data_lines[nzchar(trimws(data_lines))]
  if (length(data_lines) == 0) return(integer(0))

  values <- character(length(data_lines))
  for (i in seq_along(data_lines)) {
    fields <- strsplit(data_lines[i], "\t", fixed = TRUE)[[1]]
    if (length(fields) < tree_pk_idx) {
      stop(sprintf("Malformed tree_files.tsv row %d: missing tree_pk field.", i + 1))
    }
    values[i] <- trimws(fields[tree_pk_idx])
  }

  invalid <- !grepl("^[0-9]+$", values)
  if (any(invalid)) {
    stop(sprintf(
      "tree_files.tsv contains invalid/non-integer tree_pk value(s): %s",
      paste(unique(values[invalid]), collapse = ", ")
    ))
  }

  as.integer(values)
}

copy_or_stop <- function(from, to, overwrite = TRUE) {
  ok <- file.copy(from, to, overwrite = overwrite, copy.mode = TRUE, copy.date = TRUE)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to copy '%s' -> '%s'", from, to))
  }
}

# Commit a set of prepared files with rollback if any replacement fails.
commit_with_rollback <- function(staged_paths, destination_paths, transaction_dir) {
  if (length(staged_paths) != length(destination_paths)) {
    stop("Internal error: staged/destination file counts do not match.")
  }

  backup_dir <- file.path(transaction_dir, "backups")
  dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)

  backup_paths <- rep(NA_character_, length(destination_paths))
  existed_before <- file.exists(destination_paths)

  for (i in seq_along(destination_paths)) {
    if (existed_before[i]) {
      backup_paths[i] <- file.path(backup_dir, sprintf("backup_%02d", i))
      copy_or_stop(destination_paths[i], backup_paths[i], overwrite = TRUE)
    }
  }

  committed <- logical(length(destination_paths))

  rollback <- function() {
    for (j in which(committed)) {
      if (existed_before[j] && !is.na(backup_paths[j]) && file.exists(backup_paths[j])) {
        file.copy(backup_paths[j], destination_paths[j], overwrite = TRUE)
      } else if (!existed_before[j] && file.exists(destination_paths[j])) {
        unlink(destination_paths[j])
      }
    }
  }

  for (i in seq_along(destination_paths)) {
    ok <- file.copy(
      staged_paths[i],
      destination_paths[i],
      overwrite = TRUE,
      copy.mode = TRUE,
      copy.date = TRUE
    )
    if (!isTRUE(ok)) {
      committed[i] <- TRUE
      rollback()
      stop(sprintf(
        "Promotion transaction failed while replacing '%s'. Previous files were restored where possible.",
        destination_paths[i]
      ))
    }
    committed[i] <- TRUE
  }

  invisible(TRUE)
}

# ── 2. Parse command-line arguments ──────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)

if (any(args %in% c("-h", "--help"))) {
  print_usage()
  quit(status = 0)
}

option_values <- list(
  promoted_tree_dir = NULL,
  tree_files_tsv = NULL,
  nhphyml_run_pk = NULL
)
positional <- character(0)

i <- 1
while (i <= length(args)) {
  arg <- args[i]

  if (startsWith(arg, "--promoted-tree-dir=")) {
    option_values$promoted_tree_dir <- sub("^--promoted-tree-dir=", "", arg)
  } else if (startsWith(arg, "--tree-files-tsv=")) {
    option_values$tree_files_tsv <- sub("^--tree-files-tsv=", "", arg)
  } else if (startsWith(arg, "--nhphyml-run-pk=")) {
    option_values$nhphyml_run_pk <- sub("^--nhphyml-run-pk=", "", arg)
  } else if (arg %in% c("--promoted-tree-dir", "--tree-files-tsv", "--nhphyml-run-pk")) {
    if (i == length(args)) {
      stop(sprintf("Option %s requires a value.", arg))
    }
    value <- args[i + 1]
    if (arg == "--promoted-tree-dir") option_values$promoted_tree_dir <- value
    if (arg == "--tree-files-tsv") option_values$tree_files_tsv <- value
    if (arg == "--nhphyml-run-pk") option_values$nhphyml_run_pk <- value
    i <- i + 1
  } else if (startsWith(arg, "--")) {
    stop(sprintf("Unknown option: %s", arg))
  } else {
    positional <- c(positional, arg)
  }

  i <- i + 1
}

if (length(positional) != 4) {
  print_usage()
  stop(sprintf("Expected 4 positional arguments; received %d.", length(positional)))
}

genomes_dir <- clean_path(positional[1])
outgroup_file <- clean_path(positional[2])
newick_file <- clean_path(positional[3])
manifest_file <- clean_path(positional[4])

promoted_tree_dir <- prompt_if_missing(
  option_values$promoted_tree_dir,
  "Enter promoted tree-set directory: "
)
promoted_tree_dir <- clean_path(promoted_tree_dir)

tree_files_tsv <- prompt_if_missing(
  option_values$tree_files_tsv,
  "Enter tree_files.tsv path: "
)
tree_files_tsv <- clean_path(tree_files_tsv)

nhphyml_run_pk_raw <- prompt_if_missing(
  option_values$nhphyml_run_pk,
  "Enter nhphyml_run_pk: "
)
if (!grepl("^[1-9][0-9]*$", nhphyml_run_pk_raw)) {
  stop("nhphyml_run_pk must be a positive integer.")
}
nhphyml_run_pk <- as.integer(nhphyml_run_pk_raw)

base_name <- tools::file_path_sans_ext(basename(newick_file))
alignments_dir <- file.path(genomes_dir, "records", "compleasm", "alignments")
rooted_dir <- file.path(alignments_dir, "06_rooted_tree")

cat("==============================================================================\n")
cat("08a_root_a_tree_v6.R\n")
cat("==============================================================================\n")
cat(sprintf("Genomes directory:      %s\n", genomes_dir))
cat(sprintf("Outgroup file:          %s\n", outgroup_file))
cat(sprintf("Input tree:             %s\n", newick_file))
cat(sprintf("Project manifest:       %s\n", manifest_file))
cat(sprintf("Promoted tree dir:      %s\n", promoted_tree_dir))
cat(sprintf("tree_files.tsv:         %s\n", tree_files_tsv))
cat(sprintf("nhphyml_run_pk:         %d\n\n", nhphyml_run_pk))

# ── 3. Validate required paths ───────────────────────────────────────────────

if (!dir.exists(genomes_dir)) {
  stop(sprintf("Genomes directory not found: %s", genomes_dir))
}
if (!dir.exists(alignments_dir)) {
  stop(sprintf("Alignments directory not found: %s", alignments_dir))
}
if (!file.exists(outgroup_file)) {
  stop(sprintf("Outgroup list not found: %s", outgroup_file))
}
if (!file.exists(newick_file)) {
  stop(sprintf("Newick file not found: %s", newick_file))
}
if (!file.exists(manifest_file)) {
  stop(sprintf("Project manifest not found: %s", manifest_file))
}
if (!file.exists(tree_files_tsv)) {
  stop(sprintf("tree_files.tsv not found: %s", tree_files_tsv))
}

if (!dir.exists(rooted_dir)) {
  dir.create(rooted_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", rooted_dir))
}

# ── 4. Load project manifest and build authoritative label map ───────────────

required_manifest_columns <- c(
  "project_name",
  "project_policy",
  "manifest_version",
  "manifest_date",
  "species_key",
  "organism_name",
  "accession",
  "accession_root",
  "assembly_name",
  "project_role",
  "repository",
  "source_database",
  "source_record_url",
  "is_preferred_for_species_at_freeze",
  "is_current_for_species_at_freeze",
  "analysis_notes"
)

manifest <- read.csv(
  manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  strip.white = TRUE
)

missing_manifest_columns <- setdiff(required_manifest_columns, names(manifest))
if (length(missing_manifest_columns) > 0) {
  stop(sprintf(
    "Project manifest is missing required column(s): %s",
    paste(missing_manifest_columns, collapse = ", ")
  ))
}

manifest$accession <- trimws(manifest$accession)
manifest$organism_name <- trimws(manifest$organism_name)

if (any(is.na(manifest$accession) | !nzchar(manifest$accession))) {
  stop("Project manifest contains blank/NA accession value(s).")
}
if (any(is.na(manifest$organism_name) | !nzchar(manifest$organism_name))) {
  stop("Project manifest contains blank/NA organism_name value(s).")
}

manifest_duplicate_accessions <- unique(manifest$accession[duplicated(manifest$accession)])
if (length(manifest_duplicate_accessions) > 0) {
  cat("\nERROR: Duplicate accession values in project manifest:\n")
  cat(paste(" -", manifest_duplicate_accessions), sep = "\n")
  cat("\n")
  stop("Project manifest accession values must be unique.")
}

label_map <- setNames(manifest$organism_name, manifest$accession)
cat(sprintf(
  "Loaded %d accession -> species name mappings from %s\n",
  length(label_map),
  basename(manifest_file)
))

# ── 5. Read outgroup and input tree ──────────────────────────────────────────

outgroup_taxa <- readLines(outgroup_file, warn = FALSE)
outgroup_taxa <- trimws(outgroup_taxa)
outgroup_taxa <- outgroup_taxa[
  nchar(outgroup_taxa) > 0 & !startsWith(outgroup_taxa, "#")
]
cat(sprintf(
  "Loaded %d outgroup taxa from %s\n",
  length(outgroup_taxa),
  basename(outgroup_file)
))

tree <- read.tree(newick_file)

if (inherits(tree, "multiPhylo")) {
  cat(sprintf("Multiple trees detected - using tree 1 of %d\n", length(tree)))
  tree <- tree[[1]]
}

cat(sprintf("Tree has %d tips\n", length(tree$tip.label)))

if (anyDuplicated(tree$tip.label)) {
  duplicated_tree_tips <- unique(tree$tip.label[duplicated(tree$tip.label)])
  stop(sprintf(
    "Input tree contains duplicate tip labels: %s",
    paste(duplicated_tree_tips, collapse = ", ")
  ))
}

# ── 6. Strict manifest/tree consistency ──────────────────────────────────────

tree_tips <- tree$tip.label
tree_missing_from_manifest <- setdiff(tree_tips, names(label_map))
manifest_missing_from_tree <- setdiff(names(label_map), tree_tips)

if (length(tree_missing_from_manifest) > 0) {
  cat("\nERROR: The following tree accessions are absent from the project manifest:\n")
  cat(paste(" -", tree_missing_from_manifest), sep = "\n")
  cat("\n")
  stop("Every input-tree tip must occur in the supplied project manifest.")
}

cat(sprintf(
  "Manifest/tree consistency: all %d tree tips are present in the project manifest\n",
  length(tree_tips)
))

if (length(manifest_missing_from_tree) > 0) {
  cat(sprintf(
    "Manifest accessions absent from this tree: %d (reported only; not fatal)\n",
    length(manifest_missing_from_tree)
  ))
  cat(paste(" -", manifest_missing_from_tree), sep = "\n")
  cat("\n")
} else {
  cat("All manifest accessions are represented in the input tree\n")
}

# ── 7. Plot original IQ-TREE tree with full manifest organism names ──────────

original_display <- rename_tips_full(tree, label_map)

original_pdf <- file.path(rooted_dir, paste0(base_name, ".original_iqtree.pdf"))
plot_tree_pdf(
  tree = original_display,
  pdf_path = original_pdf,
  title_str = sprintf(
    "Original IQ-TREE topology (unrooted) - %d tips",
    length(original_display$tip.label)
  )
)

# ── 8. Reconcile outgroup names against accession-labelled tree tips ─────────

matched_outgroup <- outgroup_taxa[outgroup_taxa %in% tree_tips]
missing_outgroup <- outgroup_taxa[!outgroup_taxa %in% tree_tips]

if (length(missing_outgroup) > 0) {
  cat("\nWARNING: The following outgroup taxa were NOT found in the tree:\n")
  cat(paste(" -", missing_outgroup), sep = "\n")
  cat("\n")
}

if (length(matched_outgroup) == 0) {
  stop("No outgroup taxa match any tip in the tree.")
}

cat(sprintf("Using %d matched outgroup taxa\n", length(matched_outgroup)))

# ── 9. Root the tree ─────────────────────────────────────────────────────────

ingroup_tips <- tree$tip.label[!tree$tip.label %in% matched_outgroup]
cat(sprintf("Ingroup has %d taxa\n", length(ingroup_tips)))

if (length(matched_outgroup) == 1) {
  cat("Single outgroup taxon - rooting on its branch\n")
  rooted_tree <- root(
    tree,
    outgroup = matched_outgroup,
    resolve.root = TRUE
  )
} else {
  cat("Multiple outgroup taxa - rooting with resolve.root\n")
  rooted_tree <- root(
    tree,
    outgroup = matched_outgroup,
    resolve.root = TRUE
  )

  cat(sprintf(
    "Outgroup MRCA node after rerooting: %d\n",
    findMRCA(rooted_tree, matched_outgroup)
  ))
}

if (!is.rooted(rooted_tree)) {
  stop("Tree is still unrooted after rerooting attempt.")
}
cat("Tree successfully rooted\n")

# ── 10. Verify outgroup monophyly ────────────────────────────────────────────

if (length(matched_outgroup) > 1) {
  mono_check <- is.monophyletic(rooted_tree, matched_outgroup)
  if (!mono_check) {
    cat("\nWARNING: Outgroup is NOT monophyletic in this tree!\n")
    cat("Inspect the original IQ-TREE PDF before proceeding.\n\n")
  } else {
    cat("Outgroup monophyly confirmed\n")
  }
} else {
  mono_check <- NA
}

# ── 11. Write rooted accession-labelled tree ─────────────────────────────────

rooted_file <- file.path(rooted_dir, paste0(base_name, ".rooted.nwk"))
write.tree(rooted_tree, file = rooted_file, digits = 10)
cat(sprintf("Rooted tree written to: %s\n", rooted_file))

# ── 12. Write Genus_species rooted tree and rooted-check PDF ─────────────────

species_tree <- rename_tips_genus_species(rooted_tree, label_map)

renamed_file <- file.path(rooted_dir, paste0(base_name, ".rooted_renamed.nwk"))
write.tree(species_tree, file = renamed_file, digits = 10)
cat(sprintf("Renamed tree written to: %s\n", renamed_file))

rooted_check_pdf <- file.path(rooted_dir, paste0(base_name, ".rooted_check.pdf"))
plot_tree_pdf(
  tree = species_tree,
  pdf_path = rooted_check_pdf,
  title_str = sprintf(
    "Rooted tree - %d tips | outgroup (%d taxa)",
    length(species_tree$tip.label),
    length(matched_outgroup)
  )
)

# ── 13. Write NHPhyML topology-only tree ─────────────────────────────────────

nhphyml_tree <- species_tree

# Remove substitution-per-site branch lengths.
nhphyml_tree$edge.length <- NULL

# Remove all internal-node labels, including IQ-TREE support/confidence values
# and any root label such as "Root".
nhphyml_tree$node.label <- NULL

nhphyml_file <- file.path(
  rooted_dir,
  paste0(base_name, ".rooted_renamed_topology_only_nhphyml_format.nwk")
)
write.tree(nhphyml_tree, file = nhphyml_file)
cat(sprintf(
  "NHPhyML-formatted topology-only tree written to: %s\n",
  nhphyml_file
))

# ── 14. Promotion preflight ──────────────────────────────────────────────────

if (!dir.exists(promoted_tree_dir)) {
  dir.create(promoted_tree_dir, recursive = TRUE)
  cat(sprintf("Created promoted tree-set directory: %s\n", promoted_tree_dir))
}

if (!dir.exists(promoted_tree_dir)) {
  stop(sprintf("Unable to create promoted tree-set directory: %s", promoted_tree_dir))
}

expected_tsv_columns <- c(
  "tree_pk",
  "nhphyml_run_pk",
  "tree_type",
  "tree_text",
  "branch_length_type",
  "notes"
)

tsv_lines <- readLines(tree_files_tsv, warn = FALSE)
if (length(tsv_lines) == 0) {
  stop("tree_files.tsv is empty and does not contain a header.")
}

tsv_header <- strsplit(tsv_lines[1], "\t", fixed = TRUE)[[1]]
missing_tsv_columns <- setdiff(expected_tsv_columns, tsv_header)
if (length(missing_tsv_columns) > 0) {
  stop(sprintf(
    "tree_files.tsv is missing required column(s): %s",
    paste(missing_tsv_columns, collapse = ", ")
  ))
}

existing_tree_pks <- read_tree_pk_values(tree_files_tsv, tsv_header)
next_tree_pk <- if (length(existing_tree_pks) == 0) 1L else max(existing_tree_pks) + 1L
new_tree_pks <- next_tree_pk + 0:2

promoted_names <- c(
  "1_unrooted.treefile",
  "2_rooted_pgls_input.nwk",
  "3_rooted_nhphyml_input.nwk"
)

source_paths <- c(
  newick_file,
  renamed_file,
  nhphyml_file
)

promoted_paths <- file.path(promoted_tree_dir, promoted_names)

# Ask about all conflicts before modifying any promoted-tree/README/TSV file.
promotion_cancelled <- FALSE
existing_promoted <- file.exists(promoted_paths)

if (any(existing_promoted)) {
  cat("\nExisting promoted tree files detected. Overwrite approval is required.\n")
  for (idx in which(existing_promoted)) {
    response <- tolower(trimws(readline(sprintf(
      "%s already exists. Overwrite? [y/N]: ",
      promoted_names[idx]
    ))))
    if (!response %in% c("y", "yes")) {
      promotion_cancelled <- TRUE
    }
  }
}

if (promotion_cancelled) {
  cat("\nPromotion cancelled because at least one existing promoted file was not approved for overwrite.\n")
  cat("No promoted files, README entries, or tree_files.tsv rows were changed.\n")
  cat("The normal rooted-tree outputs in 06_rooted_tree were retained.\n\n")
  cat("==============================================================================\n")
  cat("08a_root_a_tree_v6.R COMPLETE (promotion cancelled)\n")
  cat("==============================================================================\n")
  quit(status = 0)
}

# ── 15. Stage promoted tree files and metadata updates ───────────────────────

transaction_dir <- tempfile(
  pattern = ".08a_root_a_tree_v6_transaction_",
  tmpdir = promoted_tree_dir
)
dir.create(transaction_dir, recursive = TRUE)

staged_promoted_paths <- file.path(transaction_dir, promoted_names)
for (idx in seq_along(source_paths)) {
  copy_or_stop(source_paths[idx], staged_promoted_paths[idx], overwrite = TRUE)
}

# Validate staged tree contents before preparing metadata updates.
staged_tree_texts <- vapply(
  staged_promoted_paths,
  read_tree_text_one_line,
  character(1),
  USE.NAMES = FALSE
)

run_date <- format(Sys.Date(), "%Y-%m-%d")
run_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

# Prepare README.md as a full replacement that preserves prior history.
readme_path <- file.path(promoted_tree_dir, "README.md")
existing_readme <- if (file.exists(readme_path)) {
  readLines(readme_path, warn = FALSE)
} else {
  "# Promoted phylogenetic tree sets"
}

readme_section <- c(
  "",
  sprintf("## Promotion run — %s", run_timestamp),
  "",
  sprintf("nhphyml_run_pk: %d", nhphyml_run_pk),
  "",
  "| file | og_file | og_path |",
  "|---|---|---|",
  sprintf("| %s | %s | %s |", promoted_names[1], basename(source_paths[1]), source_paths[1]),
  sprintf("| %s | %s | %s |", promoted_names[2], basename(source_paths[2]), source_paths[2]),
  sprintf("| %s | %s | %s |", promoted_names[3], basename(source_paths[3]), source_paths[3])
)

staged_readme <- file.path(transaction_dir, "README.md")
writeLines(c(existing_readme, readme_section), staged_readme, useBytes = TRUE)

# Prepare three tree_files.tsv rows in memory.
new_records <- list(
  list(
    tree_pk = new_tree_pks[1],
    nhphyml_run_pk = nhphyml_run_pk,
    tree_type = "unrooted",
    tree_text = staged_tree_texts[1],
    branch_length_type = "substitution_per_site",
    notes = sprintf(
      "IQ-TREE 3; promoted by 08a_root_a_tree_v6.R; created %s",
      run_date
    )
  ),
  list(
    tree_pk = new_tree_pks[2],
    nhphyml_run_pk = nhphyml_run_pk,
    tree_type = "rooted",
    tree_text = staged_tree_texts[2],
    branch_length_type = "substitution_per_site",
    notes = sprintf(
      "IQ-TREE 3; rooted/formatted by 08a_root_a_tree_v6.R; PGLS input; created %s",
      run_date
    )
  ),
  list(
    tree_pk = new_tree_pks[3],
    nhphyml_run_pk = nhphyml_run_pk,
    tree_type = "rooted",
    tree_text = staged_tree_texts[3],
    branch_length_type = "none",
    notes = sprintf(
      "IQ-TREE 3; rooted/formatted by 08a_root_a_tree_v6.R; NHPhyML topology-only input; created %s",
      run_date
    )
  )
)

new_tsv_lines <- vapply(
  new_records,
  make_tsv_line,
  character(1),
  header_fields = tsv_header,
  USE.NAMES = FALSE
)

# Verify that TSV serialization preserves every value exactly, especially
# species names in tree_text and the literals "unrooted", "rooted", and
# "substitution_per_site". This would have caught the v5 corruption.
for (record_idx in seq_along(new_records)) {
  expected_fields <- vapply(
    tsv_header,
    function(field) {
      if (field %in% names(new_records[[record_idx]])) {
        as.character(new_records[[record_idx]][[field]])
      } else {
        ""
      }
    },
    character(1),
    USE.NAMES = FALSE
  )
  # Compare serialized lines directly: this also handles optional extra
  # TSV columns whose value is blank (including trailing empty columns).
  if (!identical(new_tsv_lines[record_idx], paste(expected_fields, collapse = "\t"))) {
    stop(sprintf(
      "TSV serialization changed at least one field in new tree row %d; refusing to write.",
      record_idx
    ))
  }
}

if (any(!nzchar(vapply(new_records, function(x) x$tree_text, character(1))))) {
  stop("At least one staged tree_text value is empty; refusing to update tree_files.tsv.")
}

staged_tsv <- file.path(transaction_dir, "tree_files.tsv")
writeLines(c(tsv_lines, new_tsv_lines), staged_tsv, useBytes = TRUE)

# Final staged-file validation.
all_staged <- c(staged_promoted_paths, staged_readme, staged_tsv)
if (!all(file.exists(all_staged))) {
  stop("Promotion staging failed: one or more staged files are missing.")
}
if (any(file.info(all_staged)$size <= 0)) {
  stop("Promotion staging failed: one or more staged files are empty.")
}

# ── 16. Commit promotion + README + tree_files.tsv together ─────────────────

commit_destinations <- c(
  promoted_paths,
  readme_path,
  tree_files_tsv
)
commit_staged <- c(
  staged_promoted_paths,
  staged_readme,
  staged_tsv
)

commit_with_rollback(
  staged_paths = commit_staged,
  destination_paths = commit_destinations,
  transaction_dir = transaction_dir
)

# Verify committed promoted trees and metadata files.
if (!all(file.exists(promoted_paths))) {
  stop("Post-commit validation failed: one or more promoted tree files are missing.")
}
if (!file.exists(readme_path)) {
  stop("Post-commit validation failed: README.md is missing.")
}
if (!file.exists(tree_files_tsv)) {
  stop("Post-commit validation failed: tree_files.tsv is missing.")
}

# ── 17. Final report ─────────────────────────────────────────────────────────

cat("\n==============================================================================\n")
cat("08a_root_a_tree_v6.R COMPLETE\n")
cat("==============================================================================\n\n")

cat("Manifest:\n")
cat(sprintf("  %s\n", manifest_file))
cat(sprintf("  Loaded %d accession -> species name mappings\n\n", length(label_map)))

cat("Outgroup:\n")
cat(sprintf("  %s\n", basename(outgroup_file)))
cat(sprintf("  Loaded %d outgroup taxa\n\n", length(outgroup_taxa)))

cat("Manifest/tree consistency:\n")
cat(sprintf("  %d/%d tree tips found in manifest\n", length(tree_tips), length(tree_tips)))
cat(sprintf("  %d manifest accessions absent from tree\n\n", length(manifest_missing_from_tree)))

cat("Rooting:\n")
cat("  Tree successfully rooted\n")
if (length(matched_outgroup) > 1) {
  cat(sprintf("  Outgroup monophyly: %s\n\n", ifelse(isTRUE(mono_check), "confirmed", "NOT confirmed")))
} else {
  cat("  Outgroup monophyly: not applicable (single outgroup taxon)\n\n")
}

cat("Rooted-tree outputs:\n")
cat(sprintf("  %s\n", original_pdf))
cat(sprintf("  %s\n", rooted_check_pdf))
cat(sprintf("  %s\n", rooted_file))
cat(sprintf("  %s\n", renamed_file))
cat(sprintf("  %s\n\n", nhphyml_file))

cat("Promoted tree set:\n")
cat(sprintf("  Directory: %s\n", promoted_tree_dir))
for (idx in seq_along(promoted_names)) {
  overwrite_note <- if (existing_promoted[idx]) " (overwritten)" else ""
  cat(sprintf("  %s%s\n", promoted_names[idx], overwrite_note))
}
cat("\n")

cat("tree_files.tsv:\n")
cat(sprintf("  %s\n", tree_files_tsv))
cat(sprintf("  nhphyml_run_pk = %d\n", nhphyml_run_pk))
cat(sprintf("  Added tree_pk = %s\n\n", paste(new_tree_pks, collapse = ", ")))

cat("README:\n")
cat(sprintf("  %s\n", readme_path))
cat("==============================================================================\n")

unlink(transaction_dir, recursive = TRUE, force = TRUE)
