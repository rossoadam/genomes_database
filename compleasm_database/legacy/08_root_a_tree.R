#!/usr/bin/env Rscript
# root_for_nhphyml.R
# Usage: Rscript root_for_nhphyml.R <genomes_dir> <outgroup_list> <newick_file>

suppressPackageStartupMessages({
  library(phytools)
  library(ape)
})

# ── 1. Parse and validate arguments ──────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  cat("Usage: Rscript root_for_nhphyml.R <genomes_dir> <outgroup_list> <newick_file>\n")
  quit(status = 1)
}

genomes_dir   <- args[1]
outgroup_file <- args[2]
newick_file   <- args[3]

base_name <- tools::file_path_sans_ext(basename(newick_file))

alignments_dir <- file.path(genomes_dir, "records", "compleasm", "alignments")
rooted_dir     <- file.path(alignments_dir, "06_rooted_tree")

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

if (!dir.exists(rooted_dir)) {
  dir.create(rooted_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", rooted_dir))
}

# ── 2. Load metadata and build label map ─────────────────────────────────────

metadata_file <- file.path(genomes_dir, "records", "compleasm", "records", "metadata.csv")

label_map <- NULL

if (!file.exists(metadata_file)) {
  cat(sprintf("WARNING: metadata.csv not found at %s - labels will be raw accessions\n",
              metadata_file))
} else {
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE,
                       header = TRUE, strip.white = TRUE, quote = "")
  label_map <- setNames(metadata$organism_name, metadata$accession)
  cat(sprintf("Loaded %d accession -> species name mappings\n", nrow(metadata)))
}

rename_tips <- function(tree, label_map) {
  if (is.null(label_map)) return(tree)
  current_tips  <- tree$tip.label
  tree$tip.label <- ifelse(
    current_tips %in% names(label_map),
    label_map[current_tips],
    current_tips
  )
  n_renamed <- sum(current_tips %in% names(label_map))
  n_missing <- sum(!current_tips %in% names(label_map))
  cat(sprintf("Renamed %d tips | %d tips had no match (kept as accession)\n",
              n_renamed, n_missing))
  return(tree)
}

plot_tree_pdf <- function(tree, pdf_path, title_str) {
  n_tips   <- length(tree$tip.label)
  fig_h    <- max(8, n_tips * 0.18)
  fnt_size <- min(0.7, 25 / n_tips)
  pdf(pdf_path, width = 12, height = fig_h)
  plotTree(tree,
           fsize     = fnt_size,
           lwd       = 0.8,
           mar       = c(1, 0, 2, 1),
           direction = "rightwards")
  title(main = title_str, cex.main = 0.9)
  dev.off()
  cat(sprintf("PDF written: %s\n", pdf_path))
}

# ── 3. Read inputs ────────────────────────────────────────────────────────────

outgroup_taxa <- readLines(outgroup_file)
outgroup_taxa <- trimws(outgroup_taxa)
outgroup_taxa <- outgroup_taxa[nchar(outgroup_taxa) > 0 & !startsWith(outgroup_taxa, "#")]

cat(sprintf("Loaded %d outgroup taxa (geckos)\n", length(outgroup_taxa)))

tree <- read.tree(newick_file)

if (inherits(tree, "multiPhylo")) {
  cat(sprintf("Multiple trees detected - using tree 1 of %d\n", length(tree)))
  tree <- tree[[1]]
}

cat(sprintf("Tree has %d tips\n", length(tree$tip.label)))

# ── 4. Plot original IQ-TREE tree with species names ─────────────────────────

original_display <- rename_tips(tree, label_map)

plot_tree_pdf(
  tree      = original_display,
  pdf_path  = file.path(rooted_dir, paste0(base_name, ".original_iqtree.pdf")),
  title_str = sprintf("Original IQ-TREE topology (unrooted) - %d tips", length(original_display$tip.label))
)

# ── 5. Reconcile outgroup names against tree tips ─────────────────────────────

tree_tips        <- tree$tip.label
matched_outgroup <- outgroup_taxa[outgroup_taxa %in% tree_tips]
missing_outgroup <- outgroup_taxa[!outgroup_taxa %in% tree_tips]

if (length(missing_outgroup) > 0) {
  cat("\nWARNING: The following outgroup taxa were NOT found in the tree:\n")
  cat(paste(" -", missing_outgroup, collapse = "\n"), "\n\n")
}

if (length(matched_outgroup) == 0) {
  stop("No outgroup taxa match any tip in the tree.")
}

cat(sprintf("Using %d matched gecko outgroup taxa\n", length(matched_outgroup)))

# ── 6. Root the tree ──────────────────────────────────────────────────────────

ingroup_tips <- tree$tip.label[!tree$tip.label %in% matched_outgroup]
cat(sprintf("Ingroup has %d taxa\n", length(ingroup_tips)))

if (length(matched_outgroup) == 1) {
  
  cat("Single outgroup taxon - rooting on its branch\n")
  rooted_tree <- root(tree,
                      outgroup      = matched_outgroup,
                      resolve.root  = TRUE)
  
} else {
  
  cat("Multiple outgroup taxa - rooting with resolve.root\n")
  
  # Use all matched outgroup taxa — ape::root() with a multi-tip outgroup
  # finds their MRCA branch and roots there, resolve.root = TRUE forces
  # a clean bifurcation instead of leaving a trifurcation at the root
  rooted_tree <- root(tree,
                      outgroup      = matched_outgroup,
                      resolve.root  = TRUE)
  
  cat(sprintf("Outgroup MRCA node after rerooting: %d\n",
              findMRCA(rooted_tree, matched_outgroup)))
}

if (!is.rooted(rooted_tree)) {
  stop("Tree is still unrooted after rerooting attempt")
}

cat("Tree successfully rooted\n")

# ── 7. Verify outgroup monophyly ──────────────────────────────────────────────

if (length(matched_outgroup) > 1) {
  mono_check <- is.monophyletic(rooted_tree, matched_outgroup)
  if (!mono_check) {
    cat("\nWARNING: Gecko outgroup is NOT monophyletic in this tree!\n")
    cat("Inspect the original IQ-TREE PDF before proceeding.\n\n")
  } else {
    cat("Outgroup monophyly confirmed\n")
  }
}

# ── 8. Write rooted newick ────────────────────────────────────────────────────

output_file <- file.path(rooted_dir, paste0(base_name, ".rooted.nwk"))
write.tree(rooted_tree, file = output_file, digits = 10)
cat(sprintf("\nRooted tree written to: %s\n", output_file))

# ── 9. Write renamed rooted newick and PDF ────────────────────────────────────

display_tree  <- rename_tips(rooted_tree, label_map)

renamed_file <- file.path(rooted_dir, paste0(base_name, ".rooted_renamed.nwk"))
write.tree(display_tree, file = renamed_file, digits = 10)
cat(sprintf("Renamed tree written to: %s\n", renamed_file))

plot_tree_pdf(
  tree      = display_tree,
  pdf_path  = file.path(rooted_dir, paste0(base_name, ".rooted_check.pdf")),
  title_str = sprintf("Rooted tree - %d tips | gecko outgroup (%d taxa)",
                      length(display_tree$tip.label), length(matched_outgroup))
)

cat("\nDone. Feed the .rooted.nwk to nhphyml with -u flag.\n")