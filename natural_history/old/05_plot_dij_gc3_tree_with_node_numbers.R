#!/usr/bin/env Rscript
# 05_plot_dij_gc3_tree.R
#
# Plot a tree using d_ij values as branch lengths and GC3 values as branch colors.
#
# Required arguments:
#   --tree /path/to/tree.nwk
#   --dij /path/to/dij_results.json OR /path/to/dij_results.tsv
#   --gc3-summary /path/to/gc3_summary.tsv
#
# Optional but strongly recommended:
#   --node-numbers /path/to/node_numbers_rebuilt.csv
#
# Optional:
#   --outdir /path/to/output_dir
#   --prefix output_prefix
#
# Outputs:
#   <prefix>.dij_branch_lengths.nwk
#   <prefix>.dij_branch_lengths_gc3_gradient.pdf
#   <prefix>.edge_plot_table.tsv
#   <prefix>.parse_diagnostics.log
#   <prefix>.ape_to_pipeline_node_map.tsv
#
# Visual defaults:
#   - branch line width = 3
#   - GC3 color gradient = black for low GC3, light blue for high GC3
#   - missing GC3 edges are colored grey80 but not included in the legend

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

# ── 1. Argument parsing ──────────────────────────────────────────────────────

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected positional argument: %s", key))
    }
    if (i == length(args)) {
      stop(sprintf("Missing value after argument: %s", key))
    }
    value <- args[i + 1]
    clean_key <- sub("^--", "", key)
    clean_key <- gsub("-", "_", clean_key)
    out[[clean_key]] <- value
    i <- i + 2
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

required_args <- c("tree", "dij", "gc3_summary")
missing_args <- required_args[!required_args %in% names(args)]

if (length(missing_args) > 0) {
  cat("Usage:\n")
  cat("  Rscript 05_plot_dij_gc3_tree.R \\\n")
  cat("    --tree /path/to/tree.nwk \\\n")
  cat("    --dij /path/to/dij_results.tsv \\\n")
  cat("    --gc3-summary /path/to/gc3_summary.tsv \\\n")
  cat("    --node-numbers /path/to/node_numbers_rebuilt.csv \\\n")
  cat("    --outdir /path/to/output_dir \\\n")
  cat("    [--prefix output_prefix]\n\n")
  stop(sprintf("Missing required argument(s): %s", paste(missing_args, collapse = ", ")))
}

tree_file <- args$tree
dij_file <- args$dij
gc3_summary_file <- args$gc3_summary
node_numbers_file <- if ("node_numbers" %in% names(args)) args$node_numbers else NA_character_

if (!file.exists(tree_file)) stop(sprintf("Tree file not found: %s", tree_file))
if (!file.exists(dij_file)) stop(sprintf("dij file not found: %s", dij_file))
if (!file.exists(gc3_summary_file)) stop(sprintf("gc3 summary file not found: %s", gc3_summary_file))
if (!is.na(node_numbers_file) && !file.exists(node_numbers_file)) {
  stop(sprintf("node_numbers file not found: %s", node_numbers_file))
}

base_name <- tools::file_path_sans_ext(basename(tree_file))
outdir <- if ("outdir" %in% names(args)) args$outdir else file.path(dirname(tree_file), "gc3_dij_tree_plots")
outdir <- path.expand(outdir)
prefix <- if ("prefix" %in% names(args)) args$prefix else paste0(base_name, ".dij_gc3")

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(outdir, paste0(prefix, ".parse_diagnostics.log"))
if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = log_file, append = TRUE)
}

log_section <- function(title) {
  line <- paste(rep("=", nchar(title)), collapse = "")
  log_msg("")
  log_msg(line)
  log_msg(title)
  log_msg(line)
}

log_msg("Log file: %s", log_file)
log_msg("Output directory: %s", normalizePath(outdir, mustWork = FALSE))
log_msg("Tree file: %s", tree_file)
log_msg("dij file: %s", dij_file)
log_msg("GC3 summary file: %s", gc3_summary_file)
if (!is.na(node_numbers_file)) {
  log_msg("node_numbers file: %s", node_numbers_file)
}

# ── 2. Helpers ───────────────────────────────────────────────────────────────

normalize_species_name <- function(x) {
  x <- trimws(as.character(x))
  if (grepl("_", x)) {
    parts <- unlist(strsplit(x, "_"))
    parts <- parts[nchar(parts) > 0]
  } else {
    parts <- unlist(strsplit(x, "\\s+"))
    parts <- parts[nchar(parts) > 0]
  }

  if (length(parts) >= 2) {
    return(tolower(paste(parts[1], parts[2], sep = "_")))
  } else if (length(parts) == 1) {
    return(tolower(parts[1]))
  }
  return("unknown_species")
}

coerce_node_name <- function(x) {
  x <- trimws(as.character(x))
  if (is.na(x) || nchar(x) == 0) return("")

  low <- tolower(x)
  if (low %in% c("ancestor", "ancestral", "root")) return("ancestral")

  if (grepl("^(node[ _-]*)?[0-9]+$", low)) {
    num <- sub("^node[ _-]*", "", low)
    return(paste0("node_", as.integer(num)))
  }

  if (grepl(" ", x) && !grepl("_", x)) return(normalize_species_name(x))

  if (!startsWith(low, "node_")) return(normalize_species_name(x))

  return(low)
}

read_table_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) {
    return(read.delim(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (ext == "csv") {
    return(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  # Fallback to tab-delimited
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

descendant_tips_from_tree_node <- function(tree, node_number) {
  desc <- phytools::getDescendants(tree, node_number)
  tip_desc <- desc[desc <= length(tree$tip.label)]
  sort(vapply(tree$tip.label[tip_desc], normalize_species_name, character(1)))
}

signature_from_tips <- function(tips) {
  paste(sort(unique(tips)), collapse = "|")
}

read_node_relationships <- function(path) {
  df <- read_table_auto(path)

  lower_names <- tolower(trimws(names(df)))

  child_candidates <- c("daughter", "daughter_node", "child", "child_node", "node",
                        "node_id", "descendant", "descendant_node", "son", "son_node")
  parent_candidates <- c("parent", "parent_node", "ancestor", "ancestor_node", "mother", "parent_id")

  child_col_idx <- match(TRUE, lower_names %in% child_candidates)
  parent_col_idx <- match(TRUE, lower_names %in% parent_candidates)

  edges <- data.frame(child = character(), parent = character(), stringsAsFactors = FALSE)

  if (!is.na(child_col_idx) && !is.na(parent_col_idx) && child_col_idx != parent_col_idx) {
    child_col <- names(df)[child_col_idx]
    parent_col <- names(df)[parent_col_idx]
    edges <- data.frame(
      child = vapply(df[[child_col]], coerce_node_name, character(1)),
      parent = vapply(df[[parent_col]], coerce_node_name, character(1)),
      stringsAsFactors = FALSE
    )
    edges <- edges[nchar(edges$child) > 0 & nchar(edges$parent) > 0, ]
  } else {
    # Wide format: each column name is child, first non-empty row is parent.
    rows <- list()
    idx <- 1
    for (col in names(df)) {
      child <- coerce_node_name(col)
      parent_values <- df[[col]]
      parent_values <- parent_values[!is.na(parent_values) & nchar(trimws(as.character(parent_values))) > 0]
      if (length(parent_values) == 0) next

      for (p in parent_values) {
        parent <- coerce_node_name(p)
        if (nchar(child) > 0 && nchar(parent) > 0) {
          rows[[idx]] <- data.frame(child = child, parent = parent, stringsAsFactors = FALSE)
          idx <- idx + 1
        }
      }
    }
    if (length(rows) > 0) edges <- do.call(rbind, rows)
  }

  edges <- unique(edges)
  if (nrow(edges) == 0) stop(sprintf("No node relationships parsed from %s", path))
  edges
}

build_pipeline_descendant_signatures <- function(node_edges) {
  children_by_parent <- split(node_edges$child, node_edges$parent)
  memo <- new.env(parent = emptyenv())

  get_desc_tips <- function(node) {
    node <- coerce_node_name(node)
    if (exists(node, envir = memo, inherits = FALSE)) {
      return(get(node, envir = memo, inherits = FALSE))
    }

    if (!startsWith(node, "node_") && node != "ancestral") {
      tips <- normalize_species_name(node)
      assign(node, tips, envir = memo)
      return(tips)
    }

    kids <- children_by_parent[[node]]
    if (is.null(kids) || length(kids) == 0) {
      tips <- character()
      assign(node, tips, envir = memo)
      return(tips)
    }

    tips <- sort(unique(unlist(lapply(kids, get_desc_tips), use.names = FALSE)))
    assign(node, tips, envir = memo)
    tips
  }

  internal_nodes <- sort(unique(c(
    node_edges$parent[startsWith(node_edges$parent, "node_") | node_edges$parent == "ancestral"],
    node_edges$child[startsWith(node_edges$child, "node_")]
  )))

  out <- data.frame(
    pipeline_node = internal_nodes,
    descendant_signature = character(length(internal_nodes)),
    n_descendant_tips = integer(length(internal_nodes)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(internal_nodes)) {
    tips <- get_desc_tips(internal_nodes[i])
    out$descendant_signature[i] <- signature_from_tips(tips)
    out$n_descendant_tips[i] <- length(tips)
  }

  out
}

build_ape_to_pipeline_node_map <- function(tree, node_numbers_path, outdir, prefix) {
  node_edges <- read_node_relationships(node_numbers_path)
  pipeline_sig <- build_pipeline_descendant_signatures(node_edges)

  n_tips <- length(tree$tip.label)
  ape_internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)

  ape_sig <- data.frame(
    ape_node = ape_internal_nodes,
    descendant_signature = character(length(ape_internal_nodes)),
    n_descendant_tips = integer(length(ape_internal_nodes)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(ape_internal_nodes)) {
    tips <- descendant_tips_from_tree_node(tree, ape_internal_nodes[i])
    ape_sig$descendant_signature[i] <- signature_from_tips(tips)
    ape_sig$n_descendant_tips[i] <- length(tips)
  }

  map <- merge(
    ape_sig,
    pipeline_sig,
    by = "descendant_signature",
    all.x = TRUE,
    suffixes = c("_ape", "_pipeline")
  )

  # Remove ambiguous duplicate pipeline signatures if any.
  map <- map[order(map$ape_node), ]

  out_path <- file.path(outdir, paste0(prefix, ".ape_to_pipeline_node_map.tsv"))
  write.table(map, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)

  mapped <- sum(!is.na(map$pipeline_node))
  log_section("node_numbers mapping")
  log_msg("Parsed %d node relationship rows from node_numbers", nrow(node_edges))
  log_msg("APE internal nodes mapped to pipeline node IDs: %d/%d", mapped, nrow(map))
  log_msg("Mapping table written: %s", out_path)

  if (mapped < nrow(map)) {
    missing <- map[is.na(map$pipeline_node), c("ape_node", "n_descendant_tips_ape")]
    log_msg("WARNING: %d APE internal nodes did not map to a pipeline node.", nrow(missing))
    capture.output(print(head(missing, 50)), file = log_file, append = TRUE)
  }

  lookup <- map$pipeline_node
  names(lookup) <- as.character(map$ape_node)
  lookup
}

make_edge_child_keys <- function(tree, ape_to_pipeline = NULL) {
  keys <- character(nrow(tree$edge))
  n_tips <- length(tree$tip.label)

  for (i in seq_len(nrow(tree$edge))) {
    child <- tree$edge[i, 2]

    if (child <= n_tips) {
      keys[i] <- normalize_species_name(tree$tip.label[child])
    } else if (!is.null(ape_to_pipeline) && as.character(child) %in% names(ape_to_pipeline) &&
               !is.na(ape_to_pipeline[[as.character(child)]])) {
      keys[i] <- ape_to_pipeline[[as.character(child)]]
    } else {
      internal_index <- child - n_tips

      if (!is.null(tree$node.label) &&
          length(tree$node.label) >= internal_index &&
          !is.na(tree$node.label[internal_index]) &&
          nchar(trimws(tree$node.label[internal_index])) > 0) {
        keys[i] <- coerce_node_name(tree$node.label[internal_index])
      } else {
        keys[i] <- paste0("node_", child)
      }
    }
  }

  keys
}

read_dij <- function(path) {
  suffix <- tolower(tools::file_ext(path))

  if (suffix %in% c("tsv", "txt")) {
    dij <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)

    child_col <- intersect(c("node_or_species", "child", "daughter", "node"), names(dij))[1]
    parent_col <- intersect(c("parent", "parent_node"), names(dij))[1]
    dij_col <- intersect(c("dij", "d_ij", "distance", "branch_length"), names(dij))[1]

    if (is.na(child_col) || is.na(parent_col) || is.na(dij_col)) {
      stop(sprintf(
        "Could not identify child/parent/dij columns in %s. Columns found: %s",
        path, paste(names(dij), collapse = ", ")
      ))
    }

    out <- data.frame(
      child = vapply(dij[[child_col]], coerce_node_name, character(1)),
      parent = vapply(dij[[parent_col]], coerce_node_name, character(1)),
      dij = as.numeric(dij[[dij_col]]),
      stringsAsFactors = FALSE
    )
    out <- out[!is.na(out$dij) & nchar(out$child) > 0, ]
    return(out)
  }

  if (suffix == "json") {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Reading JSON dij files requires jsonlite. Install with install.packages('jsonlite'), or pass dij_results.tsv instead.")
    }

    obj <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    rows <- list()
    idx <- 1

    for (child_raw in names(obj)) {
      parent_obj <- obj[[child_raw]]
      if (length(parent_obj) == 0) next

      for (parent_raw in names(parent_obj)) {
        rows[[idx]] <- data.frame(
          child = coerce_node_name(child_raw),
          parent = coerce_node_name(parent_raw),
          dij = as.numeric(parent_obj[[parent_raw]]),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1
      }
    }

    if (length(rows) == 0) stop(sprintf("No dij rows parsed from JSON file: %s", path))
    out <- do.call(rbind, rows)
    out <- out[!is.na(out$dij) & nchar(out$child) > 0, ]
    return(out)
  }

  stop(sprintf("Unsupported dij file extension: %s. Use .tsv, .txt, or .json", suffix))
}

read_gc3_summary <- function(path) {
  gc3 <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)

  node_col <- intersect(c("node_or_species", "genus_species", "node", "species"), names(gc3))[1]
  gc3_col <- intersect(c("mean_gc3", "gc3", "GC3", "mean_gc"), names(gc3))[1]

  if (is.na(node_col) || is.na(gc3_col)) {
    stop(sprintf(
      "Could not identify node/species and GC3 columns in %s. Columns found: %s",
      path, paste(names(gc3), collapse = ", ")
    ))
  }

  keys <- vapply(gc3[[node_col]], coerce_node_name, character(1))
  values <- as.numeric(gc3[[gc3_col]])
  names(values) <- keys

  values <- values[!is.na(values) & nchar(names(values)) > 0]
  return(values)
}

color_from_gc3 <- function(gc3_values, keys) {
  pal <- colorRampPalette(c("black", "#00B8C2"))(100)

  vals <- as.numeric(gc3_values[keys])
  valid <- !is.na(vals)

  colors <- rep("grey80", length(keys))
  if (!any(valid)) return(colors)

  min_gc <- min(vals[valid], na.rm = TRUE)
  max_gc <- max(vals[valid], na.rm = TRUE)

  if (isTRUE(all.equal(min_gc, max_gc))) {
    colors[valid] <- pal[50]
    return(colors)
  }

  scaled <- (vals[valid] - min_gc) / (max_gc - min_gc)
  color_index <- pmax(1, pmin(100, round(scaled * 99) + 1))
  colors[valid] <- pal[color_index]
  colors
}

plot_tree_pdf <- function(tree, edge_colors, pdf_path, title_str, gc3_values) {
  n_tips <- length(tree$tip.label)
  fig_h <- max(8, n_tips * 0.18)
  fnt_size <- min(0.7, 25 / n_tips)

  valid_gc3 <- gc3_values[!is.na(gc3_values)]

  pdf(pdf_path, width = 15, height = fig_h)
  layout(matrix(c(1, 2), nrow = 1), widths = c(4.8, 1.2))

  par(mar = c(1, 1, 3, 0.5))
  plot.phylo(
    tree,
    type = "phylogram",
    direction = "rightwards",
    show.tip.label = TRUE,
    cex = fnt_size,
    edge.width = 3,
    edge.color = edge_colors,
    no.margin = FALSE
  )
  title(main = title_str, cex.main = 0.9)

  par(mar = c(3, 0.5, 3, 2))
  plot.new()

  if (length(valid_gc3) > 0) {
    pal <- colorRampPalette(c("black", "#00B8C2"))(100)
    y <- seq(0.15, 0.85, length.out = 100)

    for (i in seq_len(99)) {
      rect(
        xleft = 0.25,
        ybottom = y[i],
        xright = 0.5,
        ytop = y[i + 1],
        col = pal[i],
        border = NA
      )
    }

    text(0.62, 0.85, labels = sprintf("%.2f", max(valid_gc3)), adj = 0, cex = 0.8)
    text(0.62, 0.15, labels = sprintf("%.2f", min(valid_gc3)), adj = 0, cex = 0.8)
    text(0.25, 0.93, labels = "GC3", adj = 0, cex = 0.9, font = 2)
    text(0.25, 0.08, labels = "low = black, high = light blue", adj = 0, cex = 0.65)
  }

  dev.off()
  log_msg("PDF written: %s", pdf_path)
}

# ── 3. Read tree, dij, GC3 summary, and optional node map ────────────────────

tree <- read.tree(tree_file)

if (inherits(tree, "multiPhylo")) {
  log_msg("Multiple trees detected - using tree 1 of %d", length(tree))
  tree <- tree[[1]]
}

log_section("Tree parsing")
log_msg("Tree has %d tips and %d edges", length(tree$tip.label), nrow(tree$edge))
log_msg("First 10 tip labels: %s", paste(head(tree$tip.label, 10), collapse = ", "))
if (!is.null(tree$node.label)) {
  log_msg("Tree has %d node labels. First 20 node labels: %s",
          length(tree$node.label), paste(head(tree$node.label, 20), collapse = ", "))
} else {
  log_msg("Tree has no explicit node labels.")
}

dij <- read_dij(dij_file)
gc3_values <- read_gc3_summary(gc3_summary_file)

log_section("Input tables")
log_msg("Read %d dij rows", nrow(dij))
log_msg("Read %d GC3 node/species values", length(gc3_values))
log_msg("First 10 dij children: %s", paste(head(dij$child, 10), collapse = ", "))
log_msg("First 10 GC3 keys: %s", paste(head(names(gc3_values), 10), collapse = ", "))

ape_to_pipeline <- NULL
if (!is.na(node_numbers_file)) {
  ape_to_pipeline <- build_ape_to_pipeline_node_map(
    tree = tree,
    node_numbers_path = node_numbers_file,
    outdir = outdir,
    prefix = prefix
  )
} else {
  log_section("node_numbers mapping")
  log_msg("No --node-numbers file supplied. Internal node mapping will use tree node labels if present, otherwise ape internal node IDs.")
}

# ── 4. Replace branch lengths with d_ij ──────────────────────────────────────

edge_keys <- make_edge_child_keys(tree, ape_to_pipeline = ape_to_pipeline)

log_section("Edge key construction")
preview_n <- min(20, nrow(tree$edge))
edge_key_table_preview <- data.frame(
  edge_index = seq_len(preview_n),
  parent_ape_node = tree$edge[seq_len(preview_n), 1],
  child_ape_node = tree$edge[seq_len(preview_n), 2],
  child_key = edge_keys[seq_len(preview_n)],
  stringsAsFactors = FALSE
)
capture.output(print(edge_key_table_preview), file = log_file, append = TRUE)

log_msg("Unique child keys from tree edges: %d", length(unique(edge_keys)))
log_msg("Tree edge child keys matching dij children: %d/%d", sum(edge_keys %in% dij$child), length(edge_keys))
log_msg("Tree edge child keys matching GC3 keys: %d/%d", sum(edge_keys %in% names(gc3_values)), length(edge_keys))

dij_by_child <- dij[!duplicated(dij$child), ]
dij_lookup <- dij_by_child$dij
names(dij_lookup) <- dij_by_child$child

dup_dij_children <- dij$child[duplicated(dij$child)]
if (length(dup_dij_children) > 0) {
  log_msg("WARNING: %d duplicated dij child keys detected. First duplicates: %s",
          length(unique(dup_dij_children)), paste(head(unique(dup_dij_children), 20), collapse = ", "))
}

new_edge_lengths <- as.numeric(dij_lookup[edge_keys])
missing_dij <- is.na(new_edge_lengths)

log_section("dij branch-length matching")
log_msg("dij matches: %d/%d edges", sum(!missing_dij), length(missing_dij))

if (any(missing_dij)) {
  log_msg(
    "WARNING: %d/%d edges had no matching dij value. They will be assigned a small plotting length.",
    sum(missing_dij), length(missing_dij)
  )

  missing_dij_table <- data.frame(
    edge_index = which(missing_dij),
    parent_ape_node = tree$edge[missing_dij, 1],
    child_ape_node = tree$edge[missing_dij, 2],
    child_key = edge_keys[missing_dij],
    stringsAsFactors = FALSE
  )
  capture.output(print(head(missing_dij_table, 100)), file = log_file, append = TRUE)

  fallback <- suppressWarnings(min(new_edge_lengths[!is.na(new_edge_lengths)], na.rm = TRUE))
  if (!is.finite(fallback)) fallback <- 1
  new_edge_lengths[missing_dij] <- fallback * 0.05
}

tree$edge.length <- new_edge_lengths

# ── 5. Create branch colors from child-node GC3 values ───────────────────────

edge_colors <- color_from_gc3(gc3_values, edge_keys)

missing_gc3_vec <- is.na(gc3_values[edge_keys])
missing_gc3 <- sum(missing_gc3_vec)
matched_gc3 <- length(edge_keys) - missing_gc3

log_section("GC3 branch-color matching")
log_msg("GC3 color matches: %d/%d edges", matched_gc3, length(edge_keys))

if (missing_gc3 > 0) {
  log_msg("WARNING: %d/%d edges had no matching GC3 value and were colored grey80.", missing_gc3, length(edge_keys))

  missing_gc3_table <- data.frame(
    edge_index = which(missing_gc3_vec),
    parent_ape_node = tree$edge[missing_gc3_vec, 1],
    child_ape_node = tree$edge[missing_gc3_vec, 2],
    child_key = edge_keys[missing_gc3_vec],
    has_dij = edge_keys[missing_gc3_vec] %in% dij$child,
    stringsAsFactors = FALSE
  )
  capture.output(print(head(missing_gc3_table, 100)), file = log_file, append = TRUE)
}

log_section("Targeted diagnostic patterns")
interesting_patterns <- c("acritoscincus", "carinascincus", "anolis", "node_108", "node_109")
for (pat in interesting_patterns) {
  matching_keys <- unique(edge_keys[grepl(pat, edge_keys, ignore.case = TRUE)])
  if (length(matching_keys) > 0) {
    log_msg("Keys matching pattern '%s': %s", pat, paste(matching_keys, collapse = ", "))
    for (key in matching_keys) {
      log_msg(
        "  key=%s | in_dij=%s | in_gc3=%s | dij=%s | gc3=%s",
        key,
        key %in% dij$child,
        key %in% names(gc3_values),
        ifelse(key %in% names(dij_lookup), as.character(dij_lookup[[key]]), "NA"),
        ifelse(key %in% names(gc3_values), as.character(gc3_values[[key]]), "NA")
      )
    }
  } else {
    log_msg("No edge keys matched pattern '%s'", pat)
  }
}

# ── 6. Write outputs ─────────────────────────────────────────────────────────

newick_out <- file.path(outdir, paste0(prefix, ".dij_branch_lengths.nwk"))
pdf_out <- file.path(outdir, paste0(prefix, ".dij_branch_lengths_gc3_gradient.pdf"))
edge_table_out <- file.path(outdir, paste0(prefix, ".edge_plot_table.tsv"))

write.tree(tree, file = newick_out, digits = 10)
log_msg("Newick with d_ij branch lengths written: %s", newick_out)

edge_table <- data.frame(
  edge_index = seq_len(nrow(tree$edge)),
  parent_ape_node = tree$edge[, 1],
  child_ape_node = tree$edge[, 2],
  child_key = edge_keys,
  dij_branch_length = tree$edge.length,
  gc3_for_color = as.numeric(gc3_values[edge_keys]),
  branch_color = edge_colors,
  stringsAsFactors = FALSE
)
write.table(edge_table, file = edge_table_out, sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Edge plotting table written: %s", edge_table_out)

plot_tree_pdf(
  tree = tree,
  edge_colors = edge_colors,
  pdf_path = pdf_out,
  title_str = sprintf("Tree with d_ij branch lengths and GC3 branch colors - %d tips", length(tree$tip.label)),
  gc3_values = gc3_values
)

log_section("Output files")
log_msg("Newick: %s", newick_out)
log_msg("PDF: %s", pdf_out)
log_msg("Edge table: %s", edge_table_out)
log_msg("Log: %s", log_file)

cat("\nDone.\n")
