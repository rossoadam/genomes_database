#!/usr/bin/env Rscript
# 16_plot_dij_gc3_tree_v2.R
#
# Plot a rooted phylogram with d_ij branch lengths and GC3 percentage colors.
# Compatible with 15_calculate_gc_divergence_v6.py:
#   gc3_summary.tsv: node_or_species_key_short, mean_gc3_nhphyml_unweighted
#   dij_results.tsv: node_or_species_key_short, parent, dij
#   dij_results.json: { child: { parent: dij } }
#
# Input GC3 means are proportions (0..1), converted to percentages ONLY for
# plotting. d_ij is already in percentage-point units and is NEVER rescaled.
# Missing or ambiguous branches/nodes are errors; no data are fabricated.
#
# Required:
#   --tree /path/to/tree.nwk
#   --dij /path/to/dij_results.tsv  (or .json)
#   --gc3-summary /path/to/gc3_summary.tsv
#
# Recommended (required unless Newick contains validated explicit node_N IDs):
#   --node-numbers /path/to/node_numbers_rebuilt.csv
#
# Optional:
#   --outdir /path/to/output_dir
#   --prefix output_prefix
#   --gradient-steps 40
#   --tree-x-padding 1.12
#   --pdf-width 13
#   --nhphyml-dir /path/to/nhphyml_output  (deprecated; ignored)
#
# Visual defaults: rightwards tree; line width 6 for colored branches;
# black (#000000) to turquoise (#00B8C2); per-figure GC3 percentage scale.

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

usage <- function() {
  paste(
    "Usage:",
    "  Rscript 16_plot_dij_gc3_tree_v2.R \\",
    "    --tree /path/to/tree.nwk \\",
    "    --dij /path/to/dij_results.tsv \\",
    "    --gc3-summary /path/to/gc3_summary.tsv \\",
    "    --node-numbers /path/to/node_numbers_rebuilt.csv \\",
    "    --outdir /path/to/output_dir \\",
    "    [--prefix output_prefix] \\",
    "    [--gradient-steps 40] [--tree-x-padding 1.12] [--pdf-width 13]",
    "",
    "The --node-numbers option may be omitted only if every internal node",
    "has an explicit, uniquely validated node_N Newick label.",
    sep = "\n"
  )
}

parse_args <- function(argv) {
  if (length(argv) == 1L && argv[1L] %in% c("--help", "-h")) {
    cat(usage(), "\n")
    quit(save = "no", status = 0L)
  }
  allowed <- c("tree", "dij", "gc3_summary", "node_numbers", "nhphyml_dir",
               "outdir", "prefix", "gradient_steps", "tree_x_padding", "pdf_width")
  result <- list()
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[i]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key), call. = FALSE)
    name <- gsub("-", "_", sub("^--", "", key))
    if (!(name %in% allowed)) stop(sprintf("Unrecognized option: %s", key), call. = FALSE)
    if (i == length(argv) || startsWith(argv[i + 1L], "--")) {
      stop(sprintf("Missing value for %s", key), call. = FALSE)
    }
    if (!is.null(result[[name]])) stop(sprintf("Option supplied more than once: %s", key), call. = FALSE)
    result[[name]] <- argv[i + 1L]
    i <- i + 2L
  }
  needed <- c("tree", "dij", "gc3_summary")
  missing <- setdiff(needed, names(result))
  if (length(missing)) {
    stop(sprintf("Missing required argument(s): %s\n\n%s",
                 paste(paste0("--", gsub("_", "-", missing)), collapse = ", "), usage()),
         call. = FALSE)
  }
  result
}

positive_option <- function(value, default, name, min_value) {
  if (is.null(value)) return(default)
  val <- suppressWarnings(as.numeric(value))
  if (length(val) != 1L || !is.finite(val) || val < min_value) {
    stop(sprintf("--%s requires a finite number >= %s", name, min_value), call. = FALSE)
  }
  val
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
tree_file <- path.expand(args$tree)
dij_file <- path.expand(args$dij)
gc3_summary_file <- path.expand(args$gc3_summary)
node_numbers_file <- if (!is.null(args$node_numbers)) path.expand(args$node_numbers) else NULL
base_name <- tools::file_path_sans_ext(basename(tree_file))
outdir <- if (!is.null(args$outdir)) path.expand(args$outdir) else {
  file.path(dirname(tree_file), "gc3_dij_tree_plots")
}
prefix <- if (!is.null(args$prefix)) args$prefix else paste0(base_name, ".dij_gc3")
if (!nzchar(prefix) || grepl("[/\\\\]", prefix)) {
  stop("--prefix must be a nonempty filename prefix without path separators", call. = FALSE)
}
gradient_steps <- positive_option(args$gradient_steps, 40, "gradient-steps", 2)
if (gradient_steps != floor(gradient_steps)) {
  stop("--gradient-steps must be a whole number", call. = FALSE)
}
gradient_steps <- as.integer(gradient_steps)
tree_x_padding <- positive_option(args$tree_x_padding, 1.12, "tree-x-padding", 1.02)
pdf_width <- positive_option(args$pdf_width, 13, "pdf-width", 8)

if (!file.exists(tree_file)) stop(sprintf("Tree file not found: %s", tree_file), call. = FALSE)
if (!file.exists(dij_file)) stop(sprintf("d_ij file not found: %s", dij_file), call. = FALSE)
if (!file.exists(gc3_summary_file)) {
  stop(sprintf("GC3 summary file not found: %s", gc3_summary_file), call. = FALSE)
}
if (!is.null(node_numbers_file) && !file.exists(node_numbers_file)) {
  stop(sprintf("Node relationships file not found: %s", node_numbers_file), call. = FALSE)
}
if (!dir.exists(outdir)) {
  if (!dir.create(outdir, recursive = TRUE, showWarnings = FALSE)) {
    stop(sprintf("Could not create output directory: %s", outdir), call. = FALSE)
  }
}

log_file <- file.path(outdir, paste0(prefix, ".parse_diagnostics.log"))
cat("", file = log_file)  # Start a fresh log for every run.
log_msg <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = log_file, append = TRUE)
  invisible(msg)
}
log_section <- function(title) {
  log_msg("")
  log_msg("%s", paste(rep("=", nchar(title)), collapse = ""))
  log_msg("%s", title)
  log_msg("%s", paste(rep("=", nchar(title)), collapse = ""))
}
log_table <- function(x, maximum = 30L) {
  capture.output(print(utils::head(x, maximum), row.names = FALSE),
                 file = log_file, append = TRUE)
}

# Species identifiers are exact strings: no case folding, genus reconstruction,
# removal of subspecies, or silent character replacement.
coerce_node_name <- function(value) {
  if (length(value) != 1L || is.na(value)) return("")
  x <- trimws(as.character(value))
  if (is_missing(x)) return("")
  if (tolower(x) %in% c("ancestral", "ancestor", "root")) return("ancestral")
  if (grepl("^(node[ _-]*)?[0-9]+$", x, ignore.case = TRUE)) {
    digits <- sub("^node[ _-]*", "", x, ignore.case = TRUE)
    num <- suppressWarnings(as.integer(digits))
    if (is.na(num)) stop(sprintf("Invalid internal-node identifier: %s", x), call. = FALSE)
    return(paste0("node_", num))
  }
  x
}

is_internal_node <- function(x) grepl("^node_[0-9]+$", x)

is_missing <- function(x) {
  is.na(x) | tolower(trimws(as.character(x))) %in%
    c("", "na", "nan", "none", "null", "<na>", "\\n")
}

read_table_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (!(ext %in% c("tsv", "txt", "csv"))) {
    stop(sprintf("Unsupported table extension for %s (expected .tsv or .csv)", path), call. = FALSE)
  }
  sep <- if (ext == "csv") "," else "\t"
  df <- utils::read.table(path, header = TRUE, sep = sep, colClasses = "character",
                          stringsAsFactors = FALSE, check.names = FALSE,
                          na.strings = character(0), quote = '"', comment.char = "",
                          fileEncoding = "UTF-8-BOM")
  if (nrow(df) == 0L) stop(sprintf("Input table has no rows: %s", path), call. = FALSE)
  df
}

require_columns <- function(df, fields, path) {
  absent <- setdiff(fields, names(df))
  if (length(absent)) {
    stop(sprintf("Missing column(s) %s in %s; found: %s",
                 paste(absent, collapse = ", "), path, paste(names(df), collapse = ", ")),
         call. = FALSE)
  }
}

parse_number <- function(x, description) {
  if (is_missing(x)) stop(sprintf("Missing number for %s", description), call. = FALSE)
  result <- suppressWarnings(as.numeric(x))
  if (length(result) != 1L || !is.finite(result)) {
    stop(sprintf("Invalid or nonfinite number %s for %s", shQuote(as.character(x)), description),
         call. = FALSE)
  }
  result
}

read_gc3_summary <- function(path) {
  df <- read_table_auto(path)
  require_columns(df, c("node_or_species_key_short", "mean_gc3_nhphyml_unweighted"), path)
  keys <- vapply(df$node_or_species_key_short, coerce_node_name, character(1))
  if (any(!nzchar(keys))) stop(sprintf("Blank GC3 node/species identifier in %s", path), call. = FALSE)
  if (anyDuplicated(keys)) {
    stop(sprintf("Duplicate GC3 identifier(s) after node-alias handling: %s",
                 paste(unique(keys[duplicated(keys)]), collapse = ", ")), call. = FALSE)
  }
  props <- vapply(seq_along(keys), function(i) {
    parse_number(df$mean_gc3_nhphyml_unweighted[i],
                 sprintf("GC3 for %s in %s", keys[i], path))
  }, numeric(1))
  if (any(props < 0 | props > 1)) {
    problem <- which(props < 0 | props > 1)
    stop(sprintf("GC3 must be a proportion (0..1), but %s has %s in %s",
                 keys[problem[1]], format(props[problem[1]]), path), call. = FALSE)
  }
  if (!("ancestral" %in% keys)) {
    stop(sprintf("Missing required 'ancestral' GC3 row in %s", path), call. = FALSE)
  }
  result <- props * 100
  names(result) <- keys
  log_msg("Read %d GC3 records; converted proportions to percentages (x100)", length(result))
  result
}

read_dij <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) {
    df <- read_table_auto(path)
    require_columns(df, c("node_or_species_key_short", "parent", "dij"), path)
    raw_child <- df$node_or_species_key_short
    raw_parent <- df$parent
    raw_value <- df$dij
  } else if (ext == "json") {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Reading d_ij JSON requires the R package jsonlite (or provide dij_results.tsv)",
           call. = FALSE)
    }
    obj <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    if (!is.list(obj) || is.null(names(obj))) {
      stop(sprintf("Expected nested child-to-parent object in JSON: %s", path), call. = FALSE)
    }
    raw_child <- character()
    raw_parent <- character()
    raw_value <- character()
    for (child in names(obj)) {
      parents <- obj[[child]]
      if (is.list(parents) && length(parents) == 0L) next
      if (!is.list(parents) || is.null(names(parents))) {
        stop(sprintf("Malformed d_ij JSON object for child %s", child), call. = FALSE)
      }
      for (parent in names(parents)) {
        val <- parents[[parent]]
        raw_child <- c(raw_child, child)
        raw_parent <- c(raw_parent, parent)
        raw_value <- c(raw_value, if (is.null(val)) "\\N" else as.character(val))
      }
    }
  } else {
    stop(sprintf("Unsupported d_ij file extension: %s (.tsv, .txt, or .json expected)", path),
         call. = FALSE)
  }
  if (!length(raw_child)) stop(sprintf("No d_ij relationships found in %s", path), call. = FALSE)
  child <- vapply(raw_child, coerce_node_name, character(1))
  parent <- vapply(raw_parent, coerce_node_name, character(1))
  if (any(!nzchar(child) | !nzchar(parent))) {
    bad <- which(!nzchar(child) | !nzchar(parent))[1]
    stop(sprintf("Empty d_ij child/parent in record %d of %s", bad, path), call. = FALSE)
  }
  pair_key <- paste(parent, child, sep = "\r")
  if (anyDuplicated(pair_key)) {
    bad <- which(duplicated(pair_key))[1]
    stop(sprintf("Duplicate d_ij parent-child pair %s -> %s in %s",
                 parent[bad], child[bad], path), call. = FALSE)
  }
  values <- vapply(seq_along(child), function(i) {
    parse_number(raw_value[i], sprintf("d_ij %s -> %s in %s", parent[i], child[i], path))
  }, numeric(1))
  if (any(values < 0)) {
    bad <- which(values < 0)[1]
    stop(sprintf("Negative d_ij for %s -> %s: %.12g", parent[bad], child[bad], values[bad]),
         call. = FALSE)
  }
  data.frame(child = child, parent = parent, dij = values, stringsAsFactors = FALSE)
}

read_node_relationships <- function(path) {
  df <- read_table_auto(path)
  column_names <- tolower(trimws(names(df)))
  child_names <- c("daughter", "daughter_node", "child", "child_node", "node",
                   "node_id", "descendant", "descendant_node", "son", "son_node")
  parent_names <- c("parent", "parent_node", "ancestor", "ancestor_node",
                    "mother", "parent_id")
  child_idx <- match(TRUE, column_names %in% child_names)
  parent_idx <- match(TRUE, column_names %in% parent_names)
  child <- character()
  parent <- character()
  if (!is.na(child_idx) && !is.na(parent_idx) && child_idx != parent_idx) {
    child <- vapply(df[[child_idx]], coerce_node_name, character(1))
    parent <- vapply(df[[parent_idx]], coerce_node_name, character(1))
    keep <- nzchar(child) & nzchar(parent)
    if (any(!keep)) {
      log_msg("NOTE: node_numbers contains %d incomplete parent-child row(s), ignored (possible root entries)",
              sum(!keep))
    }
    child <- child[keep]
    parent <- parent[keep]
  } else {
    # Wide format: each column names one child; nonempty cells give its parents.
    for (name in names(df)) {
      this_child <- coerce_node_name(name)
      if (!nzchar(this_child)) next
      these_parents <- vapply(df[[name]], coerce_node_name, character(1))
      these_parents <- these_parents[nzchar(these_parents)]
      child <- c(child, rep(this_child, length(these_parents)))
      parent <- c(parent, these_parents)
    }
  }
  edges <- unique(data.frame(child = child, parent = parent, stringsAsFactors = FALSE))
  if (!nrow(edges)) stop(sprintf("No parent-child relationships in %s", path), call. = FALSE)
  if (any(edges$child == edges$parent)) {
    bad <- which(edges$child == edges$parent)[1]
    stop(sprintf("Self-parent relationship for %s in %s", edges$child[bad], path), call. = FALSE)
  }
  edges
}

signature_from_tips <- function(tips) paste(sort(unique(tips)), collapse = "|")

descendant_tips_from_tree_node <- function(tree, ape_node) {
  n_tips <- length(tree$tip.label)
  descendants <- phytools::getDescendants(tree, ape_node)
  tips <- descendants[descendants <= n_tips]
  # No name rewriting: these are the exact Newick labels.
  sort(unique(tree$tip.label[tips]))
}

build_pipeline_descendant_signatures <- function(node_edges) {
  children_by_parent <- split(node_edges$child, node_edges$parent)
  memo <- new.env(parent = emptyenv())
  visiting <- new.env(parent = emptyenv())

  collect_tips <- function(node) {
    if (exists(node, envir = memo, inherits = FALSE)) {
      return(get(node, envir = memo, inherits = FALSE))
    }
    if (exists(node, envir = visiting, inherits = FALSE)) {
      stop(sprintf("Cycle detected in node relationships involving %s", node), call. = FALSE)
    }
    if (!is_internal_node(node) && node != "ancestral") {
      tips <- node
    } else {
      assign(node, TRUE, envir = visiting)
      kids <- children_by_parent[[node]]
      if (is.null(kids) || !length(kids)) {
        tips <- character()
      } else {
        tips <- sort(unique(unlist(lapply(kids, collect_tips), use.names = FALSE)))
      }
      rm(list = node, envir = visiting)
    }
    assign(node, tips, envir = memo)
    tips
  }

  # 'ancestral' can be an extra parent of the true root. It is an annotation,
  # NOT a candidate for the plotted APE root, even with the same tip set.
  internal_nodes <- sort(unique(c(node_edges$parent[is_internal_node(node_edges$parent)],
                                  node_edges$child[is_internal_node(node_edges$child)])))
  if (!length(internal_nodes)) stop("node_numbers contains no node_N internal nodes", call. = FALSE)
  tips <- lapply(internal_nodes, collect_tips)
  data.frame(pipeline_node = internal_nodes,
             descendant_signature = vapply(tips, signature_from_tips, character(1)),
             n_descendant_tips_pipeline = lengths(tips), stringsAsFactors = FALSE)
}

write_node_map <- function(map, outdir, prefix) {
  path <- file.path(outdir, paste0(prefix, ".ape_to_pipeline_node_map.tsv"))
  utils::write.table(map, file = path, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "\\N")
  log_msg("Internal-node mapping table written: %s", path)
  invisible(path)
}

build_ape_to_pipeline_node_map <- function(tree, node_edges, outdir, prefix) {
  n_tips <- length(tree$tip.label)
  ape_nodes <- seq.int(n_tips + 1L, n_tips + tree$Nnode)
  ape_sigs <- vapply(ape_nodes, function(node) {
    signature_from_tips(descendant_tips_from_tree_node(tree, node))
  }, character(1))
  ape_n_tips <- vapply(ape_nodes, function(node) {
    length(descendant_tips_from_tree_node(tree, node))
  }, integer(1))

  if (!is.null(node_edges)) {
    pipeline <- build_pipeline_descendant_signatures(node_edges)
    pipeline <- pipeline[pipeline$n_descendant_tips_pipeline > 0L, , drop = FALSE]
    if (anyDuplicated(pipeline$descendant_signature)) {
      dup <- unique(pipeline$descendant_signature[duplicated(pipeline$descendant_signature)])
      offending <- pipeline$pipeline_node[pipeline$descendant_signature %in% dup]
      stop(sprintf("Ambiguous node_numbers: internal nodes have identical descendant-tip sets: %s",
                   paste(offending, collapse = ", ")), call. = FALSE)
    }
    matched <- match(ape_sigs, pipeline$descendant_signature)
    map <- data.frame(
      ape_node = ape_nodes, descendant_signature = ape_sigs,
      n_descendant_tips_ape = ape_n_tips,
      pipeline_node = pipeline$pipeline_node[matched],
      n_descendant_tips_pipeline = pipeline$n_descendant_tips_pipeline[matched],
      stringsAsFactors = FALSE
    )
  } else {
    labels <- tree$node.label
    if (is.null(labels) || length(labels) != tree$Nnode) {
      stop("Without --node-numbers, all Newick internal nodes need explicit node_N labels",
           call. = FALSE)
    }
    # Numeric-only labels (e.g., bootstrap 100) are NOT reliable node IDs.
    valid_labels <- !is.na(labels) & grepl("^node[ _-]*[0-9]+$", labels, ignore.case = TRUE)
    if (!all(valid_labels)) {
      bad <- which(!valid_labels)[1]
      stop(sprintf("Cannot use unlabeled/bootstrap-labeled internal node %d (label %s); supply --node-numbers",
                   ape_nodes[bad], shQuote(labels[bad])), call. = FALSE)
    }
    pipeline_nodes <- vapply(labels, coerce_node_name, character(1))
    map <- data.frame(
      ape_node = ape_nodes, descendant_signature = ape_sigs,
      n_descendant_tips_ape = ape_n_tips,
      pipeline_node = pipeline_nodes, n_descendant_tips_pipeline = ape_n_tips,
      stringsAsFactors = FALSE
    )
  }
  map$mapping_status <- ifelse(is.na(map$pipeline_node), "unmatched", "matched")
  write_node_map(map, outdir, prefix)
  log_section("Internal-node mapping")
  log_msg("Mapped internal nodes: %d/%d", sum(!is.na(map$pipeline_node)), nrow(map))
  if (anyNA(map$pipeline_node)) {
    log_table(map[is.na(map$pipeline_node), , drop = FALSE])
    stop(sprintf("%d internal Newick nodes could not be matched to node_numbers by descendant species",
                 sum(is.na(map$pipeline_node))), call. = FALSE)
  }
  if (anyDuplicated(map$pipeline_node)) {
    log_table(map[duplicated(map$pipeline_node) | duplicated(map$pipeline_node, fromLast = TRUE), , drop = FALSE])
    stop("Ambiguous internal node mapping (multiple APE nodes map to the same pipeline node)",
         call. = FALSE)
  }
  lookup <- map$pipeline_node
  names(lookup) <- as.character(map$ape_node)
  lookup
}

make_node_key <- function(tree, ape_node, ape_to_pipeline) {
  n_tips <- length(tree$tip.label)
  if (ape_node <= n_tips) return(tree$tip.label[ape_node])
  key <- ape_to_pipeline[as.character(ape_node)]
  if (length(key) != 1L || is.na(key) || !nzchar(key)) {
    stop(sprintf("Missing mapped pipeline ID for APE node %d", ape_node), call. = FALSE)
  }
  unname(key)
}

make_edge_keys <- function(tree, ape_to_pipeline) {
  data.frame(
    parent_key = vapply(tree$edge[, 1], function(n) make_node_key(tree, n, ape_to_pipeline), character(1)),
    child_key = vapply(tree$edge[, 2], function(n) make_node_key(tree, n, ape_to_pipeline), character(1)),
    stringsAsFactors = FALSE
  )
}

validate_edges <- function(tree, edge_keys, dij, gc3_values, node_edges) {
  log_section("Exact parent-child edge matching")
  plotted_pairs <- paste(edge_keys$parent_key, edge_keys$child_key, sep = "\r")
  source_pairs <- paste(dij$parent, dij$child, sep = "\r")
  log_msg("Tree edges: %d; d_ij relationships: %d", nrow(tree$edge), nrow(dij))
  if (anyDuplicated(plotted_pairs)) {
    bad <- which(duplicated(plotted_pairs))[1]
    stop(sprintf("Duplicate tree edge %s -> %s", edge_keys$parent_key[bad],
                 edge_keys$child_key[bad]), call. = FALSE)
  }
  if (!is.null(node_edges)) {
    expected_pairs <- paste(node_edges$parent, node_edges$child, sep = "\r")
    absent <- which(!(plotted_pairs %in% expected_pairs))
    log_msg("Tree edges consistent with node_numbers: %d/%d", nrow(tree$edge) - length(absent), nrow(tree$edge))
    if (length(absent)) {
      problem <- data.frame(edge_index = absent, parent = edge_keys$parent_key[absent],
                            child = edge_keys$child_key[absent])
      log_table(problem)
      stop(sprintf("Tree edge %s -> %s (edge %d) is absent from node_numbers relationships",
                   problem$parent[1], problem$child[1], problem$edge_index[1]), call. = FALSE)
    }
  }
  matches <- match(plotted_pairs, source_pairs)
  missing <- which(is.na(matches))
  log_msg("Exact d_ij matches: %d/%d", length(matches) - length(missing), length(matches))
  if (length(missing)) {
    problem <- data.frame(edge_index = missing,
                          parent_ape_node = tree$edge[missing, 1],
                          child_ape_node = tree$edge[missing, 2],
                          parent_key = edge_keys$parent_key[missing],
                          child_key = edge_keys$child_key[missing])
    log_table(problem)
    stop(sprintf("Missing d_ij for edge %d: %s -> %s (no fallback branch lengths allowed)",
                 missing[1], problem$parent_key[1], problem$child_key[1]), call. = FALSE)
  }
  extras <- which(!(source_pairs %in% plotted_pairs))
  log_msg("Unused d_ij relationships: %d", length(extras))
  if (length(extras)) {
    log_msg("NOTE: Extra relationships do not add branches to the supplied Newick tree.")
    log_table(dij[extras, , drop = FALSE])
  }

  # GC3 is required at BOTH endpoints, not only the child.
  required_nodes <- unique(c(edge_keys$parent_key, edge_keys$child_key))
  missing_gc3 <- setdiff(required_nodes, names(gc3_values))
  log_msg("Required GC3 node/species values matched: %d/%d",
          length(required_nodes) - length(missing_gc3), length(required_nodes))
  if (length(missing_gc3)) {
    log_msg("Missing GC3 identifiers: %s", paste(missing_gc3, collapse = ", "))
    stop(sprintf("Missing required GC3 for node/species: %s", missing_gc3[1]), call. = FALSE)
  }
  log_msg("Child GC3 matches: %d/%d", nrow(tree$edge), nrow(tree$edge))
  log_msg("Parent GC3 matches: %d/%d", nrow(tree$edge), nrow(tree$edge))
  unused_gc3 <- setdiff(names(gc3_values), c(required_nodes, "ancestral"))
  log_msg("GC3 node/species records unused by plotted edges (excluding ancestral): %d", length(unused_gc3))
  if (length(unused_gc3)) {
    log_msg("Unused GC3 IDs (first 30): %s", paste(head(unused_gc3, 30), collapse = ", "))
  }
  lengths <- dij$dij[matches]
  if (any(!is.finite(lengths) | lengths < 0)) {
    stop("A required d_ij branch length is nonfinite or negative", call. = FALSE)
  }
  log_msg("Zero-length branches (preserved): %d", sum(lengths == 0))
  list(branch_lengths = lengths, unused_dij = extras, plotted_gc3 = gc3_values[required_nodes])
}


# Preserve the original script's focused edge diagnostics. These are reports
# only: matching and plot calculations never depend on these taxa/node names.
log_targeted_patterns <- function(edge_keys, dij, gc3_values) {
  log_section("Targeted diagnostic patterns (informational)")
  patterns <- c("acritoscincus", "carinascincus", "anolis", "node_108", "node_109")
  pair_keys <- paste(dij$parent, dij$child, sep = "\r")
  for (pat in patterns) {
    idx <- which(grepl(pat, edge_keys$child_key, ignore.case = TRUE) |
                 grepl(pat, edge_keys$parent_key, ignore.case = TRUE))
    if (!length(idx)) {
      log_msg("No plotted edges matched pattern '%s'", pat)
      next
    }
    log_msg("Edges matching '%s': %s", pat, paste(head(idx, 20), collapse = ", "))
    for (i in head(idx, 20)) {
      parent <- edge_keys$parent_key[i]
      child <- edge_keys$child_key[i]
      dij_index <- match(paste(parent, child, sep = "\r"), pair_keys)
      log_msg("  edge=%d | parent=%s | child=%s | parent_GC3=%.4f%% | child_GC3=%.4f%% | d_ij=%.8g",
              i, parent, child, unname(gc3_values[parent]), unname(gc3_values[child]),
              dij$dij[dij_index])
    }
  }
}

make_color_function <- function(min_gc, max_gc) {
  pal <- grDevices::colorRampPalette(c("black", "#00B8C2"))(100)
  function(x) {
    if (any(!is.finite(x))) stop("Invalid GC3 provided to color interpolation", call. = FALSE)
    if (isTRUE(all.equal(min_gc, max_gc))) return(rep(pal[50], length(x)))
    scaled <- pmin(1, pmax(0, (x - min_gc) / (max_gc - min_gc)))
    pal[pmax(1L, pmin(100L, round(scaled * 99) + 1L))]
  }
}

draw_gradient_branch <- function(x0, x1, y, parent_gc3, child_gc3,
                                 color_fun, steps = 40L, lwd = 6) {
  xs <- seq(x0, x1, length.out = steps + 1L)
  values <- seq(parent_gc3, child_gc3, length.out = steps)
  cols <- color_fun(values)
  for (i in seq_len(steps)) {
    graphics::segments(xs[i], y, xs[i + 1L], y, col = cols[i],
                       lwd = lwd, lend = "round")
  }
}

draw_vertical_connectors <- function(tree, last_plot, edge_keys,
                                     gc3_values, color_fun, lwd = 6) {
  n_tips <- length(tree$tip.label)
  nodes <- unique(tree$edge[, 1])
  nodes <- nodes[nodes > n_tips]
  for (node in nodes) {
    children <- which(tree$edge[, 1] == node)
    if (length(children) < 2L) next
    x <- last_plot$xx[node]
    ys <- last_plot$yy[tree$edge[children, 2]]
    key <- edge_keys$parent_key[children[1]]
    col <- color_fun(unname(gc3_values[key]))
    graphics::segments(x, min(ys), x, max(ys), col = col, lwd = lwd,
                       lend = "round")
  }
}

plot_tree_pdf <- function(tree, edge_keys, path, gc3_values, plotted_gc3,
                          ancestral_gc_percent, steps = 40L,
                          tree_x_padding = 1.12, pdf_width = 13) {
  n_tips <- length(tree$tip.label)
  fig_h <- max(8, n_tips * 0.19)
  fnt_size <- min(1.1, 200 / n_tips)
  min_gc <- min(plotted_gc3)
  max_gc <- max(plotted_gc3)
  color_fun <- make_color_function(min_gc, max_gc)
  log_msg("Color range for plotted nodes: %.6f%% to %.6f%% (independent per figure)",
          min_gc, max_gc)
  log_msg("Ancestral GC3 displayed: %.6f%%", ancestral_gc_percent)

  grDevices::pdf(path, width = pdf_width, height = fig_h)
  device_open <- TRUE
  on.exit(if (device_open) grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(c(1, 2), nrow = 1L), widths = c(5.4, 1))
  graphics::par(mar = c(1, 0.3, 3, 0.05), xpd = NA)
  tree_depth <- max(ape::node.depth.edgelength(tree))
  # An all-zero tree stays all-zero: provide nonzero plotting extent without
  # adding any artificial evolutionary distance to the edges.
  extent <- if (tree_depth > 0) tree_depth else 1
  if (tree_depth == 0) log_msg("WARNING: All tree edges have zero d_ij; using display-only x-axis extent 1.")
  ape::plot.phylo(tree, type = "phylogram", direction = "rightwards",
                  show.tip.label = TRUE, cex = fnt_size,
                  edge.width = 3, edge.color = "white", no.margin = FALSE,
                  label.offset = extent * 0.002,
                  x.lim = c(0, extent * tree_x_padding))
  last_plot <- get("last_plot.phylo", envir = getFromNamespace(".PlotPhyloEnv", "ape"))
  draw_vertical_connectors(tree, last_plot, edge_keys, gc3_values, color_fun, lwd = 6)
  for (i in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    draw_gradient_branch(
      x0 = last_plot$xx[parent], x1 = last_plot$xx[child],
      y = last_plot$yy[child],
      parent_gc3 = unname(gc3_values[edge_keys$parent_key[i]]),
      child_gc3 = unname(gc3_values[edge_keys$child_key[i]]),
      color_fun = color_fun, steps = steps, lwd = 6
    )
  }

  graphics::par(mar = c(3, 0, 3, 0.3), xpd = NA)
  graphics::plot.new()
  pal <- if (isTRUE(all.equal(min_gc, max_gc))) {
    rep(color_fun(min_gc), 100)
  } else {
    grDevices::colorRampPalette(c("black", "#00B8C2"))(100)
  }
  ys <- seq(0.15, 0.85, length.out = 100)
  bar_x0 <- 0.37
  bar_x1 <- 0.68
  text_x <- 0.37
  for (i in seq_len(99)) {
    graphics::rect(bar_x0, ys[i], bar_x1, ys[i + 1L], col = pal[i], border = NA)
  }
  graphics::text(text_x, 0.86, labels = sprintf("%.2f%%", max_gc), adj = 0, cex = 1.5)
  graphics::text(text_x, 0.14, labels = sprintf("%.2f%%", min_gc), adj = 0, cex = 1.5)
  graphics::text(bar_x0, 0.9, labels = "GC3 (%)", adj = 0, cex = 2, font = 2)
  # Retain original annotation placement relative to the legend panel.
  graphics::text(-5.75, 0.98,
                 labels = sprintf("Ancestral GC3 = %.2f%%", ancestral_gc_percent),
                 adj = 0, cex = 1.5, font = 2)
  grDevices::dev.off()
  device_open <- FALSE
  log_msg("PDF written: %s", path)
}

write_outputs <- function(tree, edge_keys, gc3_values, plotted_gc3, ancestral_gc,
                          gradient_steps, tree_x_padding, pdf_width) {
  newick_out <- file.path(outdir, paste0(prefix, ".dij_branch_lengths.nwk"))
  pdf_out <- file.path(outdir, paste0(prefix, ".dij_branch_lengths_gc3_gradient.pdf"))
  edge_table_out <- file.path(outdir, paste0(prefix, ".edge_plot_table.tsv"))
  # More precision than v1 (digits=10). Never alter the source d_ij lengths.
  ape::write.tree(tree, file = newick_out, digits = 16)
  written <- ape::read.tree(newick_out)
  if (inherits(written, "multiPhylo")) stop("Serialized Newick unexpectedly contains multiple trees", call. = FALSE)
  if (!setequal(written$tip.label, tree$tip.label) ||
      length(written$edge.length) != length(tree$edge.length) ||
      any(abs(sort(written$edge.length) - sort(tree$edge.length)) > 1e-10)) {
    stop(sprintf("Newick branch-length serialization verification failed: %s", newick_out),
         call. = FALSE)
  }
  log_msg("Newick with verified, unscaled d_ij branch lengths: %s", newick_out)

  edge_table <- data.frame(
    edge_index = seq_len(nrow(tree$edge)),
    parent_ape_node = tree$edge[, 1],
    child_ape_node = tree$edge[, 2],
    parent_key = edge_keys$parent_key,
    child_key = edge_keys$child_key,
    dij_branch_length = tree$edge.length,
    parent_gc3_for_gradient_start = unname(gc3_values[edge_keys$parent_key]),
    child_gc3_for_gradient_end = unname(gc3_values[edge_keys$child_key]),
    stringsAsFactors = FALSE
  )
  utils::write.table(edge_table, file = edge_table_out, sep = "\t",
                     quote = FALSE, row.names = FALSE)
  log_msg("Edge plotting table written: %s", edge_table_out)
  log_msg("Both edge-table GC3 gradient endpoint columns are PERCENTAGES (0..100).")

  plot_tree_pdf(tree, edge_keys, pdf_out, gc3_values, plotted_gc3,
                ancestral_gc_percent = ancestral_gc,
                steps = gradient_steps, tree_x_padding = tree_x_padding,
                pdf_width = pdf_width)
  log_section("Output files")
  log_msg("Newick: %s", newick_out)
  log_msg("PDF: %s", pdf_out)
  log_msg("Edge table: %s", edge_table_out)
  log_msg("Node map: %s", file.path(outdir, paste0(prefix, ".ape_to_pipeline_node_map.tsv")))
  log_msg("Log: %s", log_file)
}

run_main <- function() {
  log_msg("Script: 16_plot_dij_gc3_tree_v2.R")
  log_msg("Output directory: %s", normalizePath(outdir, mustWork = FALSE))
  log_msg("Tree input: %s", tree_file)
  log_msg("d_ij input: %s", dij_file)
  log_msg("GC3 summary input: %s", gc3_summary_file)
  log_msg("node_numbers input: %s", if (is.null(node_numbers_file)) "not supplied" else node_numbers_file)
  log_msg("Gradient steps per branch: %d", gradient_steps)
  log_msg("Tree x padding: %.3f; PDF width: %.3f", tree_x_padding, pdf_width)
  if (!is.null(args$nhphyml_dir)) {
    log_msg("DEPRECATED: --nhphyml-dir is ignored; ancestral GC3 comes from gc3_summary.tsv")
  }
  tree <- ape::read.tree(tree_file)
  if (inherits(tree, "multiPhylo")) {
    log_msg("WARNING: Multiple trees supplied; using first of %d (original script behavior)", length(tree))
    tree <- tree[[1]]
  }
  if (!inherits(tree, "phylo") || is.null(tree$edge) || !nrow(tree$edge)) {
    stop(sprintf("No valid phylogenetic tree read from %s", tree_file), call. = FALSE)
  }
  if (anyNA(tree$tip.label) || any(!nzchar(tree$tip.label)) || anyDuplicated(tree$tip.label)) {
    stop("The Newick tree must contain unique, nonempty terminal species labels", call. = FALSE)
  }
  if (!ape::is.rooted(tree)) {
    stop("The Newick tree must be rooted to interpret parent-child d_ij relationships", call. = FALSE)
  }
  if (!is.null(tree$root.edge)) {
    stop("Input Newick contains a root stem edge (root.edge) with no corresponding plotted parent-child d_ij; remove or resolve this upstream rather than retaining an unverified length",
         call. = FALSE)
  }
  log_section("Tree parsing")
  log_msg("Terminal species: %d; internal nodes: %d; edges: %d",
          length(tree$tip.label), tree$Nnode, nrow(tree$edge))
  log_msg("First 10 original tip labels: %s", paste(head(tree$tip.label, 10), collapse = ", "))
  log_msg("Input Newick branch lengths are replaced, not reused")

  dij <- read_dij(dij_file)
  gc3_values <- read_gc3_summary(gc3_summary_file)
  log_section("Input tables")
  log_msg("d_ij records: %d", nrow(dij))
  log_msg("GC3 node/species records: %d", length(gc3_values))
  log_msg("GC3 unit conversion: raw proportions x 100; d_ij unchanged")
  missing_tips <- setdiff(tree$tip.label, names(gc3_values))
  log_msg("Exact tree-tip GC3 matches: %d/%d",
          length(tree$tip.label) - length(missing_tips), length(tree$tip.label))
  if (length(missing_tips)) {
    log_msg("Missing exact GC3 species labels: %s", paste(missing_tips, collapse = ", "))
    stop(sprintf("Species-name mismatch: Newick tip %s not in GC3 summary", missing_tips[1]),
         call. = FALSE)
  }

  node_edges <- if (!is.null(node_numbers_file)) read_node_relationships(node_numbers_file) else NULL
  if (!is.null(node_edges)) log_msg("Read %d unique node_numbers relationships", nrow(node_edges))
  ape_to_pipeline <- build_ape_to_pipeline_node_map(tree, node_edges, outdir, prefix)
  edge_keys <- make_edge_keys(tree, ape_to_pipeline)
  validated <- validate_edges(tree, edge_keys, dij, gc3_values, node_edges)
  tree$edge.length <- validated$branch_lengths
  ancestral_gc <- unname(gc3_values["ancestral"])
  log_targeted_patterns(edge_keys, dij, gc3_values)

  log_section("Plotting inputs")
  log_msg("GC3 percentages on plotted edges: %.6f%% to %.6f%%",
          min(validated$plotted_gc3), max(validated$plotted_gc3))
  log_msg("Ancestral GC3 from gc3_summary.tsv: %.6f%%", ancestral_gc)
  log_msg("d_ij branch lengths: min %.12g, max %.12g (not rescaled)",
          min(tree$edge.length), max(tree$edge.length))
  write_outputs(tree, edge_keys, gc3_values, validated$plotted_gc3,
                ancestral_gc, gradient_steps, tree_x_padding, pdf_width)
  log_msg("DONE: all plotted d_ij branches and GC3 endpoints validated; no fallback data used")
}

tryCatch(run_main(), error = function(e) {
  log_msg("ERROR: %s", conditionMessage(e))
  stop(e)
})
