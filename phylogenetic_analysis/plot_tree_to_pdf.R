#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(phytools)
})

option_list <- list(
  make_option(c("--tree"), type = "character", help = "Path to input tree file in Newick/Nexus format."),
  make_option(c("--output"), type = "character", help = "Path to output PDF file."),
  make_option(c("--format"), type = "character", default = "newick",
              help = "Tree format: newick or nexus [default: %default]"),
  make_option(c("--layout"), type = "character", default = "phylogram",
              help = "Layout: phylogram, cladogram, fan, radial [default: %default]"),
  make_option(c("--width"), type = "double", default = 12,
              help = "PDF width in inches [default: %default]"),
  make_option(c("--height"), type = "double", default = 16,
              help = "PDF height in inches [default: %default]"),
  make_option(c("--cex"), type = "double", default = 0.5,
              help = "Tip label size [default: %default]"),
  make_option(c("--label_offset"), type = "double", default = 0,
              help = "Offset for tip labels [default: %default]"),
  make_option(c("--line_width"), type = "double", default = 1,
              help = "Branch line width [default: %default]"),
  make_option(c("--no_labels"), action = "store_true", default = FALSE,
              help = "Suppress tip labels."),
  make_option(c("--show_node_labels"), action = "store_true", default = FALSE,
              help = "Plot node labels if present."),
  make_option(c("--ladderize"), action = "store_true", default = FALSE,
              help = "Ladderize the tree before plotting."),
  make_option(c("--rightwards"), action = "store_true", default = TRUE,
              help = "Plot rightwards instead of leftwards for rectangular layouts."),
  make_option(c("--main"), type = "character", default = NULL,
              help = "Optional plot title."),
  make_option(c("--use_phytools"), action = "store_true", default = FALSE,
              help = "Use phytools::plotTree instead of ape::plot.phylo for rectangular trees.")
)

parser <- OptionParser(
  usage = "%prog --tree input.nwk --output tree.pdf [options]",
  option_list = option_list,
  description = "Plot a phylogenetic tree to PDF using ape/phytools."
)
opt <- parse_args(parser)

if (is.null(opt$tree) || is.null(opt$output)) {
  print_help(parser)
  stop("Both --tree and --output are required.", call. = FALSE)
}

read_tree_auto <- function(path, fmt) {
  fmt <- tolower(fmt)
  if (fmt == "newick") {
    return(read.tree(path))
  }
  if (fmt == "nexus") {
    return(read.nexus(path))
  }
  stop("Unsupported --format. Use 'newick' or 'nexus'.", call. = FALSE)
}

validate_layout <- function(layout) {
  allowed <- c("phylogram", "cladogram", "fan", "radial")
  if (!(layout %in% allowed)) {
    stop(sprintf("Unsupported --layout '%s'. Choose one of: %s",
                 layout, paste(allowed, collapse = ", ")), call. = FALSE)
  }
}

layout <- tolower(opt$layout)
validate_layout(layout)

tree <- read_tree_auto(opt$tree, opt$format)

if (!inherits(tree, "phylo")) {
  stop("The input file did not produce a valid 'phylo' object.", call. = FALSE)
}

if (opt$ladderize) {
  tree <- ladderize(tree)
}

n_tips <- length(tree$tip.label)
n_nodes <- tree$Nnode

message("Loaded tree: ", normalizePath(opt$tree, mustWork = FALSE))
message("Tips: ", n_tips)
message("Internal nodes: ", n_nodes)
message("Writing PDF: ", normalizePath(opt$output, mustWork = FALSE))

pdf(opt$output, width = opt$width, height = opt$height, onefile = TRUE)
on.exit(dev.off(), add = TRUE)

par(mar = c(1, 1, 3, 1))

show_tip_labels <- !isTRUE(opt$no_labels)
show_node_labels <- isTRUE(opt$show_node_labels)
plot_direction <- if (isTRUE(opt$rightwards)) "rightwards" else "leftwards"

if (isTRUE(opt$use_phytools) && layout %in% c("phylogram", "cladogram")) {
  phytools::plotTree(
    tree,
    type = layout,
    ftype = if (show_tip_labels) "i" else "off",
    fsize = opt$cex,
    lwd = opt$line_width,
    direction = plot_direction,
    offset = opt$label_offset,
    main = opt$main
  )
} else {
  ape::plot.phylo(
    tree,
    type = layout,
    use.edge.length = TRUE,
    show.tip.label = show_tip_labels,
    cex = opt$cex,
    label.offset = opt$label_offset,
    direction = plot_direction,
    no.margin = FALSE,
    main = opt$main,
    edge.width = opt$line_width
  )
}

if (show_node_labels && !is.null(tree$node.label)) {
  ape::nodelabels(tree$node.label, frame = "none", cex = max(0.4, opt$cex * 0.8))
}

message("Done.")
