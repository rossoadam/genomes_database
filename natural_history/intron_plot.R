#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript plot_intron_gc_panels.R input.tsv output.pdf")
}

input_file <- args[1]
output_file <- args[2]

data <- read.delim(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

required_cols <- c("gc", "intron_length", "species_normalized")
missing_cols <- setdiff(required_cols, names(data))

if (length(missing_cols) > 0) {
  stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
}

data <- data[complete.cases(data[, required_cols]), ]

data$gc <- as.numeric(data$gc)
data$intron_length <- as.numeric(data$intron_length)

species <- c(
  "phrynosoma_platyrhinos",
  "cyclura_pinguis",
  "varanus_acanthurus"
)

species <- species[species %in% unique(data$species_normalized)]

if (length(species) != 3) {
  warning(paste("Expected 3 species but found", length(species)))
}

xlim <- range(data$intron_length, na.rm = TRUE)
ylim <- range(data$gc, na.rm = TRUE)

pretty_species_name <- function(sp) {
  if (sp == "phrynosoma_platyrhinos") {
    return(expression(italic("Phrynosoma platyrhinos")))
  } else if (sp == "cyclura_nubila") {
    return(expression(italic("Cyclura nubila")))
  } else if (sp == "varanus_komodoensis") {
    return(expression(italic("Varanus komodoensis")))
  } else {
    return(gsub("_", " ", sp))
  }
}

pdf(output_file, width = 12, height = 4)

par(
  mfrow = c(1, 3),
  mar = c(5, 5, 4, 1),
  oma = c(0, 0, 2, 0),
  cex.axis = 1.2,
  cex.lab = 1.3,
  cex.main = 1.2
)

results <- data.frame()

for (sp in species) {
  sp_data <- subset(data, species_normalized == sp)
  
  model <- lm(gc ~ intron_length, data = sp_data)
  s <- summary(model)
  
  slope <- coef(model)[2]
  intercept <- coef(model)[1]
  pval <- s$coefficients["intron_length", "Pr(>|t|)"]
  r2 <- s$r.squared
  
  plot(
    sp_data$intron_length,
    sp_data$gc,
    pch = 16,
    cex = 0.6,
    xlim = xlim,
    ylim = ylim,
    xlab = "Intron length (bp)",
    ylab = "GC content",
    main = pretty_species_name(sp)
  )
  
  abline(model, lwd = 3)
  
  legend(
    "topright",
    legend = c(
      paste0("n = ", nrow(sp_data)),
      paste0("Slope = ", signif(slope, 3)),
      paste0("R² = ", signif(r2, 3)),
      paste0("p = ", signif(pval, 3))
    ),
    bty = "n",
    cex = 0.85
  )
  
  results <- rbind(
    results,
    data.frame(
      species_normalized = sp,
      intercept = intercept,
      slope = slope,
      r_squared = r2,
      p_value = pval,
      n = nrow(sp_data)
    )
  )
}

dev.off()

results_file <- sub("\\.pdf$", "_lm_results.tsv", output_file)

write.table(
  results,
  results_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Plot written to:", output_file, "\n")
cat("Model summary written to:", results_file, "\n")