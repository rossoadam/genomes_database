# Load necessary libraries
# install.packages('geiger')
# install.packages('nlme')
# install.packages('phytools')
# install.packages('ape')
# install.packages("MuMIn")  # only if not installed
library('MuMIn')
library('phytools')
library('ape')
library('geiger')
library('nlme')
library('rcompanion')

#X combine_dir = "~/projects/mass_predicts_dna_dynamics/"
tree_path = "~/projects/genomes/records/compleasm/alignments/t1.0_e1_o0_with_s_punctatus/06b_rooted_tree_for_pgls/rev_dna_renamed_keep_lengths.rooted_renamed_pgls_altered.nwk"
master_output_path = "~/projects/genomes/records/compleasm/alignments/t1.0_e1_o0_with_s_punctatus/combined_natural_history_ortholog_intron.tsv"
outdir= "~/projects/genomes/records/compleasm/alignments/t1.0_e1_o0_with_s_punctatus/figures"
#gc4_results="~/projects/genomes/records/compleasm/alignments/t1.0_e1_o26_with_s_punctatus/08_final_outputs/gc_metrics_outputs_wtih_thermal/gc4_species_summary.tsv"
setwd(outdir)
list.files(combine_dir)

# Load my phylogenetic tree
tree <- read.tree(tree_path)
class(tree)
print(tree)

# Load my data
#X data_with_anc <- read.csv("combined_data.csv")
data <- read.csv(master_output_path, sep="\t")
#gc4_results <- read.csv(gc4_results, sep="\t")
# omit the ancestor row - not necessary for the new data
#X data <- data_with_anc[c(1:43),]
# Prepare the data for PGLS, by setting the rownames to genus_species
#X rownames(data) <- data$genus_species
rownames(data) <- data$species_normalized
#rownames(gc4_results) <- gc4_results$sequence_id
# this will specify the order of taxa in my data, this important for matching tree tips to data rows
order_ <- rownames(data)
#order_2 <- rownames(gc4_results)
#data$gc4 <- gc4_results$gc4
class(data$mass_preferred)
class(data$genome_size)
class(data$mean_intron_gc)
mass_vector <- setNames(data$mass_preferred,order_)
class(mass_vector)
data$log_mass = log(as.numeric(mass_vector), base = 10)
data$ct_max <- as.numeric(data$ct_max)
data$genome_size <- as.numeric(data$genome_size)
class(data$log_mass)
print(data$log_mass)
# Ensure the species names in the data match the tip labels in the tree
plot(tree)
name.check(tree,data) # this function checks the names of tips of class phylo with rownames of data
plot(data[,c("log_mass","mean_gc3")])
print(colnames(data))
plot(data[,c("mean_gc3","sd_gc3")])
plot(data[,c("ct_max","mean_gc3")])
#plot(data[,c("ct_min","mean_gc3")])
plot(density(data$mean_gc4, na.rm=TRUE), main="GC4 Spread"); rug(data$gc4)
# convert phylogenetic tree into a special type of R object called a correlation structure
corBM <- corBrownian(phy=tree, form= ~order_)
corP <- corPagel(value = 1, phy = tree, fixed = FALSE, form = ~species_normalized)

# Perform PGLS; the gls model is very similar to lm, but I need to include the correlation argument 
model_gls_default <- gls(mean_gc3 ~ log_mass + genome_size + ct_max, data=data, correlation = corBM, method="ML", na.action=na.omit)
model_gls_complex <- gls(mean_gc3 ~ log_mass * genome_size, data=data, correlation = corBM, method="ML", na.action = na.omit)
model_lm_mass_gc3 <- lm(mean_gc3 ~ log_mass, data=data)
model_lm_mass_gc4 <- lm(mean_gc4 ~ log_mass, data=data)
model_lm_mass_intron <- lm(mean_intron_gc ~ log_mass, data=data)
model_gls_mass_intron <- gls(mean_intron_gc ~ log_mass, data=data, correlation = corBM, method="ML", na.action=na.omit)
model_gls_mass_gc3 <- gls(mean_gc3 ~ log_mass, data=data, correlation = corBM, method="ML")
model_gls_mass_gc4 <- gls(mean_gc4 ~ log_mass, data=data, correlation = corBM, method="ML")
model_gls_ctmax_gc3 <- gls(mean_gc3 ~ ct_max, data = data, correlation = corBM, method="ML", na.action=na.omit)
model_gls_genome_size_gc3 <- gls(mean_gc3 ~ genome_size, data = data, correlation = corBM, method="ML", na.action=na.omit)
model_gls_genome_size_gc4 <- gls(mean_gc4 ~ genome_size, data = data, correlation = corBM, method="ML", na.action=na.omit)
summary(model_gls_mass_gc3)
anova(model_gls_mass_gc3)
summary(model_gls_mass_gc4)
summary(model_gls_mass_intron)
summary(model_lm_mass_gc4)
summary(model_lm_mass_gc3)
summary(model_lm_mass_intron)
summary(model_gls_ctmax_gc3)
summary(model_gls_default)
summary(model_gls_genome_size_gc3)
summary(model_gls_genome_size_gc4)

# model_gls_ctmin <- gls(mean_gc3 ~ ctmin, data = data, correlation = corBM, method="ML", na.action=na.omit)
# summary(model_gls_ctmin)
# anova(model_gls_default, model_gls_mass)
# check residuals
qqnorm(residuals(model_gls_gc4)); qqline(residuals(model_gls_gc4))

#model_gls_p_default <- gls(mean_gc3 ~ log_mass + genome_size + ct_max, data=data, correlation = corPagel(1, phy = tree, fixed = FALSE, form = ~species_normalized), method="ML", na.action=na.omit)
model_gls_p_complex <- gls(mean_gc3 ~ log_mass * genome_size, data=data, correlation = corPagel(1, phy = tree, fixed = FALSE, form = ~species_normalized), method="ML", na.action = na.omit)
model_gls_p_mass_gc4 <- gls(mean_gc4 ~ log_mass, data=data, correlation = corPagel(1, phy = tree, fixed = FALSE, form = ~species_normalized), method="ML")
model_gls_p_mass_gc3 <- gls(mean_gc3 ~ log_mass, data=data, correlation = corPagel(1, phy = tree, fixed = FALSE, form = ~species_normalized), method="ML")
model_gls_p_ctmax_gc3 <- gls(mean_gc3 ~ ct_max, data = data, correlation = corPagel(1, phy = tree, fixed = FALSE, form = ~species_normalized), method="ML", na.action=na.omit)
model_gls_p_genome_size_gc3 <- gls(mean_gc3 ~ genome_size, data = data, correlation = corPagel(1, phy = tree, fixed = FALSE, form = ~species_normalized), method="ML", na.action=na.omit)
model_gls_p_genome_size_gc4 <- gls(mean_gc4 ~ genome_size, data = data, correlation = corPagel(1, phy = tree, fixed = FALSE, form = ~species_normalized), method="ML", na.action=na.omit)
summary(model_gls_p_mass_gc4)
summary(model_gls_p_mass_gc3)
summary(model_gls_p_ctmax_gc3)
summary(model_gls_p_genome_size_gc3)
summary(model_gls_p_genome_size_gc4)

# test for different types of evolution
m_mass_bm <- fitContinuous(tree, log_mass_vector, model="BM")
m_mass_ou <- fitContinuous(tree, log_mass_vector, model="OU")
m_mass_eb <- fitContinuous(tree, log_mass_vector, model="EB")
m_mass_wn <- fitContinuous(tree, log_mass_vector, model="white")
m_mass_lam <- fitContinuous(tree, log_mass_vector, model="lambda")
m_mass_eblb <- fitContinuous(tree, log_mass_vector, model="EB", bounds=list(a=c(-2,100)))

m_mass_lam$opt$lambda

AIC(m_mass_bm,m_mass_ou,m_mass_wn,m_mass_lam,m_mass_eblb)

summary(model_gls_mass)
nagelkerke(model_gls_mass)


model_gls_stdev_gc3 <- gls(sd_gc3 ~ mean_gc3, data = data,
                   correlation = corBrownian(1, phy = tree, form = ~species_normalized),
                   method = "ML")
model_gls_stdev_gc4 <- gls(sd_gc4 ~ mean_gc4, data = data,
                           correlation = corBrownian(1, phy = tree, form = ~species_normalized),
                           method = "ML")

model_gls_p_stdev_gc3 <- gls(sd_gc3 ~ mean_gc3, data = data,
                           correlation = corPagel(1, phy = tree, fixed = FALSE, form = ~species_normalized),
                           method = "ML")
model_gls_p_stdev_gc4 <- gls(sd_gc4 ~ mean_gc4, data = data,
                             correlation = corPagel(1, phy = tree, fixed = FALSE, form = ~species_normalized),
                             method = "ML")

# plotting

plot_gls_model_base <- function(model, data, x, y, model_label = "model",
                            xlab = x, ylab = y, main = NULL) {
  d <- data[, c(x, y)]
  d <- d[complete.cases(d), ]
  
  slope <- coef(model)[x]
  pval <- summary(model)$tTable[x, "p-value"]  # works for nlme::gls
  
  if (is.na(pval)) {
    pval <- summary(model)$coefficients[x, "Pr(>|t|)"]  # fallback for lm
  }
  
  plot(d[[x]], d[[y]],
       pch = 21, bg = "gray80", col = "black",
       xlab = xlab, ylab = ylab,
       main = ifelse(is.null(main), paste(y, "~", x), main))
  
  newdat <- data.frame(xseq = seq(min(d[[x]]), max(d[[x]]), length.out = 200))
  names(newdat) <- x
  
  pred <- predict(model, newdata = newdat)
  
  lines(newdat[[x]], pred, lwd = 2)
  
  legend("topright",
         legend = c(
           paste0("Model: ", model_label),
           paste0("Slope = ", round(slope, 4)),
           paste0("p = ", signif(pval, 3))
         ),
         bty = "n")
}

plot_model_base <- function(model,
                            data,
                            x,
                            y,
                            model_label = "model",
                            outdir = ".",
                            filename = NULL,
                            device = "pdf",
                            width = 7,
                            height = 6,
                            res = 300,
                            xlab = x,
                            ylab = y,
                            main = NULL) {
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  if (is.null(filename)) {
    filename <- paste0(model_label, "_", y, "_vs_", x)
    filename <- gsub(" ", "_", filename)
  }
  
  outfile <- file.path(outdir,
                       paste0(filename, ".", device))
  
  if (device == "pdf") {
    pdf(outfile, width = width, height = height)
  }
  
  if (device == "png") {
    png(outfile,
        width = width,
        height = height,
        units = "in",
        res = res)
  }
  
  d <- data[, c(x, y)]
  d <- d[complete.cases(d), ]
  
  coefs <- coef(summary(model))
  
  slope <- coef(model)[x]
  
  if ("p-value" %in% colnames(coefs)) {
    pval <- coefs[x, "p-value"]
  } else if ("Pr(>|t|)" %in% colnames(coefs)) {
    pval <- coefs[x, "Pr(>|t|)"]
  } else {
    pval <- NA
  }
  
  plot(d[[x]], d[[y]],
       pch = 21,
       bg = "gray80",
       col = "black",
       xlab = xlab,
       ylab = ylab,
       main = ifelse(is.null(main),
                     paste(y, "~", x),
                     main))
  
  newdat <- data.frame(seq(min(d[[x]]),
                           max(d[[x]]),
                           length.out = 200))
  
  names(newdat) <- x
  
  pred <- predict(model, newdata = newdat)
  
  lines(newdat[[x]], pred, lwd = 2)
  
  legend("topright",
         legend = c(
           paste0("Model: ", model_label),
           paste0("Slope = ", round(slope, 4)),
           paste0("p = ", signif(pval, 3))
         ),
         bty = "n")
  
  dev.off()
  
  cat("Saved plot to:\n", outfile, "\n")
}

plot_model_base(
  model = model_lm_mass_gc3,
  data = data,
  x = "log_mass",
  y = "mean_gc3",
  model_label = "Linear",
  outdir = "plots",
  device = "pdf"
)

plot_model_base(
  model = model_gls_mass_gc3,
  data = data,
  x = "log_mass",
  y = "mean_gc3",
  model_label = "GLS_Brownian",
  outdir = "plots",
  device = "pdf"
)


plot_model_base(
  model = model_gls_mass_gc4,
  data = data,
  x = "log_mass",
  y = "mean_gc4",
  model_label = "GLS_Brownian",
  outdir = "plots",
  device = "pdf"
)

plot_model_base(
  model = model_gls_ctmax_gc3,
  data = data,
  x = "ct_max",
  y = "mean_gc3",
  model_label = "GLS Brownian",
  outdir = "plots",
  device = "pdf"
)

plot_model_base(
  model = model_gls_stdev_gc3,
  data = data,
  x = "mean_gc3",
  y = "sd_gc3",
  model_label = "GLS Brownian",
  outdir = "plots",
  device = "pdf"
)

plot_model_base(
  model = model_gls_stdev_gc4,
  data = data,
  x = "mean_gc4",
  y = "sd_gc4",
  model_label = "GLS Brownian",
  outdir = "plots",
  device = "pdf"
)

plot_model_base(
  model = model_gls_genome_size_gc3,
  data = data,
  x = "genome_size",
  y = "mean_gc3",
  model_label = "GLS Brownian",
  outdir = "plots",
  device = "pdf"
)

plot_model_base(
  model = model_gls_mass_intron,
  data = data,
  x = "log_mass",
  y = "mean_intron_gc",
  model_label = "GLS Brownian",
  outdir = "plots",
  device = "pdf"
)

plot_model_base(
  model = model_lm_mass_intron,
  data = data,
  x = "log_mass",
  y = "mean_intron_gc",
  model_label = "Linear",
  outdir = "plots",
  device = "pdf"
)

