# Load necessary libraries
# install.packages('geiger')
# install.packages('nlme')
# install.packages('phytools')
# install.packages('ape')
install.packages("MuMIn")  # only if not installed
library('MuMIn')
library('phytools')
library('ape')
library('geiger')
library('nlme')
library('rcompanion')

#X combine_dir = "~/projects/mass_predicts_dna_dynamics/"
tree_path = "~/projects/genomes/records/compleasm/alignments/t1.0_e1_o7_with_s_punctatus/06c_rooted_tree_for_pgls/rev_dna_renamed_keep_lengths_pgls_altered.rooted_renamed.nwk"
master_ouput_path = "~/projects/genomes/records/compleasm/alignments/t1.0_e1_o26_with_s_punctatus/08_final_outputs/gc_metrics_outputs_wtih_thermal/master_output.tsv"
outdir="~/projects/genomes/records/compleasm/alignments/t1.0_e1_o26_with_s_punctatus/08_final_outputs/figures"
gc4_results="~/projects/genomes/records/compleasm/alignments/t1.0_e1_o26_with_s_punctatus/08_final_outputs/gc_metrics_outputs_wtih_thermal/gc4_species_summary.tsv"
setwd(outdir)
list.files(combine_dir)

# Load my phylogenetic tree
tree <- read.tree(tree_path)
class(tree)
print(tree)

# Load my data
#X data_with_anc <- read.csv("combined_data.csv")
data <- read.csv(master_ouput_path, sep="\t")
gc4_results <- read.csv(gc4_results, sep="\t")
# omit the ancestor row - not necessary for the new data
#X data <- data_with_anc[c(1:43),]
# Prepare the data for PGLS, by setting the rownames to genus_species
#X rownames(data) <- data$genus_species
rownames(data) <- data$genus_species
rownames(gc4_results) <- gc4_results$sequence_id
# this will specify the order of taxa in my data, this important for matching tree tips to data rows
order_ <- rownames(data)
order_2 <- rownames(gc4_results)
data$gc4 <- gc4_results$gc4
mass_vector <- setNames(master_output$mass_preferred,order_)
data$log_mass = log(mass_vector, base = 10)
log_mass_vector <- log(mass_vector, base = 10)
class(log_mass_vector)
print(log_mass_vector)
# Ensure the species names in the data match the tip labels in the tree
plot(tree)
name.check(tree,data) # this function checks the names of tips of class phylo with rownames of data
plot(data[,c("log_mass","mean_gc3")])
print(colnames(data))
plot(data[,c("mean_gc3","sd_gc3")])
plot(data[,c("ctmax","mean_gc3")])
plot(data[,c("ctmin","mean_gc3")])
plot(density(data$gc4, na.rm=TRUE), main="GC4 Spread"); rug(data$gc4)
# convert phylogenetic tree into a special type of R object called a correlation structure
corBM <- corBrownian(phy=tree, form= ~order_)

# Perform PGLS; the gls model is very similar to lm, but I need to include the correlation argument 
model_gls_default <- gls(mean_gc3 ~ log_mass + genome_size + ctmax, data=data, correlation = corBM, method="ML", na.action=na.omit)
model_lm_mass <- lm(mean_gc3 ~ log_mass, data=data)
model_gls_complex <- gls(mean_gc3 ~ log_mass * genome_size, data=data, correlation = corBM, method="ML", na.action = na.omit)
model_gls_gc4 <- gls(gc4 ~ log_mass, data=data, correlation = corBM, method="ML")
qqnorm(residuals(model_gls_gc4)); qqline(residuals(model_gls_gc4))
model_lm_gc4 <- lm(gc4 ~ log_mass, data=data)

model_gls_mass <- gls(mean_gc3 ~ log_mass, data=data, correlation = corBM, method="ML")
model_gls_ctmax <- gls(mean_gc3 ~ ctmax, data = data, correlation = corBM, method="ML", na.action=na.omit)
model_gls_ctmin <- gls(mean_gc3 ~ ctmin, data = data, correlation = corBM, method="ML", na.action=na.omit)
anova(model_gls_default, model_gls_mass)

model_gls_genome_size <- gls(mean_gc3 ~ genome_size, data = data, correlation = corBM, method="ML", na.action=na.omit)

summary(model_gls_mass)
anova(model_gls_mass)
summary(model_gls_gc4)
summary(model_lm_gc4)
summary(model_lm_mass)
summary(model_gls_ctmax)
summary(model_gls_ctmin)
summary(model_gls_default)
summary(model_gls_genome_size)

anova(model_mass_genome, model_mass)

# this is my best model
model_mass <- gls(gc3_nhphyml~log.mass., data=data, correlation = corBM)

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


model_gls_stdev <- gls(sd_gc3 ~ mean_gc3, data = data,
                   correlation = corBrownian(1, phy = tree, form = ~genus_species),
                   method = "ML")

# Print the model summary
summary(model_gls_mass)
anova(model_gls_mass)
summary(model_gls_stdev)
anova(model_gls_stdev)

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
  model = model_gls_mass,
  data = data,
  x = "log_mass",
  y = "mean_gc3",
  model_label = "GLS_Brownian",
  outdir = "plots",
  device = "pdf"
)


plot_model_base(
  model = model_lm_mass,
  data = data,
  x = "log_mass",
  y = "mean_gc3",
  model_label = "Linear",
  outdir = "plots",
  device = "pdf"
)

plot_model_base(
  model = model_gls_ctmax,
  data = data,
  x = "ctmax",
  y = "mean_gc3",
  model_label = "GLS Brownian",
  outdir = "plots",
  device = "pdf"
)

plot_model_base(
  model = model_gls_stdev,
  data = data,
  x = "mean_gc3",
  y = "sd_gc3",
  model_label = "GLS Brownian",
  outdir = "plots",
  device = "pdf"
)

plot_model_base(
  model = model_gls_genome_size,
  data = data,
  x = "genome_size",
  y = "mean_gc3",
  model_label = "GLS Brownian",
  outdir = "plots",
  device = "pdf"
)

plot_model_base(
  model = model_gls_gc4,
  data = data,
  x = "log_mass",
  y = "gc4",
  model_label = "GLS Brownian",
  outdir = "plots",
  device = "pdf"
)
