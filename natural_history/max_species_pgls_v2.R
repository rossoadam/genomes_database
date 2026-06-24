# Load necessary libraries
# install.packages('geiger')
# install.packages('nlme')
# install.packages('phytools')
# install.packages('ape')
# install.packages("MuMIn")  # only if not installed
# install.packages("AED")  # only if not installed
library('MuMIn')
library('phytools')
library('ape')
library('geiger')
library('nlme')
library('rcompanion')
library('AED') # couldn't get this working so I am just going to use the base r

tree_path = "~/projects/genomes/records/compleasm/alignments/t1.0_e1_o7_with_s_punctatus/06c_rooted_tree_for_pgls/rev_dna_renamed_keep_lengths_pgls_altered.rooted_renamed.nwk"
master_ouput_path = "/Users/rossoaa/projects/genomes/records/compleasm/alignments/t1.0_e1_o7_with_s_punctatus/11_final_outputs/gc_metrics_outputs/master_output.tsv"
outdir="/Users/rossoaa/projects/genomes/records/compleasm/alignments/t1.0_e1_o7_with_s_punctatus/11_final_outputs/figures"
# most recently gc4 results were included in the master out. old results were moved to an archive file
setwd(outdir)

# Load my phylogenetic tree
tree <- read.tree(tree_path)
class(tree)
print(tree)

# Load my data
data <- read.csv(master_ouput_path, sep="\t")
# Prepare the data for PGLS, by setting the rownames to genus_species
rownames(data) <- data$genus_species
# this will specify the order of taxa in my data, this important for matching tree tips to data rows
order_ <- rownames(data)
# how does the distribution of the data appear?
hist(data[,c("annual_number_of_eggs")])   # log
hist(data[,c("svl_mm_title")])            # log
hist(data[,c("mass_preferred")])          # log
hist(data[,c("range_size")])              # log
hist(data[,c("max_longevity_years")])     # log
hist(data[,c("genome_size")])             # appears normal
hist(data[,c("ctmax")])                   # appears normal

# do the transformations suggested above
mass_vector <- setNames(data$mass_preferred, order_)
log_mass_vector <- log(mass_vector, base = 10)
data$log_mass = log_mass_vector

svl_vector <- setNames(data$max_female_length_svl_mm_title, order_)
log_svl_vector <- log(svl_vector, base = 10)
data$log_svl_title <- log_svl_vector

annual_egg_vector <- setNames(data$annual_number_of_eggs,order_)
log_annual_eggs <- log(annual_egg_vector, base = 10)
data$log_annual_eggs <- log_annual_eggs

longevity_vector <- setNames(data$max_longevity_years, order_)
log_longevity <- log(longevity_vector, base = 10)
data$log_longevity <- log_longevity

# check that these appear as they should
class(log_mass_vector)
print(log_mass_vector)

# check distribution of svl compared to gc4
plot(density(log_svl_vector, na.rm = TRUE), main="SVL Title Spread"); rug(data$svl_mm_title)
plot(density(data$gc4, na.rm=TRUE), main="GC4 Spread"); rug(data$gc4)

# look for collinearity
variables_explanatory <- c("log_mass","log_svl_title","log_longevity")
x_longevity <- data[,variables_explanatory]
x_longevity <- na.omit(x_longevity)

# pairwise correlation
cor(x_longevity, method = "pearson")

# visual check
pairs(x_longevity)

# function for variance inflation factors
vif_base <- function(df) {
  out <- numeric(ncol(df))
  names(out) <- names(df)
  for (v in names(df)) {
    others <- setdiff(names(df), v)
    form <- as.formula(paste(v, "~", paste(others, collapse = "+")))
    r2 <- summary(lm(form, data = df))$r.squared
    out[v] <- 1 / (1 - r2)
  }
  out
}

vif_base(x_longevity)
# following zuur logic

# VIF > 3: worth inspecting
# VIF > 5: concerning
# VIF > 10: serious collinearity

# mass + svl = NO
# svl + longevity = okay

# I am rechecking for collinearity among more variables
# omitting mass and and eggs
# omitting mass allows me to increase sample size
variables_explanatory <- c("log_svl_title","log_longevity","log_annual_eggs")
x_eggs <- data[,variables_explanatory]
x_eggs <- na.omit(x_eggs)

# pairwise correlation
cor(x_eggs, method = "pearson")

# visual check
pairs(x_eggs)

vif_base(x_eggs)

# svl + longevity + eggs = okay?

# Ensure the species names in the data match the tip labels in the tree
plot(tree)
name.check(tree,data) # this function checks the names of tips of class phylo with rownames of data
print(colnames(data))

# plot the relationships
# gc vs svl
plot(data[,c("mean_gc3","log_svl_title")])
plot(data[,c("gc4","log_svl_title")])
# gc vs eggs
plot(data[,c("mean_gc3","log_annual_eggs")])
plot(data[,c("gc4","log_annual_eggs")])
# gc vs longevity
plot(data[,c("mean_gc3","log_longevity")])
plot(data[,c("gc4","log_longevity")])
# gc vs ctmax
plot(data[,c("mean_gc3","ctmax")])
plot(data[,c("gc4","ctmax")])
# gc vs sd
plot(data[,c("mean_gc3","sd_gc3")])

vars_natural_history_gc3 <- c("mean_gc3","log_svl_title","log_longevity","log_annual_eggs")
vars_natural_history_gc4 <- c("gc4","log_svl_title","log_longevity","log_annual_eggs")
data_nh_gc3 <- na.omit(data[,vars_natural_history_gc3])
data_nh_gc4 <- na.omit(data[,vars_natural_history_gc4])
vars2 <- c("mean_gc3","log_svl_title","log_longevity","log_annual_eggs","genome_size","ctmax")

sum(complete.cases(data[, vars_natural_history_gc3]))
sum(complete.cases(data[, vars_natural_history_gc4]))

# new cor variable:
data_nh_gc3 <- na.omit(data[, vars_natural_history_gc3])
rownames(data_nh_gc3) <- data_nh_gc3$species_normalized
tree_nh_gc3 <- drop.tip(
  tree,
  setdiff(tree$tip.label, data_nh_gc3$species_normalized)
)
data_nh_gc3 <- data_nh_gc3[tree_nh_gc3$tip.label, ]


# convert phylogenetic tree into a special type of R object called a correlation structure
corBM <- corBrownian(phy=tree, form= ~order_)

# Perform PGLS; the gls model is very similar to lm, but I need to include the correlation argument 
# gc3 ~ explanatory variables
model_gls_svl_3_gc3 <- gls(mean_gc3 ~ log_svl_vector + genome_size + ctmax, data=data, correlation = corBM, method="ML", na.action=na.omit)
summary(model_gls_svl_3_gc3)
model_gls_svl_2_gc3 <- gls(mean_gc3 ~ log_svl_vector + genome_size, data=data, correlation = corBM, method="ML", na.action=na.omit)
summary(model_gls_svl_2_gc3)
model_gls_svl_1_gc3 <- gls(mean_gc3 ~ log_svl_vector, data=data, correlation = corBM, method="ML", na.action=na.omit)
summary(model_gls_svl_1_gc3)

model_gls_svl_3_gc4 <- gls(gc4 ~ log_svl_vector + genome_size + ctmax, data=data, correlation = corBM, method="ML", na.action=na.omit)
summary(model_gls_svl_3_gc4)
model_gls_svl_2_gc4 <- gls(gc4 ~ log_svl_vector + genome_size, data=data, correlation = corBM, method="ML", na.action=na.omit)
summary(model_gls_svl_2_gc4)
model_gls_svl_1_gc4 <- gls(gc4 ~ log_svl_vector , data=data, correlation = corBM, method="ML", na.action=na.omit)
summary(model_gls_svl_1_gc4)
model_gls_egg_3_gc3 <- gls(mean_gc3 ~ annual_egg_vector + genome_size + ctmax, data=data, correlation = corBM, method="ML", na.action=na.omit)

model_gls_natural_history_3_gc3 <- gls(mean_gc3 ~ log_svl_title + log_annual_eggs + log_longevity, data=data_nh_gc3, correlation = corBM, method="ML", na.action=na.omit)
summary(model_gls_natural_history_3_gc3)

model_gls_natural_history_2_gc3 <- gls(mean_gc3 ~ log_svl_title + log_longevity, data=data_nh_gc3, correlation = corBM, method="ML", na.action=na.omit)
summary(model_gls_natural_history_2_gc3)

model_gls_natural_history_1_gc3 <- gls(mean_gc3 ~ log_svl_vector , data=data_nh_gc3, correlation = corBM, method="ML", na.action=na.omit)
summary(model_gls_natural_history_1_gc3)






model_gls_egg_3_gc4 <- gls(gc4 ~ annual_egg_vector + genome_size + ctmax, data=data, correlation = corBM, method="ML", na.action=na.omit)
summary(model_gls_egg_3_gc4)
model_gls_egg_2_gc4 <- gls(gc4 ~ annual_egg_vector + genome_size , data=data, correlation = corBM, method="ML", na.action=na.omit)
summary(model_gls_egg_2_gc4)
model_gls_egg_1_gc4 <- gls(gc4 ~ genome_size , data=data, correlation = corBM, method="ML", na.action=na.omit)
summary(model_gls_egg_1_gc4)

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

summary(model_gls_egg_gc3)
summary(model_gls_egg_gc4)
summary(model_gls_mass_gc4)
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
