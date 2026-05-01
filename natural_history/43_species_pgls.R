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

combine_dir = "~/projects/mass_predicts_dna_dynamics/"
setwd(combine_dir)
list.files(combine_dir)

# Load my phylogenetic tree
tree <- read.tree("~/projects/mass_predicts_dna_dynamics/tree_rooted_with_geckos_gc_final_oct_29_treeonly.tre")
class(tree)
print(tree)

# Load my data
data_with_anc <- read.csv("combined_data.csv")
# omit the ancestor row
data <- data_with_anc[c(1:43),]
# Prepare the data for PGLS, by setting the rownames to genus_species
rownames(data) <- data$genus_species
# this will specify the order of taxa in my data, this important for matching tree tips to data rows
order_ <- rownames(data)
mass_vector <- setNames(data$log.mass.,order_)
class(mass_vector)
print(mass_vector)
# Ensure the species names in the data match the tip labels in the tree
plot(tree)
name.check(tree,data) # this function checks the names of tips of class phylo with rownames of data
plot(data[,c("log.mass.","gc3_nhphyml")])
print(colnames(data))
plot(data[,c("gc3_nhphyml","stdev")])

# convert phylogenetic tree into a special type of R object called a correlation structure
corBM <- corBrownian(phy=tree, form= ~order_)

# Perform PGLS; the gls model is very similar to lm, but I need to include the correlation argument 
model_default <- gls(gc3_nhphyml~log.mass.*Genome.Size, data=data, correlation = corBM, method="ML")
model_mass_genome <- gls(gc3_nhphyml~log.mass.+Genome.Size, data=data, correlation = corBM, method="ML")
anova(model_default, model_mass_genome)
model_mass <- gls(gc3_nhphyml~log.mass., data=data, correlation = corBM, method="ML")
model <- lm(gc3_nhphyml~log.mass., data=data)
summary(model)
summary(model_mass)
anova(model_mass_genome, model_mass)

# this is my best model
model_mass <- gls(gc3_nhphyml~log.mass., data=data, correlation = corBM)

# test for different types of evolution
m_mass_bm <- fitContinuous(tree, mass_vector, model="BM")
m_mass_ou <- fitContinuous(tree, mass_vector, model="OU")
m_mass_eb <- fitContinuous(tree, mass_vector, model="EB")
m_mass_wn <- fitContinuous(tree, mass_vector, model="white")
m_mass_lam <- fitContinuous(tree, mass_vector, model="lambda")
m_mass_eblb <- fitContinuous(tree, mass_vector, model="EB", bounds=list(a=c(-2,100)))

AIC(m_mass_bm,m_mass_ou,m_mass_wn,m_mass_lam,m_mass_eblb)



summary(model_mass)
nagelkerke(model_mass)


model_stdev <- gls(stdev ~ gc3_nhphyml, data = data,
                   correlation = corBrownian(1, phy = tree, form = ~genus_species),
                   method = "ML")

# Print the model summary
summary(model_mass)
anova(model_mass)
summary(model_stdev)
anova(model_stdev)


