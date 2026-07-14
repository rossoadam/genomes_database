# install.packages(c(
#   "tidyverse",
#   # "data.table",
#   "janitor",
#   "GGally",
#   "car",
#   "broom",
#   "patchwork",
#   # "nlme"
# ))

install.packages('R.utils')

library(tidyverse)
library(data.table)
library(janitor)
library(GGally)
library(car)
library(broom)
library(patchwork)
library(nlme)
library(R.utils)

# set paths

sql_tsv_dir <- "/Users/rossoaa/projects/genomes/records/sql_tsvs"

genome_summary_dir <- file.path(sql_tsv_dir, "genome_summary")
sequence_summary_dir <- file.path(sql_tsv_dir, "sequence_summary")

out_dir <- file.path(sql_tsv_dir, "zuur_gc_decay_exploration")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# read tsvs

read_sql_tsv <- function(path) {
  readr::read_tsv(
    path,
    na = c("\\N", "NA", "NaN", "", "NULL", "null"),
    show_col_types = FALSE
  ) %>%
    clean_names()
}

natural_history <- read_sql_tsv(file.path(sql_tsv_dir, "natural_history.tsv"))
genomes <- read_sql_tsv(file.path(sql_tsv_dir, "genomes.tsv"))
window_set <- read_sql_tsv(file.path(sql_tsv_dir, "window_set.tsv"))
analysis_run <- read_sql_tsv(file.path(sql_tsv_dir, "analysis_run.tsv"))

# check that the key fields are present

names(natural_history)
names(genomes)
names(window_set)
names(analysis_run)
names(genome_summary)

# read compressed genome_summary chunks

genome_summary_files <- list.files(
  genome_summary_dir,
  pattern = "\\.tsv\\.gz$",
  full.names = TRUE
)

length(genome_summary_files)
head(genome_summary_files)

# read all chunks

genome_summary <- data.table::rbindlist(
  lapply(genome_summary_files, fread),
  fill = TRUE
) %>%
  as_tibble() %>%
  clean_names() %>%
  mutate(across(where(is.character), ~na_if(.x, "\\N")))

# inspect

glimpse(genome_summary)
summary(genome_summary$sd_gc)
table(genome_summary$standard_window_size_bp, useNA = "ifany")

# join genome_summary to genomes and natural history -
# use the genome_pk and species_pk keys

dat0 <- genome_summary %>%
  left_join(
    genomes %>% select(genome_pk, accession_id, is_current),
    by = "genome_pk"
  ) %>%
  left_join(natural_history, by = "species_pk") %>%
  mutate(
    accession_id = as.factor(accession_id),
    species_normalized = as.factor(species_normalized),
    window_kb = standard_window_size_bp / 1000,
    window_mb = standard_window_size_bp / 1000000,
    
    mass_g = coalesce(mass_meiri, mass_title, mass_ji),
    log10_mass_g = log10(mass_g),
    
    log_var_gc = log(var_gc),
    log10_var_gc = log10(var_gc),
    log_sd_gc = log(sd_gc),
    log10_sd_gc = log10(sd_gc),
    log10_window_bp = log10(standard_window_size_bp),
    log10_window_kb = log10(window_kb)
  )

# sanity checks

dat0 %>%
  count(accession_id, species_normalized)

# added sd_gc from walk through
dat0 %>%
  select(accession_id, species_normalized, standard_window_size_bp, sd_gc, var_gc, mass_g, genome_size, ct_min, ct_max) %>%
  arrange(accession_id, standard_window_size_bp) %>%
  print(n = 50)

# remove rows that cannot be used to in log-scale models
# if the sd_gc is NA or less than 0 remove the row
# the standard_window_size must be larger than 0
# the standard_window_size must not be NA
# the mass cannot be NA

dat <- dat0 %>%
  filter(
    !is.na(sd_gc),
    sd_gc > 0,
    !is.na(standard_window_size_bp),
    standard_window_size_bp > 0,
    !is.na(mass_g),
    mass_g > 0
  )

# look at dropped rows

anti_join(
  dat0 %>% mutate(row_id = row_number()),
  dat %>% mutate(row_id = row_number()),
  by = "row_id"
) %>%
  select(accession_id.x, species_normalized.x, standard_window_size_bp.x, sd_gc.x, mass_g.x)

anti_join(
  dat0 %>% mutate(row_id = row_number()),
  dat %>% mutate(row_id = row_number()),
  by = "row_id"
)

# looks like no rows were dropped

nrow(dat0)
nrow(dat)

# Zuur style first exploration
# are there extreme values in the response or explanatory variables?

# possible response-like variables

# var_gc
# sd_gc
# iqr_gc
# q05_gc
# q95_gc

# possible explanatory variables

# window_size
# mass_g
# genome_size
# ct_min
# ct_max
# calllable fraction
# largest_seq_fraction
# seq_count_used


pdf(file.path(out_dir, "01_cleveland_dotplots_genome_summary.pdf"), width = 10, height = 8)

op <- par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

dotchart(dat$sd_gc, main = "sd_gc", xlab = "GC standard deviation")
dotchart(dat$log10_sd_gc, main = "log10(sd_gc)", xlab = "log10 GC standard deviation")
dotchart(dat$standard_window_size_bp, main = "window size", xlab = "bp")
dotchart(dat$log10_window_bp, main = "log10 window size", xlab = "log10 bp")
dotchart(dat$mass_g, main = "mass_g", xlab = "g")
dotchart(dat$log10_mass_g, main = "log10 mass_g", xlab = "log10 g")
dotchart(dat$genome_size, main = "genome_size", xlab = "genome size")
dotchart(dat$mean_callable_fraction, main = "mean callable fraction", xlab = "fraction")
dotchart(dat$largest_seq_fraction, main = "largest seq fraction", xlab = "fraction")

par(op)
dev.off()



