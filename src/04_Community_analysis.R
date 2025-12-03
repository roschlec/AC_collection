# Bacterial community structure of natural Arabidopsis thaliana from New Zealand
# Statistical analysis
# Author: Rudolf Schlechter

# Libraries ---------------------------------------------------------------

library(here)
library(tidyverse)
library(vegan)
library(rstatix)
library(phyloseq)
set.seed(20250725)

# Data Input --------------------------------------------------------------

pseq <- readRDS(here('data_input', 'phyloseqData.rds'))

# Community Matrix --------------------------------------------------------

taxa_abundance <-
  pseq |> 
  psmelt() |> 
  select(OTU, Sample, Abundance) |> 
  pivot_wider(id_cols = OTU, names_from = Sample, values_from = Abundance, values_fill = 0) |> 
  column_to_rownames(var = "OTU") |> 
  mutate(across(all_of(everything()), ~ replace_na(.x, 0))) |> 
  as.matrix()

# Distance matrix ---------------------------------------------------------

taxa_rel <- decostand(taxa_abundance |> t(), method = "total")
dist <- vegdist(taxa_rel, method = "bray")

# Alpha diversity ---------------------------------------------------------

alpha_diversity <- 
  taxa_rel |> 
  as_tibble() |> 
  reframe(
    Shannon = diversity(across(), index = "shannon"),
    InvSimpson = 1/(diversity(across(), index = "simpson")),
    sample = rownames(taxa_rel)) |> 
  mutate(Compartment = ifelse(grepl("epi", sample), "Epi", "Endo"),
         Repl = case_when(grepl("1", sample) ~ "rep1",
                          grepl("2", sample) ~ "rep2",
                          grepl("3", sample) ~ "rep3",
                          grepl("4", sample) ~ "rep4",
                          grepl("5", sample) ~ "rep5")) |> 
  arrange(Compartment, Repl)

# Stats
# Mean Inverse Simpson Index
alpha_diversity |> 
  group_by(Compartment) |> 
  summarise(meanInvSimpson = mean(InvSimpson))

# Shapiro-Wilk Normality Test
alpha_diversity |> 
  group_by(Compartment) |> 
  shapiro_test(InvSimpson)

# Levene's Test for homogeneity of variances
alpha_diversity |> 
  levene_test(InvSimpson ~ as.factor(Compartment))

# Student's t test
alpha_diversity_stat <-
  alpha_diversity |> 
  t_test(InvSimpson ~ Compartment, paired = TRUE, alternative = "two.sided") |> 
  add_xy_position(x = "Compartment")

# Beta diversity ----------------------------------------------------------

# NMDS
nmds <- metaMDS(dist)
stressplot(nmds)

nmds.df <- 
  scores(nmds) |> 
  as_tibble(rownames = "sample") |> 
  mutate(Compartment = factor(ifelse(grepl("epi", sample), "Epi", "Endo"))) |> 
  column_to_rownames(var = "sample")

nmds.df <- nmds.df[match(attr(dist, "Labels"), rownames(nmds.df)), ]

#  Beta dispersion
mod <- betadisper(dist, group = nmds.df$Compartment)
permutest(mod, permutations = 9999)

#  PERMANOVA
perm <- adonis2(dist ~ Compartment, data = nmds.df, permutations = 999) 
(perm.df <- as.data.frame(perm))
