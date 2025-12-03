# Colour Palette - AC collection paper
# Author: Rudolf Schlechter

library(tidyverse)
library(colorspace)
library(phyloseq)
set.seed(12133)

# Data Input --------------------------------------------------------------

# Compartment
compartment <- c("Endo", "Epi")

# AC collection Data
ac_taxa <- 
  readRDS(here("data_input", "taxonomy_ac.rds"))

ac_phylum <- ac_taxa |> pull(phylum) |> unique()
ac_class <- ac_taxa |> pull(class) |> unique()
ac_family <- ac_taxa |> pull(family) |> unique()

# 16S Data (Class level)
# Set a threshold to define "rare" taxa
threshold <- 0.01  # i.e., 1%

taxa_16s <- 
  readRDS(here('data_input', 'phyloseqData.rds')) |> 
  tax_glom(taxrank = "Class") |> 
  transform_sample_counts(function(x) x / sum(x)) |> 
  psmelt() |> 
  mutate(Class = ifelse(Abundance < threshold, "Other", Class)) |> 
  pull(Class) |> 
  unique()

taxa_unique_16s <- setdiff(taxa_16s, ac_class) # 16S taxa without AC


# Colours -----------------------------------------------------------------

colors_compartment <- divergingx_hcl(length(compartment), palette = "ArmyRose")
names_colors_compartment <- setNames(colors_compartment, compartment)

colors_ac_phylum <- qualitative_hcl(n = 10, palette = "Dynamic")[c(4, 7, 1, 9)]
names_colors_ac_phylum <- setNames(colors_ac_phylum, ac_phylum)

colors_ac_family <- qualitative_hcl(length(ac_family), palette = "Dark 3")
names_colors_ac_family <- setNames(colors_ac_family, ac_family)

colors_ac_class <- divergingx_hcl(length(ac_class), palette = "Zissou 1")
names_colors_ac_class <- setNames(colors_ac_class, ac_class)

colors_16s <- qualitative_hcl(length(taxa_unique_16s), palette = "Set 2") 
names_colors_16S <- setNames(colors_16s, taxa_unique_16s)
names_colors_16S["Other"] <- "gray50"

# Combination
name_colors <- c(names_colors_compartment, names_colors_ac_phylum, names_colors_ac_family, names_colors_ac_class, names_colors_16S)


# Save --------------------------------------------------------------------

saveRDS(name_colors, here('data_input', 'colour_palette.rds'))
