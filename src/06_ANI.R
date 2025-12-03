# FastANI -----------------------------------------------------------------
# Author: Rudolf Schlechter

# Libraries ---------------------------------------------------------------

library(here)
library(tidyverse)
library(ggtree)
library(ape)
library(reshape2)

# Data Input --------------------------------------------------------------

# AC Taxa
ac_taxa <- 
  readRDS(here("data_input", "taxonomy_ac.rds")) |> 
  select(strain, Compartment)

# CheckM for completeness and contamination
checkm <-
  read_csv(here("data_input", "checkM.csv"), show_col_types = FALSE) |> 
  rename(strain = `Bin Name`, lineage = `Marker Lineage`)

# FastANI
# Data Frame
ani.df <- read.csv(here("data_input", "fastani.csv"))

# Matrix
ani.matrix <- 
  ani.df |> 
  pivot_wider(id_cols = query, names_from = reference, values_from = ANI, values_fill = 0) |> 
  column_to_rownames(var = "query") |> 
  as.matrix()
ani.matrix <- ani.matrix[, rownames(ani.matrix)]

# Hierarchical clustering
stopifnot(is.matrix(ani.matrix), all(rownames(ani.matrix) == colnames(ani.matrix)))

ani.dist <- 100 - ani.matrix

ani.clust <- hclust(as.dist(ani.dist), method = "average")

cl <- cutree(ani.clust, h = 5)

# ANI clades --------------------------------------------------------------

# Members of the AC collection
members <- tibble(strain = names(cl), clade = cl)

# Keep only clades with >= 2 members (true “>95% ANI clades”)
valid_clades <- 
  members |> 
  count(clade, name = "size")  |> 
  filter(size >= 2)

members_valid <- 
  members |> 
  semi_join(valid_clades, by = "clade") |> 
  left_join(checkm, by = "strain") |> 
  select(-`# Genomes`:-`5+`)

clade_representative <- 
  members_valid |> 
  left_join(ac_taxa, by = "strain") |> 
  group_by(Compartment, clade)  |> 
  arrange(desc(Completeness), Contamination, strain)  |> 
  slice(1) |> 
  ungroup()

reps <-
  clade_representative |> 
  select(strain, clade) |> 
  mutate(type = "clade_representative")

singletons <- 
  members |> 
  count(clade) |>
  filter(n == 1) |>
  left_join(members, by = "clade") |>
  select(-n) |> 
  mutate(type = "singleton")

strain_list <- 
  rbind(singletons, reps) |> 
  left_join(ac_taxa, by = "strain") |> 
  arrange(clade)

# Save --------------------------------------------------------------------

write_rds(members, here("data_input", "ac_aniclades.rds"))
write_rds(strain_list, here("data_input", "ac_filter.rds"))
