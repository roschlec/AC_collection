# Epiphytic and endophytic bacterial community analysis - ASV Classification
# Author: Rudolf Schlechter

# Load libraries ----------------------------------------------------------

library(here)
library(tidyverse)
library(dada2)
library(DECIPHER)
library(phyloseq)

# Input Data --------------------------------------------------------------

# Sequence table
seqtab <- 
  readRDS(here("data_input", "seqtab.rds"))

seq <-
  read_csv(here("data_input", "ASV_sequences.csv"), show_col_types = FALSE) |> 
  column_to_rownames(var = "ASV")

# Classification ----------------------------------------------------------
taxa <- 
  assignTaxonomy(
    seq$Sequence, 
    # Use here() to add SILVA database trainset
    multithread = TRUE)

rownames(taxa) <- rownames(seq)

# Analysis ----------------------------------------------------------------

samples.out <- rownames(seqtab)

compartment <- str_extract(samples.out, "^[a-zA-Z]+")

samdf <- data.frame(label = compartment) |> 
  mutate(compartment = case_when(
    label == "endo" ~ "Endophytic",
    label == "epi" ~ "Epiphytic"))

rownames(samdf) <- samples.out

# Full data set
ps <- 
  phyloseq(
    otu_table(seqtab, taxa_are_rows=FALSE), 
    sample_data(samdf), 
    tax_table(taxa))

# Filter out Eukaryota, chloroplast and mitochondrial reads
ps_sub <- subset_taxa(ps, 
                      Kingdom != "Eukaryota" & 
                        !Order %in% c("Chloroplast", "Rickettsiales"))

# Save --------------------------------------------------------------------
saveRDS(taxa, here('data_input', 'taxa_silva.rds'))
saveRDS(ps, here('data_input', 'fullphyloseqData.rds'))
saveRDS(ps_sub, here('data_input', 'phyloseqData.rds'))