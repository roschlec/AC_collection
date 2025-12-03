# AC collection - Taxonomy
# Author: Rudolf Schlechter

# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)

# Data input --------------------------------------------------------------

# GTDBTK Classification
ac_list <- 
  read_csv(here("data_input", "AC_samples.csv"))

ac_gtdbtk <-
  read_csv(here("data_input", "gtdbtk_classification.csv"))

ac_tax <- 
  ac_gtdbtk |> 
  separate(classification,
           into = c("domain","phylum","class","order","family","genus","species"),
           sep = ";") |> 
  mutate(across(everything(), ~ sub(".*__", "", .))) |> 
  mutate(across(everything(), ~ sub("_[A-z]", "", .))) |> 
  mutate(
    class = case_when(
      order == "Burkholderiales" ~ "Betaproteobacteria",
      order == "Cytophagales" ~ "Cytophagia",
      order == "Flavobacteriales" ~ "Flavobacteriia",
      order == "Sphingobacteriales" ~ "Sphingobacteriia",
      class == "Actinomycetia" ~ "Actinomycetes", 
      TRUE ~ class),
    order = case_when(
      family %in% c("Dermatophilaceae", "Microbacteriaceae", "Micrococcaceae") ~ "Micrococcales",
      family == "Kineococcaceae" ~ "Kineosporiaceae", TRUE ~ order)) |> 
  left_join(ac_list, by = "strain")

# Save --------------------------------------------------------------------

saveRDS(ac_tax, here("data_input", "taxonomy_ac.rds"))
