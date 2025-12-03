# Genome Structure Analysis -----------------------------------------------
# Author: Rudolf Schlechter

# Libraries ---------------------------------------------------------------

library(here)
library(tidyverse)
library(ape)
library(caper)
library(phytools)
library(phylolm)
library(patchwork)

# Functions ---------------------------------------------------------------
source(here("src", "helper_functions.R"))

# Load Data ---------------------------------------------------------------

# Global options for plotting
# Colour palette
palette <- readRDS(here('data_input', 'colour_palette.rds'))

# Compartment labels
compartment <- c("Endo", "Epi")
compartment_labels <- setNames(c("Endophytic", "Epiphytic"), compartment)

# AC LIST & TAXONOMY
ac_list <- read_rds(here("data_input", "ac_filter.rds"))
ac_taxa <- read_rds(here("data_input", "taxonomy_ac.rds"))

# Phylogenetic information
tree <- read.tree(here("data_input", "SpeciesTree_rooted_node_labels.txt"))
keep <- unique(ac_list$strain) # Vector of strains to KEEP
keep_in_tree <- intersect(keep, tree$tip.label) # only keep tips that exist in the tree
tr_trim <- keep.tip(tree, keep_in_tree) # prune
tr_trim <- force.ultrametric(tr_trim, method = "extend")

# Genome Features
genome <- read_rds(here("data_input", "genome.rds"))
genome.category <- data.frame(category = "Genome", feature = names(genome))

# Phage
phage <- read_rds(here("data_input", "phage.rds"))
phage.category <- data.frame(category = "Phage", feature = names(phage))

# Plasmids
plasmid <- readRDS(here("data_input", "plasmids.rds")) |> rename("Plasmid" = n, "strain" = Strain)
plasmid.category <- data.frame(category = "Plasmid", feature = names(plasmid))

# Combine datasets --------------------------------------------------------
categories.features <- 
  list(genome.category, phage.category, plasmid.category) |> 
  bind_rows() |> 
  filter(feature != "strain")
  
list.features <- list(ac_taxa, genome, phage, plasmid)

feature.df <-
  reduce(list.features, left_join, by = "strain") |>
  filter(strain %in% ac_list$strain) |> 
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

# Phylogenetic Linear Model -----------------------------------------------

ac_df <- 
  as.data.frame(feature.df) |> 
  mutate(across(where(is.numeric), ~ scale(.)[,1]),
         sp = strain)

rownames(ac_df) <- ac_df$strain

compdat <- caper::comparative.data(tr_trim, ac_df, names.col = "strain", vcv = TRUE)

num_feats <- 
  ac_df |> 
  dplyr::select(where(is.numeric)) |> 
  dplyr::select(-any_of(c("some_meta_to_exclude_if_numeric"))) |> 
  names()

# Phylogeny LM
phylo_lm <- 
  run_phylolm_features(ac_df, tr_trim, num_feats, response = "Compartment")

# Pagan's Lambda
sig_tab <- 
  map_df(num_feats, ~ sig_one(tr_trim, ac_df, .x)) |> 
  mutate(FDR_lambda = p.adjust(p_value, method = "fdr")) |> 
  rename(raw_lambda = lambda, p_lambda = p_value)

# Combine
trait.df <- 
  list(categories.features, phylo_lm, sig_tab) |> 
  reduce(left_join, "feature") |> 
  mutate(feature = case_when(
    feature %in% c("GC", "tRNAs", "tmRNAs", "rRNAs", "ncRNAs", "CDSs", "sORFs") ~ feature,
    TRUE ~ str_to_title(feature))) |> 
  mutate(feature = case_when(
    feature == "Size_kb" ~ "Genome Size",
    feature == "Density" ~ "Coding Density",
    feature == "Regions" ~ "Cis-regulatory ncRNAs",
    feature == "Arrays" ~ "CRISPR Arrays",
    feature == "N_phages" ~ "Number of Phages",
    feature == "Phage_gc" ~ "Phage %GC",
    feature == "Phage_percentage" ~ "Phage-like Percentage",
    TRUE ~ feature
  ))


# Plot --------------------------------------------------------------------

plt_figure4a <-
  trait.df |> 
  mutate(signif = ifelse(p_value < 0.05, "sig", "ns")) |> 
  ggplot(aes(x = estimate, y = reorder(feature, estimate), colour = signif)) +
  facet_grid(rows = vars(category), scales = "free", space = "free") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(aes(x = conf_low, xend = conf_high, yend = feature)) +
  geom_point() +
  theme_bw() +
  guides(colour = "none") +
  labs(
    y = "Feature",
    x = expression(paste("Effect size (", beta  %+-% CI, ")")))

plt_figure4b <-
  trait.df |> 
  ggplot(aes(y = reorder(feature, estimate))) +
  facet_grid(rows = vars(category), scales = "free", space = "free") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(aes(x = raw_lambda, size = -log10(FDR_lambda))) +
  scale_size_continuous(name = expression(-log[10]~"(p-value)"), range = c(0.5, 3.5)) +
  theme_bw() +
  labs(
    y = "Feature",
    x = expression(paste("Pagel's ", lambda)))

# Phylogenetic PCA --------------------------------------------------------
X <- ac_df |> dplyr::select(-strain:-Compartment, -sp)
ppca <- phytools::phyl.pca(tr_trim, X, method = "lambda")

# Loadings & scores
ppca$L
scores <- ppca$S

scores_df <- as.data.frame(scores)
scores_df$strain <- rownames(scores_df)
scores_df <- left_join(scores_df, ac_taxa, by = "strain")

# Plot

plt_figure4c <-
  ggplot(scores_df, aes(x = PC1, y = PC2, color = Compartment)) +
  geom_point(size = 2) +
  stat_ellipse(level = 0.95) +
  labs(x = paste0("PC1 (", round(ppca$Eval[1,1]/sum(ppca$Eval)*100,1), "%)"),
       y = paste0("PC2 (", round(ppca$Eval[2,2]/sum(ppca$Eval)*100,1), "%)")) +
  scale_colour_manual(name = "Compartment", values = palette, labels = compartment_labels) +
  theme_bw()

# Final Plot --------------------------------------------------------------

plt_figure4 <-
  plt_figure4a + plt_figure4b + plt_figure4c +
  plot_annotation(tag_levels = "A") &
  theme(
    axis.text.y = element_text(size = 8, colour = "black"),
    strip.background = element_blank(),
    strip.text.y.right = element_text(angle = 0),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

ggsave(plot = plt_figure4,
  here("output", "Figure4.tiff"), dpi = 300, width = 12, height = 5)
