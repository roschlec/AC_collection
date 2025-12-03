# Genome Structure Analysis -----------------------------------------------
# Libraries ---------------------------------------------------------------

library(here)
library(tidyverse)
library(rstatix)
library(vegan)
library(ape)
library(caper)
library(phytools)
library(phylolm)
library(RRPP)
library(patchwork)
library(colorspace)

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

# KEGG
kegg_ref <-
  read_rds(here("data_input", "kegg_annotation.rds")) |> 
  mutate(feature = str_remove(pathway, " \\[PATH:[a-z]+[0-9]+\\]")) |> 
  mutate(feature = str_replace_all(feature, ",|\\s|\\-|/|\\(|\\)", "_")) |>
  distinct(broad, subcategory, feature) |> 
  rename(category = subcategory)

# Biosynthetic Gene Cluster
bgc_ref <- 
  read_csv(here("data_input", "bcg.csv"), show_col_types = FALSE) |>
  mutate(feature = str_replace_all(product, "-", "_")) |> 
  mutate(feature = ifelse(grepl("^[A-Z]", feature), feature, str_to_title(feature))) |> 
  distinct(category, feature) |> 
  mutate(broad = "BGC")

# TSS
tss_ref <- 
  read_rds(here("data_input", "tss.rds")) |> 
  ungroup() |> 
  distinct(tss_category, gene_name, gene)

tss_ref2 <- 
  tss_ref |> 
  dplyr::select(category = tss_category, feature = gene_name) |> 
  mutate(broad = "Secretion System")

# CAZyme
cazyme_ref <-
  read_rds(here("data_input", "cazyme.df")) |> 
  distinct(class, type) |> 
  rename(category = class, feature = type) |> 
  mutate(broad = "CAZyme")

# Functional Datasets -----------------------------------------------------

# KEGG
KEGG <- 
  readRDS(here("data_input", "kegg_pathway.rds")) |> 
  as_tibble(rownames = "label") |> 
  separate(label, into = c("Compartment", "strain")) |> 
  dplyr::select(-Compartment)
colnames(KEGG) <- str_replace_all(colnames(KEGG), ",|\\s|\\-|/|\\(|\\)", "_")

# Biosynthetic Gene Cluster
BGC <- 
  read_csv(here("data_input", "bcg.csv"), show_col_types = FALSE) |>
  mutate(product = str_replace_all(product, "-", "_")) |> 
  rename(strain = Strain) |> 
  group_by(strain, product) |>
  tally() |>
  mutate(product = ifelse(grepl("^[A-Z]", product), product, str_to_title(product))) |> 
  pivot_wider(id_cols = "strain", names_from = product, 
              values_from = n, values_fill = 0)

# TSS
TSS <-
  read_csv(here("data_input", "tss.csv"), show_col_types = FALSE) |> 
  dplyr::select(strain = Strain, gene_name) |>
  left_join(tss_ref, by = "gene_name") |>
  mutate(gene_name = str_replace_all(gene_name, "\\-", "_")) |> 
  group_by(strain, gene_name) |> 
  tally() |> 
  pivot_wider(id_cols = "strain", names_from = gene_name, 
              values_from = n, values_fill = 0)

# CAZyme
CAZyme <-
  read_rds(here("data_input", "cazymes_type.rds"))

# Combine datasets --------------------------------------------------------
# References
categories.features <- 
  list(kegg_ref, bgc_ref, tss_ref2, cazyme_ref) |> 
  bind_rows() |> 
  mutate(
    category = case_when(
      category == "other" ~ "Other",
      category == "terpene" ~ "Terpene",
      TRUE ~ category),
    label = str_replace_all(category, ",|\\s|\\-|/|\\(|\\)", "_"))

# Matrices
list.features <- list(ac_taxa, KEGG, BGC, TSS, CAZyme)

feature.df <-
  reduce(list.features, left_join, by = "strain") |>
  filter(strain %in% ac_list$strain) |> 
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

# Phylogenetic Linear Model -----------------------------------------------

ac_df <- 
  as.data.frame(feature.df) |> 
  mutate(sp = strain)

rownames(ac_df) <- ac_df$strain

num_feats <- 
  ac_df |> 
  dplyr::select(where(is.numeric)) |> 
  dplyr::select(-any_of(c("some_meta_to_exclude_if_numeric"))) |> 
  names()

good_traits <- num_feats[ sapply(ac_df[num_feats], function(x) length(unique(na.omit(x))) > 1) ]

ac_df2 <- ac_df[c("strain", "sp", "Compartment", good_traits)]

compdat <- caper::comparative.data(tr_trim, ac_df2, names.col = "strain", vcv = TRUE)

# Phylogeny LM
phylo_lm <- 
  run_phylolm_features(ac_df2, tr_trim, good_traits, response = "Compartment")

# Pagel's Lambda
sig_tab <- 
  map_df(good_traits, ~ sig_one(tr_trim, ac_df2, .x)) |> 
  mutate(FDR_lambda = p.adjust(p_value, method = "fdr")) |> 
  rename(raw_lambda = lambda, p_lambda = p_value)

# Combine
trait.df <- 
  list(categories.features, phylo_lm, sig_tab) |> 
  reduce(left_join, "feature") |> 
  na.exclude() |> 
  mutate(
    feature = case_when(
      feature == "Biofilm_formation___Vibrio_cholerae" ~ "Biofilm formation",
      feature == "Arginine_biosynthesis" ~ "Arginine biosynthesis",
      feature == "GH50" ~ "Glycoside Hydrolase family 50",
      feature == "PL27" ~ "Glucuronate Lyase",
      feature == "ComM_comGD" ~ "comGD",
      feature == "ComM_comGE" ~ "comGE",
      feature == "MSH_mshP" ~ "mshP",
      TRUE ~ feature
    ))

# Phylogenetic PCA --------------------------------------------------------
X <- ac_df2 |> dplyr::select(-strain:-Compartment)

# Remove highly correlated covariates
cor_mat <- cor(X, use = "pairwise.complete.obs")
to_drop <- caret::findCorrelation(cor_mat, cutoff = 0.90)
X <- X[ , -to_drop]
X <- scale(X)

pc <- prcomp(X, scale. = FALSE)
scores <- 
  as.data.frame(pc$x[,1:2]) |> 
  rownames_to_column(var = "strain") |> 
  left_join(ac.df, by = "strain")

apply(pc$x[,1:2], 2, function(pcaxis)
  phytools::phylosig(tr_trim, pcaxis, method = "lambda", test = TRUE)$lambda)

fit <- lm.rrpp(X ~ ac_df2$Compartment, Cov = vcv.phylo(tr_trim), iter = 999)
anova(fit)

# Plot --------------------------------------------------------------------

figure5a <-
  trait.df |> 
  filter(p_value < 0.05) |> 
  ggplot(aes(x = estimate, y = reorder(feature, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(aes(x = conf_low, xend = conf_high, yend = feature)) +
  geom_point() +
  theme_bw() +
  scale_y_discrete(position = "right") +
  guides(colour = "none") +
  labs(
    y = "",
    x = expression(paste("Effect size (", beta  %+-% CI, ")")))

figure5b <-
  trait.df |> 
  filter(FDR_lambda < 0.05) |>
  ggplot(aes(x = raw_lambda, y = fct_reorder(category, raw_lambda))) +
  facet_grid(rows = vars(broad), scales = "free", space = "free") +
  geom_vline(xintercept = 1.0, linetype = "dashed") +
  geom_boxplot(fill = "gray95", outlier.alpha = 0, width = 0.4, size = 0.2) + 
  geom_point(size = 0.25) +
  scale_x_continuous(limits = c(0, 1.06), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  theme_bw() +
  labs(
    y = "",
    x = expression(paste("Pagel's ", lambda)))

figure5c <-
  ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Compartment), size = 2, pch = 21, alpha = 0.8) +
  theme_bw() +
  labs(x = paste0("PC1 (", round(pc$sdev[1]/sum(pc$sdev)*100,1), "%)"),
       y = paste0("PC2 (", round(pc$sdev[2]/sum(pc$sdev)*100,1), "%)")) +
  scale_fill_manual(values = palette, labels = compartment_labels)


# Compile plots -----------------------------------------------------------

figure5 <-
  free(figure5a) / free(figure5b) / free(figure5c) +
  plot_annotation(tag_levels = "A") &
  theme(
    axis.title.y = element_text(margin = margin(r = 2)),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.ticks.y = element_line(linewidth = 0.2),
    axis.ticks.x = element_line(linewidth = 0.2),
    strip.background = element_blank(),
    strip.text.y.right = element_text(angle = 0, hjust = 0),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# Output ------------------------------------------------------------------

ggsave(plot = figure5, 
       here("output", "Figure5.tiff"), dpi = 300, width = 12, height = 6)
