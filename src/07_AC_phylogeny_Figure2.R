# Orthogroup_tree ---------------------------------------------------------

# Libraries ---------------------------------------------------------------

library(here)
library(tidyverse)
library(ggtree)
library(ggnewscale)
library(ggtext)
library(ape)
library(picante)

# Function ----------------------------------------------------------------

safe_getMRCA <- function(tr, tips) {
  tips <- intersect(unique(tips), tr$tip.label)
  if (length(tips) == 0) return(NA_integer_)
  if (length(tips) == 1) return(which(tr$tip.label == tips))
  getMRCA(tr, tips)
}

# Data Input --------------------------------------------------------------

# Global options for plotting
# Colour palette
palette <- readRDS(here('data_input', 'colour_palette.rds'))

# Compartment labels
compartment <- c("Endo", "Epi")
compartment_labels <- setNames(c("Endophytic", "Epiphytic"), compartment)

# ANI clades
ani_clades <- read_rds(here("data_input", "ac_aniclades.rds"))

# List of representative AC members
ac_list <- read_rds(here("data_input", "ac_filter.rds"))

# Full AC taxonomy data set
ac_taxonomy <- 
  read_rds(here("data_input", "taxonomy_ac.rds")) |> 
  mutate(genus = ifelse(genus == "", "n.d.", genus)) |> 
  mutate(mutate(across(everything(), ~ sub("_[A-z]", "", .))))

# Leaf compartments
df_compartment <-
  ac_taxonomy |> 
  column_to_rownames(var = "strain") |> 
  dplyr::select(Compartment)


# Phylogenetic tree data --------------------------------------------------

# Phylogenetic tree
tree <- read.tree(here("data_input", "SpeciesTree_rooted_node_labels.txt"))

# Tree metadata
meta_tree <-
  data.frame(strain = tree$tip.label) |> 
  left_join(ac_taxonomy, by = "strain") |> 
  left_join(ani_clades, by = "strain")

# Clade nodes
nodes_clades <- 
  meta_tree |> 
  group_by(clade) |> 
  summarise(
    n = length(clade),
    node = safe_getMRCA(tree, strain), .groups = "drop") |> 
  mutate(ani_clade = factor(clade))

# Family grouping annotation
nodes_taxa <-
  meta_tree |> 
  group_by(family) |> 
  summarise(node = safe_getMRCA(tree, strain), .groups = "drop")

# Genus annotations
nodes_genus <-
  meta_tree |> 
  group_by(genus) |> 
  summarise(node = safe_getMRCA(tree, strain), .groups = "drop")


# Plot --------------------------------------------------------------------

p <- ggtree(tree, layout = "circular") +
  
  geom_cladelab(data = filter(nodes_clades, !is.na(node) & node > 1), 
                mapping = aes(node = node, label = ani_clade), 
                fontsize = 0, 
                barcolour = "darkgray",
                barsize = 0.8,
                align = TRUE, 
                offset = 0.04) +
  
  geom_tiplab(size = 2.5, align = TRUE, linesize = 0.1, 
              linetype = "dotted", offset = 0.08) +
  
  geom_hilight(data = nodes_taxa,
               mapping = aes(node = node, fill = family),
               align="right", alpha = 0.5, colour = "black", linewidth = 0.2) +
  
  geom_cladelab(data = nodes_genus,
                mapping = aes(node = node, label = genus, fontface = "italic"), 
                angle="auto", vjust = 0.5, hjust = 0, 
                offset.text = 0.02, 
                offset = 0.15,
                align = TRUE, size = 2, lineheight = 1) +
  
  geom_tree(linewidth = 0.5) +
  theme_tree() +
  theme(legend.position = "right") +
  scale_fill_manual(name = "Family",
                      values = palette,
                      na.translate = FALSE)  +
  
  guides(colour = "none",
         fill = guide_legend(ncol = 5, order = 2, position = "bottom",
                               override.aes = list(size = 2, shape = 16))) + 
  new_scale_fill()

plt_figure2 <- 
  gheatmap(p, df_comp, offset = -0.02, width = 0.05, hjust = 1,
           colnames_angle = 95, colnames = FALSE) +
  scale_fill_manual(name = "Compartment",
                    values = palette, labels = compartment_labels,
                    na.translate = FALSE) +
  guides(fill = guide_legend(ncol = 1, order = 2, position = "bottom",
                             override.aes = list(size = 2, shape = 16)))


# Output ------------------------------------------------------------------

ggsave(plot = plt_figure2,
       here("output", "Figure2.tiff"), dpi = 600, width = 12, height = 15)
