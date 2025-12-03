# FastANI -----------------------------------------------------------------
# Author: Rudolf Schlechter

# Libraries ---------------------------------------------------------------

library(here)
library(tidyverse)
library(ggtree)
library(ape)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

# Data Input --------------------------------------------------------------
# Colour palette
palette <- readRDS(here('data_input', 'colour_palette.rds'))

# Compartment labels
compartment <- c("Endo", "Epi")
compartment_labels <- setNames(c("Endophytic", "Epiphytic"), compartment)

# Data
source(here("src", "06_ANI.R"))

# Complex Heatmap ---------------------------------------------------------

# Gradient for 80-100
col_fun <- colorRamp2(
  c(0, 79.999, 80, 90, 100),
  c("gray85", "gray85", "#EDDD53", "#57C785", "#08306b")
)

# Legends
legend_breaks <- seq(80, 100, by = 5)

# Heatmap
ani.ht <-
  Heatmap(ani.matrix, 
        name = "ANI",
        col = col_fun, 
        border = TRUE,
        cluster_rows = ani.clust,
        cluster_columns = ani.clust,
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7),
        heatmap_legend_param = list(
          at = seq(80, 100, by = 5),
          title = "ANI (%)")
        )

# Dendrogram --------------------------------------------------------------

ani.tree <- as.phylo(ani.clust)

ac.df <-
  ac_taxa |> 
  left_join(strain_list, by = c("strain", "Compartment")) |> 
  mutate(select = ifelse(is.na(clade), "No", "Yes"))

p <- 
  ggtree(ani.tree, layout = "dendrogram", branch.length = 0.2) +
  theme(plot.margin = margin(2, 0, 3, 0, "cm"))

p2 <- p %<+% ac.df

p3 <-
  p2 + 
  geom_tiplab(size = 3, offset = -0.5, aes(colour = Compartment)) +
  geom_tippoint(aes(shape = select)) +
  scale_colour_manual(values = palette) +
  scale_shape_manual(values = c(NA, 16)) +
  guides(shape = 'none')


# Output ------------------------------------------------------------------

# Figure S3A
tiff(here("output", "FigureS3a.tiff"), res = 300, height = 3200, width = 3800)
draw(ani.ht)
dev.off()

# Figure S3B
ggsave(plot = p3,
       here("output", "FigureS3b.tiff"), dpi = 300, width = 14, height = 5)
