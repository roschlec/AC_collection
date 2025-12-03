# Create a barplot representing the taxonomic groups of AC strains
# Author: Rudolf Schlechter

# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(ggh4x)

# Data input --------------------------------------------------------------

# Global options for plotting
# Colour palette
palette <- readRDS(here('data_input', 'colour_palette.rds'))

# Compartment labels
compartment <- c("Endo", "Epi")
compartment_labels <- setNames(c("Endophytic", "Epiphytic"), compartment)

# AC taxonomy data set
ac_taxa <- readRDS(here("data_input", "taxonomy_ac.rds"))

# Reconstruct taxonomy of AC ----------------------------------------------

family_count <-
  ac_taxa |> 
  group_by(Compartment, phylum, class, order, family) |> 
  tally() |> 
  group_by(Compartment) |> 
  mutate(Relative_abundance = n/sum(n),
         signed_count = if_else(Compartment == "Endo", 
                                -Relative_abundance, 
                                Relative_abundance))

# Family Level ------------------------------------------------------------

plt_figure3 <- 
  family_count |> 
  ggplot(aes(x = signed_count * 100, 
             y = reorder(family, abs(signed_count)))) +
  
  facet_nested_wrap(~phylum+class+order, 
                    strip.position = "left", 
                    scales = "free_y", 
                    ncol = 1, 
                    axes = FALSE) +
  
  geom_segment(aes(xend = 0, color = Compartment)) +
  geom_point(aes(color = Compartment), size = 2) +
  geom_vline(xintercept = 0) +
  
  geom_text_repel(
    data =  family_count |> filter(signed_count < 0),
    aes(label = n),
    force_pull   = 0, # do not pull toward data points
    nudge_x = -0.01,
    direction = "x",
    segment.size = 0.1,
    size = 3) +
  
  geom_text_repel(
    data =  family_count |> filter(signed_count > 0),
    aes(label = n),
    force_pull   = 0, # do not pull toward data points
    nudge_x = 0,
    direction = "x",
    segment.size = 0.1,
    size = 3) +

  guides(colour = guide_legend(title = "Compartment")) +
  scale_x_continuous(name = "Percentage of strains per compartment [%]", 
                     expand = c(0.01,0), labels = abs,
                     breaks = seq(-50, 30, 10), limits = c(-50, 30)) +
  scale_colour_manual(values = palette[1:2], labels = compartment_labels) +
  theme_classic() +
  theme(
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),          
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    strip.placement = "outside",
    strip.switch.pad.wrap = unit(1, "mm"),
    strip.background = element_rect(linewidth = 0.2),
    strip.text.y.left = element_text(size = 8, angle = 0, color = 'black'),
    axis.text = element_text(size = 8, color = 'black'),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = "gray"),
    panel.spacing.y = unit(0.5, "mm"))

# Save --------------------------------------------------------------------

pdf(here("output", "Figure3.pdf"), height = 6, width = 12)
plt_figure3
dev.off()
