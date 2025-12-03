# Bacterial community structure of natural Arabidopsis thaliana from New Zealand
# Statistical analysis
# Author: Rudolf Schlechter

# Libraries ---------------------------------------------------------------

library(here)
library(tidyverse)
library(vegan)
library(rstatix)
library(ggpubr)
library(phyloseq)
library(patchwork)
library(gridtext)
library(grid)

# Function ----------------------------------------------------------------

extract_ellipse_coords <- 
  function(ordi_ellipses, segments = 60) {
    ellipse_df <- do.call(rbind, lapply(names(ordi_ellipses), function(grp) {
      v <- vegan:::veganCovEllipse(ordi_ellipses[[grp]]$cov,
                                   ordi_ellipses[[grp]]$center,
                                   ordi_ellipses[[grp]]$scale)
      v <- as.data.frame(v)
      v$Compartment <- grp
      v
    }))
    colnames(ellipse_df)[1:2] <- c("NMDS1", "NMDS2")
    ellipse_df
  }

# Data Input --------------------------------------------------------------

# Global options for plotting
# Colour palette
palette <- readRDS(here('data_input', 'colour_palette.rds'))

# Compartment labels
compartment <- c("Endo", "Epi")
compartment_labels <- setNames(c("Endophytic", "Epiphytic"), compartment)

# Data

source(here("src/04_Community_analysis.R"))

# Relative abundance
# Set a threshold to define "rare" taxa
threshold <- 0.01  # 1%

ps_lvl <- 
  pseq |> 
  tax_glom(taxrank = "Class") |> 
  transform_sample_counts(function(x) x / sum(x)) |> 
  psmelt() |> 
  mutate(Class = factor(ifelse(Abundance < threshold, "Other", Class))) |> 
  group_by(Sample, compartment, label, Abundance, Class) |> 
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Beta diversity
#   Ellipse
ordi <- ordiellipse(nmds, nmds.df$Compartment, kind = "sd", conf = 0.95, draw = "polygon")
ellipse_df <- extract_ellipse_coords(ordi, segments = 60)

#   Extract R2 and p-value
R2 <- round(perm.df[1, "R2"], 3)
pval <- signif(perm.df[1, "Pr(>F)"], 3)

#   Annotation text
permanova_label <-
  paste(
    "PERMANOVA<br>",
    paste("*R*<sup>2</sup> = ", R2, "<br>"),
    paste("*p* = ", pval))

info_box <- grobTree(
  textbox_grob(
    permanova_label, x = 0.05, y = 0.95, hjust = 0, vjust = 1,
    gp = gpar(fontsize = 8)))


# Plots -------------------------------------------------------------------

# Relative abundance

plt_fig1a <-
  ps_lvl |> 
  ggplot(aes(x = Sample, y = Abundance, fill = fct_reorder(Class, -Abundance))) +
  geom_bar(stat = "identity") +
  facet_wrap(~compartment, scales = "free_x") +
  theme_bw() +
  ylab("Relative Abundance") +
  xlab(" ") +
  scale_fill_manual(name = "Class", values = palette) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text.x = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank())

# Alpha diversity

plt_fig1b <-
  alpha_diversity |> 
  ggplot(aes(x = Compartment, y = InvSimpson)) +
  stat_summary(aes(colour = Compartment), fun.data = "mean_cl_boot", linewidth = 0.5, size = 0.5) +
  geom_bracket(data = alpha_diversity_stat, 
               label = paste("italic(p) ==", alpha_diversity_stat$p), type = "expression",
               tip.length = c(0.4, 0.02),
               label.size = 3)+
  geom_jitter(colour = "black", width = 0.08, size = 2) +
  theme_bw() +
  scale_y_continuous(name = "Inverse Simpson", limits = c(1.5, 2.9)) +
  scale_x_discrete(labels = compartment_labels) +
  scale_color_manual(values = palette) +   
  guides(color = "none")

# Beta diversity

plt_fig1c <-
  nmds.df |> 
  ggplot(aes(x = NMDS1, y = NMDS2, fill = Compartment)) +
  geom_path(data = ellipse_df, aes(x = NMDS1, y = NMDS2, color = Compartment)) +
  geom_point(pch = 21, size = 3) +
  annotation_custom(info_box) +
  theme_bw() +
  scale_fill_manual(values = palette, labels = compartment_labels) +
  scale_colour_manual(values = palette, labels = compartment_labels) +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  scale_y_continuous(limits = c(-1.5, 2))


# Consolidate plots -------------------------------------------------------

plt_fig1 <-
  plt_fig1a + plt_fig1b + plt_fig1c +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1.5, 1, 1)) &
  theme(text = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, 'cm'),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.2, linetype = "dashed"))

# Output ------------------------------------------------------------------

ggsave(plot = plt_fig1,
  here("output", "Figure1.tiff"), dpi = 600, width = 12, height = 3)
