#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Script: 05_plot_abundance.R
# Purpose: Generate a clustered bar plot with species on the X-axis and 
#          3 independent bars per species representing NonIBD, CD, and UC.
# ------------------------------------------------------------------------------

library(tidyverse)

# 1. Load the merged species abundance matrix
matrix_path <- "../data/results/merged_species_abundance_matrix.tsv"
if (!file.exists(matrix_path)) {
  stop("Error: Merged abundance matrix not found.")
}

merged_wide <- read_tsv(matrix_path, show_col_types = FALSE)

# 2. Match columns containing species taxonomy strings
species_cols <- colnames(merged_wide)[grep("\\|s__", colnames(merged_wide))]

# 3. Reshape table from Wide to Long format
merged_long <- merged_wide %>%
  pivot_longer(
    cols = all_of(species_cols),
    names_to = "Full_Taxonomy",
    values_to = "Abundance"
  ) %>%
  # Clean up long prefixes to isolate just the clear genus + species text
  mutate(
    Species = str_extract(Full_Taxonomy, "s__[a-zA-Z0-9_.-]+$") %>% 
              str_replace("^s__", "") %>% 
              str_replace_all("_", " ")
  )

# 4. Collapse biological replicates by calculating the mean for each condition
cohort_averaged <- merged_long %>%
  group_by(Condition, Species) %>%
  summarise(Mean_Abundance = mean(Abundance), .groups = "drop")

# 5. Filter for the top most abundant species to keep the X-axis clear and clean
# (Adjust n = 8 or 10 based on how wide you want your final graph layout)
top_species <- cohort_averaged %>%
  group_by(Species) %>%
  summarise(total = sum(Mean_Abundance)) %>%
  slice_max(order_by = total, n = 8, with_ties = FALSE) %>%
  pull(Species)

plot_data <- cohort_averaged %>%
  filter(Species %in% top_species)

# 6. Set factor structures to lock down ordering rules
# Groups the clustered side-by-side bars Left-to-Right: NonIBD -> CD -> UC
plot_data <- plot_data %>%
  mutate(Condition = factor(Condition, levels = c("NonIBD", "CD", "UC")))

# 7. Generate the Clustered Bar Plot
clustered_plot <- ggplot(plot_data, aes(x = reorder(Species, -Mean_Abundance), y = Mean_Abundance, fill = Condition)) +
  # position = position_dodge() splits groups into 3 distinct side-by-side bars
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", size = 0.2) +
  # Custom, high-contrast color scheme for your 3 specific cohorts
  scale_fill_manual(values = c("NonIBD" = "#66c2a5", "CD" = "#fc8d62", "UC" = "#8da0cb")) +
  labs(
    x = "Microbial Taxa (Species Level)",
    y = "Mean Relative Abundance (%)",
    fill = "Cohort / Condition"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(), # Cleans background lines between columns
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", color = "black"), # Slants scientific names elegantly
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.background = element_rect(fill = "white", color = NA)
  )

# 8. Save output figure
output_image <- "../data/results/taxonomic_relative_abundance_plot.png"
ggsave(output_image, plot = clustered_plot, width = 11, height = 7, dpi = 300)

message("--- Plotting Success! ---")
message(paste("Saved clustered 3-bar profile plot to:", output_image))