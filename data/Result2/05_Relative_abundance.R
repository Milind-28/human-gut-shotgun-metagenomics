# ------------------------------------------------------------------------------
# Purpose: Generate a clustered bar plot with species on the X-axis and 
#          independent cohort bars adorned with Kruskal-Wallis p-values.
# ------------------------------------------------------------------------------

library(tidyverse)

# 1. Load Data
message("Reading cleaned abundance data...")
df <- read_tsv("Species_abundance.tsv", show_col_types = FALSE) %>% filter(!is.na(Cohort))

# 2. Reshape table from Wide to Long format for statistical modeling
long_df <- df %>%
  pivot_longer(
    cols = -c(Sample, Cohort),
    names_to = "Species",
    values_to = "Abundance"
  ) %>%
  mutate(Cohort = factor(Cohort, levels = c("NonIBD", "CD", "UC")))

# 3. Mass Screening: Run Non-Parametric Kruskal-Wallis test per individual species
message("Calculating Kruskal-Wallis significance profiles...")
stat_results <- long_df %>%
  group_by(Species) %>%
  summarise(
    p_value = kruskal.test(Abundance ~ Cohort)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    # Map p-values to publication-grade significance notation
    Significance = case_when(
      p_value <= 0.0001 ~ "****",
      p_value <= 0.001  ~ "***",
      p_value <= 0.01   ~ "**",
      p_value <= 0.05   ~ "*",
      TRUE              ~ "ns"
    )
  )

# Save statistical summary table to .csv file
write_csv(stat_results, "Species_Abundance_Kruskal_Stats.csv")
message("Saved statistical matrix logs to 'Species_Abundance_Kruskal_Stats.csv'")

# 4. Isolate the top 15 most abundant species overall (Matches your reference)
top_species <- long_df %>%
  group_by(Species) %>%
  summarise(Global_Mean = mean(Abundance), .groups = "drop") %>%
  slice_max(order_by = Global_Mean, n = 15) %>%
  pull(Species)

# Filter long dataset and collapse biological replicates to mean abundance
plot_data <- long_df %>%
  filter(Species %in% top_species) %>%
  group_by(Cohort, Species) %>%
  summarise(Mean_Abundance = mean(Abundance), .groups = "drop")

# Identify the peak abundance height per species to float labels cleanly
max_heights <- plot_data %>%
  group_by(Species) %>%
  summarise(Max_Height = max(Mean_Abundance), .groups = "drop")

# Merge mean values, maximum heights, and significance labels together
final_plot_df <- plot_data %>%
  left_join(stat_results, by = "Species") %>%
  left_join(max_heights, by = "Species")

# 5. Generate Publication-Grade Clustered Bar Plot with Statistics
message("Rendering statistical abundance visualization...")
significance_legend <- "Significance (Global Kruskal-Wallis): ns: p > 0.05, *: p <= 0.05, **: p <= 0.01, ***: p <= 0.001, ****: p <= 0.0001"

# Order the species on the X-axis by descending baseline abundance values
final_plot_df <- final_plot_df %>%
  mutate(Species = reorder(Species, -Mean_Abundance))

stat_bar_plot <- ggplot(final_plot_df, aes(x = Species, y = Mean_Abundance, fill = Cohort)) +
  # Create side-by-side clustered bars
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.25) +
  # Add the significance labels directly centered above each triplet cluster
  geom_text(
    # THE FIX: Explicitly mapped 'x = Species' inside the local aesthetic block
    aes(x = Species, y = Max_Height + (max(final_plot_df$Mean_Abundance) * 0.02), label = Significance),
    check_overlap = TRUE,
    size = 4,
    fontface = "bold",
    color = "black",
    inherit.aes = FALSE,
    data = . %>% distinct(Species, .keep_all = TRUE) # Ensures single text label per cluster
  ) +
  # Match target reference PCoA/Bar color hex codes precisely
  scale_fill_manual(values = c("NonIBD" = "#66c2a5", "CD" = "#fc8d62", "UC" = "#8da0cb")) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, size = 9.5, face = "italic", color = "grey30", margin = margin(t = 12))
  ) +
  labs(
    x = "Microbial Taxa (Species Level)",
    y = "Mean Abundance",
    fill = "Cohort / Condition",
    title = "Differential Taxa Abundance Profiles Across Cohorts",
    caption = significance_legend
  )

# 6. Save High-Resolution Figures
ggsave("Taxonomic_Abundance_Stats_Plot.png", plot = stat_bar_plot, width = 11, height = 6.5, dpi = 300)
message("Success! Image generated and saved as 'Taxonomic_Abundance_Stats_Plot.png'.")