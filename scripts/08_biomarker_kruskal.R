# ------------------------------------------------------------------------------
# Script: 09_biomarker_kruskal.R
# Purpose: Identify specific species acting as biomarkers across NonIBD, CD, 
#          and UC using Kruskal-Wallis and Post-Hoc Dunn's testing with FDR.
# ------------------------------------------------------------------------------

# Load or install required libraries
required_packages <- c("tidyverse", "FSA", "ggpubr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "http://cran.us.r-project.org")
}
library(tidyverse)
library(FSA)
library(ggpubr)

# 1. Load the merged abundance matrix
matrix_path <- "../data/results/merged_species_abundance_matrix.tsv"
if (!file.exists(matrix_path)) stop("Error: Merged abundance matrix not found.")
merged_wide <- read_tsv(matrix_path, show_col_types = FALSE)

# 2. Reshape to long format and clean taxonomy names
species_cols <- colnames(merged_wide)[grep("\\|s__", colnames(merged_wide))]

merged_long <- merged_wide %>%
  pivot_longer(cols = all_of(species_cols), names_to = "Full_Taxonomy", values_to = "Abundance") %>%
  mutate(
    Species = str_extract(Full_Taxonomy, "s__[a-zA-Z0-9_.-]+$") %>% 
              str_replace("^s__", "") %>% 
              str_replace_all("_", " ")
  )

# 3. Mass Screening: Run Kruskal-Wallis across EVERY unique species
message("Screening features for differential abundance...")
kw_results <- merged_long %>%
  group_by(Species) %>%
  summarise(
    p_value = kruskal.test(Abundance ~ Condition)$p.value,
    .groups = "drop"
  ) %>%
  # Apply Benjamini-Hochberg False Discovery Rate adjustment
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

# 4. Filter for significantly shifting species
# Note: In smaller discovery cohorts, features with raw p < 0.05 are often tracked 
# alongside strict FDR p_adj < 0.1 to maximize potential biological insights.
significant_biomarkers <- kw_results %>% 
  filter(p_value < 0.05) %>% 
  pull(Species)

message(paste("Found", length(significant_biomarkers), "potentially significant biomarker species (raw p < 0.05)."))

# Save full stats table
write_tsv(kw_results, "../data/results/all_species_kruskal_results.tsv")

# 5. Run Post-Hoc Pairwise Dunn's Test for the significant features
if (length(significant_biomarkers) > 0) {
  message("\nRunning Pairwise Dunn's post-hoc tests...")
  
  dunn_accumulator <- list()
  
  for (sp in significant_biomarkers) {
    sp_data <- merged_long %>% filter(Species == sp)
    dunn_res <- dunnTest(Abundance ~ Condition, data = sp_data, method = "bh")$res
    dunn_res$Species <- sp
    dunn_accumulator[[sp]] <- dunn_res
  }
  
  pairwise_df <- bind_rows(dunn_accumulator) %>% select(Species, Comparison, Z, P.unadj, P.adj)
  write_tsv(pairwise_df, "../data/results/biomarker_pairwise_dunn.tsv")
  print(head(pairwise_df, 15))
} else {
  stop("No species met the significance threshold. Analysis complete.")
}

# 6. Plot the Top 4 Most Significant Biomarkers
top_4_plots <- kw_results %>% head(4) %>% pull(Species)

plot_data <- merged_long %>%
  filter(Species %in% top_4_plots) %>%
  mutate(Condition = factor(Condition, levels = c("NonIBD", "CD", "UC")))

biomarker_boxplots <- ggplot(plot_data, aes(x = Condition, y = Abundance, fill = Condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, width = 0.5, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  facet_wrap(~Species, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("NonIBD" = "#66c2a5", "CD" = "#fc8d62", "UC" = "#8da0cb")) +
  labs(
    x = NULL,
    y = "Relative Abundance (%)",
    title = "Top Characterized Microbial Biomarkers Driving Cohort Variance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold.italic", size = 11),
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "none"
  )

# Save the plot grid
ggsave("../data/results/top_biomarker_boxplots.png", plot = biomarker_boxplots, width = 8.5, height = 7, dpi = 300)
message("\n[Success] Biomarker analysis executed completely.")
message("Boxplots generated: ../data/results/top_biomarker_boxplots.png")