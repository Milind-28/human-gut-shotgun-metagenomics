#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Script: 07_complete_alpha_diversity.R
# Purpose: Comprehensive Alpha Diversity Testing (Richness, Shannon, Simpson, Evenness)
#          with automated distribution checks, statistics, and multi-panel plotting.
# ------------------------------------------------------------------------------

# 1. Setup and Package Initialization
required_packages <- c("tidyverse", "vegan", "FSA", "ggpubr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "http://cran.us.r-project.org")
}
library(tidyverse)
library(vegan)
library(FSA)
library(ggpubr)

# 2. Load and Isolate Taxonomy Abundance Matrix
matrix_path <- "../data/results/merged_species_abundance_matrix.tsv"
if (!file.exists(matrix_path)) stop("Error: Abundance matrix missing.")

merged_wide <- read_tsv(matrix_path, show_col_types = FALSE)
species_cols <- colnames(merged_wide)[grep("\\|s__", colnames(merged_wide))]

# Separate Data
metadata <- merged_wide %>% select(SampleID, Condition)
abundance <- merged_wide %>% select(all_of(species_cols))

# 3. Calculate Core Metrics Parallelly
message("Computing multiple alpha diversity metrics...")
alpha_df <- metadata %>%
  mutate(
    Observed_Richness = specnumber(abundance),
    Shannon           = diversity(abundance, index = "shannon"),
    Simpson           = diversity(abundance, index = "simpson"),
    Pielou_Evenness   = Shannon / log(Observed_Richness)
  )

# Fix infinite/NaN values if an empty sample exists safely
alpha_df[is.na(alpha_df)] <- 0

# 4. Automated Statistical Significance Profiling Function
run_alpha_stats <- function(metric_name, df) {
  message(paste("\n--- Statistical Profiling for:", metric_name, "---"))
  
  # Test for normality
  shapiro_p <- shapiro.test(df[[metric_name]])$p.value
  
  if (shapiro_p > 0.05) {
    message("Data is normally distributed (Parametric). Running ANOVA...")
    fit <- aov(as.formula(paste(metric_name, "~ Condition")), data = df)
    print(summary(fit))
    anova_p <- summary(fit)[[1]]["Condition", "Pr(>F)"]
    
    if(anova_p < 0.05) {
      message("Significant differences detected. Running Tukey HSD post-hoc:")
      print(TukeyHSD(fit))
    }
  } else {
    message("Data is non-normally distributed (Non-Parametric). Running Kruskal-Wallis...")
    kw_fit <- kruskal.test(as.formula(paste(metric_name, "~ Condition")), data = df)
    print(kw_fit)
    
    if(kw_fit$p.value < 0.05) {
      message("Significant differences detected. Running Dunn's post-hoc (BH Adjusted):")
      print(dunnTest(as.formula(paste(metric_name, "~ Condition")), data = df, method = "bh")$res)
    }
  }
}

# Run the stats for all four pillars
metrics <- c("Observed_Richness", "Shannon", "Simpson", "Pielou_Evenness")
walk(metrics, ~run_alpha_stats(.x, alpha_df))

# 5. Build Comprehensive Multi-Panel Visualization
alpha_df <- alpha_df %>%
  mutate(Condition = factor(Condition, levels = c("NonIBD", "CD", "UC")))

# Helper function to generate standardized box plots
make_box_plot <- function(metric, label, color_palette) {
  ggplot(alpha_df, aes_string(x = "Condition", y = metric, fill = "Condition")) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, width = 0.4, color = "black") +
    geom_jitter(width = 0.12, size = 1.8, alpha = 0.5) +
    scale_fill_manual(values = color_palette) +
    labs(x = NULL, y = label) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(face = "bold")
    )
}

palette <- c("NonIBD" = "#66c2a5", "CD" = "#fc8d62", "UC" = "#8da0cb")

p1 <- make_box_plot("Observed_Richness", "Species Count (S)", palette)
p2 <- make_box_plot("Shannon", "Shannon Index (H')", palette)
p3 <- make_box_plot("Simpson", "Simpson Index (D)", palette)
p4 <- make_box_plot("Pielou_Evenness", "Pielou's Evenness (J')", palette)

# Arrange into a single 4-panel multi-figure grid
final_grid <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))

# 6. Save Artifacts Securely
output_dir <- "../data/results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(output_dir, "complete_alpha_diversity_grid.png"), plot = final_grid, width = 9, height = 8, dpi = 300)
write_tsv(alpha_df, file.path(output_dir, "complete_alpha_diversity_metrics.tsv"))

message("\n[Success] Complete Alpha Diversity pipeline completed. Matrix and figures exported successfully.")