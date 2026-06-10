# ------------------------------------------------------------------------------
# Script: 08_beta_diversity_permanova.R
# Purpose: Perform Beta Diversity assessment using Bray-Curtis distance, 
#          PERMANOVA testing, and generate a PCoA ordination plot.
# ------------------------------------------------------------------------------

# Load or install required libraries
required_packages <- c("tidyverse", "vegan", "pairwiseAdonis")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "pairwiseAdonis") {
      # Install pairwiseAdonis from GitHub if missing
      if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "http://cran.us.r-project.org")
      remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
    } else {
      install.packages(pkg, repos = "http://cran.us.r-project.org")
    }
  }
}
library(tidyverse)
library(vegan)
library(pairwiseAdonis)

# 1. Load the merged abundance matrix
matrix_path <- "../data/results/merged_species_abundance_matrix.tsv"
if (!file.exists(matrix_path)) {
  stop("Error: Merged abundance matrix file not found. Run your upstream merger script first.")
}
merged_wide <- read_tsv(matrix_path, show_col_types = FALSE)

# 2. Separate metadata and taxonomy values
species_cols <- colnames(merged_wide)[grep("\\|s__", colnames(merged_wide))]
metadata <- merged_wide %>% select(SampleID, Condition)
abundance_matrix <- merged_wide %>% select(all_of(species_cols))

# 3. Compute Bray-Curtis Distance Matrix
message("Calculating Bray-Curtis distance matrix...")
bray_dist <- vegdist(abundance_matrix, method = "bray")

# 4. Execute Global PERMANOVA Testing
message("\n=============================================")
message("          GLOBAL PERMANOVA RESULTS           ")
message("=============================================")
set.seed(42) # Lock randomization seed for reproducibility
permanova_global <- adonis2(bray_dist ~ Condition, data = metadata, permutations = 999)
print(permanova_global)

global_p <- permanova_global$`Pr(>F)`[1]
r_squared <- permanova_global$R2[1]

# 5. Check Multivariate Homogeneity of Dispersions (PERMDISP)
# This checks if the variance (spread) inside the groups is significantly different.
message("\n=============================================")
message("          PERMDISP DISPERSION RESULTS        ")
message("=============================================")
dispersion <- betadisper(bray_dist, metadata$Condition)
permdisp_test <- permutest(dispersion, permutations = 999)
print(permdisp_test)

# 6. Pairwise Post-Hoc PERMANOVA
# Pinpoints exactly which pairs of cohorts have distinct community dynamics
message("\n=============================================")
message("       PAIRWISE POST-HOC PERMANOVA           ")
message("=============================================")
pairwise_results <- pairwise.adonis2(bray_dist ~ Condition, data = metadata, nperm = 999)

# Clean and display pairwise data neatly
for (comparison in names(pairwise_results)) {
  if (comparison != "parent_call") {
    cat("\nComparison Group:", comparison, "\n")
    print(pairwise_results[[comparison]])
  }
}

# 7. Generate Coordinates for PCoA Plotting
message("\nGenerating PCoA Ordination Plot data...")
pcoa_coor <- cmdscale(bray_dist, k = 2, eig = TRUE)

# Calculate percentage variance explained by Axis 1 and Axis 2
var_explained <- round(100 * (pcoa_coor$eig / sum(pcoa_coor$eig)), 1)

pcoa_df <- data.frame(
  SampleID = metadata$SampleID,
  Condition = factor(metadata$Condition, levels = c("NonIBD", "CD", "UC")),
  PCoA1 = pcoa_coor$points[,1],
  PCoA2 = pcoa_coor$points[,2]
)

# 8. Plot with ggplot2
beta_pcoa_plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Condition, fill = Condition)) +
  # Draw a 95% confidence ellipse around each group cluster
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(group = Condition), linetype = "dashed", size = 0.4) +
  geom_point(size = 3.5, alpha = 0.8, shape = 21, color = "black") +
  scale_color_manual(values = c("NonIBD" = "#66c2a5", "CD" = "#fc8d62", "UC" = "#8da0cb")) +
  scale_fill_manual(values = c("NonIBD" = "#66c2a5", "CD" = "#fc8d62", "UC" = "#8da0cb")) +
  labs(
    x = paste0("PCoA Coordinate 1 (", var_explained[1], "%)"),
    y = paste0("PCoA Coordinate 2 (", var_explained[2], "%)"),
    title = "Beta Diversity PCoA Ordination Matrix",
    subtitle = paste0("Global PERMANOVA: R² = ", round(r_squared, 3), " (p = ", global_p, ")")
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

# 9. Export Results safely
output_dir <- "../data/results"
ggsave(file.path(output_dir, "beta_diversity_pcoa.png"), plot = beta_pcoa_plot, width = 7.5, height = 6.5, dpi = 300)

message("\n[Success] Beta diversity profiling successfully completed.")
message("PCoA Ordination Plot exported to: ", file.path(output_dir, "beta_diversity_pcoa.png"))