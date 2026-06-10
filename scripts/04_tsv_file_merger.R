# ------------------------------------------------------------------------------
# Script: 04_tsv_file_merger.R
# Purpose: Merge individual MetaPhlAn TSV profiles into a unified abundance matrix
#          optimized for downstream comparative diversity analysis.
# ------------------------------------------------------------------------------

# Load required libraries safely
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse", repos = "http://cran.us.r-project.org")
}
library(tidyverse)

# Define relative path to the taxonomic profiles directory
# (Assuming script is run from the 'scripts' directory)
base_dir <- "../data/results/taxonomic_profiles"
categories <- c("CD", "UC", "NonIBD")

message("--- Starting MetaPhlAn Profile Merger ---")

# 1. Loop through directories and find all target files
file_queue <- list.files(
  path = file.path(base_dir, categories),
  pattern = ".*_taxonomic_profile_3\\.tsv$",
  full.names = TRUE
)

if (length(file_queue) == 0) {
  stop("Error: No taxonomic profile files found. Check your relative directory paths.")
}

message(paste("Found", length(file_queue), "profile files across categories."))

# 2. Function to parse an individual MetaPhlAn file cleanly
parse_metaphlan <- function(file_path) {
  
  # Extract Sample ID from filename (e.g., "SRR5936131")
  sample_id <- str_extract(basename(file_path), "^SRR[0-9]+")
  
  # Identify the condition based on the parent folder name
  condition <- basename(dirname(file_path))
  
  # Read table, skipping comment lines (MetaPhlAn headers are prefixed with #)
  # Keeping only columns 1 and 3 (clade name and relative abundance)
  data <- read_tsv(
    file_path, 
    comment = "#", 
    col_names = FALSE, 
    show_col_types = FALSE
  ) %>%
    select(clade_name = X1, relative_abundance = X3)
  
  # Filter for Species-level annotations only (avoids double-counting abundances)
  # Ensures we don't accidentally pull strain-level ('t__') if present
  species_data <- data %>%
    filter(str_detect(clade_name, "s__") & !str_detect(clade_name, "t__")) %>%
    mutate(
      SampleID = sample_id,
      Condition = condition
    )
  
  return(species_data)
}

# 3. Read and combine all files into a long-format data frame
message("Processing and filtering files...")
merged_long <- map_dfr(file_queue, parse_metaphlan)

# 4. Pivot into a wide matrix (Samples as rows, Species as columns)
# Fills structural zeros with 0.0 for species not detected in a sample
message("Pivoting data into abundance matrix...")
merged_wide <- merged_long %>%
  pivot_wider(
    names_from = clade_name, 
    values_from = relative_abundance,
    values_fill = 0.0
  )

# 5. Export the merged datasets
output_dir <- "../data/results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Export the clean matrix
write_tsv(merged_wide, file.path(output_dir, "merged_species_abundance_matrix.tsv"))

# Optional: Export metadata tracking mapping file
metadata_map <- merged_wide %>% 
  select(SampleID, Condition)
write_tsv(metadata_map, file.path(output_dir, "sample_metadata_map.tsv"))

message("--- Merger Complete! ---")
message(paste("Abundance matrix saved to:", file.path(output_dir, "merged_species_abundance_matrix.tsv")))