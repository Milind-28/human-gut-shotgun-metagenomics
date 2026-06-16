# Ensure required string manipulating library is installed
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
library(stringr)

# 1. Define Paths
input_file  <- "humann3_function_profile.tsv"
output_file <- "humann3_functional_profile_species_wise.csv"

# 2. Fast-read the raw dataset (keeps column names original)
data <- read.delim(input_file, sep = "\t", check.names = FALSE)

all_cols <- colnames(data)
feature_col_name <- all_cols[1]
sample_cols <- all_cols[-1]

# Extract pure sample base IDs (e.g., "CSM67U9H")
sample_bases <- gsub("(_P|_TR)?_pathabundance_cpm$", "", sample_cols)
unique_bases <- unique(sample_bases)

# 3. Apply Industry Prioritization Filter (1 column per unique subject ID)
selected_sample_cols <- c()

for (base in unique_bases) {
  # Get all raw columns associated with this specific sample base
  matches <- sample_cols[sample_bases == base]
  
  if (length(matches) == 1) {
    # No duplicate exists, keep the only column present
    selected_sample_cols <- c(selected_sample_cols, matches)
  } else {
    # Duplicate exists: Pick according to formal sequence priority
    p_run  <- paste0(base, "_P_pathabundance_cpm")
    std_run <- paste0(base, "_pathabundance_cpm")
    tr_run  <- paste0(base, "_TR_pathabundance_cpm")
    
    if (p_run %in% matches) {
      selected_sample_cols <- c(selected_sample_cols, p_run)       # Priority 1
    } else if (std_run %in% matches) {
      selected_sample_cols <- c(selected_sample_cols, std_run)     # Priority 2
    } else if (tr_run %in% matches) {
      selected_sample_cols <- c(selected_sample_cols, tr_run)      # Priority 3
    } else {
      selected_sample_cols <- c(selected_sample_cols, matches[1])  # Absolute fallback
    }
  }
}

# Sub-select only the prioritized columns from your large data frame
data_filtered <- data[, c(feature_col_name, selected_sample_cols)]

# Clean up column headers so they match the clean, un-suffixed Sample IDs
colnames(data_filtered)[-1] <- unique_bases

# 4. Extract Species-stratified rows (Drop total community rows)
stratified_data <- data_filtered[grepl("\\|", data_filtered[[1]]), ]

# 5. Fast Matrix-style Text Parsing (Vectorized C-level parsing)
feature_vectors <- stratified_data[[1]]
parts <- str_split_fixed(feature_vectors, "\\|", 2)
id_desc <- str_split_fixed(parts[, 1], ": ", 2)

# Format species strings cleanly (e.g., "g__Bacteroides.s__Bacteroides_fragilis" -> "Bacteroides fragilis")
species_clean <- gsub(".*s__", "", parts[, 2])
species_clean <- gsub("_", " ", species_clean)

# 6. Generate final database-structured table
final_profile <- data.frame(
  Feature_ID          = id_desc[, 1],
  Feature_Description = id_desc[, 2],
  Species             = species_clean,
  stratified_data[, -1],
  check.names         = FALSE
)

# 7. Write out to structured CSV
write.csv(final_profile, output_file, row.names = FALSE)
print("SUCCESS: Species-wise standard functional profile generated using industry-standard deduplication mapping!")