# Comparative Metagenomic Profiling of the Human Gut Microbiome Across Inflammatory Bowel Disease (IBD) Phenotypes

A comprehensive, end-to-end solo computational biology framework designed to isolate species-level structural signatures of the human gut microbiota across three clinical groups: Non-Inflammatory Bowel Disease controls (**NonIBD**), Crohn's Disease (**CD**), and Ulcerative Colitis (**UC**).

---

## Project Overview & Technical Adaptation

This project maps out a complete, multi-stage metagenomics pipeline demonstrating proficiency in two distinct phases of computational biology engineering:

**Phase 1: Raw Sequence Processing Infrastructure (Bash/Python Engine)**
Developed a robust raw read cleaning, automated adapter trimming, base-quality filtration (`fastp`), and host genome subtraction (`Bowtie2` against the GRCh38 human reference index) workflow to handle raw paired-end sequence archives (N=15) from the iHMP database.

**Phase 2: Downstream Ecological Biostatistics & Biomarker Discovery (R Engine)**
Due to local hardware constraints — specifically the high RAM, storage, and CPU overhead required to locally index massive databases like CHOCOPhlAn — the processing pipeline was pivotally adapted. Instead of generating profile tables locally, pre-computed, peer-reviewed MetaPhlAn v3.0 abundance files were fetched directly from the official **iHMP Inflammatory Bowel Disease Multi-omics Database (IBDMDB)** portal. A full downstream comparative biostatistical suite was engineered from scratch in R to handle abundance matrix curation, alpha/beta diversity modeling, non-parametric multi-variable permutation tests (PERMANOVA/PERMDISP), and mass feature screening for clinical biomarker discovery.

---

## Repository Directory Architecture

```
human-gut-shotgun-metagenomics/
├── README.md
├── scripts/
│   ├── 01_download_data.py            # Automated sequential SRA toolkit raw fetch module
│   ├── 02_qc_and_profile.py           # Core fastp quality trimming & MultiQC dashboard generation
│   ├── 03_host_removal.py             # Bowtie2 human genome (GRCh38) exclusion engine
│   ├── download_result_files.py       # API downloader fetching consortium MetaPhlAn outputs
│   ├── 04_tsv_file_merger.R           # MetaPhlAn profile regex parsing and wide-matrix synthesis
│   ├── 05_plot_abundance.R            # Cluster-averaged relative abundance bar plot generation
│   ├── 06_alpha_diversity.R           # Normality-audited within-sample complexity calculator
│   ├── 07_beta_diversity_permanova.R  # Bray-Curtis distance, PCoA, PERMANOVA, and PERMDISP scripts
│   └── 08_biomarker_kruskal.R         # Mass feature screening via BH-FDR corrected Dunn's tests
└── data/
    ├── metadata/                      # Study flat files (cd_acc.txt, nonibd_acc.txt, uc_acc.txt, SraRunTable.csv)
    ├── raw/                           # [Local Cache] Temporary directory for raw compressed sequence footprints
    ├── fastp_trimmed/                 # [Local Cache] Intermediary quality-trimmed read pairs
    └── results/
        ├── merged_species_abundance_matrix.tsv
        ├── complete_alpha_diversity_metrics.tsv
        ├── all_species_kruskal_results.tsv
        ├── biomarker_pairwise_dunn.tsv
        ├── taxonomic_relative_abundance_plot.png
        ├── complete_alpha_diversity_grid.png
        ├── beta_diversity_pcoa.png
        ├── top_biomarker_boxplots.png
        ├── multiqc_report.html
        ├── bowtie2_host_removed/
        └── qc/
```

---

## Methodological Workflow & Execution Pipeline

### Step 1: Upstream Raw Read Clean-up & Quality Control

**Scripts:** `01_download_data.py`, `02_qc_and_profile.py`

The pipeline uses `prefetch` and `fasterq-dump` (SRA Toolkit) to fetch raw FASTQ data. `fastp` then isolates high-confidence paired reads by removing low-quality bases (Q < 20), correcting base mismatches in overlapping regions, and stripping residual sequencing adapters. MultiQC compiles individual run telemetry into a clean interactive dashboard (`multiqc_report.html`).

### Step 2: Host Contamination Subtraction

**Script:** `03_host_removal.py`

Reads are aligned against the human reference genome index (GRCh38) using `Bowtie2`. By streaming files and writing output summaries to `/dev/null`, the module avoids saving heavy intermediary SAM/BAM alignments — accurately logging sample quality while preserving local storage capacity. Achieves a ~0.00% human sequence validation rate.

### Step 3: Consortium Matrix Retrieval & Curation

**Scripts:** `download_result_files.py`, `04_tsv_file_merger.R`

Pre-calculated MetaPhlAn v3.0 files are retrieved from the IBDMDB repository. `04_tsv_file_merger.R` reads long taxonomic lineage strings and uses regular expressions to selectively isolate species-level assignments (`s__`), while stripping strain-level properties (`t__`). Individual profiles are merged into a unified wide-format data matrix (15 × S) where missing values are assigned an abundance of 0.0%.

### Step 4: Core Taxonomic Distribution Analysis

**Script:** `05_plot_abundance.R`

Computes cohort-averaged abundances to build a clear overview of taxonomic changes across the three phenotypes, visualizing the core differences between healthy baselines and dysbiotic states.

<img width="3300" height="2100" alt="image" src="https://github.com/user-attachments/assets/42b871f5-3163-44c4-b93e-57df291c92cf" />

**Fig:** Comparative Species-Level Abundance Profile Across IBD Phenotypes. Clustered bar chart showing the mean relative abundance (%) of the eight most dominant gut microbial species across study cohorts: Non-Inflammatory Bowel Disease (NonIBD, green), Crohn’s Disease (CD, orange), and Ulcerative Colitis (UC, blue). The horizontal axis grids taxa sequentially by overall dataset abundance. Individual bars within each cluster represent the cohort-wide mean abundance values calculated from underlying biological replicates. Distinct dysbiotic signatures are identifiable, characterized by a sharp expansion of specific Bacteroides members in active Crohn's disease and an inverse dominance shift in Prevotella copri across Ulcerative Colitis states.

### Step 5: Intrasample Micro-Complexity (Alpha Diversity)

**Script:** `06_alpha_diversity.R`

Tracks alpha diversity across three mathematical perspectives:

- **Richness** — Observed Count
- **Evenness** — Pielou's Evenness (J′)
- **Compositional Balance** — Shannon (H′) and Simpson (D)

The script automatically applies a Shapiro-Wilk normality test to select the appropriate statistical method (Parametric ANOVA vs. Non-Parametric Kruskal-Wallis).

### Step 6: Intersample Partitioning & Distance Calculations (Beta Diversity)

**Script:** `07_beta_diversity_permanova.R`

Computes a pairwise **Bray-Curtis dissimilarity matrix** across all samples and projects coordinates into a 2D space using **Principal Coordinate Analysis (PCoA)** with 95% confidence ellipses. Runs a global **PERMANOVA** (`adonis2`, 999 permutations) to assess compositional differences, accompanied by a **PERMDISP** homogeneity check (`betadisper`) to confirm that grouping is driven by true taxonomic shifts rather than unequal group variances.

### Step 7: Clinical Biomarker Discovery via Mass Screening

**Script:** `08_biomarker_kruskal.R`

Loops through every species in the dataset and performs a mass Kruskal-Wallis screening. Applies a **Benjamini-Hochberg (BH) False Discovery Rate (FDR)** adjustment to control for false positives across the large taxonomic profile. Features that pass the significance threshold are evaluated using a post-hoc pairwise **Dunn's test** to discover the key species driving differences between phenotypes.

---

## Key Biostatistical Insights & Results

### 1. Alpha Diversity Exhibits Structural Stability (p > 0.05)

The normality check identified normal distributions for Observed Richness (p = 0.847), Shannon (p = 0.419), and Pielou's Evenness (p = 0.364), analyzed using a global one-way ANOVA. The Simpson index was analyzed via Kruskal-Wallis (p = 0.4025). No metric crossed the significance threshold (α = 0.05).

**Biological Takeaway:** The absolute count and volume distribution of bacterial species remain stable across these samples. Dysbiosis in this cohort is driven by a targeted replacement of specific microbes (taxonomic turnover) rather than a wholesale collapse of the gut's total niche capacity.

<img width="2700" height="2400" alt="image" src="https://github.com/user-attachments/assets/c35aa618-448c-4f67-8975-43310dcc5f05" />

**Fig:** Comprehensive Alpha Diversity Profiling Across IBD Phenotypes. Four-panel multi-metric evaluation of within-sample microbial community complexity across Non-Inflammatory Bowel Disease (NonIBD, green), Crohn’s Disease (CD, orange), and Ulcerative Colitis (UC, blue) cohorts. Individual points represent distinct biological profiling replicates overlaying box-and-whisker distributions. (A) Species Count (S): Measures absolute observed species richness. (B) Shannon Index (H'): Assesses community entropy sensitive to rare taxa dropout. (C) Simpson Index (D): Illustrates community dominance patterns driven by major taxonomic features. (D) Pielou's Evenness (J'): Tracks the structural equitability of abundance spreads across detected community components. Global statistical verification tags note variations via automated distribution-appropriate significance testing.

### 2. Beta Diversity Confirms Phenotypic Clustering (p = 0.013)

<img width="2250" height="1950" alt="image" src="https://github.com/user-attachments/assets/805e9acf-de58-4f3f-a391-6151aa4ad877" />

**Fig:** Beta Diversity Principal Coordinate Analysis (PCoA) of Microbial Communities. Principal Coordinate Analysis ordination plot derived from a pairwise Bray-Curtis dissimilarity matrix tracking species-level composition across study cohorts: Non-Inflammatory Bowel Disease (NonIBD, green), Crohn’s Disease (CD, orange), and Ulcerative Colitis (UC, blue). Individual points map single biological profiling replicates, with dashed lines illustrating the 95% confidence ellipse boundaries for each specific clinical group. The horizontal (Axis 1) and vertical (Axis 2) coordinates explain the top dimensions of community variance within the ecosystem. The embedded statistical annotation summarizes multivariate significance metrics derived from global permutation profiling.

PCoA ordination revealed distinct sample grouping based on disease status. Clinical group membership accounted for **23.2%** of total community-wide variation (R² = 0.232, global PERMANOVA p = 0.058). The PERMDISP check was non-significant (p = 0.422), confirming uniform group variances and validating the spatial patterns in the PCoA plot.

**Pairwise Breakthrough:** Post-hoc pairwise comparisons revealed a highly significant compositional split specifically between the **Crohn's Disease and Ulcerative Colitis** groups (R² = 0.258, p = 0.013), indicating that CD and UC drive the gut microbiota toward fundamentally distinct dysbiotic states.

### 3. Highly Specific Differential Biomarkers Discovered

Mass feature screening and post-hoc pairwise Dunn's tests identified the key organisms driving separation on the PCoA plot (P.adj < 0.05):

<img width="2550" height="2100" alt="image" src="https://github.com/user-attachments/assets/6c2c1277-1dc5-4baf-a80c-ab3c1b46beb2" />

**Fig:** Relative Abundance Fluctuations of Key Identified Microbial Biomarkers. Faceted box-and-whisker plots tracking the distribution of the most significantly altered species identified through mass Kruskal-Wallis feature screening and post-hoc pairwise Dunn's testing. Individual points reflect independent biological replicates categorized into Non-Inflammatory Bowel Disease (NonIBD, green), Crohn’s Disease (CD, orange), and Ulcerative Colitis (UC, blue) cohorts. The vertical y-axes dynamically scale to track relative abundance percentages (%) per organism. Statistically verified contrasts reveal a clear division in dysbiosis styles: a pronounced expansion of specific Bacteroides species characterizes the CD landscape, whereas a distinct depletion of protective Coprococcus comes marks the UC cohort.

**`Bacteroides uniformis` & `Bacteroides fragilis` (CD Drivers)**
*B. uniformis* was significantly enriched in Crohn's Disease compared to both NonIBD controls (P.adj = 0.024) and UC (P.adj = 0.021). *B. fragilis* further highlighted the split between phenotypes, showing significant elevation exclusively in CD compared to UC (P.adj = 0.025).

**`Coprococcus comes` (UC Attrition Signature)**
This beneficial, butyrate-producing commensal was significantly depleted in the UC cohort compared to both NonIBD controls (P.adj = 0.041) and CD patients (P.adj = 0.020).

**`Alistipes finegoldii` (Health Guard)**
This organism served as an indicator of a healthy microbiome baseline, showing significant depletion in the active Crohn's Disease cohort (P.adj = 0.009, Z = -2.95).

---

## Environment Requirements & Replication Quickstart

Requires an active R environment (v4.3+) and a Python environment with the listed dependencies.

### 1. Install Dependencies

```bash
# Python packages
pip install tqdm multiqc requests

# R packages
R -e "install.packages(c('tidyverse', 'vegan', 'FSA', 'ggpubr', 'remotes'), repos='http://cran.us.r-project.org')"
R -e "remotes::install_github('pmartinezarbizu/pairwiseAdonis/pairwiseAdonis')"
```

### 2. Execute Downstream Pipeline

```bash
cd scripts/

# Step A: Parse raw profiles into a clean wide matrix
Rscript 04_tsv_file_merger.R

# Step B: Generate abundance and alpha diversity charts
Rscript 05_plot_abundance.R
Rscript 06_alpha_diversity.R

# Step C: Beta diversity calculations and PERMANOVA tests
Rscript 07_beta_diversity_permanova.R

# Step D: Mass feature screening for biomarker identification
Rscript 08_biomarker_kruskal.R
```

All statistical tables, matrix files, and publication-ready figures are saved automatically to `data/results/`.

---

## About the Author

**Milind Shrivastava**
M.Sc. Biotechnology — Thapar Institute of Engineering and Technology (TIET)
CGPA: 9.15 | GATE 2026 Qualified

Research focus: Shotgun metagenomics, microbiome comparative ecology, computational biology pipelines, and clinical biomarker discovery.

Feel free to connect or open an issue for collaboration or pipeline inquiries.
