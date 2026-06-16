# Comparative Metagenomic & Functional Profiling of the Human Gut Microbiome Across IBD Phenotypes

> A high-throughput, end-to-end computational biology framework for isolating species-level taxonomic signatures and functional metabolic profiles of the human gut microbiota across three clinical cohorts — <strong>NonIBD</strong>, <strong>Crohn's Disease (CD)</strong>, and <strong>Ulcerative Colitis (UC)</strong> — using the HMP2 (iHMP) IBD Multi-omics Database (IBDMDB).

**Author:** Milind Shrivastava | M.Sc. Biotechnology, Thapar Institute of Engineering and Technology (TIET) | CGPA: 9.15 | GATE 2026 Qualified

---

## Table of Contents

- [Overview](#overview)
- [Repository Architecture](#repository-architecture)
- [Methodological Workflow](#methodological-workflow)
- [Execution Pipeline](#execution-pipeline)
- [Key Results & Biostatistical Insights](#key-results--biostatistical-insights)
  - [Alpha Diversity](#1-alpha-diversity)
  - [Beta Diversity](#2-beta-diversity--the-variance-paradox)
  - [Biomarker Discovery](#3-high-confidence-biomarker-discovery)
  - [Functional Pathway Profiling](#4-functional-metabolic-pathway-dysregulation)
- [Environment & Dependencies](#environment--dependencies)
- [Notes on Path Alignment](#notes-on-path-alignment)

---

## Overview

This pipeline scales from raw paired-end FASTQ pre-processing to downstream mass biostatistical modeling. It integrates:

- **fastp** for quality trimming and adapter removal
- **Bowtie2** for host-genome (GRCh38) contamination subtraction
- **MetaPhlAn v3.0** for species-level taxonomic profiling
- **HUMAnN3** for functional metabolic pathway quantification
- A full **R-based analytics suite** for alpha/beta diversity, biomarker discovery (MaAsLin2), and pathway dysregulation analysis

The complete cohort spans **N = 1,605** individual biological and longitudinal microbiome profiles across three clinical phenotypes.

---

## Repository Architecture

```
human-gut-shotgun-metagenomics/
├── scripts/
│   ├── 01_download_data.py               # Automated SRA Toolkit raw sequence fetch
│   ├── 02_qc_and_profile.py              # fastp quality trimming & MultiQC reporting
│   ├── 03_host_removal.py                # Bowtie2 GRCh38 host subtraction
│   ├── 04_tsv_file_merger.R              # Upstream cohort profile aggregation
│   ├── 05_plot_abundance.R               # Structural abundance plot generation
│   ├── 06_alpha_diversity.R              # Intrasample ecological complexity
│   ├── 07_beta_diversity_permanova.R     # Intersample community distance modeling
│   ├── 08_biomarker_kruskal.R            # Non-parametric feature screening
│   └── download_result_files.py          # Validation manifest fetch helper
└── data/
    ├── metadata/                         # Cohort accessions: cd_acc.txt, uc_acc.txt, nonibd_acc.txt
    ├── raw/                              # Raw paired-end FASTQ streams
    ├── fastp_trimmed/                    # Quality-filtered, adapter-trimmed sequences
    ├── bowtie2_host_removed/             # Host-subtracted reads by phenotype (CD, UC, NonIBD)
    ├── results/                          # MultiQC logs, merged profiles, master metrics
    └── Result2/                          # Downstream Multi-Omics Biostatistical Suite (R Engine)
        ├── 01_Data_manipulation.R        # Taxonomic regex parsing & matrix transposition
        ├── 02_Alpha_Diversity.R          # Shannon / Simpson / Richness calculation
        ├── 03_Beta_Diversity.R           # Bray-Curtis matrix & PERMANOVA modeling
        ├── 04_Biomarker_analysis.R       # Multivariable linear screening via MaAsLin2
        ├── 05_Relative_abundance.R       # Univariate non-parametric abundance screening
        ├── 06_functional_profile_.R      # HUMAnN3 CPM parsing & species stratification
        ├── 07_data_cohorting.R           # Metadata splitting into cohort-specific tables
        ├── 08_Pathway_Comparison_Stats.R # Kruskal-Wallis & Wilcoxon pathway modeling (FDR)
        ├── 09_Pathway_comparision_boxplot.R  # Relative % transformation & boxplot rendering
        ├── Alpha_diversity/              # Shaded boxplot grids and alpha statistics
        ├── Beta_diversity/               # PCoA ordinations and pairwise PERMANOVA tables
        ├── maaslin2_biomarkers/          # GLM statistics, logs, and coefficient barplots
        ├── Taxonomic_Abundance_analysis/ # Clustered relative abundance profiles
        └── Functional_Profile_analysis/  # Dysregulated pathways and multi-panel boxplots
```

---

## Methodological Workflow

```
[Raw FASTQ Stream]
        │
        ▼
[fastp: QC & Adapter Trimming]
        │
        ▼
[Bowtie2: GRCh38 Host Subtraction]
        │
        ├─────────────────────────────────────┐
        ▼                                     ▼
[MetaPhlAn v3.0: Taxonomic Profiling]   [HUMAnN3: Functional Profiling]
        │                                     │
   ┌────┴────┬─────────────┐             ┌────┴────┐
   ▼         ▼             ▼             ▼         ▼
[Alpha     [Beta       [MaAsLin2     [Pathway   [Pathway
Diversity] Diversity]  Biomarkers]   Stats]     Boxplots]
```

---

## Execution Pipeline

> **Important:** All downstream R scripts must be executed from within `data/Result2/`, or file paths must be prefixed with `data/Result2/` when running from the project root.

### Step 1 — Raw Data Acquisition & QC

```bash
python scripts/01_download_data.py      # Fetch SRA FASTQ streams (cd_acc.txt, uc_acc.txt, nonibd_acc.txt)
python scripts/02_qc_and_profile.py     # fastp trimming (Q ≥ 20) + MultiQC dashboard
```

Adapter trimming, base-quality filtration, and overlapping base-mismatch correction are applied via `fastp`. Global performance telemetry is consolidated via MultiQC into `multiqc_report.html`.

### Step 2 — Host Contamination Subtraction

```bash
python scripts/03_host_removal.py       # Bowtie2 alignment against GRCh38; unmapped reads retained
```

Human background reads are aligned and discarded (`/dev/null`); only microbial outflow streams are preserved.

### Step 3 — Downstream Analytics Suite

Navigate into the R engine workspace first:

```bash
cd data/Result2/
```

**Taxonomic Modules:**

```bash
Rscript 01_Data_manipulation.R          # Species-level regex parsing; matrix transposition
Rscript 02_Alpha_Diversity.R            # Shannon / Simpson / Richness + Kruskal-Wallis + Wilcoxon (BH)
Rscript 03_Beta_Diversity.R             # Bray-Curtis PCoA + PERMANOVA (adonis2, 999 perms) + PERMDISP
Rscript 04_Biomarker_analysis.R         # MaAsLin2: TSS + LOG, prevalence ≥ 10%, q < 0.25
Rscript 05_Relative_abundance.R         # Parallelized Kruskal-Wallis loops, significance overlays
```

**Functional Pathway Modules:**

```bash
Rscript 06_functional_profile_.R        # HUMAnN3 CPM ingestion; unmapped/ungrouped read filtering
Rscript 07_data_cohorting.R             # Metadata split into cohort-specific tables
Rscript 08_Pathway_Comparison_Stats.R   # Mass Kruskal-Wallis + Wilcoxon, BH FDR (q < 0.05)
Rscript 09_Pathway_comparision_boxplot.R # CPM → Relative % transformation; boxplot rendering
```

All statistical tables, filtered matrices, and publication-ready figures are saved automatically to their respective subdirectories within `Result2/`.

---

## Key Results & Biostatistical Insights

### 1. Alpha Diversity

Early pilot trials (N = 15) yielded non-significant alpha diversity variations (p > 0.05), incorrectly suggesting uniform within-sample ecological complexity. Scaling to the complete HMP2 cohort (N = 1,605) overturned this inference entirely — a direct demonstration of the critical role of statistical power (1 − β) in high-variance human microbiome data.

<figure>
<img src="https://github.com/user-attachments/assets/11366294-3d5c-4424-9b67-749bb5080f99" alt="Four-panel alpha diversity boxplots across IBD cohorts" width="100%"/>
<figcaption><strong>Fig. 1.1 — Comprehensive Alpha Diversity Profiling Across IBD Phenotypes.</strong> Four-panel multi-metric evaluation of within-sample microbial community complexity across NonIBD (green), CD (orange), and UC (blue) cohorts. Individual points represent distinct biological profiling replicates overlaying box-and-whisker distributions. <strong>(A)</strong> Species Count (S): absolute observed species richness. <strong>(B)</strong> Shannon Index (H′): community entropy sensitive to rare taxa dropout. <strong>(C)</strong> Simpson Index (D): community dominance patterns driven by major taxonomic features. <strong>(D)</strong> Pielou's Evenness (J′): structural equitability of abundance spreads across detected community components. Global statistical significance assessed via distribution-appropriate testing.</figcaption>
</figure>

<br>

<figure>
<img src="https://github.com/user-attachments/assets/488ccc9e-6576-4593-92c6-2f4980708b2d" alt="Enhanced alpha diversity plot comparing IBD phenotypes" width="100%"/>
<figcaption><strong>Fig. 1.2 — Impact of IBD Phenotypes on Gut Microbial Alpha Diversity</strong> (N = 1,605; NonIBD = 725, CD = 604, UC = 276). Comparative ecological profiling across three clinical cohorts. <strong>(A)</strong> Observed Species Richness, <strong>(B)</strong> Shannon Diversity Index (H′), and <strong>(C)</strong> Simpson's Index (D). Each boxplot shows the median (center line), IQR (box bounds), and 1.5×IQR (whiskers), with individual sample outliers as black circles. Global community variance assessed via Kruskal-Wallis test; post-hoc pairwise comparisons via Wilcoxon rank-sum tests with Benjamini-Hochberg (BH) FDR correction. Significance brackets: ns = not significant; ∗ p ≤ 0.05; ∗∗∗ p ≤ 0.001; ∗∗∗∗ p ≤ 0.0001.</figcaption>
</figure>

**Alpha Diversity Statistics** (`Alpha_Diversity_Stats.csv`, N = 1,605):

| Metric | Global Kruskal p | CD vs. NonIBD | UC vs. NonIBD | UC vs. CD |
|:---|:---|:---|:---|:---|
| Observed Richness | 1.04 × 10⁻¹⁸ | 3.44 × 10⁻¹⁸ | 4.51 × 10⁻⁷ | 0.386 |
| Shannon Index (H') | 1.03 × 10⁻³ | 9.98 × 10⁻⁴ | 0.0479 | 0.977 |
| Simpson Index (D) | 0.1054 | 0.0883 | 0.5146 | 0.5146 |

**Key Takeaways:**

- **Total Species Collapse:** Observed Richness drops significantly in both active disease phenotypes (global p < 2.2 × 10⁻¹⁶), reflecting a widespread collapse in unique niche availability during active inflammation.
- **Simpson Dominance Resilience:** The global Simpson index remains statistically invariant (p = 0.1054). While intermediate and rare-abundance bacterial brackets collapse (captured by Shannon/Richness), the core high-abundance species retain their dominance positions.

---

### 2. Beta Diversity — The Variance Paradox

Pairwise Bray-Curtis dissimilarity matrices evaluated via PERMANOVA (`adonis2`, 999 permutations) revealed highly significant community separation (p = 0.001) alongside a compressed effect size (R² = 0.0094).

<figure>
<img src="https://github.com/user-attachments/assets/c7ed2770-67be-4c10-8c14-c885d9c828d0" alt="Early pilot PCoA ordination plot" width="100%"/>
<figcaption><strong>Fig. 2.1 — Pilot-Scale Beta Diversity PCoA of Microbial Communities.</strong> PCoA ordination derived from a pairwise Bray-Curtis dissimilarity matrix tracking species-level composition across NonIBD (green), CD (orange), and UC (blue). Individual points map single biological profiling replicates; dashed ellipses mark 95% confidence boundaries per cohort. Axes indicate the proportion of total community variance captured by the first two principal coordinates. Embedded annotation summarizes global permutation-based multivariate significance metrics.</figcaption>
</figure>

<br>

<figure>
<img src="https://github.com/user-attachments/assets/17169e03-7255-4f8e-8a7e-42ddda3ba306" alt="Full-cohort PCoA Beta Diversity ordination across IBD phenotypes" width="100%"/>
<figcaption><strong>Fig. 2.2 — Gut Microbial Community Dissimilarity (Beta Diversity) Across IBD Phenotypes</strong> (N = 1,605; NonIBD = 725, CD = 604, UC = 276). PCoA based on pairwise Bray-Curtis dissimilarity at the species level. Samples are color-coded by clinical cohort; shaded dashed ellipses outline 95% data confidence regions. Global community separation assessed via PERMANOVA (999 permutations): R² = 0.0094, F = 7.62, p = 0.001. Axis 1 and Axis 2 indicate the percentage of total variance explained by the first two principal coordinates.</figcaption>
</figure>

**Global PERMANOVA Results (N = 1,605):**

| Term | R² | F-statistic | p-value |
|:---|:---|:---|:---|
| Clinical Group | 0.0094 | 7.62 | 0.001 |
| Residual | 0.9906 | — | — |

**Post-Hoc Pairwise Splits (FDR-Corrected):**

| Comparison | R² | p-adjusted |
|:---|:---|:---|
| NonIBD vs. Crohn's Disease | 0.0080 | 0.001 |
| NonIBD vs. Ulcerative Colitis | 0.0054 | 0.001 |
| Crohn's Disease vs. Ulcerative Colitis | 0.0056 | 0.001 |

**Interpretation:** The R² drop from the pilot framework (23.2%) to the full cohort (0.94%) is a classic hallmark of population-scale human cohort studies. Massive unmeasured confounders — dietary fluctuations, antibiotic histories, age, genetics — introduce overwhelming residual noise. However, with residual df = 1,602, the F-statistic remains robustly elevated (F = 7.62, p = 0.001), confirming that structural separation across healthy, CD, and UC phenotypes is systematic and reproducible.

---

### 3. High-Confidence Biomarker Discovery

MaAsLin2 multivariable linear models (TSS normalization, LOG transformation, prevalence ≥ 10%, q < 0.25) isolated species-level biomarkers significantly enriched or depleted relative to the NonIBD baseline.

<figure>
<img src="https://github.com/user-attachments/assets/6980ce5b-003f-4209-808f-63db32a25bb5" alt="Relative abundance boxplots of key microbial biomarkers" width="100%"/>
<figcaption><strong>Fig. 3.1 — Relative Abundance Fluctuations of Key Identified Microbial Biomarkers.</strong> Faceted box-and-whisker plots tracking the most significantly altered species identified through mass Kruskal-Wallis feature screening and post-hoc pairwise Dunn's testing. Individual points reflect independent biological replicates across NonIBD (green), CD (orange), and UC (blue) cohorts. Y-axes scale dynamically to relative abundance (%). A pronounced expansion of specific <em>Bacteroides</em> species characterises the CD landscape, whereas a distinct depletion of protective <em>Coprococcus comes</em> marks the UC cohort.</figcaption>
</figure>

<br>

<figure>
<img src="https://github.com/user-attachments/assets/e5cd09b9-2e15-4536-afec-3884332563f4" alt="MaAsLin2 coefficient barplots for CD and UC biomarkers" width="100%"/>
<figcaption><strong>Fig. 3.2 — Multivariable Linear Association Discovery of Taxon-Specific Biomarkers Across IBD Phenotypes</strong> (N = 1,616). Differential abundance analysis via MaAsLin2 (LM; TSS normalization; LOG transformation; minimum prevalence 10%). Bars show the regression β-coefficient of each species relative to the NonIBD baseline. <strong>(Left)</strong> Taxa significantly shifted in CD. <strong>(Right)</strong> Taxa significantly shifted in UC. Positive values (red) denote enrichment; negative values (blue) denote depletion. Features ordered by effect size. FDR controlled via Benjamini-Hochberg; significance threshold q &lt; 0.25.</figcaption>
</figure>

**Key Biomarkers** (`Biomarker_Differential_Abundance_Results.csv`):

| Species | Phenotype | Coefficient | q-value | Clinical Relevance |
|:---|:---|:---|:---|:---|
| _Blautia wexlerae_ | CD Enriched | +1.17 | 2.40 × 10⁻⁸ | Core commensal structural shift |
| _Sutterella parvirubra_ | CD Enriched | +1.54 | 3.74 × 10⁻⁸ | Mucosal adhesion & epithelial wall degradation |
| _Collinsella intestinalis_ | UC Enriched | +3.21 | 2.40 × 10⁻⁸ | Tight-junction disruption, increased intestinal permeability, pro-inflammatory cascade activation |

---

### 4. Functional Metabolic Pathway Dysregulation

HUMAnN3 functional profiles (CPM) were transformed into relative percentage abundances per sample to eliminate sampling depth bias. Mass Kruskal-Wallis loops with post-hoc Wilcoxon rank-sum tests (BH FDR, q < 0.05) identified the following pathway-level disruptions:

<figure>
<img src="https://github.com/user-attachments/assets/0aa3a44e-3f88-4e34-abb2-4793fd4381e3" alt="Multi-panel boxplots of dysregulated metabolic pathways across IBD phenotypes" width="100%"/>
<figcaption><strong>Fig. 4.1 — Metagenomic Functional Pathway Reconfiguration Across IBD Phenotypes</strong> (N = 1,605; NonIBD = 725, CD = 604, UC = 276). Relative percentage abundance profiles of the top six dysregulated metabolic pathways derived from HUMAnN3 CPM data (unmapped/ungrouped reads excluded). Boxplots show the median (center line), IQR (box bounds), and 1.5×IQR (whiskers) across NonIBD (green), CD (red), and UC (blue). Global functional variance assessed via omnibus Kruskal-Wallis test (italicised annotation, top-left per facet). Pairwise comparisons against NonIBD via post-hoc Wilcoxon rank-sum tests with BH FDR correction. Significance brackets: ns = not significant; ∗ p ≤ 0.05; ∗∗ p ≤ 0.01; ∗∗∗ p ≤ 0.001; ∗∗∗∗ p ≤ 0.0001.</figcaption>
</figure>

| Pathway | Direction | CD FDR | UC FDR | Interpretation |
|:---|:---|:---|:---|:---|
| `TRNA-CHARGING-PWY` (aminoacyl-tRNA biosynthesis) | Depleted | 1.42 × 10⁻²⁵ | 5.16 × 10⁻¹⁰ | Universal collapse in core translation and protein synthesis capacity |
| `PWY-7234` (IMP biosynthesis III) | Depleted | — | 4.53 × 10⁻²³ | Purine precursor starvation (AMP/GMP); major reduction in DNA/RNA salvage capacity (mean CPM: 88.13 → 49.40) |
| `ANAEROFRUCAT-PWY` (homolactic fermentation) | Enriched | 0.0358 | — | Bloom of anaerobic pathobionts during mucosal inflammation, driven by oxygen gradient shifts from epithelial micro-bleeding |

---

## Environment & Dependencies

**Requirements:** R v4.3+ | Python 3.8+

### Install Python Dependencies

```bash
pip install tqdm multiqc requests pandas
```

### Install R Dependencies

```r
install.packages(
  c('tidyverse', 'vegan', 'ggpubr', 'broom', 'BiocManager', 'remotes', 'data.table'),
  repos = 'http://cran.us.r-project.org'
)
BiocManager::install('Maaslin2')
remotes::install_github('pmartinezarbizu/pairwiseAdonis/pairwiseAdonis')
```

---

## Notes on Path Alignment

Two alignment issues to verify before running the pipeline:

**1. Working Directory:** Downstream R scripts (`Result2/`) resolve data files relative to their own location. Either `cd data/Result2/` before running, or prepend `data/Result2/` to all file paths when executing from the project root.

**2. Output Filename Spelling:** `Result2/01_Data_manipulation.R` line 5 should read:

```r
output_path <- "Species_abundance.tsv"
```

Confirm this matches the actual output filename in `Taxonomic_Abundance_analysis/` to avoid manual renaming overhead.

---

*Research focus: Shotgun metagenomics · Multi-omics workflows · Microbiome comparative ecology · Computational pipeline design for clinical biomarker discovery*
