# IGSF8-PartA

Code and data for reproducing the analyses and figures in the IGSF8 Phase I (Part A) publication.

## Repository Structure

```
├── data/
│   ├── Clinical_Biopsy_RNAseq/       # Bulk RNA-seq from clinical biopsy samples
│   │   ├── gene_TPM_of_all_samples.txt
│   │   ├── GV20-0251-100_Biopsy_RNAseq_meta_*.xlsx
│   │   ├── rnaseq_analysis_paired_melanoma_gsea_kegg.csv
│   │   └── signature_marker_genes_list.txt
│   └── PanMyeloid_scRNAseq/          # Pan-cancer myeloid scRNA-seq (from Pan-Myeloid Atlas)
│       └── PanMyeloid_*_DataTablePlot_*.csv
├── src/
│   ├── Clinical_Biopsy_RNAseq.R              # Main analysis script for clinical biopsy RNA-seq
│   ├── Biopsy_RNAseq_analysis_functions.R    # Helper functions for biopsy RNA-seq analysis
│   └── PanMyeloid_scRNASeq.R                 # Pan-myeloid scRNA-seq IGSF8 expression analysis
└── results/                           # Output figures (600 dpi PNG)
```

## Analyses

### 1. Clinical Biopsy RNA-seq (`src/Clinical_Biopsy_RNAseq.R`)

Analyzes bulk RNA-seq data from tumor biopsy samples in the GV20-0251-100 Phase I trial:

- **IGSF8 expression** by clinical response group (all baseline; melanoma baseline)
- **KEGG pathway enrichment** (GSEA) in paired pre/on-treatment melanoma biopsies
- **Leading-edge gene heatmap** for enriched pathway signatures

### 2. Pan-Myeloid scRNA-seq (`src/PanMyeloid_scRNASeq.R`)

Quantifies IGSF8 expression across myeloid cell sub-clusters in tumor tissues from multiple cancer types, using data from the [Pan-Myeloid Atlas](http://panmyeloid.cancer-pku.cn/) (Cheng S, Li Z, ..., Zhang Z. *Cell*, 2021).

## Requirements

R (>= 4.0) with the following packages:

```r
# CRAN
install.packages(c("tidyverse", "readxl", "reshape2", "RColorBrewer",
                    "rstatix", "scales", "extrafont"))

# Bioconductor
BiocManager::install(c("ComplexHeatmap", "clusterProfiler", "org.Hs.eg.db"))
```

## Usage

All scripts use **relative paths** and should be run from the project root:

```bash
cd IGSF8-PartA

# Pan-myeloid scRNA-seq analysis
Rscript src/PanMyeloid_scRNASeq.R

# Clinical biopsy RNA-seq analysis
Rscript src/Clinical_Biopsy_RNAseq.R
```

Output figures are saved to `results/`.
