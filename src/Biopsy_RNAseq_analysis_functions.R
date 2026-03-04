###############################################################################
# Helper Functions for Clinical Biopsy RNA-seq Analysis
#
# Description:
#   Shared utility functions used by Clinical_Biopsy_RNAseq.R, including:
#     - Gene ID conversion (Ensembl -> Symbol) with protein-coding filtering
#     - Expression data preprocessing (log2 TPM)
#     - Boxplot visualization for gene expression by clinical response group
#     - GSEA dotplot visualization
#
# Usage:
#   source("src/Biopsy_RNAseq_analysis_functions.R")
###############################################################################

# --- Dependencies -----------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)
library(grid)
library(rstatix)
library(clusterProfiler)
library(org.Hs.eg.db)
library(extrafont)

# --- Data file paths (relative to project root) -----------------------------

TPM_FILE                   <- "gene_TPM_of_all_samples.txt"
META_FILE_BASELINE         <- "GV20-0251-100_Biopsy_RNAseq_meta_baseline_2025Oct.xlsx"
META_FILE_MELANOMA_BL      <- "GV20-0251-100_Biopsy_RNAseq_meta_Melanoma_baseline_2025Oct.xlsx"
META_FILE_PAIRED_MELANOMA  <- "GV20-0251-100_Biopsy_RNAseq_meta_paired_Melanoma_2025Oct.xlsx"

# --- Gene ID conversion -----------------------------------------------------

#' Convert Ensembl gene IDs to HGNC symbols and aggregate by symbol.
#' Retains protein-coding and "other" biotype genes.
#'
#' @param expression_data Matrix or data.frame with Ensembl IDs as rownames.
#' @return Data.frame with gene symbols as rownames (duplicates summed).
convert_gene_ids <- function(expression_data,
                             from_type = "ENSEMBL",
                             to_type   = c("SYMBOL", "GENETYPE")) {

  gene_mapping <- bitr(
    rownames(expression_data),
    fromType = from_type,
    toType   = to_type,
    OrgDb    = org.Hs.eg.db
  )

  keep_ids <- gene_mapping$ENSEMBL[
    gene_mapping$GENETYPE %in% c("protein-coding", "other")
  ]

  expression_data <- expression_data[rownames(expression_data) %in% keep_ids, ]
  expression_data$symbol <- gene_mapping$SYMBOL[
    match(rownames(expression_data), gene_mapping[[from_type]])
  ]

  df_summed <- expression_data %>%
    group_by(symbol) %>%
    summarise(across(everything(), sum)) %>%
    column_to_rownames(var = "symbol")

  return(df_summed)
}

# --- Expression data preprocessing ------------------------------------------

#' Convert raw TPM matrix to log2-transformed, gene-symbol-indexed data.
#'
#' @param expression_data TPM matrix with Ensembl IDs as rownames.
#' @return log2(TPM + 1) matrix with gene symbols as rownames.
prepare_expression_data <- function(expression_data) {
  expression_data <- convert_gene_ids(expression_data)
  expression_data <- log2(expression_data + 1)
  return(expression_data)
}

# --- Boxplot: gene expression by clinical group -----------------------------

#' Boxplot comparing gene/signature expression between clinical response groups.
#'
#' @param markers Character vector of gene symbol(s) to average.
#' @param tpm_df  log2(TPM+1) matrix (genes x samples).
#' @param meta    Sample metadata data.frame (samples as rownames).
#' @param group   Column name in meta for grouping ("BOR" or "Clinical_benefit").
#' @return ggplot object.
plot_expression_boxplot <- function(markers, tpm_df, meta, group = "BOR") {

  marker_tpm <- na.omit(tpm_df[markers, , drop = FALSE])
  marker_aggr_df <- data.frame(
    Sample          = colnames(marker_tpm),
    Signature_score = as.numeric(colMeans(marker_tpm))
  ) %>%
    merge(meta %>% rownames_to_column(var = "Sample"), by = "Sample")

  # Wilcoxon rank-sum test
  stat_test <- marker_aggr_df %>%
    wilcox_test(as.formula(paste("Signature_score ~", group))) %>%
    add_significance()
  stat_test$y.position <- max(marker_aggr_df$Signature_score) + 0.5
  #message("Wilcoxon test p-value: ", stat_test$p)

  # Configure group aesthetics
  if (group == "BOR") {
    fill_order <- c("PD", "PR")
    bor_order  <- fill_order
  } else if (group == "Clinical_benefit") {
    marker_aggr_df <- marker_aggr_df %>%
      mutate(Clinical_benefit = case_when(
        Clinical_benefit == "Benefit"     ~ "PR + SD",
        Clinical_benefit == "Non-Benefit" ~ "PD",
        .default = Clinical_benefit
      ))
    fill_order <- c("PD", "PR + SD")
    bor_order  <- c("PR", "SD", "PD")
  }

  marker_aggr_df[[group]] <- factor(marker_aggr_df[[group]], levels = fill_order)
  marker_aggr_df[["BOR"]] <- factor(marker_aggr_df[["BOR"]], levels = bor_order)

  bor_colors <- c("PD" = "#D95F02", "SD" = "blue", "PR" = "#33A02C")
  bor_shapes <- c("PD" = 21, "PR" = 21, "SD" = 21)

  p <- ggplot(marker_aggr_df, aes(x = .data[[group]], y = Signature_score)) +
    stat_summary(
      fun = median, fun.min = median, fun.max = median,
      geom = "crossbar", width = 0.5, linewidth = 0.6, color = "black"
    ) +
    geom_jitter(
      aes(fill = BOR, shape = BOR),
      width = 0.15, size = 3.5, stroke = 0.5, color = "black"
    ) +
    scale_fill_manual(values = bor_colors) +
    scale_shape_manual(values = bor_shapes) +
    labs(
      y = expression(paste("Gene expression (", log[2], " TPM)")),
      x = NULL
    ) +
    scale_y_continuous(limits = c(0, 8), expand = c(0, 0)) +
    theme_classic(base_size = 14, base_family = "Arial") +
    theme(
      axis.line        = element_line(color = "black", linewidth = 1),
      axis.ticks       = element_line(color = "black", linewidth = 1),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text        = element_text(color = "black", size = 12),
      axis.title.y     = element_text(color = "black", size = 14,
                                      margin = margin(r = 10)),
      legend.title      = element_text(color = "black", size = 11, face = "bold"),
      legend.background = element_blank(),
      legend.key.height = unit(0.6, "cm"),
      legend.text       = element_text(color = "black", size = 11),
      aspect.ratio      = 1.3
    )

  return(p)
}

# --- Dotplot: GSEA pathway enrichment results --------------------------------

#' Dotplot of GSEA results showing top enriched/depleted pathways.
#'
#' @param gsea_results Data.frame with columns: pathway, NES, padj, size.
#' @param n_pathways   Max number of pathways per direction (up/down).
#' @param p_cut        Adjusted p-value cutoff.
#' @return ggplot object.
plot_gsea_dotplot <- function(gsea_results, n_pathways = 10, p_cut = 0.05) {

  top_pathways <- bind_rows(
    gsea_results %>% filter(padj <= p_cut, NES > 0) %>% slice_max(NES, n = n_pathways),
    gsea_results %>% filter(padj <= p_cut, NES < 0) %>% slice_min(NES, n = n_pathways)
  ) %>%
    mutate(pathway = gsub("KEGG_|HALLMARK_|hsa[0-9]+ |mmu[0-9]+ ", "", pathway))

  p <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway, NES))) +
    geom_point(aes(size = size, color = padj)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(
      x     = "Normalized Enrichment Score (NES)",
      y     = NULL,
      size  = "Gene Set Size",
      color = "Adjusted P-value"
    ) +
    theme_classic(base_size = 20, base_family = "Arial") +
    theme(
      axis.line    = element_line(linewidth = 0.5, color = "black"),
      axis.text    = element_text(color = "black", size = 14),
      axis.title.x = element_text(face = "bold", size = 14,
                                  margin = margin(t = 5)),
      legend.position   = "right",
      legend.direction   = "vertical",
      legend.title       = element_text(size = 12, face = "bold"),
      legend.text        = element_text(size = 10),
      legend.key.size    = unit(0.6, "cm"),
      legend.spacing.y   = unit(0.2, "cm"),
      legend.margin      = margin(l = 5),
      legend.box.margin  = margin(l = 10),
      plot.margin        = margin(10, 10, 10, 10)
    )

  return(p)
}
