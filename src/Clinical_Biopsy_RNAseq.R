###############################################################################
# Clinical Biopsy RNA-seq Analysis: IGSF8 Expression and Pathway Enrichment
#
# Description:
#   Analyzes bulk RNA-seq data from clinical biopsy samples (GV20-0251-100
#   Phase I trial) to characterize IGSF8 expression and treatment-associated
#   transcriptomic changes in tumor biopsies.
#
#   Analyses:
#     1. IGSF8 expression by clinical response (all baseline; melanoma subset)
#     2. KEGG pathway enrichment (GSEA) in paired pre/on-treatment melanoma
#     3. Leading-edge gene expression heatmap for enriched pathway signatures
#
# Input:
#   data/Clinical_Biopsy_RNAseq/gene_TPM_of_all_samples.txt
#   data/Clinical_Biopsy_RNAseq/GV20-0251-100_Biopsy_RNAseq_meta_*.xlsx
#   data/Clinical_Biopsy_RNAseq/rnaseq_analysis_paired_melanoma_gsea_kegg.csv
#   data/Clinical_Biopsy_RNAseq/signature_marker_genes_list.txt
#
# Output:
#   results/  (PNG figures at 600 dpi)
###############################################################################

# --- Setup -------------------------------------------------------------------

source("src/Biopsy_RNAseq_analysis_functions.R")

INPUT_DIR  <- file.path("data", "Clinical_Biopsy_RNAseq")
OUTPUT_DIR <- "results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

GENE_OF_INTEREST <- "IGSF8"

# --- Helper: load and annotate sample metadata -------------------------------

load_metadata <- function(meta_filename) {
  readxl::read_xlsx(file.path(INPUT_DIR, meta_filename)) %>%
    column_to_rownames(var = "Sample") %>%
    mutate(Cancer = sub(" \\(.*", "", `Tumor types`))
}

# =============================================================================
# 0. Load and preprocess TPM expression data
# =============================================================================

tpm_raw <- read.table(
  file.path(INPUT_DIR, TPM_FILE), header = TRUE, check.names = FALSE
)

tpm_summarized <- tpm_raw %>%
  dplyr::select(-c(Gene_Name, Chromosome, Strand, Start, End)) %>%
  group_by(Gene_ID) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames(var = "Gene_ID")

tpm_log2 <- prepare_expression_data(expression_data = tpm_summarized)

# =============================================================================
# 1. IGSF8 expression by clinical response group
# =============================================================================

# --- 1a. All baseline samples (3 PR + 1 SD + 6 PD) --------------------------

meta_baseline <- load_metadata(META_FILE_BASELINE)

p_baseline <- plot_expression_boxplot(
  markers = GENE_OF_INTEREST,
  tpm_df  = tpm_log2,
  meta    = meta_baseline,
  group   = "Clinical_benefit"
)

ggsave(
  file.path(OUTPUT_DIR, "baseline_all_IGSF8_expression_boxplot.png"),
  plot = p_baseline, device = png,
  width = 3.5, height = 4, units = "in", dpi = 600
)

# --- 1b. Melanoma baseline samples (3 PR + 4 PD) ----------------------------

meta_melanoma_bl <- load_metadata(META_FILE_MELANOMA_BL)

p_melanoma <- plot_expression_boxplot(
  markers = GENE_OF_INTEREST,
  tpm_df  = tpm_log2,
  meta    = meta_melanoma_bl,
  group   = "BOR"
)

ggsave(
  file.path(OUTPUT_DIR, "baseline_melanoma_IGSF8_expression_boxplot.png"),
  plot = p_melanoma, device = png,
  width = 3.5, height = 4, units = "in", dpi = 600
)

# =============================================================================
# 2. KEGG pathway enrichment (GSEA) — paired melanoma biopsies
# =============================================================================

gsea_kegg <- read.csv(
  file.path(INPUT_DIR, "rnaseq_analysis_paired_melanoma_gsea_kegg.csv")
)

p_gsea <- plot_gsea_dotplot(gsea_kegg, n_pathways = 15, p_cut = 0.05)

ggsave(
  file.path(OUTPUT_DIR, "paired_melanoma_gsea_kegg_dotplot.png"),
  plot = p_gsea, device = png,
  width = 10, height = 4, units = "in", dpi = 600, bg = "white"
)

# =============================================================================
# 3. Leading-edge gene expression heatmap (paired melanoma)
# =============================================================================

meta_paired <- load_metadata(META_FILE_PAIRED_MELANOMA)

# Paired sample IDs (matched by patient: pre[i] <-> post[i])
pre_treatment_ids  <- c("25431XR-01-04", "25431XR-01-11",
                         "25431XR-01-02", "25431XR-01-09")
post_treatment_ids <- c("25431XR-01-03", "25431XR-01-10",
                         "25431XR-01-01", "25431XR-01-08")

# --- 3a. Load signature gene list and rank by group difference ---------------

markers_df <- read.table(
  file.path(INPUT_DIR, "signature_marker_genes_list.txt"),
  sep = "\t", header = TRUE
)

# Z-score genes across paired samples, then compute per-group means
tpm_signature <- tpm_log2[unique(markers_df$Gene), rownames(meta_paired)]
tpm_signature_z <- t(apply(tpm_signature, 1, scale))
colnames(tpm_signature_z) <- rownames(meta_paired)

group_means <- tpm_signature_z %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("Gene") %>%
  reshape2::melt(id.vars = "Gene", variable.name = "Sample") %>%
  mutate(Timepoint = meta_paired$Timepoint[match(Sample, rownames(meta_paired))]) %>%
  group_by(Gene, Timepoint) %>%
  summarise(average = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = Timepoint, values_from = average)

# Order genes by signature group, then by on-treatment vs baseline difference
gene_order_df <- group_means %>%
  mutate(
    Signature = markers_df$Signature[match(Gene, markers_df$Gene)],
    diff      = `On-treatment` - Baseline
  ) %>%
  arrange(Signature, desc(diff))

ordered_genes      <- gene_order_df$Gene[gene_order_df$Gene %in% markers_df$Gene]
ordered_signatures <- gene_order_df$Signature[gene_order_df$Gene %in% markers_df$Gene]

signature_order <- c("Antigen_Presentation", "LAMP3_DC",
                     "NK_Cytoxicity", "PPAR_Signaling")

# --- 3b. Prepare z-scored expression matrix ----------------------------------

expr_matrix <- tpm_log2[ordered_genes, c(pre_treatment_ids, post_treatment_ids)]
expr_matrix_z <- t(apply(expr_matrix, 1, scale))
colnames(expr_matrix_z) <- c(pre_treatment_ids, post_treatment_ids)

# --- 3c. Build ComplexHeatmap annotations ------------------------------------

# Top annotation: Baseline vs On-treatment
timepoint_labels <- c(
  rep("Baseline",     length(pre_treatment_ids)),
  rep("On-treatment", length(post_treatment_ids))
)
timepoint_colors <- c("Baseline" = "#999999", "On-treatment" = "#E69F00")

top_anno <- HeatmapAnnotation(
  Group = timepoint_labels,
  col   = list(Group = timepoint_colors),
  annotation_name_side = "left",
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Group = list(labels_gp = gpar(fontsize = 12, fontfamily = "Arial"))
  )
)

# Left annotation: pathway/signature membership
n_signatures <- length(unique(ordered_signatures))
signature_colors <- setNames(
  brewer.pal(n = max(3, n_signatures), name = "Set1")[seq_len(n_signatures)],
  unique(ordered_signatures)
)

left_anno <- rowAnnotation(
  Signature = factor(ordered_signatures, levels = signature_order),
  col = list(Signature = signature_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Signature = list(labels_gp = gpar(fontsize = 12, fontfamily = "Arial"))
  )
)

# --- 3d. Draw heatmap -------------------------------------------------------

# Dashed line separating pre- and on-treatment columns
separator_layer <- function(j, i, x, y, width, height, fill, ...) {
  boundary_col <- length(pre_treatment_ids)
  if (boundary_col %in% j) {
    idx    <- which(j == boundary_col)
    x_edge <- x[idx] + 0.5 * width[idx]
    grid.lines(
      x  = unit(c(x_edge, x_edge), "native"),
      y  = unit(c(0, 1), "npc"),
      gp = gpar(col = "black", lwd = 2, lty = "dashed")
    )
  }
}

ht <- Heatmap(
  expr_matrix_z,
  name             = "Z-scored\nlog2TPM",
  top_annotation   = top_anno,
  left_annotation  = left_anno,
  row_split        = factor(ordered_signatures, levels = signature_order),
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  show_row_names   = TRUE,
  show_column_names = TRUE,
  column_labels    = rep("", ncol(expr_matrix_z)),
  row_names_gp     = gpar(fontfamily = "Arial", fontsize = 12),
  column_names_gp  = gpar(fontfamily = "Arial", fontsize = 12),
  row_title        = NULL,
  show_row_dend    = FALSE,
  gap              = unit(1, "mm"),
  layer_fun        = separator_layer,
  heatmap_legend_param = list(
    direction = "vertical",
    labels_gp = gpar(fontsize = 12, fontfamily = "Arial")
  )
)

png(
  file.path(OUTPUT_DIR, "paired_melanoma_leading_edge_heatmap.png"),
  width = 6, height = 8, res = 600, units = "in", bg = "white"
)
draw(ht, merge_legends = TRUE, heatmap_legend_side = "right")
dev.off()

message("All figures saved to: ", OUTPUT_DIR)
