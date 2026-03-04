###############################################################################
# Pan-Myeloid scRNA-seq Analysis: IGSF8 Expression in myeloid cells across Cancer Types
#
# Description:
#   Quantifies IGSF8 expression across myeloid cell sub-clusters in tumor
#   tissues from multiple cancer types using single-cell RNA-seq data from
#   the TISCH2 database (Tumor Immune Single-cell Hub).
#
#   Workflow:
#     1. Load per-cell expression data exported from Pan-Myeloid Atlas
#     2. Filter for tumor-infiltrating myeloid cells (Tissue == "T")
#     3. Compute mean IGSF8 expression per myeloid sub-cluster per cancer type
#     4. Z-score normalize within each cancer type
#     5. Generate a heatmap visualization
#
# Data Source:
#   Pan-Myeloid Atlas (http://panmyeloid.cancer-pku.cn/)
#   Cheng S, Li Z, ..., Zhang Z. A pan-cancer single-cell transcriptional
#   atlas of tumor infiltrating myeloid cells. Cell, 2021.
#   CSV files: per-cell expression with sub-cluster annotations
#
# Input:  data/PanMyeloid_scRNAseq/PanMyeloid_*_DataTablePlot_*.csv
# Output: results/IGSF8_Pan-Myeloid_across_cancers.png
###############################################################################

# --- Dependencies -----------------------------------------------------------

library(tidyverse)
library(scales)

# --- Configuration -----------------------------------------------------------

GENE_OF_INTEREST <- "IGSF8"

# Paths relative to the project root
DATA_DIR    <- file.path("data", "PanMyeloid_scRNAseq")
OUTPUT_DIR  <- "results"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load and aggregate per-cell expression data -------------------------

pan_myeloid_files <- list.files(
  path       = DATA_DIR,
  pattern    = "^PanMyeloid_.*\\.csv$",
  full.names = TRUE
)

cluster_means <- map_dfr(pan_myeloid_files, function(filepath) {
  cancer_type <- str_split(basename(filepath), "_", simplify = TRUE)[2]

  read.csv(filepath) %>%
    filter(Tissue == "T") %>%
    group_by(Sub_Cluster) %>%
    summarise(
      mean_expression = mean(.data[[GENE_OF_INTEREST]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Cancer = cancer_type)
})

# --- 2. Z-score normalization within each cancer type -----------------------

expression_zscore <- cluster_means %>%
  group_by(Cancer) %>%
  mutate(zscore = as.numeric(scale(mean_expression))) %>%
  ungroup()

# --- 3. Order sub-clusters by lineage ---------------------------------------

all_clusters <- sort(unique(expression_zscore$Sub_Cluster))

cluster_order <- c(
  grep("cDC|pDC",   all_clusters, value = TRUE),
  grep("Macro",     all_clusters, value = TRUE),
  grep("Mono",      all_clusters, value = TRUE),
  grep("Mast",      all_clusters, value = TRUE),
  grep("Myeloid",   all_clusters, value = TRUE)
)

expression_zscore <- expression_zscore %>%
  mutate(Sub_Cluster = factor(Sub_Cluster, levels = cluster_order))

# --- 4. Heatmap visualization -----------------------------------------------

heatmap_plot <- ggplot(expression_zscore,
                       aes(x = Cancer, y = Sub_Cluster, fill = zscore)) +
  geom_tile(color = "white", linewidth = 0.2) +
  coord_fixed() +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "#F7F7F7",
    high     = "#B2182B",
    midpoint = 0,
    name     = "Z-scored\nlog2CPM",
    limits   = c(-3, 3),
    oob      = squish,
    na.value = "white"
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x        = element_text(angle = 75, hjust = 1, face = "bold",
                                      family = "Arial"),
    axis.text.y        = element_text(face = "bold", family = "Arial"),
    axis.ticks         = element_blank(),
    panel.grid         = element_blank(),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title       = element_text(face = "italic"),
    legend.title.position = "top",
    legend.justification = "center",
    legend.key.width   = unit(1, "cm"),
    legend.key.height  = unit(0.3, "cm"),
    plot.margin        = margin(10, 5, 5, 5, unit = "mm")
  )

# --- 5. Save figure ----------------------------------------------------------

output_file <- file.path(OUTPUT_DIR,
                         paste0(GENE_OF_INTEREST, "_Pan-Myeloid_across_cancers.png"))

ggsave(
  filename = output_file,
  plot     = heatmap_plot,
  device   = png,
  width    = 7,
  height   = 8,
  units    = "in",
  dpi      = 600,
  bg       = "white"
)

message("Saved: ", output_file)
