suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(grid)
  library(png)
  library(fmsb)
  library(cowplot)
})

utils::globalVariables(c(
  "ensembl", "transcript", "refseq", "prot_refseq", "gene_name",
  "structural_category", "associated_gene", "associated_transcript",
  "isoform", "id_type", "subcategory", "ref_exons", "diff_to_TSS",
  "diff_to_TTS", "gene"
))

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
project_root <- "/Users/tianyuan/Desktop/github_dev/tusco-paper"
data_root <- file.path(project_root, "figs/data")

# Output path
plot_dir <- file.path(project_root, "figs/fig-5/plot")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
out_pdf <- file.path(plot_dir, "new-fig-5c.pdf")

# ------------------------------------------------------------
# Helpers (borrowed/adapted from fig-5c-5d.R)
# ------------------------------------------------------------

read_tsv_safe <- function(file_path, ...) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  tryCatch(read_tsv(file_path, show_col_types = FALSE, ...),
           error = function(e) stop("Unable to read ", file_path, ": ", e$message))
}

load_tusco_reference <- function(tusco_path) {
  tusco_lines <- readLines(tusco_path)
  tusco_data_lines <- tusco_lines[!grepl("^#", tusco_lines)]
  tusco_temp_file <- tempfile()
  writeLines(tusco_data_lines, tusco_temp_file)
  tusco_df <- read_delim(
    tusco_temp_file, delim = "\t",
    col_names = c("ensembl","transcript","gene_name","gene_id_num","refseq","prot_refseq"),
    col_types = cols(.default = "c"), trim_ws = TRUE, show_col_types = FALSE
  )
  unlink(tusco_temp_file)
  tusco_df %>%
    mutate(
      ensembl = str_remove(ensembl, "\\.\\d+$"),
      transcript = str_remove(transcript, "\\.\\d+$"),
      refseq = str_remove(refseq, "\\.\\d+$"),
      prot_refseq = str_remove(prot_refseq, "\\.\\d+$")
    ) %>%
    select(ensembl, refseq, gene_name) %>%
    distinct()
}

compute_tusco_bigcats <- function(classification_data, tusco_ref_df) {
  patterns <- list(
    ensembl   = "^(ENSG|ENSMUSG)\\d{11}(\\.\\d+)?$",
    refseq    = "^(NM_|NR_|NP_)\\d{6,}$",
    gene_name = "^[A-Za-z0-9][A-Za-z0-9_-]*[A-Za-z0-9]?$"
  )

  class_clean <- classification_data %>%
    filter(structural_category != "fusion") %>%
    bind_rows(
      classification_data %>% filter(structural_category == "fusion") %>%
        separate_rows(associated_gene, sep = "_")
    ) %>%
    mutate(
      associated_gene = str_remove(associated_gene, "\\.\\d+$"),
      associated_transcript = str_remove(associated_transcript, "\\.\\d+$")
    ) %>%
    distinct(isoform, associated_gene, .keep_all = TRUE) %>%
    mutate(id_type = case_when(
      str_detect(associated_gene, patterns$ensembl)   ~ "ensembl",
      str_detect(associated_gene, patterns$refseq)    ~ "refseq",
      str_detect(associated_gene, patterns$gene_name) ~ "gene_name",
      TRUE ~ "unknown"
    ))

  top_id_type <- class_clean %>% count(id_type, sort = TRUE) %>% slice_head(n = 1) %>% pull(id_type)
  if (length(top_id_type) == 0 || is.na(top_id_type)) top_id_type <- "gene_name"

  class_tusco <- class_clean %>% filter(associated_gene %in% tusco_ref_df[[top_id_type]])

  TP_df <- class_tusco %>%
    filter(
      subcategory == "reference_match" |
        (subcategory == "mono-exon" &
           ref_exons == 1 &
           abs(diff_to_TSS) < 50 &
           abs(diff_to_TTS) < 50)
    )
  PTP_df <- class_tusco %>%
    filter(structural_category %in% c("full-splice_match", "incomplete-splice_match"),
           !(associated_gene %in% TP_df$associated_gene))

  found_ids <- union(TP_df$associated_gene, PTP_df$associated_gene)
  ref_ids <- unique(tusco_ref_df[[top_id_type]])
  FN_count <- sum(!(ref_ids %in% found_ids))

  tibble(
    TP = nrow(TP_df),
    PTP = nrow(PTP_df),
    FP = class_tusco %>% filter(structural_category %in% c(
      "novel_in_catalog", "novel_not_in_catalog", "genic",
      "fusion", "antisense", "intergenic", "genic_intron"
    )) %>% nrow(),
    FN = FN_count,
    id_type = top_id_type
  )
}

# Compute radar metrics following the approach used in figure3a-mouse.R
compute_tusco_metrics <- function(classification_data, tusco_ref_df) {
  patterns <- list(
    ensembl   = "^(ENSG|ENSMUSG)\\d{11}(\\.\\d+)?$",
    refseq    = "^(NM_|NR_|NP_)\\d{6,}$",
    gene_name = "^[A-Za-z0-9][A-Za-z0-9_-]*[A-Za-z0-9]?$"
  )

  class_clean <- classification_data %>%
    filter(structural_category != "fusion") %>%
    bind_rows(
      classification_data %>% filter(structural_category == "fusion") %>%
        separate_rows(associated_gene, sep = "_")
    ) %>%
    mutate(
      associated_gene = str_remove(associated_gene, "\\.\\d+$"),
      associated_transcript = str_remove(associated_transcript, "\\.\\d+$")
    ) %>%
    distinct(isoform, associated_gene, .keep_all = TRUE) %>%
    mutate(id_type = case_when(
      str_detect(associated_gene, patterns$ensembl)   ~ "ensembl",
      str_detect(associated_gene, patterns$refseq)    ~ "refseq",
      str_detect(associated_gene, patterns$gene_name) ~ "gene_name",
      TRUE ~ "unknown"
    ))

  id_summary <- class_clean %>% count(id_type, sort = TRUE)
  if (nrow(id_summary) > 0) {
    top_id_type <- id_summary %>% slice_max(n, n = 1) %>% pull(id_type)
  } else {
    top_id_type <- "gene_name"
  }

  tusco_ids <- unique(tusco_ref_df[[top_id_type]])
  tusco_ids <- tusco_ids[!is.na(tusco_ids) & tusco_ids != ""]
  rTUSCO <- length(tusco_ids)

  tusco_transcripts <- class_clean %>% filter(associated_gene %in% tusco_ids)

  tusco_RM <- tusco_transcripts %>%
    filter(
      subcategory == "reference_match" |
        (subcategory == "mono-exon" & ref_exons == 1 & abs(diff_to_TSS) < 50 & abs(diff_to_TTS) < 50)
    )

  TP_tusco <- tusco_RM
  PTP_tusco <- tusco_transcripts %>%
    filter(structural_category %in% c("full-splice_match", "incomplete-splice_match") &
             !associated_gene %in% TP_tusco$associated_gene)

  fsm_ism_count <- tusco_transcripts %>%
    filter(structural_category %in% c("full-splice_match", "incomplete-splice_match")) %>%
    nrow()

  unique_tp_ids <- length(unique(TP_tusco$associated_gene))
  unique_tp_ptp_ids <- length(unique(c(TP_tusco$associated_gene, PTP_tusco$associated_gene)))

  total_detected <- nrow(tusco_transcripts)
  tp_count <- nrow(TP_tusco)
  ptp_count <- nrow(PTP_tusco)

  Sn <- if (rTUSCO > 0) unique_tp_ids / rTUSCO else NA_real_
  nrPre <- if (total_detected > 0) tp_count / total_detected else NA_real_
  redundancy <- if (unique_tp_ptp_ids > 0) fsm_ism_count / unique_tp_ptp_ids else NA_real_
  inv_red <- if (!is.na(redundancy) && redundancy != 0) 1 / redundancy else 0
  FDR <- if (total_detected > 0) (total_detected - nrow(tusco_RM)) / total_detected else NA_real_
  one_minus_FDR <- if (!is.na(FDR)) (1 - FDR) else NA_real_
  PDR <- if (rTUSCO > 0) unique_tp_ptp_ids / rTUSCO else NA_real_
  rPre <- if (total_detected > 0) (tp_count + ptp_count) / total_detected else NA_real_

  tibble(
    id_type = top_id_type,
    Sn = Sn * 100,
    nrPre = nrPre * 100,
    `1/red` = inv_red * 100,
    `1-FDR` = one_minus_FDR * 100,
    PDR = PDR * 100,
    rPre = rPre * 100
  )
}

infer_tissue <- function(sample_basename) {
  if (startsWith(sample_basename, "B")) return("Brain")
  if (startsWith(sample_basename, "K")) return("Kidney")
  return("Unknown")
}

# ------------------------------------------------------------
# Load inputs mirroring fig-5c
# ------------------------------------------------------------

tusco_universal_tsv <- file.path(data_root, "tusco/tusco_mouse.tsv")
tusco_brain_tsv     <- file.path(data_root, "tusco/tusco_mouse_brain.tsv")
tusco_kidney_tsv    <- file.path(data_root, "tusco/tusco_mouse_kidney.tsv")

isoseq_ar_root <- file.path(data_root, "nih/single_sample")

all_sample_dirs <- list.dirs(isoseq_ar_root, full.names = TRUE, recursive = FALSE)
all_sample_dirs <- all_sample_dirs[grepl("/(B|K)[0-9]+\\.isoforms$", all_sample_dirs)]
if (length(all_sample_dirs) == 0) stop("No sample directories found in ", isoseq_ar_root)

tusco_ref_univ   <- load_tusco_reference(tusco_universal_tsv)
tusco_ref_brain  <- load_tusco_reference(tusco_brain_tsv)
tusco_ref_kidney <- load_tusco_reference(tusco_kidney_tsv)

results <- list()
for (sd in all_sample_dirs) {
  sample_name <- basename(sd)
  class_file  <- file.path(sd, paste0(sample_name, "_classification.txt"))
  if (!file.exists(class_file)) next

  classification <- read_tsv_safe(class_file)
  sample_tissue <- infer_tissue(sample_name)

  mu <- compute_tusco_metrics(classification, tusco_ref_univ)
  results[[length(results) + 1]] <- mu %>% mutate(
    Sample = sample_name, Tissue = sample_tissue, Type = "Universal"
  )

  if (sample_tissue == "Brain") {
    mb <- compute_tusco_metrics(classification, tusco_ref_brain)
    results[[length(results) + 1]] <- mb %>% mutate(
      Sample = sample_name, Tissue = sample_tissue, Type = "Brain"
    )
  } else if (sample_tissue == "Kidney") {
    mk <- compute_tusco_metrics(classification, tusco_ref_kidney)
    results[[length(results) + 1]] <- mk %>% mutate(
      Sample = sample_name, Tissue = sample_tissue, Type = "Kidney"
    )
  }
}

metrics_df <- bind_rows(results) %>%
  select(Sample, Tissue, Type, `Sn`, `nrPre`, `1/red`, `1-FDR`, `PDR`, `rPre`)

metrics_summary_df <- metrics_df %>%
  group_by(Tissue, Type) %>%
  summarise(
    `Sn` = mean(`Sn`, na.rm = TRUE),
    `nrPre` = mean(`nrPre`, na.rm = TRUE),
    `1/red` = mean(`1/red`, na.rm = TRUE),
    `1-FDR` = mean(`1-FDR`, na.rm = TRUE),
    `PDR` = mean(`PDR`, na.rm = TRUE),
    `rPre` = mean(`rPre`, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Radar utilities (adapted from figure3a-mouse.R)
# ------------------------------------------------------------

radar_grob <- function(radar_data, color_map, var_labels, title = NULL) {
  # radar_data: data.frame with columns in metric order and rows: Max, Min, then one per Type
  tmpfile <- tempfile(fileext = ".png")
  png(tmpfile, width = 1200, height = 1200, res = 200)
  par(family = "Helvetica", mar = c(0, 0, 0, 0))
  # Derive colors in row order (skip first two rows which are max/min)
  polygon_colors <- unname(color_map[rownames(radar_data)[-(1:2)]])
  radarchart(
    radar_data,
    axistype = 0,
    cglcol = "grey", cglty = 1, cglwd = 1.5,
    plty = rep(1, length(polygon_colors)),
    axislabcol = "black",
    vlabels = var_labels,
    pcol = polygon_colors, plwd = 8, pty = 16,
    caxislabels = NULL, vlcex = 0
  )
  dev.off()
  img <- png::readPNG(tmpfile)
  file.remove(tmpfile)
  grid::rasterGrob(img)
}

build_radar_data_for_tissue <- function(sdf, tissue_name) {
  metrics_order <- c("Sn", "nrPre", "1/red", "1-FDR", "PDR", "rPre")
  types_order <- c("Universal", "Brain", "Kidney")
  base <- sdf %>% filter(Tissue == tissue_name) %>%
    select(Type, all_of(metrics_order))
  present_types <- intersect(types_order, unique(as.character(base$Type)))
  base <- base %>% filter(Type %in% present_types) %>% arrange(factor(Type, levels = present_types))
  mat <- as.matrix(base %>% select(all_of(metrics_order)))
  rownames(mat) <- base$Type
  df <- as.data.frame(rbind(
    rep(100, ncol(mat)),
    rep(0, ncol(mat)),
    mat
  ))
  colnames(df) <- metrics_order
  rownames(df) <- c("Max", "Min", present_types)
  df
}

# Colors consistent with fig-5c bars
type_colors <- c(
  "Universal" = "#1C9E77",
  "Brain"     = "#5893a4",
  "Kidney"    = "#506f68"
)

# Build radar grobs for Brain and Kidney
tissues <- c("Brain", "Kidney")

radar_grobs <- list()
metrics_labels <- c("Sn", "nrPre", "1/red", "1-FDR", "PDR", "rPre")

compute_cosine <- function(vec_a, vec_b) {
  if (any(is.na(vec_a)) || any(is.na(vec_b))) return(NA_real_)
  denom <- sqrt(sum(vec_a^2)) * sqrt(sum(vec_b^2))
  if (denom == 0) return(NA_real_)
  sum(vec_a * vec_b) / denom
}

for (t in tissues) {
  rd <- build_radar_data_for_tissue(metrics_summary_df, t)
  # Restrict color map to displayed types
  cm <- type_colors[rownames(rd)[-(1:2)]]

  # Cosine similarity between Universal and tissue-specific set
  other_type <- if (t == "Brain") "Brain" else "Kidney"
  v_univ <- metrics_summary_df %>% filter(Tissue == t, Type == "Universal") %>% select(all_of(metrics_labels)) %>% as.numeric()
  v_tiss <- metrics_summary_df %>% filter(Tissue == t, Type == other_type) %>% select(all_of(metrics_labels)) %>% as.numeric()
  cos_label <- NA_character_
  if (length(v_univ) == length(metrics_labels) && length(v_tiss) == length(metrics_labels)) {
    cs <- compute_cosine(v_univ, v_tiss)
    if (!is.na(cs)) cos_label <- sprintf("Cosine similarity: %.3f", cs)
  }

  radar_grobs[[t]] <- ggdraw() +
    draw_grob(radar_grob(rd, cm, var_labels = metrics_labels, title = NULL), 0, 0, 1, 1) +
    draw_label(t, x = 0.5, y = 1.02, hjust = 0.5, vjust = 0, size = 12) +
    if (!is.na(cos_label)) draw_label(cos_label, x = 0.5, y = -0.02, hjust = 0.5, vjust = 1, size = 10) else NULL
}

# Legend
legend_df <- data.frame(Type = factor(names(type_colors), levels = names(type_colors)),
                        x = c(1, 2, 3), y = c(1, 2, 3))
legend_plot <- ggplot(legend_df, aes(x = x, y = y, color = Type, linetype = Type, shape = Type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = type_colors) +
  scale_linetype_manual(values = c("Universal" = "solid", "Brain" = "solid", "Kidney" = "solid")) +
  scale_shape_manual(values = c("Universal" = 16, "Brain" = 16, "Kidney" = 16)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.5, "lines"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.direction = "horizontal"
  )
legend_grob <- cowplot::get_legend(legend_plot)
legend_draw <- ggdraw() + draw_grob(legend_grob, 0, 0, 1, 1)

# Arrange two radars side-by-side with legend
panel <- plot_grid(radar_grobs[["Brain"]], radar_grobs[["Kidney"]], ncol = 2, rel_widths = c(1, 1))
panel_with_legend <- plot_grid(panel, legend_draw, ncol = 1, rel_heights = c(1, 0.15))

# Export
pdf(out_pdf, width = 7.09, height = 3.2)
print(panel_with_legend)
dev.off()

message("Wrote radar panel to ", out_pdf)


