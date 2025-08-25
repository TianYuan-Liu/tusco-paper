suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(rtracklayer)
  library(cowplot)
  library(scales)
  library(grid)
  library(png)
  library(fmsb)
})

utils::globalVariables(c(
  "ensembl", "transcript", "refseq", "prot_refseq", "gene_name",
  "structural_category", "associated_gene", "associated_transcript",
  "isoform", "id_type", "subcategory", "ref_exons", "diff_to_TSS",
  "diff_to_TTS", "Origin", "SamplesLabel", "Value", "improvement",
  "sample1_value", "sample2_value", "annotation_y", "improvement_text",
  "Samples", "gene"
))

# ------------------------------------------------------------
# Paths and constants
# ------------------------------------------------------------
project_root <- "/Users/tianyuan/Desktop/github_dev/tusco-paper"
data_root <- file.path(project_root, "figs/data")
intersection_root <- file.path(data_root, "nih/intersection")
plot_dir <- file.path(project_root, "figs/fig-5/plot")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
out_pdf <- file.path(plot_dir, "new-fig-5b-5c.pdf")

# Inputs used by 5b metrics
sirv_file <- "/Users/tianyuan/Desktop/GitHub/TUSCO/reference_dataset/SIRVs.gtf"
tusco_file <- file.path(Sys.getenv("HOME"),
                        "Desktop/github_dev/tusco-paper/figs/data/tusco/tusco_mouse.tsv")

B_pipelines <- c("B31", "B31_B32", "B31_B32_B33", "B31_B32_B33_B34", "B31_B32_B33_B34_B35")
K_pipelines <- c("K31", "K31_K32", "K31_K32_K33", "K31_K32_K33_K34", "K31_K32_K33_K34_K35")
pipelines <- c(B_pipelines, K_pipelines)

read_tsv_safe <- function(file_path, ...) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  tryCatch(read_tsv(file_path, show_col_types = FALSE, ...),
           error = function(e) stop("Unable to read ", file_path, ": ", e$message))
}

get_sample_number <- function(pipeline_name) length(strsplit(pipeline_name, "_")[[1]])
get_tissue_type <- function(pipeline_name) if (startsWith(pipeline_name, "B")) "Brain" else "Kidney"

# ------------------------------------------------------------
# Figure 5b: metrics and plotting
# ------------------------------------------------------------
compute_metrics_from_classification <- function(classification_data) {
  sirv_gtf_df <- as.data.frame(rtracklayer::import(sirv_file))
  sirv_exons <- sirv_gtf_df[sirv_gtf_df$type == "exon", ]
  rSIRV <- length(unique(sirv_exons$transcript_id))
  classification_data_cleaned_sirv <- classification_data %>%
    dplyr::filter(structural_category != "fusion") %>%
    dplyr::bind_rows(classification_data %>%
                dplyr::filter(structural_category == "fusion") %>%
                tidyr::separate_rows(associated_transcript, sep = "_")) %>%
    dplyr::mutate(
      associated_gene = stringr::str_remove(associated_gene, "\\.\\d+$"),
      associated_transcript = stringr::str_remove(associated_transcript, "\\.\\d+$")
    ) %>%
    dplyr::distinct(isoform, associated_transcript, .keep_all = TRUE)
  sirv_chromosomes <- unique(sirv_gtf_df$seqnames)
  classification_data_cleaned_sirv <- classification_data_cleaned_sirv %>%
    dplyr::filter(chrom %in% sirv_chromosomes)
  SIRV_transcripts <- classification_data_cleaned_sirv %>% dplyr::filter(grepl("SIRV", chrom))
  SIRV_RM <- SIRV_transcripts %>%
    dplyr::filter(
      subcategory == "reference_match" |
        (subcategory == "mono-exon" & ref_exons == 1 & abs(diff_to_TSS) < 50 & abs(diff_to_TTS) < 50)
    )
  TP_sirv <- SIRV_RM
  PTP_sirv <- SIRV_transcripts %>%
    dplyr::filter(structural_category %in% c("full-splice_match", "incomplete-splice_match") &
             !associated_transcript %in% TP_sirv$associated_transcript)
  fsm_ism_count_sirv <- SIRV_transcripts %>% dplyr::filter(structural_category %in% c("full-splice_match", "incomplete-splice_match")) %>% nrow()
  non_redundant_sensitivity_sirv <- length(unique(TP_sirv$associated_transcript)) / rSIRV
  non_redundant_precision_sirv  <- if (nrow(SIRV_transcripts) > 0) nrow(TP_sirv) / nrow(SIRV_transcripts) else NA
  redundancy_sirv              <- if (length(unique(c(TP_sirv$associated_transcript, PTP_sirv$associated_transcript))) > 0)
                                     fsm_ism_count_sirv / length(unique(c(TP_sirv$associated_transcript, PTP_sirv$associated_transcript))) else NA
  redundant_precision_sirv     <- if (nrow(SIRV_transcripts) > 0) (nrow(TP_sirv) + nrow(PTP_sirv)) / nrow(SIRV_transcripts) else NA
  false_discovery_rate_sirv    <- if (nrow(SIRV_transcripts) > 0) (nrow(SIRV_transcripts) - nrow(SIRV_RM)) / nrow(SIRV_transcripts) else NA
  positive_detection_rate_sirv <- if (rSIRV > 0) length(unique(c(TP_sirv$associated_transcript, PTP_sirv$associated_transcript))) / rSIRV else NA

  tusco_lines <- readLines(tusco_file)
  tusco_data_lines <- tusco_lines[!grepl("^#", tusco_lines)]
  tusco_temp_file <- tempfile()
  writeLines(tusco_data_lines, tusco_temp_file)
  tusco_df <- read_delim(tusco_temp_file, delim = "\t", col_names = c("ensembl","transcript","gene_name","gene_id_num","refseq","prot_refseq"),
                         col_types = cols(.default = "c"), trim_ws = TRUE, show_col_types = FALSE)
  unlink(tusco_temp_file)
  tusco_df_clean <- tusco_df %>%
    dplyr::mutate(
      ensembl = stringr::str_remove(ensembl, "\\.\\d+$"),
      transcript = stringr::str_remove(transcript, "\\.\\d+$"),
      refseq = stringr::str_remove(refseq, "\\.\\d+$"),
      prot_refseq = stringr::str_remove(prot_refseq, "\\.\\d+$")
    )
  annotation_data_tusco <- tusco_df_clean %>% dplyr::select(ensembl, refseq, gene_name) %>% dplyr::distinct()
  rTUSCO <- nrow(annotation_data_tusco)

  patterns <- list(ensembl="^(ENSG|ENSMUSG)\\d{11}(\\.\\d+)?$", refseq="^(NM_|NR_|NP_)\\d{6,}$", gene_name="^[A-Za-z0-9][A-Za-z0-9_-]*[A-Za-z0-9]?$")

  classification_data_cleaned_tusco <- classification_data %>%
    dplyr::filter(structural_category != "fusion") %>%
    dplyr::bind_rows(classification_data %>% dplyr::filter(structural_category == "fusion") %>% tidyr::separate_rows(associated_gene, sep="_")) %>%
    dplyr::mutate(
      associated_gene = stringr::str_remove(associated_gene, "\\.\\d+$"),
      associated_transcript = stringr::str_remove(associated_transcript, "\\.\\d+$")
    ) %>%
    dplyr::mutate(id_type = dplyr::case_when(
             stringr::str_detect(associated_gene, patterns$ensembl)   ~ "ensembl",
             stringr::str_detect(associated_gene, patterns$refseq)    ~ "refseq",
             stringr::str_detect(associated_gene, patterns$gene_name) ~ "gene_name",
             TRUE                                            ~ "unknown")) %>%
    dplyr::distinct(isoform, associated_gene, .keep_all = TRUE)

  top_id_type_tusco <- classification_data_cleaned_tusco %>%
    dplyr::count(id_type, sort = TRUE) %>% dplyr::slice_head(n = 1) %>% dplyr::pull(id_type)

  TUSCO_transcripts <- classification_data_cleaned_tusco %>%
    dplyr::filter(associated_gene %in% annotation_data_tusco[[top_id_type_tusco]])

  TUSCO_RM <- TUSCO_transcripts %>%
    dplyr::filter(
      subcategory == "reference_match" |
        (subcategory == "mono-exon" & ref_exons == 1 & abs(diff_to_TSS) < 50 & abs(diff_to_TTS) < 50)
    )
  TP_tusco <- TUSCO_RM
  PTP_tusco <- TUSCO_transcripts %>%
    dplyr::filter(structural_category %in% c("full-splice_match","incomplete-splice_match") &
             !associated_gene %in% TP_tusco$associated_gene)

  fsm_ism_count_tusco <- TUSCO_transcripts %>%
    dplyr::filter(structural_category %in% c("full-splice_match","incomplete-splice_match")) %>% nrow()

  non_redundant_sensitivity_tusco <- length(unique(TP_tusco$associated_gene)) / rTUSCO
  non_redundant_precision_tusco  <- if (nrow(TUSCO_transcripts) > 0) nrow(TP_tusco) / nrow(TUSCO_transcripts) else NA
  redundancy_tusco               <- if (length(unique(c(TP_tusco$associated_gene, PTP_tusco$associated_gene))) > 0)
                                       fsm_ism_count_tusco / length(unique(c(TP_tusco$associated_gene, PTP_tusco$associated_gene))) else NA
  false_discovery_rate_tusco     <- if (nrow(TUSCO_transcripts) > 0) (nrow(TUSCO_transcripts) - nrow(TUSCO_RM)) / nrow(TUSCO_transcripts) else NA
  positive_detection_rate_tusco  <- if (rTUSCO > 0) length(unique(c(TP_tusco$associated_gene, PTP_tusco$associated_gene))) / rTUSCO else NA
  f1_tusco <- if (!is.na(non_redundant_precision_tusco) && !is.na(non_redundant_sensitivity_tusco) &&
                   (non_redundant_precision_tusco + non_redundant_sensitivity_tusco) > 0)
                (2 * non_redundant_precision_tusco * non_redundant_sensitivity_tusco) /
                (non_redundant_precision_tusco + non_redundant_sensitivity_tusco) else NA

  list(
    metrics_sirv = c(
      Sensitivity = non_redundant_sensitivity_sirv * 100,
      `Non-redundant Precision` = non_redundant_precision_sirv * 100,
      `Inv. Redundancy` = if (!is.na(redundancy_sirv) && redundancy_sirv != 0) (1/redundancy_sirv) * 100 else 0,
      `1 - FDR` = if (!is.na(false_discovery_rate_sirv)) (100 - (false_discovery_rate_sirv * 100)) else 0,
      `PDR` = if (!is.na(positive_detection_rate_sirv)) positive_detection_rate_sirv * 100 else 0,
      `Redundant Precision` = if (!is.na(redundant_precision_sirv)) redundant_precision_sirv * 100 else 0
    ),
    metrics_tusco = c(
      Sensitivity = non_redundant_sensitivity_tusco * 100,
      `Non-redundant Precision` = non_redundant_precision_tusco * 100,
      `Inv. Redundancy` = if (!is.na(redundancy_tusco) && redundancy_tusco != 0) (1/redundancy_tusco) * 100 else 0,
      `1 - FDR` = if (!is.na(false_discovery_rate_tusco)) (100 - (false_discovery_rate_tusco * 100)) else 0,
      `PDR` = if (!is.na(positive_detection_rate_tusco)) positive_detection_rate_tusco * 100 else 0,
      `F1` = if (!is.na(f1_tusco)) f1_tusco * 100 else 0
    )
  )
}

load_classification_for <- function(root_dir, pipeline_prefix, mode) {
  pipeline_dir <- file.path(root_dir, pipeline_prefix)
  if (!dir.exists(pipeline_dir)) stop("Missing pipeline dir: ", pipeline_dir)
  files <- list.files(pipeline_dir, full.names = TRUE)
  if (tolower(mode) == "union") {
    f <- files[grepl("_union_classification\\.txt$", files)]
  } else {
    f <- files[grepl("(?<!_union)_classification\\.txt$", files, perl = TRUE)]
  }
  if (length(f) == 0) stop("No classification file for ", pipeline_prefix, " in ", mode)
  read_tsv_safe(f[1])
}

process_pipeline_for_barplot <- function(root_dir, pipeline_prefix) {
  classification_data <- load_classification_for(root_dir, pipeline_prefix, mode = "intersection")
  metrics <- compute_metrics_from_classification(classification_data)
  sample_num <- get_sample_number(pipeline_prefix)
  tissue <- get_tissue_type(pipeline_prefix)
  tusco_df <- data.frame(
    Pipeline = pipeline_prefix,
    Samples = sample_num,
    Tissue = tissue,
    Dataset = "TUSCO",
    Origin = "Intersection",
    Metric = names(metrics$metrics_tusco),
    Value = as.numeric(metrics$metrics_tusco),
    stringsAsFactors = FALSE
  )
  tusco_df
}

all_metrics_list <- lapply(pipelines, function(p) process_pipeline_for_barplot(intersection_root, p))
all_metrics_df_all <- do.call(rbind, all_metrics_list)

intersection_best_df <- all_metrics_df_all %>%
  dplyr::mutate(Origin = "Intersection") %>%
  dplyr::group_by(Dataset, Tissue, Metric, Samples) %>%
  dplyr::slice_max(order_by = Value, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

# Only keep Brain and Kidney (remove Mixed bars entirely)
all_metrics_df <- intersection_best_df %>% dplyr::filter(Tissue %in% c("Brain", "Kidney"))

all_metrics_df$Tissue <- factor(all_metrics_df$Tissue, levels = c("Brain", "Kidney"))
all_metrics_df$Dataset <- factor(all_metrics_df$Dataset, levels = c("TUSCO"))
metric_order <- c("Sensitivity", "Non-redundant Precision", "Inv. Redundancy", "1 - FDR", "PDR", "F1")
all_metrics_df$Metric <- factor(all_metrics_df$Metric, levels = metric_order)

# ------------------------------------------------------------
# Export bar-point matrix (per-bar values) for Figure 5b
# ------------------------------------------------------------
points_df <- all_metrics_df %>%
  dplyr::mutate(
    Combo = Pipeline,
    SamplesLabel = dplyr::case_when(TRUE ~ as.character(Samples))
  ) %>%
  dplyr::select(Tissue, Samples, Combo, Metric, Value, SamplesLabel) %>%
  dplyr::arrange(Tissue, Samples, Metric)

points_tsv <- file.path(plot_dir, "new-fig-5b_points.tsv")
readr::write_tsv(points_df, points_tsv)
message("Wrote bar-point matrix to ", points_tsv)

create_metric_barplot_single <- function(metric_name, data_df, color_value) {
  metric_data <- data_df %>% dplyr::filter(Metric == metric_name)
  metric_data <- metric_data %>% dplyr::mutate(
    SamplesLabel = dplyr::case_when(
      TRUE ~ as.character(Samples)
    )
  )
  metric_data$SamplesLabel <- factor(metric_data$SamplesLabel,
                                     levels = c("1","2","3","4","5","2x10M","pad","3x7M"))

  improvement_base <- metric_data %>% dplyr::filter(Origin == "Intersection")
  improvement_data <- improvement_base %>%
    dplyr::filter(Samples %in% c(1, 2), Tissue %in% c("Brain", "Kidney")) %>%
    dplyr::group_by(Tissue) %>% dplyr::arrange(Samples) %>% dplyr::filter(n() == 2) %>%
    dplyr::summarise(
      improvement = Value[Samples == 2] - Value[Samples == 1],
      sample1_value = Value[Samples == 1], sample2_value = Value[Samples == 2], .groups = "drop"
    ) %>%
    dplyr::mutate(improvement_text = ifelse(improvement >= 0, paste0("+", round(improvement, 1), "pp"), paste0(round(improvement, 1), "pp")),
                  annotation_y = pmax(sample1_value, sample2_value) + 12)

  position_d <- position_dodge(width = 0.75)
  ggplot(metric_data, aes(x = SamplesLabel, y = Value, fill = Origin)) +
    geom_bar(data = metric_data %>% dplyr::filter(Tissue %in% c("Brain", "Kidney")),
             stat = "identity", position = position_d, width = 0.7,
             color = "white", linewidth = 0.3) +
    geom_bar(data = metric_data %>% dplyr::filter(Origin == "Intersection", Samples %in% c(1, 2), Tissue %in% c("Brain", "Kidney")),
             aes(x = SamplesLabel, y = Value),
             stat = "identity", position = position_d, fill = "transparent", width = 0.7,
             color = "#d62728", linewidth = 1.0) +
    geom_text(data = improvement_data,
              aes(x = 1.5, y = annotation_y, label = improvement_text),
              inherit.aes = FALSE, hjust = 0.5, vjust = 0.5, size = 7/.pt, fontface = "bold", color = "#d62728") +
    facet_grid(. ~ Tissue, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c(Intersection = color_value, Merged_10M = alpha(color_value, 0.6), Merged_7M = alpha("#1f77b4", 0.7)),
                      breaks = c("Intersection", "Merged_10M", "Merged_7M"), name = "Origin") +
    scale_y_continuous(limits = c(0, 130), labels = scales::percent_format(scale = 1), breaks = seq(0, 100, 20)) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)),
                     labels = function(lbl) { out <- ifelse(lbl == "pad", "", lbl); out[out == "2x10M"] <- "2x\n10M"; out[out == "3x7M"] <- "3x\n7M"; out }) +
    coord_cartesian(clip = "off") +
    labs(title = metric_name, x = NULL, y = "Performance (%)") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 7, hjust = 0.5, face = "bold", margin = margin(b = 6)),
      axis.title.x = element_blank(), axis.title.y = element_text(size = 7, margin = margin(r = 4)),
      axis.text = element_text(size = 7), axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5, vjust = 0.8),
      strip.text = element_text(size = 7, face = "bold", margin = margin(b = 4)),
      panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.spacing = unit(0.5, "cm"), plot.margin = margin(2, 2, 2, 2), legend.position = "none"
    )
}

metric_levels <- metric_order
tusco_plots <- lapply(metric_levels, function(metric) create_metric_barplot_single(metric, all_metrics_df, "#2ca02c"))
fig_5b_panel <- plot_grid(plotlist = tusco_plots, ncol = 3, nrow = 2, align = "hv")

# ------------------------------------------------------------
# Figure 5c: radar metrics and plotting
# ------------------------------------------------------------
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

radar_grob <- function(radar_data, color_map, var_labels, title = NULL) {
  tmpfile <- tempfile(fileext = ".png")
  png(tmpfile, width = 2000, height = 2000, res = 600)
  par(family = "Helvetica", ps = 7, cex = 1, mar = c(0, 0, 0, 0))
  present_types <- rownames(radar_data)[-(1:2)]
  polygon_colors <- unname(color_map[present_types])
  linetype_map <- c("Universal" = 1, "Brain" = 2, "Kidney" = 3)
  polygon_lty <- unname(linetype_map[present_types])
  radarchart(
    radar_data,
    axistype = 0,
    cglcol = "grey85", cglty = 1, cglwd = 0.4,
    plty = polygon_lty,
    axislabcol = "black",
    vlabels = var_labels,
    pcol = polygon_colors, plwd = 1.6, pty = 16,
    caxislabels = NULL, vlcex = 1.8
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

type_colors <- c(
  "Universal" = "#1C9E77",
  "Brain"     = "#5893a4",
  "Kidney"    = "#506f68"
)

# Load inputs for radar metrics
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

metrics_labels_internal <- c("Sn", "nrPre", "1/red", "1-FDR", "PDR", "rPre")
metrics_labels_short    <- c("Sn", "nrPre", "1/red", "1 - FDR", "PDR", "rPre")

compute_cosine <- function(vec_a, vec_b) {
  if (any(is.na(vec_a)) || any(is.na(vec_b))) return(NA_real_)
  denom <- sqrt(sum(vec_a^2)) * sqrt(sum(vec_b^2))
  if (denom == 0) return(NA_real_)
  sum(vec_a * vec_b) / denom
}

build_radar_plot_for_tissue <- function(summary_df, tissue_label) {
  rd <- build_radar_data_for_tissue(summary_df, tissue_label)
  cm <- type_colors[rownames(rd)[-(1:2)]]
  # cosine similarity between Universal and tissue-specific vector
  metrics_order <- c("Sn", "nrPre", "1/red", "1-FDR", "PDR", "rPre")
  v_univ <- summary_df %>% filter(Tissue == tissue_label, Type == "Universal") %>%
    select(all_of(metrics_order)) %>% as.numeric()
  v_tiss <- summary_df %>% filter(Tissue == tissue_label, Type == tissue_label) %>%
    select(all_of(metrics_order)) %>% as.numeric()
  cos_label <- NA_character_
  if (length(v_univ) == length(metrics_order) && length(v_tiss) == length(metrics_order)) {
    cs <- compute_cosine(v_univ, v_tiss)
    if (!is.na(cs)) cos_label <- sprintf("Cosine similarity: %.3f", cs)
  }
  ggdraw() +
    # Slightly overfill horizontally to crop internal whitespace from the raster
    draw_grob(radar_grob(rd, cm, var_labels = metrics_labels_short, title = NULL), x = -0.05, y = 0.10, width = 1.10, height = 0.86) +
    draw_label(tissue_label, x = 0.5, y = 1.0, hjust = 0.5, vjust = 0, size = 7)
}

radar_brain <- build_radar_plot_for_tissue(metrics_summary_df, "Brain")
radar_kidney <- build_radar_plot_for_tissue(metrics_summary_df, "Kidney")

# Inline legend (overlay on bottom radar to keep radar area to exactly two rows)
legend_df <- data.frame(Type = factor(names(type_colors), levels = names(type_colors)),
                        x = c(1, 2, 3), y = c(1, 2, 3))
# Custom legend labels
legend_labels <- c(
  "Universal" = "TUSCO (Mouse)",
  "Brain"     = "TUSCO (Mouse Brain)",
  "Kidney"    = "TUSCO (Mouse Kidney)"
)
legend_plot <- ggplot(legend_df, aes(x = x, y = y, color = Type, linetype = Type, shape = Type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = type_colors, labels = legend_labels) +
  scale_linetype_manual(values = c("Universal" = "solid", "Brain" = "dashed", "Kidney" = "dotted"), labels = legend_labels) +
  scale_shape_manual(values = c("Universal" = 16, "Brain" = 16, "Kidney" = 16), labels = legend_labels) +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 7))
legend_grob <- cowplot::get_legend(legend_plot)

wrap_with_legend_space <- function(radar_plot, legend = NULL) {
  if (is.null(legend)) legend <- grid::nullGrob()
  ggdraw() +
    # Slight horizontal overfill to crop left/right whitespace and reserve larger bottom space for legend
    draw_plot(radar_plot, -0.05, 0.14, 1.10, 0.82) +
    draw_grob(legend, x = 0.5, y = 0.03, width = 0.95, height = 0.10)
}

# Build a dedicated legend row to guarantee visibility
legend_only <- ggdraw() + draw_grob(legend_grob, x = 0.5, y = 0.5, width = 0.95, height = 1)

# Stack: radar brain, radar kidney, legend
fig_5c_panel <- plot_grid(
  radar_brain,
  radar_kidney,
  legend_only,
  ncol = 1,
  rel_heights = c(1, 1, 0.18)
)
fig_5c_panel <- fig_5c_panel + theme(plot.margin = margin(0, 0, 0, 0))

# ------------------------------------------------------------
# Combine: side-by-side in one row (barplot wider), same width/height as 5b
# ------------------------------------------------------------
combined <- plot_grid(
  fig_5b_panel, fig_5c_panel,
  ncol = 2,
  # Give more space to barplots and reduce radar panel width
  rel_widths = c(3.0, 1.0),
  align = "hv",
  axis = "tb",
  greedy = FALSE,
  labels = c("c", "d"),
  label_size = 10,
  label_fontface = "bold"
)

pdf(out_pdf, width = 7.09, height = 3.55)
print(combined)
dev.off()

message("Wrote combined figure to ", out_pdf)