#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(readr)
  library(stringr)
  library(rtracklayer)
  library(cowplot)
  library(scales)
})
utils::globalVariables(c(
  "ensembl", "transcript", "refseq", "prot_refseq", "gene_name",
  "structural_category", "associated_gene", "associated_transcript",
  "isoform", "id_type", "subcategory", "ref_exons", "diff_to_TSS",
  "diff_to_TTS", "Origin", "SamplesLabel", "Value", "improvement",
  "sample1_value", "sample2_value", "annotation_y", "improvement_text",
  "Samples"
))

# Minimal script to create Figure 5b (TUSCO main with merged bars)
# Writes: figs/fig-5/plot/fig-5b.pdf

project_root <- "/Users/tianyuan/Desktop/github_dev/tusco-paper"
data_root <- file.path(project_root, "figs/data")
intersection_root <- file.path(data_root, "nih/intersection")
merge10_root <- "/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/data/nih/10M_2reps_kidney/K31_K32K34_K35/sqanti3_out"
merge7_root  <- "/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/data/nih/7M_3reps_mixed/sqanti3_out"
output_dir <- file.path(project_root, "figs/fig-5/plot")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

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
  redundant_precision_tusco      <- if (nrow(TUSCO_transcripts) > 0) (nrow(TP_tusco) + nrow(PTP_tusco)) / nrow(TUSCO_transcripts) else NA
  false_discovery_rate_tusco     <- if (nrow(TUSCO_transcripts) > 0) (nrow(TUSCO_transcripts) - nrow(TUSCO_RM)) / nrow(TUSCO_transcripts) else NA
  positive_detection_rate_tusco  <- if (rTUSCO > 0) length(unique(c(TP_tusco$associated_gene, PTP_tusco$associated_gene))) / rTUSCO else NA

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
      `Redundant Precision` = if (!is.na(redundant_precision_tusco)) redundant_precision_tusco * 100 else 0
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

merge10_class <- file.path(merge10_root, "sqanti3_classification.txt")
merge7_class  <- file.path(merge7_root,  "sqanti3_classification.txt")
if (!file.exists(merge10_class)) stop("Missing: ", merge10_class)
if (!file.exists(merge7_class))  stop("Missing: ", merge7_class)

merge10_metrics <- compute_metrics_from_classification(read_tsv_safe(merge10_class))
merge7_metrics  <- compute_metrics_from_classification(read_tsv_safe(merge7_class))

merge10_tusco_df <- data.frame(
  Pipeline = "K31_K32K34_K35",
  Samples = 2,
  Tissue = "Mixed",
  Dataset = "TUSCO",
  Origin = "Merged_10M",
  Metric = names(merge10_metrics$metrics_tusco),
  Value = as.numeric(merge10_metrics$metrics_tusco),
  stringsAsFactors = FALSE
)
merge7_tusco_df <- data.frame(
  Pipeline = "7M_3reps_mixed",
  Samples = 3,
  Tissue = "Mixed",
  Dataset = "TUSCO",
  Origin = "Merged_7M",
  Metric = names(merge7_metrics$metrics_tusco),
  Value = as.numeric(merge7_metrics$metrics_tusco),
  stringsAsFactors = FALSE
)

intersection_best_df <- all_metrics_df_all %>%
  dplyr::mutate(Origin = "Intersection") %>%
  dplyr::group_by(Dataset, Tissue, Metric, Samples) %>%
  dplyr::slice_max(order_by = Value, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

all_metrics_df <- dplyr::bind_rows(intersection_best_df, merge10_tusco_df, merge7_tusco_df)

all_metrics_df$Tissue <- factor(all_metrics_df$Tissue, levels = c("Brain", "Kidney", "Mixed"))
all_metrics_df$Dataset <- factor(all_metrics_df$Dataset, levels = c("TUSCO"))
metric_order <- c("Sensitivity", "Non-redundant Precision", "Inv. Redundancy", "1 - FDR", "PDR", "Redundant Precision")
all_metrics_df$Metric <- factor(all_metrics_df$Metric, levels = metric_order)

create_metric_barplot_single <- function(metric_name, data_df, color_value) {
  metric_data <- data_df %>% dplyr::filter(Metric == metric_name)
  metric_data <- metric_data %>% dplyr::mutate(
    SamplesLabel = dplyr::case_when(
      Tissue == "Mixed" & Origin == "Merged_10M" & Samples == 2 ~ "2x10M",
      Tissue == "Mixed" & Origin == "Merged_7M"  & Samples == 3 ~ "3x7M",
      TRUE ~ as.character(Samples)
    )
  )
  metric_data$SamplesLabel <- factor(metric_data$SamplesLabel,
                                     levels = c("1","2","3","4","5","2x10M","pad","3x7M"))
  if (any(metric_data$Tissue == "Mixed")) {
    origins_mixed <- unique(metric_data$Origin[metric_data$Tissue == "Mixed"])
    pad_rows <- data.frame(
      Pipeline = NA_character_, Samples = NA_integer_, Tissue = "Mixed",
      Dataset = "TUSCO", Origin = origins_mixed, Metric = metric_name, Value = NA_real_,
      SamplesLabel = factor("pad", levels = c("1","2","3","4","5","2x10M","pad","3x7M")),
      stringsAsFactors = FALSE
    )
    metric_data$SamplesLabel <- factor(metric_data$SamplesLabel,
                                       levels = c("1","2","3","4","5","2x10M","pad","3x7M"))
    metric_data <- dplyr::bind_rows(metric_data, pad_rows)
  }
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
    geom_bar(data = metric_data %>% dplyr::filter(Tissue != "Mixed"),
             stat = "identity", position = position_d, width = 0.7,
             color = "white", linewidth = 0.3) +
    geom_bar(data = metric_data %>% dplyr::filter(Tissue == "Mixed", Origin == "Intersection"),
             stat = "identity", position = position_d, width = 0.7,
             color = "white", linewidth = 0.3) +
    geom_bar(data = metric_data %>% dplyr::filter(Tissue == "Mixed", Origin != "Intersection"),
             stat = "identity", position = position_d, width = 0.7,
             color = "white", linewidth = 0.3) +
    geom_bar(data = metric_data %>% dplyr::filter(Origin == "Intersection", Samples %in% c(1, 2), Tissue != "Mixed"),
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
      panel.spacing = unit(0.5, "cm"), plot.margin = margin(4, 4, 4, 4), legend.position = "none"
    )
}

metric_levels <- metric_order
tusco_plots <- lapply(metric_levels, function(metric) create_metric_barplot_single(metric, all_metrics_df, "#2ca02c"))
fig_tusco <- plot_grid(plotlist = tusco_plots, ncol = 3, nrow = 2, align = "hv")

pdf(file.path(output_dir, "fig-5b.pdf"), width = 7.09, height = 3.55)
print(fig_tusco)
dev.off()

message("Saved fig-5b.pdf to ", output_dir)


