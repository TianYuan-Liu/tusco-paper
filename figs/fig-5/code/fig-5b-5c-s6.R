# Determine LOCAL_ONLY early to conditionally load heavy packages
local_only <- identical(Sys.getenv("LOCAL_ONLY"), "1")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(cowplot)
  library(scales)
  library(grid)
  library(png)
  library(fmsb)
})

# rtracklayer is only required when computing from raw classifications
if (!local_only) {
  suppressPackageStartupMessages(library(rtracklayer))
}

utils::globalVariables(c(
  "ensembl", "transcript", "refseq", "prot_refseq", "gene_name",
  "structural_category", "associated_gene", "associated_transcript",
  "isoform", "id_type", "subcategory", "ref_exons", "diff_to_TSS",
  "diff_to_TTS", "Origin", "SamplesLabel", "Value", "improvement",
  "sample1_value", "sample2_value", "sample3_value", "annotation_y", "improvement_text", "bracket_y",
  "Samples", "gene", "chrom", "type", "n", "Mean", "SD", "Type", "Metric", "Tissue",
  "Sample", "Combo", "Sn", "nrPre", "1/red", "1-FDR", "PDR", "rPre",
  "Max", "Min", ".pt", "N", "SE", "CI_t", "CI_lower", "CI_upper", "x_start", "x_end", "x_mid"
))

# ------------------------------------------------------------
# Paths and constants (support local figure dir and repo figs/data)
# ------------------------------------------------------------
# Figure directory is either cwd if inside figs/fig-5, or figs/fig-5 from repo root
detect_fig_dir <- function() {
  cwd <- normalizePath(getwd(), mustWork = FALSE)
  if (basename(cwd) == "fig-5" && dir.exists(file.path(cwd, "code"))) return(cwd)
  if (dir.exists(file.path(cwd, "figs", "fig-5"))) return(normalizePath(file.path(cwd, "figs", "fig-5"), mustWork = FALSE))
  cwd
}
fig_dir <- detect_fig_dir()
plot_dir <- file.path(fig_dir, "plot")
tsv_dir  <- file.path(fig_dir, "tsv")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(tsv_dir)) dir.create(tsv_dir, recursive = TRUE)
out_pdf <- file.path(plot_dir, "fig-5b-5c.pdf")

# Repository root: two levels above figs/fig-5
repo_root <- normalizePath(file.path(fig_dir, "..", ".."), mustWork = FALSE)
# Data root under repo
data_root <- file.path(repo_root, "figs", "data")
intersection_root <- file.path(data_root, "nih", "intersection")
intersection_all_comb_root <- file.path(data_root, "nih", "intersection_all_comb")
isoseq_ar_root <- file.path(data_root, "nih", "single_sample")

# Local-only flag already set at top

# Inputs used by metrics
# Input reference files
detect_sirv_gtf <- function(data_root_path) {
  candidates <- c(
    file.path(data_root_path, "nih", "SIRVs.gtf"),
    file.path(data_root_path, "spike-ins", "lrgasp_sirvs.gtf"),
    file.path(data_root_path, "spike-ins", "lrgasp_sirvs4.gtf"),
    file.path(data_root_path, "reference", "SIRVs.gtf")
  )
  hits <- candidates[file.exists(candidates)]
  if (length(hits) > 0) return(hits[1])
  NA_character_
}
sirv_file <- if (!local_only) detect_sirv_gtf(data_root) else NA_character_
tusco_file <- if (!local_only) file.path(data_root, "tusco", "tusco_mouse.tsv") else NA_character_

read_tsv_safe <- function(file_path, ...) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  tryCatch(read_tsv(file_path, show_col_types = FALSE, ...),
           error = function(e) stop("Unable to read ", file_path, ": ", e$message))
}

get_sample_number <- function(pipeline_name) length(strsplit(pipeline_name, "_")[[1]])
get_tissue_type <- function(pipeline_name) if (startsWith(pipeline_name, "B")) "Brain" else "Kidney"

# ------------------------------------------------------------
# Metric computation (reused from new-fig-5b-5c.R)
# ------------------------------------------------------------
compute_metrics_from_classification <- function(classification_data) {
  # SIRV metrics (optional)
  if (!is.null(sirv_file) && !is.na(sirv_file) && file.exists(sirv_file)) {
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
      dplyr::mutate(mono_thresh = ifelse(!is.na(ref_length) & ref_length > 3000, 100, 50)) %>%
      dplyr::filter(
        subcategory == "reference_match" |
        (!is.na(ref_length) & ref_length > 3000 & !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
           abs(diff_to_TSS) <= 100 & abs(diff_to_TTS) <= 100) |
        (subcategory == "mono-exon" & ref_exons == 1 & !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
           abs(diff_to_TSS) <= mono_thresh & abs(diff_to_TTS) <= mono_thresh)
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
    f1_sirv <- if (!is.na(non_redundant_sensitivity_sirv) && !is.na(non_redundant_precision_sirv) &&
                    (non_redundant_sensitivity_sirv + non_redundant_precision_sirv) > 0) {
      2 * non_redundant_sensitivity_sirv * non_redundant_precision_sirv /
        (non_redundant_sensitivity_sirv + non_redundant_precision_sirv)
    } else { NA_real_ }
  } else {
    non_redundant_sensitivity_sirv <- NA_real_
    non_redundant_precision_sirv <- NA_real_
    redundancy_sirv <- NA_real_
    false_discovery_rate_sirv <- NA_real_
    positive_detection_rate_sirv <- NA_real_
    fsm_ism_count_sirv <- NA_real_
    TP_sirv <- NULL; PTP_sirv <- NULL
    f1_sirv <- NA_real_
  }

  if (is.null(tusco_file) || is.na(tusco_file) || !file.exists(tusco_file)) {
    stop("Required TUSCO reference not found: ", tusco_file)
  }
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
      (structural_category == "full-splice_match" & !is.na(ref_exons) & ref_exons == 1 &
         !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
         abs(diff_to_TSS) <= 50 & abs(diff_to_TTS) <= 50)
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

  # Compute F1 using non-redundant Sensitivity and non-redundant Precision
  f1_tusco <- if (!is.na(non_redundant_sensitivity_tusco) && !is.na(non_redundant_precision_tusco) &&
                    (non_redundant_sensitivity_tusco + non_redundant_precision_tusco) > 0) {
    2 * non_redundant_sensitivity_tusco * non_redundant_precision_tusco /
      (non_redundant_sensitivity_tusco + non_redundant_precision_tusco)
  } else { NA_real_ }

  list(
    metrics_sirv = c(
      Sensitivity = non_redundant_sensitivity_sirv * 100,
      `Non-redundant Precision` = non_redundant_precision_sirv * 100,
      `Inv. Redundancy` = if (!is.na(redundancy_sirv) && redundancy_sirv != 0) (1/redundancy_sirv) * 100 else 0,
      `1 - FDR` = if (!is.na(false_discovery_rate_sirv)) (100 - (false_discovery_rate_sirv * 100)) else 0,
      `PDR` = if (!is.na(positive_detection_rate_sirv)) positive_detection_rate_sirv * 100 else 0,
      `F1` = if (!is.na(f1_sirv)) f1_sirv * 100 else NA_real_
    ),
    metrics_tusco = c(
      Sensitivity = non_redundant_sensitivity_tusco * 100,
      `Non-redundant Precision` = non_redundant_precision_tusco * 100,
      `Inv. Redundancy` = if (!is.na(redundancy_tusco) && redundancy_tusco != 0) (1/redundancy_tusco) * 100 else 0,
      `1 - FDR` = if (!is.na(false_discovery_rate_tusco)) (100 - (false_discovery_rate_tusco * 100)) else 0,
      `PDR` = if (!is.na(positive_detection_rate_tusco)) positive_detection_rate_tusco * 100 else 0,
      `F1` = if (!is.na(f1_tusco)) f1_tusco * 100 else NA_real_
    )
  )
}

# ------------------------------------------------------------
# Helpers to find classification files
# ------------------------------------------------------------
find_classification_file <- function(dir_path) {
  files <- list.files(dir_path, full.names = TRUE)
  # Prefer non-union intersection classification
  f <- files[grepl("(?<!_union)_classification\\.txt$", files, perl = TRUE)]
  if (length(f) > 0) return(f[1])
  # Fallback to union
  f2 <- files[grepl("_union_classification\\.txt$", files)]
  if (length(f2) > 0) return(f2[1])
  # Last resort: any classification
  f3 <- files[grepl("classification\\.txt$", files)]
  if (length(f3) > 0) return(f3[1])
  stop("No classification file in ", dir_path)
}

# ------------------------------------------------------------
# Shared helpers: collect points and summarise for bar plots
# ------------------------------------------------------------
collect_points_for_metrics <- function(metrics_key) {
  # Single-sample points
  single_sample_dirs <- list.dirs(isoseq_ar_root, full.names = TRUE, recursive = FALSE)
  single_sample_dirs <- single_sample_dirs[grepl("/(B|K)[0-9]+\\.isoforms$", single_sample_dirs)]

  single_points <- list()
  for (sd in single_sample_dirs) {
    sample_name <- basename(sd)
    class_file  <- file.path(sd, paste0(sample_name, "_classification.txt"))
    if (!file.exists(class_file)) next
    tissue <- if (startsWith(sample_name, "B")) "Brain" else if (startsWith(sample_name, "K")) "Kidney" else "Unknown"
    classification <- read_tsv_safe(class_file)
    m <- compute_metrics_from_classification(classification)
    mv <- m[[metrics_key]]
    df <- tibble(
      Tissue = tissue,
      Samples = 1L,
      Combo = sample_name,
      Metric = names(mv),
      Value = as.numeric(mv)
    )
    single_points[[length(single_points) + 1]] <- df
  }
  single_points_df <- dplyr::bind_rows(single_points)

  # Multi-sample combinations (pure Brain or pure Kidney)
  combo_dirs <- list.dirs(intersection_all_comb_root, full.names = TRUE, recursive = FALSE)
  combo_dirs <- combo_dirs[grepl("/(B|K)[0-9]+(_(B|K)[0-9]+)*$", combo_dirs)]
  is_pure_tissue <- function(name) {
    parts <- strsplit(name, "_")[[1]]
    starts_with_b <- all(startsWith(parts, "B"))
    starts_with_k <- all(startsWith(parts, "K"))
    starts_with_b || starts_with_k
  }
  combo_dirs <- combo_dirs[vapply(basename(combo_dirs), is_pure_tissue, logical(1))]

  combo_points <- list()
  for (cd in combo_dirs) {
    combo_name <- basename(cd)
    parts <- strsplit(combo_name, "_")[[1]]
    tissue <- if (all(startsWith(parts, "B"))) "Brain" else if (all(startsWith(parts, "K"))) "Kidney" else "Mixed"
    if (tissue == "Mixed") next
    samples_n <- length(parts)
    cf <- find_classification_file(cd)
    classification <- read_tsv_safe(cf)
    m <- compute_metrics_from_classification(classification)
    mv <- m[[metrics_key]]
    df <- tibble(
      Tissue = tissue,
      Samples = samples_n,
      Combo = combo_name,
      Metric = names(mv),
      Value = as.numeric(mv)
    )
    combo_points[[length(combo_points) + 1]] <- df
  }
  combo_points_df <- dplyr::bind_rows(combo_points)

  bind_rows(single_points_df, combo_points_df)
}

summarize_points_df <- function(points_df) {
  points_df %>%
    group_by(Tissue, Samples, Metric) %>%
    summarise(
      N = sum(!is.na(Value)),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      SE = ifelse(N > 0, SD / sqrt(pmax(N, 1)), NA_real_),
      CI_t = ifelse(N > 1, qt(0.975, df = N - 1), NA_real_),
      CI_lower = ifelse(!is.na(CI_t), Mean - CI_t * SE, NA_real_),
      CI_upper = ifelse(!is.na(CI_t), Mean + CI_t * SE, NA_real_),
      SamplesLabel = as.character(Samples)
    )
}

# ------------------------------------------------------------
# Figure 5b: compute metrics for ALL combinations (Brain and Kidney)
# ------------------------------------------------------------
metric_order <- c("Sensitivity", "Non-redundant Precision", "Inv. Redundancy", "1 - FDR", "PDR", "F1")
metric_order_sirv <- c("Sensitivity", "Non-redundant Precision", "Inv. Redundancy", "1 - FDR", "PDR", "F1")

if (!local_only) {
  # Compute TUSCO points and summaries using shared helpers
  all_points_df <- collect_points_for_metrics("metrics_tusco")
  all_points_df <- all_points_df %>%
    filter(Tissue %in% c("Brain", "Kidney")) %>%
    mutate(
      Tissue = factor(Tissue, levels = c("Brain", "Kidney")),
      Metric = factor(Metric, levels = metric_order),
      SamplesLabel = as.character(Samples)
    )

  summary_df <- summarize_points_df(all_points_df)
} else {
  # LOCAL_ONLY mode: read precomputed data from ./plot/
  points_5b_path <- file.path(plot_dir, "fig-5b_points.tsv")
  bars_5b_path   <- file.path(plot_dir, "fig-5b_bars.tsv")
  points_5c_path <- file.path(plot_dir, "fig-5c_points.tsv")

  if (!file.exists(points_5b_path) || !file.exists(bars_5b_path) || !file.exists(points_5c_path)) {
    stop("LOCAL_ONLY set but required TSVs are missing in ./plot/. Expected: ",
         basename(points_5b_path), ", ", basename(bars_5b_path), ", ", basename(points_5c_path))
  }
  all_points_df <- readr::read_tsv(points_5b_path, show_col_types = FALSE) %>%
    mutate(
      Tissue = factor(Tissue, levels = c("Brain", "Kidney")),
      Metric = factor(Metric, levels = metric_order),
      SamplesLabel = as.character(Samples)
    )
  summary_df <- readr::read_tsv(bars_5b_path, show_col_types = FALSE) %>%
    mutate(SamplesLabel = as.character(Samples))
}

# ------------------------------------------------------------
# Figure S6: compute SIRV metrics for ALL combinations (Brain and Kidney)
# ------------------------------------------------------------
points_s6_path <- file.path(plot_dir, "fig-s6_points.tsv")
bars_s6_path   <- file.path(plot_dir, "fig-s6_bars.tsv")
points_s6_tsv  <- file.path(tsv_dir,  "fig-s6_points.tsv")
bars_s6_tsv    <- file.path(tsv_dir,  "fig-s6_bars.tsv")

if (!local_only) {
  if (is.null(sirv_file) || is.na(sirv_file) || !file.exists(sirv_file)) {
    stop("SIRV GTF not found. Checked spike-ins and reference locations under ", data_root)
  }

  # Compute SIRV points and summaries using shared helpers
  all_points_sirv_df <- collect_points_for_metrics("metrics_sirv") %>%
    filter(Tissue %in% c("Brain", "Kidney")) %>%
    mutate(
      Tissue = factor(Tissue, levels = c("Brain", "Kidney")),
      Metric = factor(Metric, levels = metric_order_sirv),
      SamplesLabel = as.character(Samples)
    )

  summary_sirv_df <- summarize_points_df(all_points_sirv_df)

  # Save per-figure TSVs to support LOCAL_ONLY reruns (in plot) and archive in tsv
  readr::write_tsv(all_points_sirv_df, points_s6_path)
  readr::write_tsv(summary_sirv_df, bars_s6_path)
  readr::write_tsv(all_points_sirv_df, points_s6_tsv)
  readr::write_tsv(summary_sirv_df, bars_s6_tsv)
} else {
  # LOCAL_ONLY: read precomputed SIRV data
  if (!file.exists(points_s6_path) || !file.exists(bars_s6_path)) {
    stop("LOCAL_ONLY set but required SIRV TSVs are missing in ./plot/. Expected: ",
         basename(points_s6_path), ", ", basename(bars_s6_path))
  }
  all_points_sirv_df <- readr::read_tsv(points_s6_path, show_col_types = FALSE) %>%
    mutate(
      Tissue = factor(Tissue, levels = c("Brain", "Kidney")),
      Metric = factor(Metric, levels = metric_order_sirv),
      SamplesLabel = as.character(Samples)
    )
  summary_sirv_df <- readr::read_tsv(bars_s6_path, show_col_types = FALSE) %>%
    mutate(SamplesLabel = as.character(Samples))
}

create_metric_barplot_all <- function(metric_name, points_df, summary_df, color_value = "#2ca02c") {
  metric_points <- points_df %>% filter(Metric == metric_name)
  metric_summary <- summary_df %>% filter(Metric == metric_name)

  # Ensure x order 1..5 even if some are missing
  x_levels <- as.character(sort(unique(points_df$Samples)))
  x_levels <- intersect(as.character(1:5), x_levels)
  metric_points$SamplesLabel <- factor(metric_points$SamplesLabel, levels = as.character(1:5))
  metric_summary$SamplesLabel <- factor(metric_summary$SamplesLabel, levels = as.character(1:5))

  # Compute improvement from 1 to 2 samples per Tissue using mean bars
  # Now compute improvement from 1 to 3 samples per Tissue
  improvement_data <- metric_summary %>%
    filter(Samples %in% c(1, 3), Tissue %in% c("Brain", "Kidney")) %>%
    arrange(Tissue, Samples) %>%
    group_by(Tissue) %>%
    filter(n() == 2) %>%
    summarise(
      improvement = Mean[Samples == 3] - Mean[Samples == 1],
      sample1_value = Mean[Samples == 1],
      sample3_value = Mean[Samples == 3],
      .groups = "drop"
    ) %>%
    mutate(
      improvement_text = ifelse(improvement >= 0, paste0("+", round(improvement, 1), "pp"), paste0(round(improvement, 1), "pp")),
      annotation_y = pmax(sample1_value, sample3_value) + 18,
      bracket_y = pmax(sample1_value, sample3_value) + 8,
      x_start = factor("1", levels = as.character(1:5)),
      x_end   = factor("3", levels = as.character(1:5)),
      x_mid   = factor("2", levels = as.character(1:5))
    )

  ggplot() +
    # Bars: mean
    geom_bar(
      data = metric_summary,
      aes(x = SamplesLabel, y = Mean, fill = Tissue),
      stat = "identity", position = position_dodge2(width = 0.75, preserve = "single"),
      width = 0.7, color = "white", linewidth = 0.3
    ) +
    # Error bars: mean +/- SD
    geom_errorbar(
      data = metric_summary,
      aes(x = SamplesLabel, ymin = pmax(0, CI_lower), ymax = pmin(130, CI_upper), group = Tissue),
      position = position_dodge2(width = 0.75, preserve = "single"),
      width = 0.5, linewidth = 0.5, color = "black"
    ) +
    # Points: each combination
    geom_point(
      data = metric_points,
      aes(x = SamplesLabel, y = Value, group = Tissue),
      position = position_jitter(width = 0.08, height = 0, seed = 1),
      size = 0.35, alpha = 0.45, color = "#1f1f1f"
    ) +
    # Improvement label (1 to 3 samples) per facet
    geom_text(
      data = improvement_data,
      aes(x = x_mid, y = annotation_y, label = improvement_text),
      inherit.aes = FALSE,
      hjust = 0.5, vjust = 0.5, size = 7/.pt, fontface = "bold", color = "#d62728"
    ) +
    # Bracket under label to indicate covered samples (1..3)
    geom_segment(
      data = improvement_data,
      aes(x = x_start, xend = x_end, y = bracket_y + 2, yend = bracket_y + 2, group = Tissue),
      inherit.aes = FALSE, color = "#d62728", linewidth = 0.4
    ) +
    geom_segment(
      data = improvement_data,
      aes(x = x_start, xend = x_start, y = bracket_y, yend = bracket_y + 2, group = Tissue),
      inherit.aes = FALSE, color = "#d62728", linewidth = 0.4
    ) +
    geom_segment(
      data = improvement_data,
      aes(x = x_end, xend = x_end, y = bracket_y, yend = bracket_y + 2, group = Tissue),
      inherit.aes = FALSE, color = "#d62728", linewidth = 0.4
    ) +
    facet_grid(. ~ Tissue, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c(Brain = alpha(color_value, 0.90), Kidney = alpha(color_value, 0.70))) +
    scale_y_continuous(limits = c(0, 130), labels = scales::percent_format(scale = 1), breaks = seq(0, 100, 20)) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +
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
plots_5b <- lapply(metric_levels, function(m) create_metric_barplot_all(m, all_points_df, summary_df, "#2ca02c"))
fig_5b_panel <- plot_grid(plotlist = plots_5b, ncol = 3, nrow = 2, align = "hv")

# Build SIRV bar panel (Fig S6)
metric_levels_sirv <- metric_order_sirv
plots_s6 <- lapply(metric_levels_sirv, function(m) create_metric_barplot_all(m, all_points_sirv_df, summary_sirv_df, "#cab2d6"))
fig_s6_panel <- plot_grid(plotlist = plots_s6, ncol = 3, nrow = 2, align = "hv")

# ------------------------------------------------------------
# Figure 5c: radar metrics and plotting (unchanged logic: single_sample single-sample)
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
      (structural_category == "full-splice_match" & !is.na(ref_exons) & ref_exons == 1 &
         !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
         abs(diff_to_TSS) <= 50 & abs(diff_to_TTS) <= 50)
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

if (!local_only) {
  # Load inputs for radar metrics
  tusco_universal_tsv <- file.path(data_root, "tusco/tusco_mouse.tsv")
  tusco_brain_tsv     <- file.path(data_root, "tusco/tusco_mouse_brain.tsv")
  tusco_kidney_tsv    <- file.path(data_root, "tusco/tusco_mouse_kidney.tsv")

  all_sample_dirs_radar <- list.dirs(isoseq_ar_root, full.names = TRUE, recursive = FALSE)
  all_sample_dirs_radar <- all_sample_dirs_radar[grepl("/(B|K)[0-9]+\\.isoforms$", all_sample_dirs_radar)]
  if (length(all_sample_dirs_radar) == 0) stop("No sample directories found in ", isoseq_ar_root)

  tusco_ref_univ   <- load_tusco_reference(tusco_universal_tsv)
  tusco_ref_brain  <- load_tusco_reference(tusco_brain_tsv)
  tusco_ref_kidney <- load_tusco_reference(tusco_kidney_tsv)

  results <- list()
  for (sd in all_sample_dirs_radar) {
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
} else {
  # LOCAL_ONLY: read precomputed radar metrics
  metrics_df <- readr::read_tsv(file.path(plot_dir, "fig-5c_points.tsv"), show_col_types = FALSE) %>%
    select(Sample, Tissue, Type, `Sn`, `nrPre`, `1/red`, `1-FDR`, `PDR`, `rPre`)
}

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

metrics_labels_short <- c("Sn", "nrPre", "1/red", "1 - FDR", "PDR", "rPre")

compute_cosine <- function(vec_a, vec_b) {
  if (any(is.na(vec_a)) || any(is.na(vec_b))) return(NA_real_)
  denom <- sqrt(sum(vec_a^2)) * sqrt(sum(vec_b^2))
  if (denom == 0) return(NA_real_)
  sum(vec_a * vec_b) / denom
}

build_radar_plot_for_tissue <- function(summary_df, tissue_label) {
  rd <- build_radar_data_for_tissue(summary_df, tissue_label)
  cm <- type_colors[rownames(rd)[-(1:2)]]
  # Print detailed cosine similarity between Universal and tissue-specific vectors
  metrics_order <- c("Sn", "nrPre", "1/red", "1-FDR", "PDR", "rPre")
  v_univ <- summary_df %>% filter(Tissue == tissue_label, Type == "Universal") %>%
    select(all_of(metrics_order)) %>% as.numeric()
  v_tiss <- summary_df %>% filter(Tissue == tissue_label, Type == tissue_label) %>%
    select(all_of(metrics_order)) %>% as.numeric()
  if (length(v_univ) == length(metrics_order) && length(v_tiss) == length(metrics_order)) {
    cs <- compute_cosine(v_univ, v_tiss)
    message(sprintf(
      paste0(
        "[Radar] Tissue=%s\n",
        "  Metrics order: %s\n",
        "  Universal: [%s]\n",
        "  %s:        [%s]\n",
        "  Cosine similarity: %.4f\n"
      ),
      tissue_label,
      paste(metrics_order, collapse = ", "),
      paste(round(v_univ, 3), collapse = ", "),
      tissue_label,
      paste(round(v_tiss, 3), collapse = ", "),
      cs
    ))
  }
  ggdraw() +
    draw_grob(radar_grob(rd, cm, var_labels = metrics_labels_short, title = NULL), x = -0.02, y = 0.05, width = 1.04, height = 0.80) +
    draw_label(tissue_label, x = 0.5, y = 0.985, hjust = 0.5, vjust = 1, size = 7, fontface = "bold")
}

radar_brain <- build_radar_plot_for_tissue(metrics_summary_df, "Brain")
radar_kidney <- build_radar_plot_for_tissue(metrics_summary_df, "Kidney")

# Inline legend (overlay on bottom radar to keep radar area to exactly two rows)
legend_df <- data.frame(Type = factor(names(type_colors), levels = names(type_colors)),
                        x = c(1, 2, 3), y = c(1, 2, 3))
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
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = c(0.5, 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.width = unit(10, "pt"),
    legend.key.height = unit(10, "pt"),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  guides(color = guide_legend(nrow = 1), linetype = "none", shape = "none")

legend_only <- NULL
spacer_bottom <- ggdraw()

# ------------------------------------------------------------
# Combine panels
# ------------------------------------------------------------
fig_5c_panel <- plot_grid(
  radar_brain,
  radar_kidney,
  spacer_bottom,
  ncol = 1,
  rel_heights = c(1, 1, 0.22)
)
fig_5c_panel <- fig_5c_panel + theme(plot.margin = margin(0, 0, 0, 0))

combined <- plot_grid(
  fig_5b_panel, fig_5c_panel,
  ncol = 2,
  rel_widths = c(3.0, 1.0),
  align = "hv",
  axis = "tb",
  greedy = FALSE,
  labels = c("b", "c"),
  label_size = 10,
  label_fontface = "bold"
)

# ------------------------------------------------------------
# Export underlying data used for this PDF as a single TSV
# ------------------------------------------------------------
underlying_out <- file.path(tsv_dir, "fig-5b-5c.tsv")

all_points_out <- all_points_df %>%
  arrange(Tissue, Samples, Combo, Metric) %>%
  mutate(figure_id = "fig-5", panel_id = "5b_points")

bar_means_out <- summary_df %>%
  arrange(Tissue, Samples, Metric) %>%
  select(Tissue, Samples, Metric, N, Mean, SD, SE, CI_lower, CI_upper) %>%
  mutate(figure_id = "fig-5", panel_id = "5b_bars")

metrics_points_out <- metrics_df %>%
  arrange(Tissue, Type, Sample) %>%
  mutate(figure_id = "fig-5", panel_id = "5c_points")

metrics_summary_out <- metrics_summary_df %>%
  mutate(figure_id = "fig-5", panel_id = "5c_summary")

export_df <- bind_rows(
  all_points_out %>% mutate(dataset = "5b_points"),
  bar_means_out %>% mutate(dataset = "5b_bars"),
  metrics_points_out %>% mutate(dataset = "5c_points"),
  metrics_summary_out %>% mutate(dataset = "5c_summary")
)
readr::write_tsv(export_df, underlying_out)
message("Wrote underlying data to ", underlying_out)

pdf(out_pdf, width = 7.09, height = 3.55)
print(combined)
dev.off()

message("Wrote combined figure to ", out_pdf)

# ------------------------------------------------------------
# Export SIRV underlying data and figure S6
# ------------------------------------------------------------
sirv_underlying_out <- file.path(tsv_dir, "fig-s6.tsv")

sirv_all_points_out <- all_points_sirv_df %>%
  arrange(Tissue, Samples, Combo, Metric) %>%
  mutate(figure_id = "fig-s6", panel_id = "s6_points")

sirv_bar_means_out <- summary_sirv_df %>%
  arrange(Tissue, Samples, Metric) %>%
  select(Tissue, Samples, Metric, N, Mean, SD, SE, CI_lower, CI_upper) %>%
  mutate(figure_id = "fig-s6", panel_id = "s6_bars")

sirv_export_df <- bind_rows(
  sirv_all_points_out %>% mutate(dataset = "s6_points"),
  sirv_bar_means_out %>% mutate(dataset = "s6_bars")
)
readr::write_tsv(sirv_export_df, sirv_underlying_out)
message("Wrote SIRV underlying data to ", sirv_underlying_out)

s6_out_pdf <- file.path(plot_dir, "fig-s6.pdf")
pdf(s6_out_pdf, width = 7.09, height = 3.55)
print(fig_s6_panel)
dev.off()
message("Wrote SIRV spike-in figure to ", s6_out_pdf)
